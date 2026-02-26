import yaml
import datetime
import os
import pandas as pd
import time
import sys
import numpy as np
import subprocess

def cur_status(jobid):

    # Run sacct command
    cmd = f"sacct -j {jobid} --format=JobID,State"
    result = subprocess.run(cmd.split(), capture_output=True, text=True)

    # Parse output and get states from lines after header
    lines = result.stdout.strip().split('\n')[2:]  # Skip header lines
    states = {line.split()[1] for line in lines if len(line.split()) >= 2}

    if states == {"COMPLETED"}:
        return "COMPLETED"

    if len(states) == 0:
        return "SUBMITTED"

    if "RUNNING" in states:
        return "RUNNING"

    return "FAILED"

def resubmit_array_job(array_num, path, gb, sim_time, n_iter, end_chr, max_concurrent=100):
   cmd = [
       "sbatch",
       f"--array={array_num}%{max_concurrent}",  # Just submit specific array number
       "--parsable",
       f"--mem={gb}gb",
       sim_time,
       "--partition=batch", 
       f"--output={path}/slurm/simulation_array%a_%j.out",
       "--wrap",
       f"bash array_job.sh {path} $SLURM_ARRAY_TASK_ID {n_iter} {end_chr}"
   ]

   # Run sbatch command
   result = subprocess.run(cmd, capture_output=True, text=True)
   
   # Get job ID and append array number
   job_id = result.stdout.strip()
   return job_id

def post_process_job(path, iter_n, gb, ibdne_time, nthreads):
   cmd = [
       "sbatch",
       "--parsable",
       "--ntasks=1",
       f"--cpus-per-task={nthreads}",
       f"--mem={2*gb}gb",
       ibdne_time,  # Already includes --time=
       "--partition=batch", 
       f"--output={path}/slurm/IBDNe_iter{iter_n}_%j.out",
       "--wrap",
       f"bash post_process.sh {path} {iter_n}"
   ]

   result = subprocess.run(cmd, capture_output=True, text=True)
   if not result.stdout:
       print("Error:", result.stderr)
   
   job_id = result.stdout.strip()
   return job_id


def plot_ibdne(path, depend):
   cmd = [
       "sbatch",
       "--parsable",
       "--mem=2gb",
       "--time=10:00",
       "--partition=batch", 
       f"--output={path}/slurm/plot_ibdne_%j.out",
       "--wrap",
       f"python plot.py {path}",
       f"--dependency=afterok:{depend}"
   ]

   result = subprocess.run(cmd, capture_output=True, text=True)
   if not result.stdout:
       print("Error:", result.stderr)
   
   job_id = result.stdout.strip()
   return job_id

def create_status_df(jobid, n_iter, end_chr):
    n_jobs = n_iter * end_chr

    df = pd.DataFrame({"array": np.arange(1, n_jobs+1)})
    df["jobid"] = df["array"].apply(lambda x: str(jobid) + f"_{x}")
    df["iteration"] = df["array"].apply(lambda x: ((x-1) // end_chr)+1)
    df["chromosome"] = df["array"].apply(lambda x: ((x-1) % end_chr)+1)
    df["status"] = df["jobid"].apply(cur_status)
    df["times_failed"] = 0

    return df

def percent_status(df, status="COMPLETED"):
    return (df["status"]==status).mean()

def get_n_from_array(array):
    n = 0
    array = array.split(",")
    for i in array:
        if "-" in i:
            i,j = i.split("-")
            n += (int(j)-int(i)+1)
        else:
            n += 1
    return n

def main():
    path = sys.argv[1]
    args = yaml.safe_load(open(f"{path}/args.yaml"))
    n_iter = args["iter"]
    end_chr = args["end_chr"]
    sim_id = sys.argv[2]
    array = sys.argv[3]

    complete_goal = get_n_from_array(array)
    
    status_df = create_status_df(sim_id, n_iter, end_chr)
    
    # For dynamic counter
    completed_iters = set()

    open(f"{path}/failed_jobs.log", "w").close()

    n_resubmitted = 0  # Track total resubmissions
    post_process_jobids = {}

    while True:
        # Update the status
        status_df["status"] = status_df["jobid"].apply(cur_status)

        # Go through the iterations
        for iter_n, iter_df in status_df.groupby("iteration"):
            if percent_status(iter_df) == 1 and iter_n not in completed_iters:
                jobid = post_process_job(path, iter_n, args["gb"], args["ibdne_time"], args["nthreads"])
                post_process_jobids[iter_n] = jobid
                completed_iters.add(iter_n)

        # Handle failed jobs
        failed_df = status_df[status_df.status=="FAILED"]
        if failed_df.shape[0] > 0:  # Only submit if there are failed jobs
            failed_jobs = failed_df["array"].apply(str).values
            new_jobid = resubmit_array_job(",".join(failed_jobs), path, args["gb"], args["sim_time"], n_iter, end_chr)

            # Log all resubmissions
            with open(f"{path}/failed_jobs.log", "a") as f:
                for index, row in failed_df.iterrows():
                    old_jobid = row["jobid"]
                    array_specific_jobid = new_jobid + f"_{row['array']}"
                    status_df.at[index,"jobid"] = array_specific_jobid
                    status_df.at[index,"times_failed"] += 1
                    status_df.at[index, "status"] = "SUBMITTED"

                    n_resubmitted += 1
                    
                    # Log: old_jobid, new_jobid, iteration, chromosome
                    f.write(f"{old_jobid}\t{array_specific_jobid}\t{row['iteration']}\t{row['chromosome']}\n")

        n_completed = len(status_df[status_df.status=="COMPLETED"])
        n_running = len(status_df[status_df.status=="RUNNING"])
        print(f"Completed: {n_completed}/{complete_goal} | Running: {n_running} | Resubmitted: {n_resubmitted} | Completed iterations: {len(completed_iters)}", end='\r')
            

        if percent_status(status_df) == 1 or max(status_df["times_failed"].values) > 4 or complete_goal == n_completed:
            break

        if percent_status(status_df, "FAILED") == 1:
            print("All jobs failed!")
            break

        time.sleep(30)

    # Clear the dynamic counter line at the end
    print()

    status_df["status"] = status_df["jobid"].apply(cur_status)
    status_df.to_csv(f"{path}/status.tsv", sep="\t", index=False)

    for iter_n, jobid in post_process_jobids.items():
        print(f"Iteration {iter_n}: jobid {jobid}")

    if args["run_ibdne"]:
        plot_id = plot_ibdne(path, ",".join([str(i) for i in post_process_jobids.values()]))
        print(f"Submitted plotting: {plot_id}")



if __name__ == "__main__":

    main()


    # if input_arg.endswith(".yaml"):

    #     path = initialize_sim(yaml_file)

    # elif os.path.exists(f"{input_arg}/args.yaml"):

    #     path = input_arg

    # args = yaml.safe_load(open(f"{path}/args.yaml"))
    # n_iter = args["n_iter"]
    # end_chr = args["end_chr"]

    # if not args["post_process_only"]:
    #     submit_2d_array(n_iter, end_chr)
    #     print("Submitted 2d slurm array!")
    #     time.sleep(5)
    # else:
    #     print("Post process only!")

    # # Open the status dataframe
    # status_df = load_status(path)
    # while True:
    #     # Update the current status of the job
    #     status_df = check_status(status_df)
    #     # Resubmit any jobs that have failed
    #     status_df = resubmit_jobs(status_df)
        
    #     # Check to see if all jobs within an iteration have completed
    #     for iter_n, iter_df in status_df.groupby("iteration"):
    #         percent_completed = (iter_df["status"]=="completed").mean()

    #         if percent_completed == 1:
    #             submit_post_process(iter_n)

    #     if (status_df["status"]=="completed").mean() == 1:
    #         break




    



    #     args = yaml.safe_load(open(yaml_file))



    #     path = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S.%f")

    #     os.makedirs(timestamp, exist_ok=True)

    #     shutil.copy2(input_arg, f'{path}/args.yaml')

    #     yaml_file = f'{path}/args.yaml'
 
    # elif os.path.exists(f"{input_arg}/args.yaml"):

    #     yaml_file = f"{input_arg}/args.yaml"
    
    # args = yaml.safe_load(open(yaml_file))


