import subprocess

def get_job_info(jobid):
    # Run the sacct command
    command = ["sacct", "-j", jobid, "--format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode,MaxRSS"]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Parse the output
    lines = result.stdout.splitlines()
    lines = [i for i in lines if "batch" in i.split()[1]]

    if len(lines) == 0:
        return "SUBMITTED", 0
    
    line = lines[0].split()

    try:
        gb = float(line[6].replace("K", ".0")) / 1000000
    except:
        gb = 0

    return line[4], gb


# Example usage
jobid = "8512722"

print(get_job_info(jobid))
