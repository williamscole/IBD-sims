#!/bin/bash

args=${1}

if [[ -f "$args/args.yaml" ]]
    then
    path=${args}
else
    path="$(date +%Y-%m-%d_%H-%M-%S.%N)"
    mkdir $path
    mkdir $path/slurm
    mkdir $path/errors
    cat $args > ${path}/args.yaml
    python yaml_tools.py $path
    fi

echo "All results will be found in ${path}"

retrieve_arg() {
    awk -v search="$1:" '$1 == search {print $2}' "${path}/args.yaml"
}


n_iter=$(retrieve_arg "iter")
end_chr=$(retrieve_arg "end_chr")
gb=$(retrieve_arg "gb")
post_process_only=$(retrieve_arg "post_process_only")
simulate_only=$(retrieve_arg "simulate_only")
sim_time=$(retrieve_arg "sim_time")
ibdne_time=$(retrieve_arg "ibdne_time")
pedigree_mode=$(retrieve_arg "pedigree_mode")


if [[ $post_process_only == "true" ]]
then
    echo "Running post-processing only!"
    # Submit IBDNe jobs as an array
    sbatch --array=1-${n_iter} $ibdne_time --mem $((2*gb))g --ntasks 1 \
        --cpus-per-task $(retrieve_arg "nthreads") \
        --output=${path}/slurm/IBDNe_iter%a_%j.out \
        --wrap "bash post_process.sh ${path} \$SLURM_ARRAY_TASK_ID"
    exit 0
fi


# Create WF pedigree if needed
if [[ $pedigree_mode == "true" ]]
    then
    ped_job=$(sbatch \
        --array=1-${n_iter} \
        --parsable \
        --mem 6gb \
        --time 20:00 \
        --partition=batch \
        --output=${path}/slurm/create_pedigree_iter%a_%j.out \
        --wrap "python simulations.py ${path} \$SLURM_ARRAY_TASK_ID 0"
    )
    depend="--dependency=afterok:${ped_job}"
else
    depend=""
fi

echo $ped_job

tot_array=$(( $end_chr * $n_iter ))

# array=$(python check_completed.py ${path} ${tot_array} | tail -1)
for array in $(python check_completed.py ${path} ${tot_array})
    do
    echo "Running array $array"
    sim_job=$(sbatch \
        --array=${array}%100 \
        --parsable \
        --mem ${gb}gb \
        ${sim_time} \
        --partition=batch \
        --output=${path}/slurm/simulation_array%a_%j.out \
        --wrap "bash array_job.sh ${path} \$SLURM_ARRAY_TASK_ID $n_iter $end_chr" \
        $depend
    )

    echo $sim_job

    python monitor_simulation.py $path $sim_job $array
    done