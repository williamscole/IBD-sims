#!/bin/bash
set -e  # Exit on any error

path=${1}
iter=${2}

retrieve_arg() {
    awk -v search="$1:" '$1 == search {print $2}' "${path}/args.yaml"
}

end_chr=$(retrieve_arg "end_chr")
gb=$(retrieve_arg "gb")
ibdne_time=$(retrieve_arg "ibdne_time")

# Function to format elapsed time
format_time() {
    local seconds=$1
    local minutes=$((seconds/60))
    local remaining_seconds=$((seconds%60))
    echo "${minutes}m ${remaining_seconds}s"
}

for chrom in $(seq 1 $end_chr); do
    echo "Starting chromosome ${chrom} simulation..."
    start_time=$(date +%s)
    
    if ! python simulations.py $path $iter $chrom; then
        elapsed=$(( $(date +%s) - start_time ))
        formatted_time=$(format_time $elapsed)
        echo "ERROR: Simulation failed for chromosome ${chrom} after ${formatted_time}"
        echo "Stopping remaining simulations"
        exit 1
    fi
    
    elapsed=$(( $(date +%s) - start_time ))
    formatted_time=$(format_time $elapsed)
    echo "Completed chromosome ${chrom} in ${formatted_time}"
done

echo "All chromosome simulations completed successfully"

bash post_sim.sh ${path} ${iter}

# Submit the post-processing
jobid=$(sbatch $ibdne_time --mem $((2*gb))g --ntasks 1 \
    --cpus-per-task $(retrieve_arg "nthreads") \
    --output=${path}/slurm/iter${iter}_postprocess_%j.out \
    --partition=batch \
    --wrap "bash post_process.sh $path $iter")

echo "Submitted post-processing! Job no. ${jobid}"