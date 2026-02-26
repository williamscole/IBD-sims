path=${1}
iter=${2}

retrieve_arg() {
    awk -v search="$1:" '$1 == search {print $2}' "${path}/args.yaml"
}

end_chr=$(retrieve_arg "end_chr")
run_ibdne=$(retrieve_arg "run_ibdne")
purple_nodes=$(retrieve_arg "purple_nodes")
dir_name=$(retrieve_arg "dir_name")
ibd_filter=$(retrieve_arg "ibd_filter")
n_samples=$(retrieve_arg "samples")

mkdir -p ${path}/${dir_name}

echo "Writing to ${path}/${dir_name}" 

cat ${path}/args.yaml > ${path}/${dir_name}/args.yaml

prefix="${path}/iter${iter}"

bash post_sim.sh $path $iter

if [[ $run_ibdne == "true" ]]
then
    echo "Running IBDNe!"
    gb=$(retrieve_arg "gb" | awk '{printf "%d", $1 * 1.8}')

    python filter_ibd.py ${prefix}.ibd.gz $n_samples $ibd_filter  | java -jar -Xmx${gb}g $(retrieve_arg "ibdne") \
    map=${prefix}.map \
    out=${path}/${dir_name}/iter${iter} \
    nthreads=$(retrieve_arg "nthreads") \
    filtersamples=$(retrieve_arg "filtersamples") \
    npairs=$(retrieve_arg "npairs") \
    nits=$(retrieve_arg "nits") \
    nboots=$(retrieve_arg "nboots") \
    mincm=$(retrieve_arg "mincm") \
    trimcm=$(retrieve_arg "trimcm") \
    gmin=$(retrieve_arg "gmin") \
    gmax=$(retrieve_arg "gmax")
fi