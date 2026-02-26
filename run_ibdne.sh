path=${1}
iter=${2}

retrieve_arg() {
    awk -v search="$1:" '$1 == search {print $2}' "${path}/args.yaml"
}

end_chr=$(retrieve_arg "end_chr")

prefix="${path}/iter${iter}"

if [[ -f ${prefix}.ibd ]]
    then
    echo "Concatenated IBD file exists!"
else
    for chrom in $(seq 1 $end_chr )
        do
        zcat "${prefix}_chr${chrom}.ibd.gz" | awk -v chrom="$chrom" '{ $5 = chrom; print }' >> ${prefix}.ibd
        cat "${prefix}_chr${chrom}.map" | awk -v chrom="$chrom" '{ $1 = chrom; print}' >> ${prefix}.map
        done
fi

gb=$(retrieve_arg "gb" | awk '{printf "%d", $1 * 0.8}')

cat ${prefix}.ibd | java -jar -Xmx${gb}g $(retrieve_arg "ibdne") \
map=${prefix}.map \
out=${prefix} \
nthreads=$(retrieve_arg "nthreads") \
filtersamples=$(retrieve_arg "filtersamples") \
npairs=$(retrieve_arg "npairs") \
nits=$(retrieve_arg "nits") \
nboots=$(retrieve_arg "nboots") \
mincm=$(retrieve_arg "mincm") \
trimcm=$(retrieve_arg "trimcm") \
gmin=$(retrieve_arg "gmin") \
gmax=$(retrieve_arg "gmax")

if [[ -f ${prefix}.ne ]]
    then
    python plot.py ${path}
    echo "Success!"
    if [[ -f "${path}/iter${iter}_chr1.ibd.gz" ]]
        then
        rm "${path}"/iter${iter}_chr*.ibd.gz
    fi
else
    echo "Job failed: ${prefix}.ne" > ${path}/errors/iter${iter}.err
    fi
