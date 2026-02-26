path=${1}
iter=${2}

retrieve_arg() {
    awk -v search="$1:" '$1 == search {print $2}' "${path}/args.yaml"
}

end_chr=$(retrieve_arg "end_chr")

prefix="${path}/iter${iter}"


# Target files:
# {path}/iter{iter}.ibd.gz
# {path}/iter{iter}.map
# {path}/iter{iter}.tmrca.gz
# {path}/iter{iter}.npy

purple_file="${prefix}.npy"
tmrca_file="${prefix}.tmrca.gz"
ibd_file="${prefix}.ibd"
map_file="${prefix}.map"

# Purple node matrix: purple_file
if [[ ! -f $purple_file ]]
then
python purple.py "${prefix}_chr1.ibd.gz" ${end_chr} $purple_file
fi

# Concatenate the TMRCAs: tmrca_file
if [[ ! -f $tmrca_file ]]
then
    # Run the concatenation script
    python concat_tmrca.py ${path} ${iter}
    
    # Check if concatenation was successful (exit code 0)
    if [ $? -eq 0 ]; then
        # Delete all intermediate chromosome files
        echo "Concatenation successful. Removing intermediate files..."
        rm ${path}/iter${iter}_chr*.tmrca.pkl
    else
        echo "Concatenation failed. Keeping intermediate files for debugging."
        exit 1
    fi
else
    echo "Final file ${tmrca_file} already exists. Skipping concatenation."
fi

# Delete and concat files: ibd_file
if [[ ! -f "${ibd_file}.gz" ]]
then
for chrom in $(seq 1 $end_chr )
    do
    zcat "${prefix}_chr${chrom}.ibd.gz" | awk -v chrom="$chrom" '{ $5 = chrom; print }' >> $ibd_file
    rm "${prefix}_chr${chrom}.ibd.gz"
    done
fi

# Delete and concat files: map_file
if [[ ! -f $map_file ]]
then
for chrom in $(seq 1 $end_chr )
    do
    cat "${prefix}_chr${chrom}.map" | awk -v chrom="$chrom" '{ $1 = chrom; print}' >> $map_file
    rm "${prefix}_chr${chrom}.map"
    done
fi


if [[ -f ${prefix}.ibd ]]
    then
    gzip ${prefix}.ibd
fi