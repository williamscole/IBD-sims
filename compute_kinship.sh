Ne=${1}
n=${2}

for mode in pedigree coalescent
    do
    for iter in {1..50}
        do
        ibdfile="${mode}/iter${iter}_Ne${Ne}_n${n}_${mode}"

        # sbatch -t 5:00 --mem 2gb --wrap "python kinship_counts.py ${ibdfile}.ibd > ${ibdfile}.kinship "
        sbatch -t 5:00 --mem 2gb --wrap "python kinship.py ${ibdfile}.ibd"

    done
done


# for chrom in {1..30}
# do
# zcat ibd/${prefix}_chr${chrom}.ibd.gz | awk -v chrom="$chrom" '{ $5 = chrom; print }' >> ibdne/${prefix}.ibd
# cat maps/${prefix}_chr${chrom}.map | awk -v chrom="$chrom" '{ $1 = chrom; print}' >> ibdne/${prefix}.map
# done

# sbatch -t 1:00:00 --mem 32gb --ntasks 1 --cpus-per-task 32 --wrap "cat ibdne/${prefix}.ibd | java -jar -Xmx28g $ibdne map=ibdne/${prefix}.map nthreads=32 out=ibdne/${prefix} mincm=2 gmax=100 gmin=1"
