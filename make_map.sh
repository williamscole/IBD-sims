prefix=${1}
chrom=${2}

bcftools query -f "%POS\n" vcfs/${prefix}_chr${chrom}.vcf.gz | head | awk -v chrom="$chrom" '{printf "%s rs%d %.6f %'\''d\n", 1, NR, $1/1000000, $1}' > maps/${prefix}_chr${chrom}.map
