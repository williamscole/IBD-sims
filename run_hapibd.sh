prefix=${1}

for i in {1..30}
do
java -jar -Xmx12g $hapibd gt=vcfs/${prefix}_chr${i}.vcf.gz map=maps/${prefix}_chr${i}.map out=ibd/${prefix}_chr${i}
done
