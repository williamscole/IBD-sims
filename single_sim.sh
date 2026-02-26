mode=${1}
Ne=${2}
n=${3}
chr=${4}
output=${5}
pedigree=${6}
dend=${7}

if [[ $mode == "coalescent" ]]
then
python msprime_scripts.py --mode $mode --constant_Ne $Ne --samples $n --chr $chr --output $output
else
python msprime_scripts.py --mode $mode --constant_Ne $Ne --samples $n --chr $chr --output $output --pedigree $pedigree --dtwf_end $dend
fi

prefix=${output}_Ne${Ne}_n${n}_${mode}_chr${chr}

echo ${prefix}.vcf.gz

java -jar -Xmx12g /users/cwilli50/hap-ibd.jar gt=${prefix}.vcf.gz map=${prefix}.map out=${prefix}
