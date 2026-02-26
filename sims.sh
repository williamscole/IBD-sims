Ne=${1}
n=${2}
monog=${3}

# mode="coalescent"

# for iter in {1..50}
#     do
#     output=${mode}/iter${iter}
#     for chrom in {1..30}
#         do

#         output=${mode}/iter${iter}

#         njobs=$(myq | grep batch | wc -l)

#         while [[ $njobs -gt  500 ]]
#             do
#             njobs=$(myq | grep batch | wc -l)
#             done

#         sbatch -t 30:00 --mem 16g --wrap "bash single_sim.sh $mode $Ne $n $chrom $output"
#         echo "bash single_sim.sh $mode $Ne $n $chrom $output submitted!"

#     done
# done

mode="pedigree"
gens=30

if [[ $monog == "monog" ]]
then
monog="--monogamous"
prefix="monogamous"
else
prefix="pedigree"
fi

for iter in 1 # {2..50}
    do
    python wf_pedigree.py --constant_Ne $Ne --samples $n --gens $gens --output ${prefix}/iter${iter} $monog

    pedigree="${prefix}/iter${iter}_Ne${Ne}_n${n}_g${gens}.pedigree"
    output=${prefix}/iter${iter}

    for chrom in {1..30}
        do

        njobs=$(myq | grep batch | wc -l)

        while [[ $njobs -gt  500 ]]
            do
            njobs=$(myq | grep batch | wc -l)
            done

        sbatch -t 30:00 --mem 16g --wrap "bash single_sim.sh $mode $Ne $n $chrom $output $pedigree $gens"
        echo "bash single_sim.sh $mode $Ne $n $chrom $output $pedigree $gens submitted!"

        done
    done






# pedigree=test/iter${iter}_Ne${Ne}_n${n}_g25.pedigree
# dend=25

# bash single_sim.sh $mode $Ne $n $chr $output $pedigree $dend

# mode=${1}
# Ne=${2}
# n=${3}
# chr=${4}
# output=${5}
# pedigree=${6}
# dend=${7}