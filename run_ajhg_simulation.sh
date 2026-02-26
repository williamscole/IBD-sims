script=${1}

for i in {1..30}
do
sbatch -t 30:00 --mem 8gb --wrap "python $script $i"
done

