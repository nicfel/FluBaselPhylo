# get all the forward raw reads
d=$1

new_vcf=${d/rawreads/vcf}
new_out=${d/rawreads/out}



mkdir $new_vcf
mkdir $new_out


for f in $d*_R1_001.fastq.gz; do

outfile_tmp=${f/rawreads/out}
outfile=${outfile_tmp/.fastq.gz/.out}

bsub -J getvcf -W 24:00 -R "rusage[mem=8000]" -o $outfile ./getVCF.sh $new_vcf $PWD $f

done


