
# copy the *.gz files to the local scratch
gz_in_R1=$3
gz_in_R2=${gz_in_R1/_R1_/_R2_}

cp $2/$gz_in_R1 $TMPDIR
cp $2/$gz_in_R2 $TMPDIR
cp reference.fasta $TMPDIR

# target directory
file_dir=$2/$1

# go to the local scratch
cd $TMPDIR


gz_loc_in_R1_tmp=*_R1_001.fastq.gz
gz_loc_in_R2_tmp=*_R2_001.fastq.gz

for f in $gz_loc_in_R1_tmp; do gz_loc_in_R1=${f};done
for f in $gz_loc_in_R2_tmp; do gz_loc_in_R2=${f};done

# unzip forward and backward reads
gunzip $gz_loc_in_R1
gunzip $gz_loc_in_R2

# outnames
fastq_R1=${gz_loc_in_R1/.fastq.gz/.fastq}
fastq_R2=${gz_loc_in_R2/.fastq.gz/.fastq}


# trim the sequences using trimmmomatic
trimmed_R1=${fastq_R1/.fastq/.trimmed.fastq}
trimmed_R2=${fastq_R2/.fastq/.trimmed.fastq}

echo $fastq_R1
echo $fastq_R2


echo $trimmed_R1
echo $trimmed_R2


# build the name for the vcf file
vcf_file=${gz_loc_in_R1/_L001_R1_001.fastq.gz/.vcf}
echo $vcf_file


ls


echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "trim sequences"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"

# ILLUMINACLIP in case adapter are used to sequence


java -jar ~/jar/trimmomatic-0.36.jar SE $fastq_R1 $trimmed_R1 ILLUMINACLIP:/cluster/home/nicmuell/jar/adapters/NexteraPE-PE.fa:2:30:10: LEADING:3 TRAILING:3  SLIDINGWINDOW:5:30 MINLEN:100
java -jar ~/jar/trimmomatic-0.36.jar SE $fastq_R2 $trimmed_R2 ILLUMINACLIP:/cluster/home/nicmuell/jar/adapters/NexteraPE-PE.fa:2:30:10: LEADING:3 TRAILING:3  SLIDINGWINDOW:5:30 MINLEN:100

echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "build reference"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"


# build the reference
bowtie2-build reference.fasta ref

echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "bowtie align"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"


# do teh mapping
bowtie2 -x ref -U $trimmed_R1,$trimmed_R2 -S sample.sam --local

echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "samtools view"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"

# conversion
samtools view -bS sample.sam > sample.bam


echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "samtools sort"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"

samtools sort sample.bam ${vcf_file/.vcf/.aligned}



echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "samtools depth"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"

samtools depth ${vcf_file/.vcf/.aligned}.bam > ${vcf_file/.vcf/.log}


echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "lofreq variant calling"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"


# filtering
~/lofreq_star-2.1.2/bin/lofreq call  --ref reference.fasta --out sample.vcf.lofreq ${vcf_file/.vcf/.aligned}.bam

echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "lofreq filtering part 2"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"


~/lofreq_star-2.1.2/bin/lofreq filter --cov-min 100 --snvqual-thresh 30 --af-min 0.01 -i sample.vcf.lofreq -o $vcf_file

echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "copy files back to local dir"
echo "------------------------------------------------"
echo "------------------------------------------------"
echo "------------------------------------------------"


# copy the beast output to the folder with the script again
cp *.vcf $file_dir/
cp *.log $file_dir/
cp *.aligned.bam $file_dir/
