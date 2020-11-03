#Specify paths
#####

#To BAM fi
GENOME=~/Genomes/Onchocerca_volvulus/ncbi/oncho_with_wolb.fna
TRIMMODIR=~/software/Trimmomatic-0.39/
THREADS=6

#####
#Prepare and index genome
#####

echo "bwa indexing genome...."

#bwa index $GENOME

echo "samtools indexing genome..."

#samtools faidx $GENOME

#start loop over sample list

echo "starting loop...."

for f in `cat list`

do

    ######
    #trim reads
    ######

    java -jar $TRIMMODIR/trimmomatic-0.39.jar PE fastq/untrim/$f"_R1_001.fastq.gz" fastq/untrim/$f"_R2_001.fastq.gz" fastq/trim/$f"_R1.trim.fq.gz" fastq/trim/$f"_R1.trim.unpaired.fq.gz" fastq/trim/$f"_R2.trim.fq.gz" fastq/trim/$f"_R2.trim.unpaired.fq.gz" ILLUMINACLIP:$TRIMMODIR/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

    ######
    #align reads
    ######

    bwa mem -t $THREADS $GENOME fastq/trim/$f"_R1.trim.fq.gz" fastq/trim/$f"_R2.trim.fq.gz" | samtools sort > bam/$f".srt.bam"

    ######
    #remove duplicates
    ######

    java -jar ~/software/picard/build/libs/picard.jar MarkDuplicates I=bam/$f".srt.bam" O=bam/$f".srt.dd.bam" M=bam/$f".M"

done
