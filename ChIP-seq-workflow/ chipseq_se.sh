#CHIP-SEQ Workflow for homo sapiens (single-end reads)

#the absolute location of the files (this will be the reference directory)
location=/home/belab/Desktop/chipseq/minimap2

#samples are stored in a variable (to be iterated) 
samples=$location/SRR/*.fastq

#this part should be changed if your data (.fastq files) do not follow order of sample_annot array
declare -i i=0
sample_annot=(bcl6_r1 bcl6_r1 bcl6_ctrl)

#This 'for' loop iterates over the fastq files and processes the data through minimap2, samtools and macs3
for sample in $samples
do
	annot=${sample_annot[$i]} #now the annotation is in annot variable
	
	echo $sample is annotated as $annot #prints the fastq file in progress
	
	set -o pipefail
	
	#minimap2 is an aligner and its sam output will be processed through samtools to mark
	#secondary reads (for that it is first 'fixed' and sorted)
	
	minimap2 -ax sr $location/human_reference/Homo_sapiens.GRCh38.dna.primary_assembly.mmi $sample --secondary no > $location/results/$annot.sam 
	
	#could not pipe minimap2 to samtool, doesn't work
	samtools fixmate -u -m -r $location/results/$annot.sam -|samtools sort -u -@10 -T $location/tmp/example_prefix|samtools markdup -@10 - $location/results/$annot.bam
	
	#this is needed for following steps (chipqc)
	samtools index $location/results/$annot.bam
	 
	echo result bam file is now in results folder
	
	#we don't need the sam file anymore
	rm $location/results/$annot.sam
	
	i=$i+1
done
	
	echo loop is finished
