#RNA-SEQ Workflow for homo sapiens (single-end reads)

#the absolute location of the files (this will be the reference directory)
location=/home/belab/Desktop/rnaseq

#samples are stored in a variable (to be iterated)
samples=$location/SRR/*.fastq

#this part should be changed if your data (.fastq files) do not follow order of sample_annot array
declare -i i=0
sample_annot=(wt1 wt2 wt3 ko1 ko2 ko3)

#This 'for' loop iterates over the fastq files and processes the data through hisat2, samtools annd stringtie pipe
for sample in $samples
do
	annot=${sample_annot[$i]} #now the annotation is in annot variable
	
        echo $sample #prints the fastq file in progress
        
        hisat2 -p 10 --dta -x $location/grch38_tran/genome_tran -U $sample --summary-file $location/SRR/$annot.txt -S $location/sam/$annot.sam
        
        #I wanted to pipe hisat2 to samtools but it gives an error saying there is not enough memory, that's why hisat2 is separate
        
        echo quality of the alignments are in the SRR file
        echo continuing to samtools and stringtie
        
        samtools sort -@ 9 -T $location/bam $location/sam/$annot.sam|stringtie -p 10 -G $location/grch38_tran/genome.gtf -o $location/gtf/$annot/$annot.gtf - -e #samtools and stringtie pipe
       
	echo samtools and stringtie part is finished
	
	#We don't need the sam file (output of hisat2) anymore
	rm $location/sam/$annot.sam
        
        i=$i+1 
done
	echo loop is finished, continuing with prepDE python script
	
	
	
#after everything is done, gtf files will be combined with a stringtie python script to give a dataframe ready to analyze in deseq2
python3 prepDE.py3 -i $location/gtf

#we don't need the gtf files anymore
rm -r $location/gtf


