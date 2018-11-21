#!/bin/sh
#PBS -l nodes=1:ppn=32
#PBS batch

module load bowtie2/2.3.1
module load samtools/0.1.18







### Variables to change
# 1.) Directory where bwt index are stored + prefix of bwt index
# 2.) Directory where the fastq are stored
# 3.) Array of chromosomes
BWT_DIR="/home/phamd/chromFa/mm10.genome"
MGS_DIR="/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/"
chromosome=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X")
cd $MGS_DIR









### Loop through each of the MGS_DO_* directory
for i in */;
do



  # Change to directory of the MGS_DO_* directory
  SAMPLE_DIR="$MGS_DIR$i"
  cd "$SAMPLE_DIR"
  
  echo $SAMPLE_DIR
  
  
  
  
  # Save the MGS_DO_* directory name without the ending '/' to a variable for bam outfile name
  sample=`echo -n $i | head -c -1`







  ### Loop through each of the week sub-directory in the MGS_DO_* directory
  for j in */;
  do

    # Save the week directory without the ending '/' as a variable for bam outfile name
    week=`echo -n $j | head -c -1`
    

    # Get the Pair End and Singleton directories
    SAMPLE_WEEK_DIR="${SAMPLE_DIR}${j}"
    SAMPLE_WEEK_PAIR_END_DIR="${SAMPLE_DIR}${j}Pair_End/"
    SAMPLE_WEEK_SINGLETON_DIR="${SAMPLE_DIR}${j}Singletons/"    
    





    # Get the singleton reads
    cd $SAMPLE_WEEK_SINGLETON_DIR
    
    #   Checking to see if 'b' version of the singleton fastq file exist
    temp=`find . -iname '*b_host.fastq.3' | wc -l`
    if [ $temp == 1 ]
    then
        singleton_1a=`ls | grep *a_host.fastq.3`
        singleton_1b=`ls | grep *b_host.fastq.3`

        singleton="${SAMPLE_WEEK_SINGLETON_DIR}${singleton_1a},${SAMPLE_WEEK_SINGLETON_DIR}${singleton_1b}"
    else

        singleton=`ls | grep *fastq.3`
        singleton="${SAMPLE_WEEK_SINGLETON_DIR}${singleton}"
    fi







    # Get the paired-end reads
    cd "$SAMPLE_WEEK_PAIR_END_DIR"
    
    #   Checking to see if 'b' version of the paired end fastq files exist
    temp=`find . -iname '*b_host.fastq*' | wc -l`
    if [ $temp == 2 ]
    then
        pair_end_1a=`ls | grep *a_host.fastq.1`
    	 pair_end_2a=`ls | grep *a_host.fastq.2`
        pair_end_1b=`ls | grep *b_host.fastq.1`
    	 pair_end_2b=`ls | grep *b_host.fastq.2`
        
	 pair_end_1="${pair_end_1a},${pair_end_1b}"
        pair_end_2="${pair_end_2a},${pair_end_2b}"
    else
     
         pair_end_1=`ls | grep *fastq.1`
        pair_end_2=`ls | grep *fastq.2` 
    fi




    # Printing the name of the read files
    echo ${pair_end_1}
    echo ${pair_end_2}
    echo ${singleton}









    ### Bowtie Alignment
    #     Save as BAM format, sort the BAM File, and finally create an BAM index file
    bowtie2 -q -p 32 -x ${BWT_DIR} -1 ${pair_end_1} -2 ${pair_end_2} -U ${singleton} | samtools view -bS - | samtools sort - "${SAMPLE_WEEK_DIR}${sample}_${week}_sorted"
    samtools index "${SAMPLE_WEEK_DIR}${sample}_${week}_sorted.bam"
    
    
    
    
    
    
    
    
  cd ${SAMPLE_DIR}; # Go back to the MGS_DO_* directory to loop through the next week
  done; # for j
  
  
cd $MGS_DIR;  # Go back to the directory where the MGS_DO directories are stored to loop through the next sample
done; # for i
