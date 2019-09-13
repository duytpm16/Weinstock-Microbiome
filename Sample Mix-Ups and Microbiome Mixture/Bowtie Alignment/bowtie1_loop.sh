#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS batch

module load bowtie/1.1.2
module load samtools/1.8







### Variables to change
# 1.) Directory where bwt index are stored
# 2.) Directory where the fastq are stored
BWT_DIR="/home/phamd/chromFa/mm10.bwt1.autosome"
MGS_DIR="/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/"
cd $MGS_DIR









### Loop through each of the MGS_DO_* directory
for (( i=$begin; i <= $end; i++ ))
do

  # Change to directory of the MGS_DO_* directory
  sample="MGS_DO_${i}"
  SAMPLE_DIR="${MGS_DIR}${sample}/"
  if [ -d "$SAMPLE_DIR" ]; then
     cd "$SAMPLE_DIR"
  
     echo $SAMPLE_DIR






     ### Loop through each of the week sub-directory in the MGS_DO_* directory
     for j in */;
     do

       # Save the week as a variable for bam outfile name
       week=`echo -n $j | head -c -1`
    

       # Change to the Pair_End directory within the  week directory 
       SAMPLE_WEEK_DIR="${SAMPLE_DIR}${j}"
       SAMPLE_WEEK_PAIR_END_DIR="${SAMPLE_DIR}${j}Pair_End/"
       SAMPLE_WEEK_SINGLETON_DIR="${SAMPLE_DIR}${j}Singletons/"    
    


       # Get the singleton reads
       cd $SAMPLE_WEEK_SINGLETON_DIR
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


       echo ${pair_end_1}
       echo ${pair_end_2}
       echo ${singleton}

    
       bam_file="${SAMPLE_WEEK_DIR}${sample}_${week}_bowtie1"
       ### Loop through each of the bwt indexes directory to perform bowtie alignment for every chromosome
       bowtie -q --quiet --sam -m 1 -p 4 ${BWT_DIR} -1 ${pair_end_1} -2 ${pair_end_2} -s ${singleton} | samtools view -bS - > "${bam_file}.bam"

       samtools sort "${bam_file}.bam" > "${bam_file}_sorted.bam"
       samtools index "${bam_file}_sorted.bam"

       cd ${SAMPLE_DIR};
       done; # for j
     fi
     cd $MGS_DIR;
done; # for i
