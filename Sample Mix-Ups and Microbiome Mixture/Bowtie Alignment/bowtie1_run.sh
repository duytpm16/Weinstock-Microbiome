#PBS -l nodes=1:ppn=2
#PBS batch

module load bowtie/1.1.2
module load samtools/1.8







### Variables to change
# 1.) Directory where bwt index are stored
# 2.) Directory where the fastq are stored
BWT_DIR="/home/phamd/chromFa/mm10.bwt1.autosome"
MGS_DIR="/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/"
cd $MGS_DIR









### Loop through each of the MGS_DO_* directory
for (( i=$start; i <= $end; i++ ))
do

  # Change to directory of the MGS_DO_* directory if it exists
  sample="MGS_DO_${i}"
  SAMPLE_DIR="${MGS_DIR}${sample}/"
  if [ -d "$SAMPLE_DIR" ]; then
     cd "$SAMPLE_DIR" 
     echo $SAMPLE_DIR




     ### Loop through each of the week sub-directory in the MGS_DO_* directory
     weeks=('Week_6' 'Week_17' 'Week_24')
     for j in ${weeks[@]};
     do
       # Getting directory names
       SAMPLE_WEEK_DIR="${SAMPLE_DIR}${j}/"
       SAMPLE_WEEK_PAIR_END_DIR="${SAMPLE_WEEK_DIR}Pair_End/"
       SAMPLE_WEEK_SINGLETON_DIR="${SAMPLE_WEEK_DIR}Singletons/"    
    
    

       # If the week exists for the sample... 
       if [ -d "$SAMPLE_WEEK_DIR" ]; 
       then

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
        
	     pair_end_1="$SAMPLE_WEEK_PAIR_END_DIR${pair_end_1a},$SAMPLE_WEEK_PAIR_END_DIR${pair_end_1b}"
             pair_end_2="$SAMPLE_WEEK_PAIR_END_DIR${pair_end_2a},$SAMPLE_WEEK_PAIR_END_DIR${pair_end_2b}"
          else

             pair_end_1=`ls | grep *fastq.1`
    	     pair_end_2=`ls | grep *fastq.2`

	     pair_end_1="$SAMPLE_WEEK_PAIR_END_DIR$pair_end_1"
	     pair_end_2="$SAMPLE_WEEK_PAIR_END_DIR$pair_end_2"
          fi

	  echo "Sample-week = $SAMPLE_WEEK_DIR"
          echo "Pair-end 1 = ${pair_end_1}"
          echo "Pair-end 2 = ${pair_end_2}"
          echo "Singleton  = ${singleton} "



	  # Creating dirctory to store bowtie1 results
          cd $SAMPLE_WEEK_DIR
          mkdir bowtie1_run_v2
          bam_file="${SAMPLE_WEEK_DIR}bowtie1_run_v2/${sample}_${week}_bowtie1"



	  # Bowtie1 begin
          bowtie -q --quiet --sam -m 1 -p 2 ${BWT_DIR} -1 ${pair_end_1} -2 ${pair_end_2} -s ${singleton} | samtools view -bS - > "${bam_file}.bam"

	  # Sort bam file
          samtools sort "${bam_file}.bam" > "${bam_file}_sorted.bam"
          samtools index "${bam_file}_sorted.bam"



          cd ${SAMPLE_DIR};
       fi
     done; # for j
   fi
  cd $MGS_DIR;
done; # for i
