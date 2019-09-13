#PBS -l nodes=1:ppn=1
#PBS batch

module load bowtie/1.1.2
module load samtools/1.8







### Variables to change
# 1.) Directory where the fastq are stored
MGS_DIR="/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/"
cd $MGS_DIR









### Loop through each of the MGS_DO_* directory
for i in */;
do

  # Change to directory of the MGS_DO_* directory
  SAMPLE_DIR="${MGS_DIR}${i}"
  cd "$SAMPLE_DIR"


  # Loop through each of the week sub-directory in the MGS_DO_* directory
  week=('Week_6' 'Week_17' 'Week_24')
  for j in ${week[@]};
  do

    WEEK_DIR="$SAMPLE_DIR${j}/"

    if [ -d "$WEEK_DIR" ];
    then

       # Change to the Pair_End directory within the  week directory 
       SAMPLE_WEEK_PAIR_END_DIR="${WEEK_DIR}Pair_End/"
       cd $SAMPLE_WEEK_PAIR_END_DIR
       gunzip *.gz

       SAMPLE_WEEK_SINGLETON_DIR="${WEEK_DIR}Singletons/"
       cd $SAMPLE_WEEK_SINGLETON_DIR
       gunzip *.gz
    fi

  done;
done; # for i
