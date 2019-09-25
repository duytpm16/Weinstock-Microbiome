#PBS -l nodes=1:ppn=1



for(( i=400; i<=900; i=i+5))
do

  qsub -v start=$i,end=$(($i+4)) bowtie1_run.sh
  
done
