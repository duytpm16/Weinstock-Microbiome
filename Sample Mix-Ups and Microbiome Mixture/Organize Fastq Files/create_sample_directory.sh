### Creating a directory for each sample after unzipping the .gz files
#   There's a better way to do this...



### Go to directory where fastq files are
cd /projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/




### Create a directory for each sample

# Samples 401-417
for i in {401..417};
do
    mkdir "MGS_DO_$i";
done;



# Samples 419-426
for i in {419..426};
do
    mkdir "MGS_DO_$i";
done; 



# Samples 428-600
for i in {428..600};
do
    mkdir "MGS_DO_$i";
done;



# Sample 675
mkdir MGS_DO_675




# Samples 678-682
for i in {678..682};
do
    mkdir "MGS_DO_$i";
done;



# Samples 701-746
for i in {701..746};
do
    mkdir "MGS_DO_$i";
done;



# Samples 801-838
for i in {801..838};
do
    mkdir "MGS_DO_$i";
done;




# Samples 840-846
for in {840..846};
do
    mkdir "MGS_DO_$i";
done;



# Sample 850
mkdir MGS_DO_850
