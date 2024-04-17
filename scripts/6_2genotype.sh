# This part will set up the qsub properties and save them in a new script

echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N NAME
#$ -o /master/kbailey/conta_BRE-LE/results/logs/gvcf.NAME.log
#$ -j y
#$ -q all.q
#$ -pe smp 10

""" >header.qsub.sh

# This is where we set up new folders and variables
mkdir /master/kbailey/conta_BRE-LE/results/genotype

REF=/master/kbailey/egg_RNA/data/reference

# This is where we write the code that will be run. Note that the main command is stored as a variable.
for interval in $(cat $REF/intervals/intervals.list); do
  NAME=$(echo ${interval} | cut -f1 -d"-")

  CMD="gatk --java-options \"-Xmx12g\" GenotypeGVCFs \
        -R $REF/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa \
        -V gendb://$master/master/kbailey/conta_BRE-LE/results/gdbconta/${NAME} \
        -O /master/kbailey/conta_BRE-LE/results/genotype/${NAME}.vcf.gz \
        --tmp-dir /master/kbailey/sm_single_gt/temp" #-L $ref/intervals/${interval}"

# Save the header as a script for each individual loop
      cat header.qsub.sh >Sm_${NAME}gvcf.sh
# Modify the header to change the name of the job and log file
      sed -i "s/NAME/Sm_${NAME}/g" Sm_${NAME}gvcf.sh
# Append the command to the new script
      echo $CMD >>Sm_${NAME}gvcf.sh
# Execute the script on the server
      qsub Sm_${NAME}gvcf.sh
done
