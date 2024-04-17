echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SAMPLE
#$ -o /master/kbailey/conta_BRE-LE/results/logs/haplo.SAMPLE_INTERVAL.log
#$ -j y
#$ -q all.q
#$ -pe smp 1

""" >header.qsub.sh

REF_DIR="/master/kbailey/egg_RNA/data/reference"
RESULTS_DIR="/master/kbailey/conta_BRE-LE/results"

mkdir /master/kbailey/conta_BRE-LE/results/haplo
mkdir /master/kbailey/conta_BRE-LE/temp

for interval in $(cat $REF_DIR/intervals/intervals.list); do
  for sample in $(cat /master/kbailey/conta_BRE-LE/data/seq_data/samples.list); do
    NUM=$(echo ${interval} | cut -f1 -d"-")
    SAM=$(echo ${sample} | cut -c 1-11)
	
    CMD="gatk --java-options \"-Xmx12g\" HaplotypeCaller \
    -R $REF_DIR/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa \
    -I \"$RESULTS_DIR/finalbam/${sample}.bam\" \
    -O \"$RESULTS_DIR/haplo/${sample}.${NUM}.g.vcf.gz\" \
    -ERC GVCF \
    -L $REF_DIR/intervals/${interval} \
    --tmp-dir /master/kbailey/conta_BRE-LE/temp"

    cat header.qsub.sh >${SAM}_$NUM.sh
    sed -i "s/SAMPLE/${SAM}/g" ${SAM}_${NUM}.sh
    sed -i "s/INTERVAL/${NUM}/g" ${SAM}_${NUM}.sh

    echo $CMD >>${SAM}_${NUM}.sh

    qsub ${SAM}_${NUM}.sh
  done
done
