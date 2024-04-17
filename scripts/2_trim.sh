echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SAMPLE
#$ -o /master/kbailey/conta_BRE-LE/results/logs/trim.SAMPLE.log
#$ -j y
#$ -q all.q
#$ -pe smp 10

""" >header.qsub.sh

mkdir /master/kbailey/conta_BRE-LE/results/trimmed_reports

REF_DIR="/master/kbailey/conta_BRE-LE/data/seq_data"
RESULTS_DIR="/master/kbailey/conta_BRE-LE/results"

for sample in $(cat /master/kbailey/conta_BRE-LE/data/seq_data/samples.list);
  do
  SAM=$(echo ${sample} | cut -f1 -d"_")
  
  CMD="trim_galore -q 28 --fastqc_args \"--outdir $RESULTS_DIR/trimmed_reports\" --illumina --max_n 1 --trim-n -o $RESULTS_DIR/trimmed_fastq --clip_R1 7 --clip_R2 7 --paired $REF_DIR/${sample}_R1_001.fastq.gz $REF_DIR/${sample}_R2_001.fastq.gz"

    cat header.qsub.sh >${SAM}.sh
    sed -i "s/SAMPLE/${SAM}/g" ${SAM}.sh

    echo $CMD >>${SAM}.sh

    qsub ${SAM}.sh
  done
