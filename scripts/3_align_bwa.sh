echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SAMPLE
#$ -o /master/kbailey/conta_BRE-LE/results/logs/align.SAMPLE.log
#$ -j y
#$ -q all.q
#$ -pe smp 10

""" >header.qsub.sh

#mkdir /master/kbailey/conta_BRE-LE/results/aligned
#mkdir /master/kbailey/conta_BRE-LE/results/nodups

REF_DIR="/master/kbailey/egg_RNA/data/reference"
RES_DIR="/master/kbailey/conta_BRE-LE/results"

for sample in $(cat /master/kbailey/conta_BRE-LE/data/seq_data/samples.list);
  do
  SAM=$(echo ${sample} | cut -c 1-11)

  CMD="bwa mem -M $REF_DIR/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa $RES_DIR/trimmed_fastq/${sample}_R1_001_val_1.fq.gz $RES_DIR/trimmed_fastq/${sample}_R2_001_val_2.fq.gz | samtools sort -o \"$RES_DIR/aligned/${sample}.bam\"
	conda activate home
	picard MarkDuplicates I=\"$RES_DIR/aligned/${sample}.bam\" O=\"$RES_DIR/nodups/${sample}.bam\" M=\"${sample}.txt\" REMOVE_DUPLICATES=true ASSUME_SORTED=true
	conda deactivate 
	samtools flagstat $RES_DIR/nodups/${sample}.bam > ${sample}.flagstat;
	  samtools index $RES_DIR/nodups/${sample}.bam"

    cat header.qsub.sh >${SAM}.sh
    sed -i "s/SAMPLE/${SAM}/g" ${SAM}.sh

    echo "$CMD" >>${SAM}.sh

    #qsub -V -cwd -S /bin/bash -N haplo_${sample}_${interval} -o /master/kbailey/ascaris_test/results/logs/haplo.o_${sample}_${interval}.log  ${sample}_${interval}.sh
    qsub ${SAM}.sh
  done
