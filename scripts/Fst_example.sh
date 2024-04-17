# Example code for Fst calculation


#!/bin/bash

# Generate a file containing two sequences from pools of interest 
samtools mpileup -B ~/conta_BRE-LE/results/finalbam/LE102621_m.bam ~/conta_BRE-LE/results/finalbam/BRE102621_m.bam > ~/conta_BRE-LE/results/fst/LE102621_m_BRE102621_m.mpileup

# Convert fo popoolation input file
java -ea -Xmx7g -jar ~/bin/popoolation/mpileup2sync.jar --input ~/conta_BRE-LE/results/fst/LE102621_m_BRE102621_m.mpileup \
--output ~/conta_BRE-LE/results/fst/LE102621_m_BRE102621_m.mpileup_java.sync --fastq-type sanger --min-qual 20 --threads 8

# Calculate Fst
perl ~/bin/popoolation/fst-sliding.pl --input ~/conta_BRE-LE/results/fst/LE102621_m_BRE102621_m.mpileup_java.sync \
--output ~/conta_BRE-LE/results/fst/LE102621_m_BRE102621_m.fst --suppress-noninformative --min-count 6 --min-coverage 50 \
--max-coverage 200 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 60:100 
