#!/bin/bash
#
# aligns XP12_lambda.fa reference with XP12_unMC_lambda_reads.fasta generated
# from poretools using marginAlign

MARGIN_ALIGN="../ont_dap/marginAlign"

# run marginAlign
bash $MARGIN_ALIGN/marginAlign Rh_20151111.fastq ../simulation/chr1.fasta Rh_20151111_mA.sam --jobTree ./jobTree 2> Rh_20151111_mA_err.txt

# convert sam to bam
samtools view -bt ../simulation/chr1.fa.fai Rh_20151111_mA.sam > Rh_20151111_mA.bam

# sort
samtools sort Rh_20151111_mA.bam Rh_20151111_mA.sorted

# index
samtools index Rh_20151111_mA.sorted.bam

# print samtools stats
samtools flagstat Rh_20151111_mA.sorted.bam

# create sam with proper header
samtools view -bt chr1.fa.fai Rh_20151111_mA.sam > Rh_20151111_mA_hdr.sam

# get alignment stats
bash $MARGIN_ALIGN/marginStats Rh_20151111_mA_hdr.sam Rh_20151111.fastq ../simulation/chr1.fasta --printAlignmentData --includeUnaligned > mA_stats_Rh_20151111.txt
