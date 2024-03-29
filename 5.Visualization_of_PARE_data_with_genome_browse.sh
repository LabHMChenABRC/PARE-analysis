#!/usr/bin/env bash
wd=$PWD

# 6.Preparation of the genomic fasta index for genome browser
samtools faidx TAIR10_chr_all.fasta

# 7.Preparation of the Tabix gff for usage in genome browser using tabix program of HTSlib [Ref4]
(grep ^"#" TAIR10_GFF3_genes.gff; grep -v ^"#" TAIR10_GFF3_genes.gff | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n) | bgzip >TAIR10_GFF3_genes.sorted.gff.gz
tabix -p gff TAIR10_GFF3_genes.sorted.gff.gz

#8.Open Jbrowse2 desktop version, and press “Open sequence” with genomic fasta and fai file [Ref.9]. 

#9.Launch Linear genome view and open track selector. Add tracks with tabix GFF and genomic bigwig files on the right panel. The demonstration browse view display PARE signals in genome coordinate. (Fig 7a)

#10.Run Jbrowse2 with transcriptome fasta and transcriptome bigwig files, it can display PARE signals in transcript coordinate (Fig7b). The transcriptome fasta can be created using gffread.

gffread -w TAIR10.transcript.fasta -g TAIR10_chr_all.fasta TAIR10_GFF3_genes.gff

