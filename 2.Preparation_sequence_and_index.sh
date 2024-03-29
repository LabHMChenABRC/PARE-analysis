#!/usr/bin/env bash
wd=$PWD
core=16
Url_GFF_file=https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
Url_Genome_file=https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
seq_dir=$PWD/Reference_sequence
GFF_file=$seq_dir/TAIR10_GFF3_genes.gff
GTF_file=${GFF_file%.gff}.gtf
Genome_file=$seq_dir/TAIR10_chr_all.fasta.gz
Genome_fasta=${Genome_file%.gz}

### Download and preparate for reference sequences and annotated file ###
mkdir -p $seq_dir
cd $seq_dir
# Download 
wget -O $GFF_file $Url_GFF_file
wget -O $Genome_file $Url_Genome_file
wget -O $miRNA_mature_file $Url_miRNA_mature_file
# Decompression
bgzip -@ $core -f -d $Genome_file >$Genome_fasta
# make fasta index (fasta.fai)
samtools faidx $Genome_fasta
# make Tabix GFF file
(grep ^"#" TAIR10_GFF3_genes.gff; grep -v ^"#" TAIR10_GFF3_genes.gff | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n) | bgzip >TAIR10_GFF3_genes.sorted.gff.gz
tabix -p gff TAIR10_GFF3_genes.sorted.gff.gz

# make non-coding RNA transcript sequences in fasta format (inculding rRNA, tRNA, snRNA and snoRNA)
# 1. make transcript fasta
gffread -w TAIR10.transcript.fasta -g $Genome_fasta $GFF_file
# 2. Fetch transcript ID of non-coding RNA from gff file
grep -P 'rRNA\t|tRNA\t|snRNA\t|snoRNA\t' $GFF_file | \
	gffread --stream --table @id >ncRNA.id.list
gffread --table @id TAIR10_ncRNAs.gff >ncRNA.id.list # Get ID list
samtools faidx --region-file ncRNA.id.list TAIR10.transcript.fasta >TAIR10.ncRNA.gff.fasta
# 3. Obtain a un-annotated rRNA fragments from NCBI
efetch -db nucleotide -id X52320.1 -format fasta >X52320.1.rRNA.fasta
# 4. Merge non-coding RNA sequences
cat TAIR10.ncRNA.gff.fasta X52320.1.rRNA.fasta >TAIR10.ncRNA.fasta
rm TAIR10_ncRNAs.gff ncRNA.id.list TAIR10.ncRNA.gff.fasta X52320.1.rRNA.fasta

# Subset mitochondria and chloroplast genomes from genome file
echo ChrC ChrM | xargs -n 1 samtools faidx $Genome_fasta >TAIR10.chrCM.fasta

# Preparation of index files for bowtie and STAR aligner
# 1. Create Bowtie Index for non-coding RNA transcripts, mitochondria and chloroplast genomes
mkdir -p $wd/BowtieIndex
bowtie-build $seq_dir/TAIR10.chrCM.fasta $wd/BowtieIndex/TAIR10.chrCM # Bowtie index of mitochondria and chloroplast genomes
bowtie-build $seq_dir/TAIR10.ncRNA.fasta $wd/BowtieIndex/TAIR10.ncRNA # Bowtie index of non-coding RNA transcripts

# 2. Create STAR Index with genome fasta and annotated file in gtf format
gffread -T -o $GTF_file $GFF_file #Covert gff to gtf format
mkdir -p $wd/STARIndex
cd $wd/STARIndex
STAR \
--runThreadN $core \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--genomeDir $wd/STARIndex/TAIR10 \
--genomeFastaFiles $Genome_fasta \
--sjdbGTFfile $GTF_file
