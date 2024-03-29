#!/usr/bin/env bash
wd=$PWD
mkdir -p $wd/Reference_sequence
Rep_info_file=$wd/Reference_sequence/TAIR10_representative_gene_models.txt
gff_file=$wd/Reference_sequence/TAIR10_GFF3_genes.gff
gtf_rep_file_URL=https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gene_lists/TAIR10_representative_gene_models
gtf_rep_file=$wd/Reference_sequence/TAIR10_representative_gene_models.gtf

# Get representative gene models information
[ ! -f $Rep_info_file ] && wget -O $Rep_info_file $gtf_rep_file_URL

# Keep representative_gene_models and discard transposable element genes
cut -f 1 $Rep_info_file | grep -P '^AT' >rep.ids.list
grep "mRNA_TE_gene\t" $gff_file | sed  's/mRNA_TE_gene/mRNA/' | gffread --table @id - >TE_ids.list
gffread --ids rep.ids.list $gff_file | gffread --nids TE_ids.list -T -o $gtf_rep_file -
rm rep.ids.list TE_ids.list
