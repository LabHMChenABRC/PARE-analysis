#!/usr/bin/env bash
wd=$PWD
Adapter=TGGAATTCTCGGGTGCCAAGG # What's 3' adapter sequence is used in libraries
core=16
two_color_platform=no # yes/no. Type yes if the dataset is sequencing from two-channel system

FqInDir=$wd/PARE/raw
OutDir_trim=$wd/PARE/trimmed
OutDir_filter=$wd/PARE/filter
OutDir_align=$wd/PARE/alignments
Bowtie_ncRNA_Index=$wd/BowtieIndex/TAIR10.ncRNA
Bowtie_chrCM_Index=$wd/BowtieIndex/TAIR10.chrCM
START_Index=$wd/STARIndex/TAIR10


mkdir -p $OutDir_trim $OutDir_filter $OutDir_align

for File in $FqInDir/*.fastq.gz
do
	ID=$(basename ${File%.fastq.gz})
	date +"%b %d %T ..... Processing $ID"
        date +"%b %d %T ..... Trim adapter, quality filtering for trimmed reads and keep 20-21 bp trimmed reads"
	if [[ $two_color_platform == "no" ]]
	then
		cutadapt \
			-o $OutDir_trim/${ID}.fastq.gz \
			--cores $core \
			-a $Adapter \
			--max-expected-errors 1 \
			--maximum-length=21 \
			--minimum-length=20 \
			$File >$OutDir_trim/$ID.trimmed.report.txt
		grep -A 99 "=== Summary" $OutDir_trim/$ID.trimmed.report.txt | grep -B 99 "Total written" --color=never
	else
		cutadapt \
                        -o $OutDir_trim/${ID}.trimmed.adapter.fastq.gz \
                        --cores $core \
                        -a $Adapter \
                        --max-expected-errors 1 \
                        --maximum-length=21 \
                        --minimum-length=20 \
			--too-long-output $OutDir_trim/$ID.overlength.fastq.gz \
                        $File >$OutDir_trim/$ID.trimmed.report.txt
		grep -A 99 "=== Summary" $OutDir_trim/$ID.trimmed.report.txt | grep -B 99 "Total written" --color=never

		date +"%b %d %T ..... Shorten over-length reads to 5' most 20 bp"
		cutadapt \
                        -o $OutDir_trim/${ID}.shorten.fastq.gz \
                        --cores $core \
			--length 20 \
                        --max-expected-errors 1 \
                        $OutDir_trim/$ID.overlength.fastq.gz >$OutDir_trim/$ID.shorten.report.txt
		grep -A 99 "=== Summary" $OutDir_trim/$ID.shorten.report.txt | grep -B 99 "Total written" --color=never
		cat $OutDir_trim/$ID.trimmed.fastq.gz $OutDir_trim/${ID}.shorten.fastq.gz >$OutDir_trim/${ID}.fastq.gz
		rm $OutDir_trim/$ID.overlength.fastq.gz $OutDir_trim/${ID}.shorten.fastq.gz
	fi

	# Removing contamination reads that are perfectly matched on the forward strand of non-coding RNAs
	date +"%b %d %T ..... Removing trimmed reads of non-coding RNAs"
	bowtie -p $core -q -k 1 -v 0 --norc \
		--un $OutDir_filter/$ID.remove.ncRNA.fastq.gz \
		-x $Bowtie_ncRNA_Index $OutDir_trim/${ID}.fastq.gz /dev/null 2>$OutDir_filter/$ID.bowtie.ncRNA.log
	cat $OutDir_filter/$ID.bowtie.ncRNA.log
	echo ""
	# Removing contamination reads that are perfectly matched on mitochondria and chloroplast genomes
	date +"%b %d %T ..... Removing trimmed reads of mitochondria and chloroplast genomes"
       	bowtie -p $core -q -k 1 -v 0 \
		--un $OutDir_filter/$ID.fastq.gz \
		-x $Bowtie_chrCM_Index $OutDir_filter/$ID.remove.ncRNA.fastq.gz /dev/null 2>$OutDir_filter/$ID.bowtie.chrCM.log
	cat $OutDir_filter/$ID.bowtie.chrCM.log
	echo ""
	rm $OutDir_filter/$ID.remove.ncRNA.fastq.gz
	
	# Map filtering reads to the genome and makes spliced alignments across known splice junctions by STAR. 
	date +"%b %d %T ..... Map remaining reads to genome and transcriptome"
	STAR \
		--runThreadN $core \
		--readFilesIn $OutDir_filter/$ID.fastq.gz \
		--outFileNamePrefix $OutDir_align/${ID}_ \
		--genomeDir $START_Index \
		--outMultimapperOrder Random \
		--alignEndsType EndToEnd \
		--outFilterMismatchNmax 0 \
		--outFilterMultimapNmax 10 \
		--alignIntronMax 10000 \
		--alignSJoverhangMin 500 \
		--alignSJDBoverhangMin 1 \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM
	# Sometime, _STARtmp folder is not removed after STAR mapping 
	if [ -d $OutDir_align/${ID}__STARtmp ]
	then
		[[ $(tail -n 1 $OutDir_align/${ID}_Log.out) == "ALL DONE!" ]] && rm -r $OutDir_align/${ID}__STARtmp
	fi

	echo ""

done
echo ######################################
echo Trimmed reads are save in $OutDir_trim
echo Filtered/reaming reads are save in $OutDir_filter
echo Genomic and transcript alignments file are save in $OutDir_align, 
echo and their filename tail with "_Aligned.sortedByCoord.out.bam" and "_Aligned.toTranscriptome.out.bam"
date +"%b %d %T ..... Finish Pre-processing and mapping"

# Maxium number of mismatches: --outFilterMismatchNmax 0
# prevent soft-cliping aligments: --alignEndsType EndToEnd
# minimum overhang for annotated junctions : --alignSJDBoverhangMin 1
# prevent novel splicing if minimum overhang for annotated junctions more than size of reads: --alignSJoverhangMin 500
# cutoff of maximum intron size and it also affect seed window size for mapping: --alignIntronMax 10000 
# Also output transcriptome bam: --quantMode TranscriptomeSAM
