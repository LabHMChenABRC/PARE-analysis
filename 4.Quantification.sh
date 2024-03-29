#!/usr/bin/env bash
wd=$PWD
core=16
STAR_alignments_dir=$wd/PARE/alignments
OutDir_BW=$STAR_alignments_dir/bigwig
mkdir -p $OutDir_BW
cd $OutDir_BW
for Genomic_Bam in $STAR_alignments_dir/*_Aligned.sortedByCoord.out.bam
do 
	ID=$(basename ${Genomic_Bam%_Aligned.sortedByCoord.out.bam})
	Transcript_Bam=${Genomic_Bam%_Aligned.sortedByCoord.out.bam}_Aligned.toTranscriptome.out.bam
	date +"%b %d %T ..... Process $ID"
	date +"%b %d %T ..... Count primary mapping reads for the scale factor"
	# 1. Obtain the total number of primary mapping reads in the genomic bam file using samtools. The count is signed to the variant “mapping_cnt”.
	mapping_cnt=$(samtools view -c -F 256 $Genomic_Bam)

	# 2. Calculate scale factor of TP40M for normalization. The value is signed to the variant “scale_factor”.
	scale_factor=$(perl -e "printf('%.2f',40*1000000/$mapping_cnt)")
	date +"%b %d %T ..... scale factor: $scale_factor primary mapping reads: $mapping_cnt"

	# 3. Create a bigwig file based genomic forward alignments thought deepTools modules bamCoverage with the scale factor [Ref.8]. 
	if [ ! -f ${Genomic_Bam}.bai ]
	then
		date +"%b %d %T ..... make index for genomic bam file"
                samtools index -@ $core $Genomic_Bam
	fi
	date +"%b %d %T ..... make genomic forward bigwig"
	bamCoverage \
	--numberOfProcessors $core \
	-b $Genomic_Bam \
	--outFileFormat bigwig \
	--Offset 1 \
	--binSize 1 \
	--samFlagExclude 256 \
	--scaleFactor $scale_factor \
	--filterRNAstrand reverse -o $ID.genomic.plus.bw
	
	# 4. Create a bigwig file based genomic reverse alignments thought deepTools modules bamCoverage with the scale factor.
	date +"%b %d %T ..... make genomic reverse bigwig"
	bamCoverage \
	--numberOfProcessors $core \
	-b $Genomic_Bam \
	--outFileFormat bigwig \
	--Offset 1 \
	--binSize 1 \
	--samFlagExclude 256 \
	--scaleFactor -$scale_factor \
	--filterRNAstrand forward -o $ID.genomic.minus.bw
	
	# 5. Transcriptome bam should be sorted and indexed. And then using bamCoverage with the scale factor to produce a bigwig file based on transcriptome forward alignment.
	if [ -f $Transcript_Bam ]
	then
		Transcript_Bam_sorted=${Transcript_Bam%.out.bam}.sorted.out.bam
		if [ ! -f $Transcript_Bam_sorted ]
		then
			date +"%b %d %T ..... sort and index transcriptome bam file"
			samtools sort -@ $core -o $Transcript_Bam_sorted $Transcript_Bam
			samtools index -@ $core $Transcript_Bam_sorted
		fi
		if [ ! -f ${Transcript_Bam_sorted}.bai ]
                then
			samtools index -@ $core $Transcript_Bam_sorted

		fi
		date +"%b %d %T ..... make transcriptome forward bigwig"
		bamCoverage \
		--numberOfProcessors $core \
		-b $Transcript_Bam_sorted \
		--outFileFormat bigwig \
		--Offset 1 \
		--binSize 1 \
		--samFlagExclude 256 \
		--scaleFactor $scale_factor \
		--filterRNAstrand reverse -o $ID.Transcriptome.bw
	else
		echo ${ID}_Aligned.toTranscriptome.out.bam is not existed!
		echo Skip make a Transcriptome bigwig file.
	fi
done
date +"%b %d %T ..... All files are processed"
echo "Forward and reverse genomic bigwig files and trancriptome bigwig file are save in $OutDir_BW"

# 6.	Preparation of the genomic fasta index for genome browser
#samtools faidx TAIR10_chr_all.fasta

# 7.	Preparation of the Tabix gff for usage in genome browser using tabix program of HTSlib [Ref4]
# (grep ^"#" TAIR10_GFF3_genes.gff; grep -v ^"#" TAIR10_GFF3_genes.gff | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n) | bgzip >TAIR10_GFF3_genes.sorted.gff.gz
# tabix -p gff TAIR10_GFF3_genes.sorted.gff.gz

#8.	Open Jbrowse2 desktop version, and press “Open sequence” with genomic fasta and fai file [Ref.9]. 

#9.	Launch Linear genome view and open track selector. Add tracks with tabix GFF and genomic bigwig files on the right panel. The demonstration browse view display PARE signals in genome coordinate. (Fig 11a)

#10.	Run Jbrowse2 with transcriptome fasta and transcriptome bigwig files, it can display PARE signals in transcript coordinate (Fig11b). The transcriptome fasta can be created using gffread.

# gffread -w TAIR10.transcript.fasta -g TAIR10_chr_all.fasta TAIR10_GFF3_genes.gff

