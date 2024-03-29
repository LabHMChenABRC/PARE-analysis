#!/usr/bin/env bash
# This script only applies to downloading single-end degradome dataset
# Only R1 file is kept if you use this script to download paired-end dataset

wd=$PWD
Dataset=$wd/PARE.Dataset.txt
# Dataset is a tab-delimited file providing this info:
# <Accession> <SampleName>
CPU_num=16
Output=$wd/PARE/raw

[ -d $Output ] || mkdir -p $Output
cd $Output
RunList=$(cat $Dataset | cut -f 1)

function SRALayout() {
    local NLine=$(fastq-dump -X 1 -Z --split-spot $1 2>/dev/null | wc -l)
    if [ $NLine -eq 4 ]
    then
	    echo "SINGLE" 
    else
	    echo "PAIRED"
    fi
}

echo Data will be saved in $Output
for Run in $RunList
do
	SampleName=$(grep ${Run} ${Dataset} | cut -f 2)
	date +"%b %d %T ..... Download ${Run} ${SampleName}"
	LibraryLayout=$( SRALayout ${Run} )
	if [ $LibraryLayout == "SINGLE" ]
	then
		date +"%b %d %T ..... Library is single-end data"
		seq_name_def='@$si'
	else
		date +"%b %d %T ..... Library is paired-end data and only R1 of mate pair is kept"
		seq_name_def='@$si/$ri'
	fi
	sleep 5 # Avoid failure caused by too fast query

	fasterq-dump --seq-defline $seq_name_def --qual-defline '+' --threads $CPU_num --temp /dev/shm --force --progress $Run
	# --temp          path of temp dir. This temporary directory will use approximately up to 10 times the size of the final output-file.
        #                 It is helpful for the speed-up, if the output-path and the scratch-path are on different file-systems.
        #                 For instance it is a good idea to point the temporary directory to a SSD if available or a RAM-disk like /dev/shm if enough RAM is available
        # --force         overwrite existing file(s)
        # --seq-defline   custom defline for sequence:
        #                 $ac=accession,  $sn=spot-name,
        #                 $sg=spot-group, $si=spot-id,
        #                 $ri=read-id,    $rl=read-length
        # --qual-defline  custom defline for qualities (see seq-defline)

	# short qname(sequence identifier) will be faster for mapping by bowtie

	if [ $LibraryLayout == "SINGLE" ]
	then
		bgzip -@ $CPU_num ${Run}.fastq
                mv ${Run}.fastq.gz ${SampleName}.fastq.gz

	else
		bgzip -@ $CPU_num ${Run}_1.fastq
                mv ${Run}_1.fastq.gz ${SampleName}.fastq.gz
                [ -f ${Run}_2.fastq ] && rm ${Run}_2.fastq
	fi
	date +"%b %d %T ..... Ouput: $Output/${SampleName}.fastq.gz"
done
date +"%b %d %T ..... Download Finish!"
