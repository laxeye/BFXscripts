#!/bin/bash
# Double index demultiplexing pipeline for Illumina reads
# 
# Based on Fadrosh et. al. "An improved dual-indexing approach for multiplexed 16S rRNA gene sequencing on the Illumina MiSeq platform. 2014"
# https://github.com/igsbma/MiSeq16S

# Requirements:
# Perl, Python
# *fq_getPairAndOrphan1.8.py (step 4) from https://github.com/igsbma/MiSeq16S 
# **Biopython >= 1.51
# *fastx-toolkit (step 1) http://hannonlab.cshl.edu/fastx_toolkit/
# sudo apt-get install fastx-toolkit
# * seqtk (step 2) from https://github.com/lh3/seqtk
# sudo apt-get install seqtk
# * seqprep from https://github.com/jstjohn/SeqPrep
# sudo apt-get install seqprep
# * fastq-mcf from ea-utils (step 0)
# sudo apt-get install ea-utils
# TagCleaner (http://tagcleaner.sourceforge.net)
# tagcleaner.pl should be in your $PATH

if ( [[ -z "$1" ]] || [[ $1 == -h* ]] || [[ $1 == -v* ]] )
then
	echo "Double index demultiplexing pipeline for Illumina based 16S amplicon sequencing"
	echo "A. Korzhenkov, N. Samarov and S. Toshchakov. 2016"
	echo "Based on"
	echo "Fadrosh et. al. An improved dual-indexing approach for multiplexed 16S rRNA gene sequencing on the Illumina MiSeq platform. 2014"
	echo "https://github.com/igsbma/MiSeq16S"
	echo "Please run as: $0 <Forward reads> <Reverse reads> <Quality score for filtering, 25 by default>"
	exit
fi

if ! ( [[ -f $1 ]] &&  [[ -f $2 ]] )
then
	echo "Input files not found!"
	exit 0
elif [[ "$1" == "$2" ]]
	then
	echo "Input files are the same!"
	exit 0
else
	R1="$1"
	R2="$2"
fi

Q=25 #Default quality score for trimming
if [[ -n "$3" ]]; then Q="$3"; fi #User-defined quality score for trimming

if ! ( which seqprep > /dev/null ); then FAIL="\nYou should install seqprep.";fi
if ! ( which fastq-mcf > /dev/null ); then	FAIL="$FAIL\nYou should install fastq-mcf from ea-utils.";fi
if ! ( which fastx_trimmer > /dev/null ); then
	FAIL="$FAIL\nYou should install fastx_trimmer from fastx-toolkit."
fi
if ! ( which seqtk > /dev/null ); then	FAIL="$FAIL\nYou should install seqtk.";fi
if ! ( which tagcleaner.pl > /dev/null ); then
	FAIL="$FAIL\nYou should have tagcleaner.pl in your PATH and make it executable."
fi
if ! ( which fq_getPairAndOrphan1.8.py > /dev/null ); then
	FAIL="$FAIL\nYou should have fq_getPairAndOrphan1.8.py in your PATH and make it executable."
fi

if [[ -n "$FAIL" ]]
	then
	echo "Requirements are not satisfied:"
	echo -e "$FAIL"
	exit 0
fi


if [[ -d tmp/ ]]; then rm -r tmp/; fi

mkdir tmp/
# 0. Quality trimming
# a. Using fastq-mcf from ea-utils
echo "Quality trimming and filtering"
if [[ $R1 == *gz ]]
then
	fastq-mcf -o tmp/R1.fq -o tmp/R2.fq n/a <(gunzip -ck $R1;) <(gunzip -ck $R2;) -q $Q -w 5 -k 0
else
	fastq-mcf -o tmp/R1.fq -o tmp/R2.fq n/a $R1 $R2 -q $Q -w 5 -k 0
fi

# 1. Generating barcodes-only fastq files
echo "Generating barcode files"
fastx_trimmer -i tmp/R1.fq -f 1 -l 12 -Q 33 | perl -pe 's/\n/\t/ if ($c+1)%4 != 0;$c++' >  tmp/R1_barcode_temp
fastx_trimmer -i tmp/R2.fq -f 1 -l 12 -Q 33 | perl -pe 's/\n/\t/ if ($c+1)%4 != 0;$c++' >  tmp/R2_barcode_temp
paste tmp/R1_barcode_temp tmp/R2_barcode_temp | awk -F"\t" '{print $5"\t"$2$6"\t"$3"\t"$4$8}' | sed 's/\t/\n/g' > tmp/Barcodes.fastq && rm tmp/R*_barcode*

# 2. Trimming barcode off the sequences
echo "Trimming barcodes off the sequences"
seqtk trimfq -b 12 tmp/R1.fq > tmp/R1_trimmed.fastq
seqtk trimfq -b 12 tmp/R2.fq > tmp/R2_trimmed.fastq

# 3. Merging paired sequences, leaving only merged reads
echo "Merging paired sequences"
seqprep -f tmp/R1_trimmed.fastq -r tmp/R2_trimmed.fastq -1 /dev/null -2 /dev/null -s tmp/Merged.fastq.gz && rm tmp/R[12]*f*q && gunzip tmp/Merged.fastq.gz


# 4. Trimming off spacers and primers
# 4.1 Mismatch evaluating (only statistics)
echo "Trimming off spacers and primers"
tagcleaner.pl -fastq tmp/Merged.fastq -verbose -stats -tag5 GGACTACHVGGGTWTCTAAT -tag3 TTACCGCGGCKGCTGVCAC
# 4.2 Trimming.
echo "Choose mismatch number for 3end:"
read MN3
echo "Choose mismatch number for 5end:"
read MN5
echo "Cutting primers with mismatch N: $MN3 and $MN5:"
tagcleaner.pl -fastq tmp/Merged.fastq -verbose -tag5 GGACTACHVGGGTWTCTAAT -tag3 TTACCGCGGCNGCTGNCAC -mm3 $MN3 -mm5 $MN5 -out tmp/Merged.clean -nomatch 3

#5. Macthing up barcodes and merged reads
echo "Macthing up barcodes and merged reads"
fq_getPairAndOrphan1.8.py tmp/Merged.clean.fastq tmp/Barcodes.fastq Reads.ready.fastq Barcodes.ready.fastq /dev/null && rm -r tmp/

echo "Workflow finished, congratulations!"
exit 1