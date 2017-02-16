#!/usr/bin/env bash

#----------------------------------------------------------------------------
# The MIT License
#
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:

#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.

#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#-----------------------------------------------------------------------------

# version 021617EH

export LANG=en_EN
StlJucLen=8
STLfile=STL96.txt
MAQ=30
READQ=0
MINREAD=1

while getopts "b:f:r:t:s:q:m:h*" opt; do
        case $opt in
                b) R1bam=$OPTARG;;
		f) R1fq=$OPTARG;;
		r) R2fq=$OPTARG;;
		t) STLfile=$OPTARG;;
		s) StlJucLen=$OPTARG;;
		q) MAQ=$OPTARG;;
		m) MINREAD=$OPTARG;;
                h) printf "Usage: dqRNASeq [options]\n-b bamfile \n-f R1_fastq \n-r R2_fastq \n-t STL_sequence_file [$STLfile] \n-s STL_size [$StlJucLen] \n-q MAQ_Cutoff [$MAQ]\n-m NREADS_Cutoff [$MINREAD]\n-h help\n"
		   exit;; 
		*) printf "Usage: dqRNASeq [options]\n-b bamfile \n-f R1_fastq \n-r R2_fastq \n-t STL_sequence_file [$STLfile] \n-s STL_size [$StlJucLen] \n-q MAQ_Cutoff [$MAQ]\n-m NREADS_Cutoff [$MINREAD]\n-h help\n"
		   exit;;      
	esac
done

OIFS=$IFS
IFS=','
genes=($genelist)
IFS=$OIFS

prefix=${R1bam%.bam}
tmpfile=$RANDOM$RANDOM$RANDOM".tmp"
Outfile1=$prefix".prefixPE"$StlJucLen".MQ"$MAQ"_MINREAD"$MINREAD".txt"
Outfile2=$prefix".STLprefixPE"$StlJucLen".MQ"$MAQ"_MINREAD"$MINREAD".txt"
Outfile3=$prefix".uniqSTLprefixPE"$StlJucLen".MQ"$MAQ"_MINREAD"$MINREAD".txt"
Outfile4=$prefix".uniqprefixPE"$StlJucLen".MQ"$MAQ"_MINREAD"$MINREAD".txt"

#for R1 in pair, make R1 bam
samtools view "$R1bam" | awk -v MAQ=$MAQ '$5>=MAQ { if ( $2==99 ){print $1"\t"$3"\t"$4"\t"$4+$9-1}}' > ${tmpfile}.R1.bam

#for R2 in pair, make R2 bam
samtools view "$R1bam" | awk -v MAQ=$MAQ '$5>=MAQ { if ( $2==163) {print $1"\t"$3"\t"$4"\t"$4+$9-1}}' > ${tmpfile}.R2.bam

#make list of headers and UMIs
zcat "$R1fq" | awk -v LEN=$StlJucLen 'NR%4==1{name=substr($1,2,length($1)-1)} NR%4==2{print name"\t"substr($1,1,LEN)}' > ${tmpfile}.R1.fq &
zcat "$R2fq" | awk -v LEN=$StlJucLen 'NR%4==1{name=substr($1,2,length($1)-1)} NR%4==2{print name"\t"substr($1,1,LEN)}' > ${tmpfile}.R2.fq &
wait;

#sort above list
sort -S16G -T ./ -b -k1,1 ${tmpfile}.R1.fq > ${tmpfile}.R1.fq.sorted &
sort -S16G -T ./ -b -k1,1 ${tmpfile}.R2.fq > ${tmpfile}.R2.fq.sorted &
wait;

if [ -s "${tmpfile}.R1.bam" ]; then
    echo -e "R1 bam exists; sorting"
    sort -S16G -T ./ -b -k1,1 ${tmpfile}.R1.bam > ${tmpfile}.R1.bam.sorted #sort R1 bam
	join -j 1 <( join -j 1 ${tmpfile}.R1.bam.sorted ${tmpfile}.R1.fq.sorted ) <( cat ${tmpfile}.R2.fq.sorted ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${tmpfile}.joined1 #join R1 bam with R1 fastq, join resulting file with R2 fastq
	sort -S16G -T ./ -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 ${tmpfile}.joined1 > ${tmpfile}.joined1.sorted
	uniq -c ${tmpfile}.joined1.sorted | awk -v MINREAD=$MINREAD '$1>=MINREAD{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1}' > $Outfile1
fi

if [ -s "${tmpfile}.R2.bam" ]; then
    echo -e "R2 bam exists; sorting"
    sort -S16G -T ./ -b -k1,1 ${tmpfile}.R2.bam > ${tmpfile}.R2.bam.sorted #sort R2 bam
	join -j 1 <( join -j 1 ${tmpfile}.R2.bam.sorted ${tmpfile}.R2.fq.sorted ) <( cat ${tmpfile}.R1.fq.sorted ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${tmpfile}.joined2 #join R2 bam with R1 fastq, join resulting file with R1 fastq
	sort -S16G -T ./ -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 ${tmpfile}.joined2 > ${tmpfile}.joined2.sorted
	uniq -c ${tmpfile}.joined2.sorted | awk -v MINREAD=$MINREAD '$1>=MINREAD{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1}' >> $Outfile1
fi


STLsize=$(head -n 1 $STLfile | wc -c)
let STLsize=STLsize-1
echo $STLsize

join -j 1 <( join -j 1 <( awk -v LEN=$STLsize '{print substr($4,1,LEN)"\t"$0}' $Outfile1 | sort -bk 1 ) <( sort -bk 1 $STLfile) | awk -v LEN=$STLsize '{print substr($6,1,LEN)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sort -bk 1 ) <( sort -bk 1 $STLfile ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'  > $Outfile2


join -j 1 <( cut -f 1,6 $Outfile2 | sort -bk 1 | awk 'BEGIN{Name="---"; Sum=0} {if (Name==$1){Sum+=$2}else{if (Name!="---") print Name"\t"Sum; Name=$1; Sum=$2}} END{print Name"\t"Sum}' ) <( join -j 1 <( cut -f 1 $Outfile2 | sort -bk 1 | uniq -c | awk '{print $2"\t"$1}' ) <( join -j 1 <( cut -f 1,2,3 $Outfile2 | sort -k1 -k2 -k3 | uniq | cut -f 1 | uniq -c | awk '{print $2"\t"$1}' | sort -bk 1 ) <( cut -f 1,4,5 $Outfile2 | sort -k1 -k2 -k3 | uniq | cut -f 1 | uniq -c | awk '{print $2"\t"$1}' | sort -bk 1 ))) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'> $Outfile3

join -j 1 <( cut -f 1,6 $Outfile1 | sort -bk 1 | awk 'BEGIN{Name="---"; Sum=0} {if (Name==$1){Sum+=$2}else{if (Name!="---")print Name"\t"Sum; Name=$1; Sum=$2}} END{print Name"\t"Sum}' ) <( join -j 1 <( cut -f 1 $Outfile1 | sort -bk 1 | uniq -c | awk '{print $2"\t"$1}' ) <( join -j 1 <( cut -f 1,2,3 $Outfile1 | sort -k1 -k2 -k3 | uniq | cut -f 1 | uniq -c | awk '{print $2"\t"$1}' | sort -bk 1 ) <( cut -f 1,4,5 $Outfile1 | sort -k1 -k2 -k3 | uniq | cut -f 1 | uniq -c | awk '{print $2"\t"$1}' | sort -bk 1 ))) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'> $Outfile4


#add headers to output files
sed -i '1itranscript_id\tstart\tstop\tR1_STL\tR2_STL\nReadPairs' ${Outfile1}
sed -i '1itranscript_id\tstart\tstop\tR1_STL\tR2_STL\nReadPairs' ${Outfile2}
sed -i '1itranscript_id\tnReadPairs\tnSTL+nUSS\tnUSS\tnSTL' ${Outfile3}
sed -i '1itranscript_id\tnReadPairs\tnSTL+nUSS\tnUSS\tnSTL' ${Outfile4}

rm *.tmp*