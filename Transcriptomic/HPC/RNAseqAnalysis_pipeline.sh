#!/bin/bash

#  Transcriptomic (RNA-seq) Analysis Pipeline version 0.1.0 -- 20180509
#  Copyright (C) 2018 	Dr Tiphaine Martin (developer, contributor) 
#        
#  This script is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this script.  If not, see <http://www.gnu.org/licenses/>.
#
#  For any bugs or problems found, please contact us at
#  tiphaine.martin@mssm.edu
# Prints a message ($1) to a file ($2)
log() {
	echo $1 >> $2
}

prefix=$1
JOBDIR=$2
WORKINGDIR_orign=$3
FASTQDIR=$4
WORKINGDIR=$5
SAMPLEFILE1=$6
SAMPLEFILE2=$7

cd $JOBDIR
#Loads the parameters
source RNAseqAnalysis_parameters.sh

#Creates a log file
LOGFILE=$WORKINGDIR/$prefix.log


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STEP 0. Prints a log message
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

log "" $LOGFILE	
log "---------------------------------------------" $LOGFILE
log "TRANCRIPTOMIC ANALYSIS PIPELINE (version 0.1.0)" $LOGFILE
log "---------------------------------------------" $LOGFILE
log "" $LOGFILE
log "Copyright (C) 2018" $LOGFILE
log "Dr. Tiphaine C. Martin <tiphaine.martin@mssm.edu>"  $LOGFILE
log "" $LOGFILE
log "This script is distributed in the hope that it will be useful"  $LOGFILE
log "but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details." $LOGFILE
log "" $LOGFILE
log "Please report comments and bugs to:" $LOGFILE
log "   - tiphaine.martin@mssm.edu"     $LOGFILE
log "" $LOGFILE
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE

sysdate=$(date)
log "Analysis starting at $sysdate." $LOGFILE
log "" $LOGFILE	
log "Analysed samples are: $SAMPLEFILE1 and $SAMPLEFILE2" $LOGFILE
log "Working directory set to $WORKINGDIR" $LOGFILE
log "Logs saved at $LOGFILE" $LOGFILE
log "New files will be saved using the '$prefix' prefix" $LOGFILE
log "" $LOGFILE	
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE	

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STEP 1. Quality assessment of the raw data file.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sysdate=$(date)
starttime=$(date +%s.%N)
log "Performing STEP 1 [Assessment of raw data quality] at $sysdate" $LOGFILE
log "" $LOGFILE	

base1=$(basename $SAMPLEFILE1 | cut -d. -f 1)
base2=$(basename $SAMPLEFILE2 | cut -d. -f 1)

### Forward reads
fastqc --quiet --noextract --format fastq --threads $threads --outdir=$WORKINGDIR $SAMPLEFILE1
mv $WORKINGDIR/${base1}_fastqc.html $WORKINGDIR/${prefix}_R1_fastqc.html
mv $WORKINGDIR/${base1}_fastqc.zip $WORKINGDIR/${prefix}_R1_fastqc.zip

### Reverse reads
fastqc --quiet --noextract --format fastq --threads $threads --outdir=$WORKINGDIR $SAMPLEFILE2
mv $WORKINGDIR/${base2}_fastqc.html $WORKINGDIR/${prefix}_R2_fastqc.html
mv $WORKINGDIR/${base2}_fastqc.zip $WORKINGDIR/${prefix}_R2_fastqc.zip

endtime=$(date +%s.%N)
exectime=$(echo "$endtime - $starttime" | bc)
sysdate=$(date)
log "           STEP 1 terminated at $sysdate ($exectime seconds)" $LOGFILE
log "" $LOGFILE	
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE	

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STEP 1a. Quantification of duplication rate. duplication reads will not be
# removed in RNA-seq cfhttps://www.nature.com/articles/srep25533
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sysdate=$(date)
starttime=$(date +%s.%N)
log "Performing STEP 1a [De-duplication] at $sysdate" $LOGFILE
log "" $LOGFILE	

clumpify.sh -Xmx$maxmemory in1=$SAMPLEFILE1 in2=$SAMPLEFILE2 out1=$WORKINGDIR/${prefix}_dedupe_R1.fq.gz out2=$WORKINGDIR/${prefix}_dedupe_R2.fq.gz qin=$qin dedupe subs=0 threads=$threads &> $WORKINGDIR/${prefix}_step1a_tmp.log

#Logs some figures about sequences passing de-duplication
log "BBduk's de-duplication stats: " $LOGFILE
log "" $LOGFILE
sed -n '/Reads In:/,/Duplicates Found:/p' $WORKINGDIR/${prefix}_step1a_tmp.log >> $LOGFILE
log "" $LOGFILE
totR=$(grep "Reads In:" $WORKINGDIR/${prefix}_step1a_tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
remR=$(grep "Duplicates Found:" $WORKINGDIR/${prefix}_step1a_tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
survivedR=$(($totR-$remR))
percentage=$(echo $survivedR $totR | awk '{print $1/$2*100}' )
log "$survivedR out of $totR paired reads survived de-duplication ($percentage%, $remR reads removed)" $LOGFILE

#rm $WORKINGDIR/${prefix}_step1a_tmp.log

endtime=$(date +%s.%N)
exectime=$(echo "$endtime - $starttime" | bc)
sysdate=$(date)
log "" $LOGFILE
log "           STEP 2 terminated at $sysdate ($exectime seconds)" $LOGFILE
log "" $LOGFILE	
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE	

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STEP 2. Trimming of low quality bases and of adapter sequences.
# Short reads are discarded. When either forward or reverse of a paired-end read 
# are discarded, the surviving end is saved on a file of unpaired reads.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sysdate=$(date)
starttime=$(date +%s.%N)
log "Performing STEP 2 [Trimming] at $sysdate" $LOGFILE
log "" $LOGFILE	

rm $WORKINGDIR/${prefix}_dedupe_R1.fq.gz $WORKINGDIR/${prefix}_dedupe_R2.fq.gz

bbduk.sh -Xmx$maxmemory in=$SAMPLEFILE1 in2=$SAMPLEFILE2 out=$WORKINGDIR/${prefix}_trimmed_R1.fq out2=$WORKINGDIR/${prefix}_trimmed_R2.fq outs=$WORKINGDIR/${prefix}_trimmed_singletons.fq ktrim=r k=$kcontaminants mink=$mink hdist=$hdist qtrim=rl trimq=$phred  minlength=$minlength ref=$adaptersPath qin=$qin threads=$threads tbo tpe ow &> $WORKINGDIR/${prefix}_step2_tmp.log

#Logs some figures about sequences passing trimming
log  "BBduk's trimming stats (trimming adapters and low quality sequences): " $LOGFILE
sed -n '/Input:/,/Result:/p' $WORKINGDIR/${prefix}_step2_tmp.log >> $LOGFILE
log "" $LOGFILE
unpairedR=$(wc -l $WORKINGDIR/${prefix}_trimmed_singletons.fq | cut -d" " -f 1)
unpairedR=$(($unpairedR/4))
log  "$unpairedR singleton reads whose mate was trimmed shorter have been preserved" $LOGFILE
log "" $LOGFILE

#$WORKINGDIR/${prefix}_step2_tmp.log
rm  $WORKINGDIR/${prefix}_dedupe_R1.fq.gz $WORKINGDIR/${prefix}_dedupe_R2.fq.gz

endtime=$(date +%s.%N)
exectime=$(echo "$endtime - $starttime" | bc)
sysdate=$(date)
log "         STEP 2 terminated at $sysdate ($exectime seconds)" $LOGFILE
log "" $LOGFILE	
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STEP 3.  Quality assessment of the trimmed data file.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sysdate=$(date)
starttime=$(date +%s.%N)
log "Performing STEP 3 [Assessment of trimmed data quality] at $sysdate" $LOGFILE
log "" $LOGFILE	

### Forward reads
fastqc --quiet --noextract --format fastq --threads $threads --outdir=$WORKINGDIR $WORKINGDIR/${prefix}_trimmed_R1.fq

### Reverse reads
fastqc --quiet --noextract --format fastq --threads $threads --outdir=$WORKINGDIR $WORKINGDIR/${prefix}_trimmed_R2.fq

#gzip files
gzip $WORKINGDIR/${prefix}_trimmed_R1.fq
gzip $WORKINGDIR/${prefix}_trimmed_R2.fq

endtime=$(date +%s.%N)
exectime=$(echo "$endtime - $starttime" | bc)
sysdate=$(date)
log "           STEP 3 terminated at $sysdate ($exectime seconds)" $LOGFILE
log "" $LOGFILE	
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE	


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STEP 4. Transcript-level Quantification. Transcript-level quantification estimates from RNA-seq data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sysdate=$(date)
starttime=$(date +%s.%N)
log "Performing STEP 4 [Transcript-level Quantification] at $sysdate" $LOGFILE
log "" $LOGFILE	

#create the forlder for output of salmon
mkdir -p $WORKINGDIR/quant

salmon quant -p $threads -i $transcriptomeIndex -l A -1 $WORKINGDIR/${prefix}_trimmed_R1.fq.gz -2 $WORKINGDIR/${prefix}_trimmed_R2.fq.gz -p 1 --validateMappings -o $WORKINGDIR/quant  &> $WORKINGDIR/${prefix}_step4_tmp.log
	
#rm $WORKINGDIR/${prefix}_step4_tmp.log 

endtime=$(date +%s.%N)
exectime=$(echo "$endtime - $starttime" | bc)
sysdate=$(date)
log "           STEP 4 terminated at $sysdate ($exectime seconds)" $LOGFILE
log "" $LOGFILE	
log "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" $LOGFILE
log "" $LOGFILE	

