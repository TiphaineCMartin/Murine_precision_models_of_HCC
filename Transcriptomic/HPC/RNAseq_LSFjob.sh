#!/bin/bash
#BSUB -J "RNAseq_data[1-60]"
#BSUB -P acc_yourAccount
#BSUB -q "yourqueue"
#BSUB -n 1
#BSUB -R rusage[mem=6000]
#BSUB -W 24:00
#BSUB -o RNAseq_data%I.stdout
#BSUB -eo RNAseq_data%I.stderr
#BSUB -L /bin/bash
#! #BSUB -x  no to put it

##not define the partition #BSUB -m mothra,manda,bondo
#  Transcirptomic Analysis Pipeline (TAP) version 0.1.0 -- 20180509
#  Copyright (C) 2017 	Dr Tiphaine Martin	(developer, contributor) 
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load Module 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module load fastqc/0.11.8 
module load bbtools/37.53
module load salmon/1.0.0

#Where the fastq data are saved
FASTQDIR="/yourmainfolder/fastq"
FASTQNAMES=`find $FASTQDIR -type f | sed 's!.*/!!' | sed  's/_R.*_001.fastq.gz//g' | sort | uniq`
set -f                      # avoid globbing (expansion of *).
FASTQNAMESARRAY=(${FASTQNAMES// / })

#Used as prefix for all the results
#array start from 0 but the jobID start from 1 !!
prefixID=`expr $LSB_JOBINDEX - 1`
prefix="${FASTQNAMESARRAY[$prefixID]}"

#Where the results should be saved
WORKINGDIR_orign="/yourmainfolder/results"

#Where the job should be saved
JOBDIR="/yourmainfolder/scripts"
jobtmp=$JOBDIR"/"${prefix}_jobtmp.sh

#Where the results should be saved
WORKINGDIR=$WORKINGDIR_orign"/"$prefix
mkdir -p $WORKINGDIR

#The paired-end compressed raw data (FASTQ):
#SAMPLEFILE1 is the forward strand, 
#SAMPLEFILE2 is the reverse strand
SAMPLEFILE1=$FASTQDIR"/"${prefix}_R1_001.fastq.gz
SAMPLEFILE2=$FASTQDIR"/"${prefix}_R2_001.fastq.gz

cd $JOBDIR

#source RNAseqAnalysis_parameters.sh
cd $JOBDIR
./RNAseqAnalysis_pipeline.sh $prefix $JOBDIR $WORKINGDIR_orign $FASTQDIR $WORKINGDIR $SAMPLEFILE1 $SAMPLEFILE2


