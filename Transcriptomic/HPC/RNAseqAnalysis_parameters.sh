#!/bin/bash

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
# General Parameters
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Number of threads to be used by the applications, and max memory allocated
threads=1
maxmemory=4G

# Input quality offset: 33 (ASCII+33) or 64 (ASCII+64)
# ASCII+33 is used by the Sanger/Illumina 1.9+ encoding, whilst ASCII+64 by the
# Illumina 1.5 encoding 
qin=33  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load Module 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module load fastqc/0.11.8 
module load bbtools/37.53
module load salmon/0.13.1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Database parameters
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Path to the BBmap FASTA file for adapter (available when installing BBmap, or at
#https://github.com/BioInfoTools/BBMap/blob/master/resources/adapters.fa)
adaptersPath="/sc/hydra/projects/parsonslab/software/BBmap/adapters.fa"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Trimming parameters (see BBduk guide for details)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Regions with average quality BELOW this will be trimmed  
phred=10 
#Reads shorter than this after trimming will be discarded
minlength=60 
#Shorter kmers at read tips to look for 
mink=11 
#Maximum Hamming distance for ref kmers   
hdist=1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Decontamination parameters (see BBduk guide for details)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#k-mer length used for finding contaminants     
kcontaminants=23
#Approximate minimum alignment identity to look for 
mind=0.95
#Longest indel to look for
maxindel=3
#Restrict alignment band to this
bwr=0.16 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Quantification parameters (see Salmon guide for details)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## reference data for TCGA 
#https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
#Where the file of transcriptome reference (fasta and fai) for quasi-mapping
## salmon index -t $transcriptomeFasta -i salmon_0.9.1/ --type quasi -k 31
transcriptomeFasta="/yourmainfolder/reference/Mus_musculus.GRCm38.cds.all.fa.gz"
transcriptomeIndex="/yourmainfolder/reference/salmon_0.12.0"

