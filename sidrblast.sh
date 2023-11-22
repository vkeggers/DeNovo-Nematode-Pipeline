#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_blastnt.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL

module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k


#makeblastdb -in /home/data/jfierst/veggers/nt_db/nt.fa -dbtype nucl -out nt

GENOME=/home/data/jfierst/veggers/Oscheius/DF5033/DF5033_hifiAssembly/DF5033ONTpb.asm.bp.p_ctg.fa
NT=/home/data/jfierst/veggers/nt_db/nt
NR=

blastx -query ${GENOME} \
    -db ${NR} \
    -outfmt 6 \
    -culling_limit 5 \
    -evalue 1e-25 \
    -out SIDRblast.txt

