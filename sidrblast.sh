#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_blastnt.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k


#makeblastdb -in /home/data/jfierst/veggers/nt_db/nt.fa -dbtype nucl -out nt

GENOME=/home/data/jfierst/veggers/Oscheius/DF5033/DF5033_hifiAssembly/DF5033ONTpb.asm.bp.p_ctg.fa
NT=/home/data/jfierst/veggers/nt_db/nt

blastn -query ${GENOME} \
    -db ${NT} \
   # -outfmt 6 \
    -culling_limit 5 \
    -evalue 1e-25 \
    -out SIDRblast.txt


################################
################################
########FIX_BLAST_OUTPUT########
################################
################################


awk '/Query=/{print; flag=1; next} flag && /^>/{print; flag=0} flag && /No hits found/{print; flag=0}' out_blastnt.log > test1.txt
awk 'NR%2==1{col1=$0} NR%2==0{print col1, $0}' test1.txt > test2.txt
awk '{print $2, $5, $6}' test2.txt > test3.txt
sed -i 's/hits found/No hits found/g' test3.txt
sed 's/ /\t/' test3.txt > blastIDs.txt
rm test*
sort -k1 blastIDs.txt > blastIDs.txt.temp
mv blastIDs.txt.temp blastIDs.txt
cp blastIDs.txt ./stats/.
echo -e "contig\tOrigin" > header.txt
cat blastIDs.txt >> header.txt
mv header.txt blastIDs.txt
