#!/bin/bash



#calculate the number of reads on the plus(forward) strand of DNA and minus(reverse) strand



#load modules
module load samtools-1.15.1-gcc-8.2.0


cd stats

###########################################################################################################################################

#pacbio hifi

#titles for columns in the output file
echo -e "contig\tPB_plus_counts" > PB_plus_counts.txt

#-F means to exclude, 16 is a flag meaning reverse complement, so we are excluding all reads which are associated with the reverse strand of DNA
#the awk script is looking at the 3rd field (contig name), and if the field matches, it increases the count. If the field does not match (new contig), then
#it prints the values for contig and count, then resets them. The END is to ensure the last line is processed and printed.
samtools view -F 16 ./../samsANDbams/PBaln_sorted.bam | awk '{if ($3 != contig) {if (contig != "") print contig, count; contig = $3; count = 0} count++} END
 {if (contig != "") print contig, count}' >> PB_plus_counts.txt

#replaces spaces with tabs, important for concatenating the table together later
sed -i 's/ /\t/' PB_plus_counts.txt

#repeat the above but for reverse strand
echo -e "contig\tPB_minus_counts" > PB_minus_counts.txt
samtools view -f 16 ./../samsANDbams/PBaln_sorted.bam | awk '{if ($3 != contig) {if (contig != "") print contig, count; contig = $3; count = 0} count++} END {if (contig != "") print contig, count}' >> PB_minus_counts.txt
sed -i 's/ /\t/' PB_minus_counts.txt

###########################################################################################################################################

#oxford nanopore

echo -e "contig\tONT_plus_counts" > ONT_plus_counts.txt
samtools view -F 16 ./../samsANDbams/ONTaln_sorted.bam | awk '{if ($3 != contig) {if (contig != "") print contig, count; contig = $3; count = 0} count++} END {if (contig != "") print contig, count}' > ONT_plus_counts.txt
sed -i 's/ /\t/' ONT_plus_counts.txt


echo -e "contig\tONT_minus_counts" > ONT_minus_counts.txt
samtools view -f 16 ./../samsANDbams/ONTaln_sorted.bam | awk '{if ($3 != contig) {if (contig != "") print contig, count; contig = $3; count = 0} count++} END {if (contig != "") print contig, count}' > ONT_minus_counts.txt
sed -i 's/ /\t/' ONT_minus_counts.txt

###########################################################################################################################################

#RNA

echo -e "contig\tONT_plus_counts" > ONT_plus_counts.txt
samtools view -F 16 ./../samsANDbams/RNAaln_sorted.bam | awk '{if ($3 != contig) {if (contig != "") print contig, count; contig = $3; count = 0} count++} END {if (contig != "") print contig, count}' > RNA_plus_counts.txt
sed -i 's/ /\t/' RNA_plus_counts.txt


echo -e "contig\tRNA_minus_counts" > RNA_minus_counts.txt
samtools view -f 16 ./../samsANDbams/RNAaln_sorted.bam | awk '{if ($3 != contig) {if (contig != "") print contig, count; contig = $3; count = 0} count++} END {if (contig != "") print contig, count}' > RNA_minus_counts.txt
sed -i 's/ /\t/' RNA_minus_counts.txt

############################################################################################################################################

cd ..
