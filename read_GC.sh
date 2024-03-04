#!/bin/bash



#calculate average GC content of the reads spanning the contig


#load modules
module load bedtools2-2.27.1-gcc-8.2.0-bxmhnwb
module load samtools-1.15.1-gcc-8.2.0


cd stats

########################################################################################################################################################

#pacbio hifi

mkdir PBgc_per_contig
cd PBgc_per_contig

#titles for columns in the output file
echo -e "contig\tPBread_GC" > PBread_GC.txt

#using bedtools, make a bed file showing the beginning and end of each read and which contig it belongs to
bedtools bamtobed -i ./../../samsANDbams/PBaln_sorted.bam > output.bed

#loop through each contig in contig_names.txt
while read -r contig;
    do

        #grep all the reads belonging to the current contig and put those in a temporary bed file. We are making a subset of the original bed file
        grep "${contig}" output.bed > "${contig}"_output.bed

        #give us back the sequences for these reads by taking the subset bed file and making a bam file of it
        #-L means to only output alignments that overlap with the region in the bed file and -h means to include headers | shuttle that sam output to a bam
        samtools view -L "${contig}"_output.bed -h ./../../samsANDbams/PBaln_sorted.bam | samtools view -b - > "${contig}"reads.bam

        #create a subset fastq from the subset bam, because I don't know how to work with sams and bams
        bedtools bamtofastq -i "${contig}"reads.bam -fq "${contig}"reads.fastq

        #use awk to print the second of four lines (2nd line is sequence) in the subset fastq | substitute any G or C or g or c to an empty character and
        #count the number of substitutions. Save this count to a variable called gc_number
        gc_number=$(awk 'NR%4 == 2 {print}' "${contig}"reads.fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')

        #use awk to print the second of four lines (2nd line is sequence) in the subset fastq | remove newline characters | count all characters and
        #save to a variable called total bases
        total_bases=$(awk 'NR%4 == 2 {print}' "${contig}"reads.fastq | tr -d '\n' | wc -c)

        #pass the variables to awk, divide gc_number by total_bases and multiply by 100
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')

        #print and append the values for the current contig to the output file
        echo -e "$contig\t$gc_content" >> PBread_GC.txt

        #remove subset files
        rm "${contig}"_output.bed
        rm "${contig}"reads.bam
        rm "${contig}"reads.fastq
    done < ./../contig_names.txt
cp PBread_GC.txt ./../.
cd ..


###################################################################################################################################################

#Oxford Nanopore

mkdir ONTgc_per_contig
cd ONTgc_per_contig

echo -e "contig\tONTread_GC" > ONTread_GC.txt

bedtools bamtobed -i ./../../samsANDbams/ONTaln_sorted.bam > output.bed

while read -r contig;
    do
        grep "${contig}" output.bed > "${contig}"_output.bed
        samtools view -L "${contig}"_output.bed -h ./../../samsANDbams/ONTaln_sorted.bam | samtools view -b - > "${contig}"reads.bam
        bedtools bamtofastq -i "${contig}"reads.bam -fq "${contig}"reads.fastq
        gc_number=$(awk 'NR%4 == 2 {print}' "${contig}"reads.fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')
        total_bases=$(awk 'NR%4 == 2 {print}' "${contig}"reads.fastq | tr -d '\n' | wc -c)
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')
        echo -e "$contig\t$gc_percent" >> ONTread_GC.txt
        rm "${contig}"_output.bed
        rm "${contig}"reads.bam
        rm "${contig}"reads.fastq
    done < ./../contig_names.txt
cp ONTread_GC.txt ./../.
cd ..


##################################################################################################################################################

#RNA

mkdir RNAgc_per_contig
cd RNAgc_per_contig

echo -e "contig\tRNAread_GC" > RNAread_GC.txt

bedtools bamtobed -i ./../../samsANDbams/RNAaln_sorted.bam > output.bed

while read -r contig;
    do
        grep "${contig}" output.bed > "${contig}"_output.bed
        samtools view -L "${contig}"_output.bed -h ./../../samsANDbams/RNAaln_sorted.bam | samtools view -b - > "${contig}"reads.bam
        bedtools bamtofastq -i "${contig}"reads.bam -fq "${contig}"reads.fastq
        gc_number=$(awk 'NR%4 == 2 {print}' "${contig}"reads.fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')
        total_bases=$(awk 'NR%4 == 2 {print}' "${contig}"reads.fastq | tr -d '\n' | wc -c)
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')
        echo -e "$contig\t$gc_percent" >> RNAread_GC.txt
        rm "${contig}"_output.bed
        rm "${contig}"reads.bam
        rm "${contig}"reads.fastq
    done < ./../contig_names.txt
cp RNAread_GC.txt ./../.
cd ..


###################################################################################################################################################

cd ..
