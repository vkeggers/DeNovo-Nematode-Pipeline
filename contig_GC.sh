#!/bin/bash


##calculate GC content of contigs, or contig_GC


#set genome variable
GENOME=/path/to/genome

#titles for columns in the output file
echo -e "contig\tlength\tgc_percent" > contig_GC.txt

#print contig names to list without the ">"
awk 'NR%2==1 {print substr($0,2) > "contig_names.txt"}' "${GENOME}"

#loop through each contig in contig_names.txt
while read -r contig;
   do
        #grep the line with the contig name, plus the following line | print only the following line |
        #substitute any G or C or g or c to an empty character, counting the number of substitutions and save as the variable gc_number.
        gc_number=$(grep -A 1 "${contig}" "${GENOME}" | awk 'NR%2==0 {print}' | awk '{count += gsub(/[gcGC]/, "");} END {print count}')

        #grep the line with the contig name, plus the following line | print only the following line |
        #remove the new line character | count all characters and save as the variable total_bases.
        total_bases=$(grep -A 1 "${contig}" "${GENOME}" | awk 'NR%2==0 {print}' | tr -d '\n' | wc -c)

        #pass the variables to awk, divide gc_number by total_bases and multiply by 100
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')

        #print and append the values for the current contig to the output file
        echo -e "$contig\t$total_bases\t$gc_percent" >> contig_GC.txt
   done < contig_names.txt
