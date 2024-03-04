#!/bin/bash


#various forms of coverage



cd stats

#########################################################################################################################################

#pacbio hifi

#titles for columns in the output file
echo -e "contig\tPB_average_coverage\tPB_coverage_percent" > PBcoverageStats.txt

#samtools depth will give you the number of reads that cover each nucleotide position in the genome
samtools depth ./../samsANDbams/PBaln_sorted.bam > PBcoverage.txt

#If column 3 of the samtools depth output (coverage at that nucleotide) is greater than zero, add 1 to covered_bases and add the value to the column to sum_
#coverage. If column 3 is not greater than zero, just add 1 to total_bases, which is basically just the length of the contig
#The END statement is a for loop that calculates and prints the average coverage over the whole contig and the coverage percent (telling you if any bases we
#re not covered by the reads) for each contig. These values are appended to the output file (PBcoverageStats.txt) which we initialized earlier
awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' PBcoverage.txt >> PBcoverageStats.txt

#replaces spacecs with tabs, will be important when concatenating the tables together at the end
sed -i 's/ /\t/g' PBcoverageStats.txt


###################################################################################################################################

#oxford nanopore

echo -e "contig\tONT_average_coverage\tONT_coverage_percent" > ONTcoverageStats.txt
samtools depth ./../samsANDbams/ONTaln_sorted.bam > ONTcoverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' ONTcoverage.txt >> ONTcoverageStats.txt
sed -i 's/ /\t/g' ONTcoverageStats.txt


####################################################################################################################################

#RNA

echo -e "contig\tRNA_average_coverage\tRNA_coverage_percent" > RNAcoverageStats.txt
samtools depth ./../samsANDbams/RNAaln_sorted.bam > RNAcoverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' RNAcoverage.txt >> RNAcoverageStats.txt
sed -i 's/ /\t/g' RNAcoverageStats.txt

####################################################################################################################################

cd ..
