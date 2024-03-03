#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH -n 8
#SBATCH --output=out_sidr.log
#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=ALL


############################
############################
########LOAD_MODULES########
############################
############################

module load miniconda3-23.5.2
module load samtools-1.15.1-gcc-8.2.0
module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k
module load minimap2-2.24
module load bwa-0.7.17-gcc-8.2.0-qgdird7
module load bedtools2-2.27.1-gcc-8.2.0-bxmhnwb


#############################
#############################
#######ARGUMENT_CENTER#######
#############################
#############################

# Initialize variables for arguments
FORWARD=""
REVERSE=""
GENOME=""
HIFI_READS=""
ONT_READS=""
NT=""
RUN_BLAST=""
BLAST_OUTPUT=""
GENUS=""

# Help message with descriptions of all arguments
print_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --forward/-f       A path to your RNA seq forward reads."
    echo "  --reverse/-r       A path to your RNA seq reverse reads."
    echo "  --genome/-g        A path to your genome assembly as a .fa file."
    echo "  --genus/-x         *REQUIRED* The genus name of your species of interest (Case- and spelling-sensitive; Ex: -g Drosophila)."
    echo "  --hifi             A path to your PacBio HiFi reads."
    echo "  --ont              A path to your Oxford Nanopore basecalled reads."
    echo "  --nt/-n            A path to your NCBI compiled NT database."
    echo "  --blast/-b         *REQUIRED* Do you need to run a blast search of your assembly? Yes/No [Yes (if not completed)/ No (previously completed)]"
    echo "  --blastoutput      Path to your blast output file."
    echo "  --help/-h          Display this help and exit."
}

# Check if no arguments were provided
if [ $# -eq 0 ]; then
    print_help
    exit 1
fi

# Loop through arguments and process them
while [ "$1" != "" ]; do
    case $1 in
        -f | --forward )     shift
                             FORWARD=$1
                             ;;
        -r | --reverse )     shift
                             REVERSE=$1
                             ;;
        -g | --genome )      shift
                             GENOME=$1
                             ;;
        -x | --genus )		 shift
                             GENUS=$1
                             ;;
        --hifi )             shift
                             HIFI_READS=$1
                             ;;
        --ont )              shift
                             ONT_READS=$1
                             ;;
        -n | --nt )          shift
                             NT=$1
                             ;;
        -b | --blast )       shift
                             RUN_BLAST=$1
                             if [ "$RUN_BLAST" != "Yes" ] && [ "$RUN_BLAST" != "No" ]; then
                                 echo "Invalid argument for --blast/-b: $RUN_BLAST"
                                 print_help
                                 exit 1
                             fi
                             ;;
        --blastoutput )      shift
                             BLAST_OUTPUT=$1
                             ;;
        -h | --help )        print_help
                             exit
                             ;;
        * )                  echo "Invalid argument: $1"
                             print_help
                             exit 1
    esac
    shift
done

# Check if --genus/-x was provided
if [ -z "$GENUS" ]; then
    echo "Error: --genus/-x argument is required."
    print_help
    exit 1
fi

# Check if --blast/-b was provided
if [ -z "$RUN_BLAST" ]; then
    echo "Error: --blast/-b argument is required."
    print_help
    exit 1
fi

# Check for blast and corresponding argument conditions
if [ "$RUN_BLAST" = "No" ] && [ -z "$NT" ]; then
    echo "Error: --nt/-n argument must be given if --blast/-b is set to 'No'."
    exit 1
elif [ "$RUN_BLAST" = "Yes" ] && [ -z "$BLAST_OUTPUT" ]; then
    echo "Error: --blastoutput argument must be given if --blast/-b is set to 'Yes'."
    exit 1
fi

########RUN_BLAST########

if [ "$RUN_BLAST" = "No" ]; then
    sbatch sidrblash.sh --genome $GENOME --nt $NT


########GENERATE_ALIGNMENTS########

mkdir samsANDbams
cd samsANDbams

#pacbio reads aligned to assembly
minimap2 -ax map-hifi \
    ${GENOME} \
    ${HIFI_READS} \
    -o PBaln.sam

#ONT reads aligned to assembly
minimap2 -ax map-ont \
    ${GENOME} \
    ${ONT_READS} \
    -o ONTaln.sam

#RNA reads aligned to assembly
bwa index ${GENOME}
bwa mem -t 4 ${GENOME} ${FORWARD} ${REVERSE} > RNAaln.sam

########MAKE_BAM_FILES########

#DNA
samtools view -Sb ./PBaln.sam -o ./PBaln.bam
samtools sort -o ./PBaln_sorted.bam ./PBaln.bam
samtools index ./PBaln_sorted.bam ./PBaln_index.bai

samtools view -Sb ./ONTaln.sam -o ./ONTaln.bam
samtools sort -o ./ONTaln_sorted.bam ./ONTaln.bam
samtools index ./ONTaln_sorted.bam ./ONTaln_index.bai


#RNA
samtools view -Sb ./RNAaln.sam -o ./RNAaln.bam
samtools sort -o ./RNAaln_sorted.bam ./RNAaln.bam
samtools index ./RNAaln_sorted.bam ./RNAaln_index.bai


########GENERATE_STATS########

cd ..
mkdir stats
cd stats

#plus and minus strand counts
samtools view -F 16 ./../samsANDbams/PBaln_sorted.bam | awk '{if ($3 != prev) {if (prev != "") print prev, count; prev = $3; count = 0} count++} END {if (prev != "") print prev, count}' > PBplus_strand_counts.txt
sed -i 's/ /\t/' PBplus_strand_counts.txt
sort -k1 PBplus_strand_counts.txt > PBplus_strand_counts.txt.temp
mv PBplus_strand_counts.txt.temp PBplus_strand_counts.txt
echo -e "contig\tPBplus_strand_counts" > header.txt
cat PBplus_strand_counts.txt >> header.txt
mv header.txt PBplus_strand_counts.txt

samtools view -f 16 ./../samsANDbams/PBaln_sorted.bam | awk '{if ($3 != prev) {if (prev != "") print prev, count; prev = $3; count = 0} count++} END {if (prev != "") print prev, count}' > PBminus_strand_counts.txt
sed -i 's/ /\t/' PBminus_strand_counts.txt
sort -k1 PBminus_strand_counts.txt > PBminus_strand_counts.txt.temp
mv PBminus_strand_counts.txt.temp PBminus_strand_counts.txt
echo -e "contig\tPBminus_strand_counts" > header.txt
cat PBminus_strand_counts.txt >> header.txt
mv header.txt PBminus_strand_counts.txt

samtools view -F 16 ./../samsANDbams/ONTaln_sorted.bam | awk '{if ($3 != prev) {if (prev != "") print prev, count; prev = $3; count = 0} count++} END {if (prev != "") print prev, count}' > ONTplus_strand_counts.txt
sed -i 's/ /\t/' ONTplus_strand_counts.txt
sort -k1 ONTplus_strand_counts.txt > ONTplus_strand_counts.txt.temp
mv ONTplus_strand_counts.txt.temp ONTplus_strand_counts.txt
echo -e "contig\tONTplus_strand_counts" > header.txt
cat ONTplus_strand_counts.txt >> header.txt
mv header.txt ONTplus_strand_counts.txt

samtools view -f 16 ./../samsANDbams/ONTaln_sorted.bam | awk '{if ($3 != prev) {if (prev != "") print prev, count; prev = $3; count = 0} count++} END {if (prev != "") print prev, count}' > ONTminus_strand_counts.txt
sed -i 's/ /\t/' ONTminus_strand_counts.txt
sort -k1 ONTminus_strand_counts.txt > ONTminus_strand_counts.txt.temp
mv ONTminus_strand_counts.txt.temp ONTminus_strand_counts.txt
echo -e "contig\tONTminus_strand_counts" > header.txt
cat ONTminus_strand_counts.txt >> header.txt
mv header.txt ONTminus_strand_counts.txt

#various forms of coverage
covered_bases[$1]++ is the number of positions (bases) that have coverage greater than zero
sum[$1] += $3 is depth at each position (takes the third column and adds it to sum)
total_bases[$1]++ is basically the length of the contig by adding to the count even if the base does not have coverage

#pacbio
samtools depth ./../samsANDbams/PBaln_sorted.bam > PBcoverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum) {
        avg_fold = sum[contig] / total_bases[contig];
        coverage_percentage = (covered_bases[contig] / total_bases[contig]) * 100;


        print contig, covered_bases[contig], avg_fold, coverage_percentage;
    }
}' PBcoverage.txt | sort -k1,1 > PBcoverageStats.txt
sed -i 's/ /\t/g' PBcoverageStats.txt
sort -k1 PBcoverageStats.txt > PBcoverageStats.txt.temp
mv PBcoverageStats.txt.temp PBcoverageStats.txt 
echo -e "contig\tPB_Covered_bases\tPB_avg_fold\tPB_coverage_percent" > header.txt
cat PBcoverageStats.txt >> header.txt
mv header.txt PBcoverageStats.txt

#ONT
samtools depth ./../samsANDbams/ONTaln_sorted.bam > ONTcoverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum) {
        avg_fold = sum[contig] / total_bases[contig];

        print contig, covered_bases[contig], avg_fold;
    }
}' ONTcoverage.txt | sort -k1,1 > ONTcoverageStats.txt
sed -i 's/ /\t/g' ONTcoverageStats.txt
sort -k1 ONTcoverageStats.txt > ONTcoverageStats.txt.temp
mv ONTcoverageStats.txt.temp ONTcoverageStats.txt
echo -e "contig\tONT_Covered_bases\tONT_avg_fold" > header.txt
cat ONTcoverageStats.txt >> header.txt
mv header.txt ONTcoverageStats.txt

#RNA
samtools depth ./../samsANDbams/RNAaln_sorted.bam > RNAcoverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum) {
        avg_fold = sum[contig] / total_bases[contig];

        print contig, covered_bases[contig], avg_fold;
    }
}' RNAcoverage.txt | sort -k1,1 > RNAcoverageStats.txt
sed -i 's/ /\t/g' RNAcoverageStats.txt
sort -k1 RNAcoverageStats.txt > RNAcoverageStats.txt.temp
mv RNAcoverageStats.txt.temp RNAcoverageStats.txt
echo -e "contig\tRNA_Covered_bases\tRNA_avg_fold" > header.txt
cat RNAcoverageStats.txt >> header.txt
mv header.txt RNAcoverageStats.txt

bash contig_GC.sh

paste contig_GC.txt RNAcoverageStats.txt | awk '{print $1, ($5 / $2) * 100}' | sed -i 's/ /\t/g' | sort -k1 > RNA_Covered_percent.txt
echo -e "contig\tRNA_Coverage_percent" > header.txt
cat RNA_Covered_percent.txt >> header.txt
mv header.txt RNA_Covered_percent.txt 

paste Ref_GC.txt ONTcoverageStats.txt | awk '{print $1, ($5 / $2) * 100}' | sed -i 's/ /\t/g' | sort -k1 > ONT_Covered_percent.txt
echo -e "contig\tONT_Coverage_percent" > header.txt
cat ONT_Covered_percent.txt >> header.txt
mv header.txt ONT_Covered_percent.txt

#calculate GC content for reads over contig region
#hifi
mkdir PBgc_per_contig
cd PBgc_per_contig


#get contig names in list
samtools view -H ./../../samsANDbams/PBaln_sorted.bam | grep '@SQ' | cut -f 2 -d ':' | cut -f 2 -d '@' > contig_names.txt

#fix names (delete LN and ?)
awk '{print $1}' contig_names.txt > fixed_names.txt
bedtools bamtobed -i ./../../samsANDbams/PBaln_sorted.bam > output.bed

#for each contig in the fixed_names.txt file, make a temporary bed, take the reads in that bed region and input to a fastq, and then calculate GC by countin#g GC and dividing by total bases
while read -r contig; do
    grep "${contig}" output.bed > "${contig}"temp_output.bed
    samtools view -L "${contig}"temp_output.bed -h ./../../samsANDbams/PBaln_sorted.bam | samtools view -b - > "${contig}"reads.bam
    rm "${contig}"temp_output.bed
    bedtools bamtofastq -i "${contig}"reads.bam -fq "${contig}"reads.fastq
    awk 'NR%4 == 2 {print $1}' "${contig}"reads.fastq | grep -o -i '[GC]' | wc -l > gc_count.txt
    total_bases=$(awk 'NR%4 == 2 {print $1}' "${contig}"reads.fastq | tr -d '\n' | wc -c)
    gc_content=$(awk '{printf "%.2f", $1 / '$total_bases'}' gc_count.txt)
    echo -e "$contig\t$gc_content" >> PBread_gc.txt
    rm "${contig}"reads.bam
    rm "${contig}"reads.fastq
done < fixed_names.txt
sort -k1 PBread_gc.txt > PBread_gc.txt.temp
mv PBread_gc.txt.temp PBread_gc.txt
cp PBread_gc.txt ./../.
cd ..
echo -e "contig\tPBread_gc" > header.txt
cat PBread_gc.txt >> header.txt
mv header.txt PBread_gc.txt

#ONT
mkdir ONTgc_per_contig
cd ONTgc_per_contig

#calculate GC content for reads over contig region
#get contig names in list
samtools view -H ./../../samsANDbams/ONTaln_sorted.bam | grep '@SQ' | cut -f 2 -d ':' | cut -f 2 -d '@' > contig_names.txt

#fix names (delete LN and ?)
awk '{print $1}' contig_names.txt > fixed_names.txt
bedtools bamtobed -i ./../../samsANDbams/ONTaln_sorted.bam > output.bed

#for each contig in the fixed_names.txt file, make a temporary bed, take the reads in that bed region and input to a fastq, and then calculate GC by countin#g GC and dividing by total bases
while read -r contig; do
    grep "${contig}" output.bed > "${contig}"temp_output.bed
    samtools view -L "${contig}"temp_output.bed -h ./../../samsANDbams/ONTaln_sorted.bam | samtools view -b - > "${contig}"reads.bam
    rm "${contig}"temp_output.bed
    bedtools bamtofastq -i "${contig}"reads.bam -fq "${contig}"reads.fastq
    awk 'NR%4 == 2 {print $1}' "${contig}"reads.fastq | grep -o -i '[GC]' | wc -l > gc_count.txt
    total_bases=$(awk 'NR%4 == 2 {print $1}' "${contig}"reads.fastq | tr -d '\n' | wc -c)
    gc_content=$(awk '{printf "%.2f", $1 / '$total_bases'}' gc_count.txt)
    echo -e "$contig\t$gc_content" >> ONTread_gc.txt
    rm "${contig}"reads.bam
    rm "${contig}"reads.fastq
done < fixed_names.txt
sort -k1 ONTread_gc.txt > ONTread_gc.txt.temp
mv ONTread_gc.txt.temp ONTread_gc.txt
cp ONTread_gc.txt ./../.
cd ..
echo -e "contig\tONTread_gc" > header.txt
cat ONTread_gc.txt >> header.txt
mv header.txt ONTread_gc.txt

#repeat of RNA
mkdir RNAgc_per_contig
cd RNAgc_per_contig

#get contig names in list
samtools view -H ./../../samsANDbams/RNAaln_sorted.bam | grep '@SQ' | cut -f 2 -d ':' | cut -f 2 -d '@' > contig_names.txt

#fix names (delete LN and ?)
awk '{print $1}' contig_names.txt > fixed_names.txt
bedtools bamtobed -i ./../../samsANDbams/RNAaln_sorted.bam > output.bed


while read -r contig; do
    grep "${contig}" output.bed > "${contig}"temp_output.bed
    samtools view -L "${contig}"temp_output.bed -h ./../../samsANDbams/RNAaln_sorted.bam | samtools view -b - > "${contig}"reads.bam
    rm "${contig}"temp_output.bed
    bedtools bamtofastq -i "${contig}"reads.bam -fq "${contig}"reads.fastq
    awk 'NR%4 == 2 {print $1}' "${contig}"reads.fastq | grep -o -i '[GC]' | wc -l > gc_count.txt
    total_bases=$(awk 'NR%4 == 2 {print $1}' "${contig}"reads.fastq | tr -d '\n' | wc -c)
    gc_content=$(awk '{printf "%.2f", $1 / '$total_bases'}' gc_count.txt)
    echo -e "$contig\t$gc_content" >> RNAread_gc.txt
    rm "${contig}"reads.bam
    rm "${contig}"reads.fastq
done < fixed_names.txt
cp RNAread_gc.txt ./../.
sort -k1 RNAread_gc.txt > RNAread_gc.txt.temp 
mv RNAread_gc.txt.temp RNAread_gc.txt
cd ..
cd ..
echo -e "contig\tRNAread_gc" > header.txt
cat RNAread_gc.txt >> header.txt
mv header.txt RNAread_gc.txt


########FIX_BLAST_OUTPUT########


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

########MAKE_TABLE########

cd ./stats/
rm *coverage.txt


ls *.txt > list
join --check-order --header -t$'\t' blastIDs.txt Ref_GC.txt > newfile
sed -i '/Ref_GC/Id' list
sed -i '/blastIDs/Id' list

while read -r file; do
    join --check-order --header -t$'\t' newfile "${file}" > newfile.temp
    mv newfile.temp newfile
done < list

rm list
mv newfile SIDRstats.tsv

########RUN_SIDR########
python xgboost.py

#make a new directory
mkdir test
# copy fasta and deny/allow lists into their own folder
cp ${GENOME} test/
cp keptContigs.csv  test/
cp contaminantContigs.csv  test/
cd test

## Splits contigs into individual files
awk '/^>/ {OUT=substr($0,2,14) ".fasta"}; OUT {print >OUT}' *.fa

##remove original fasta file
rm *.fa

##fix file
awk -F ',' '{print $2}' contaminantContigs.csv | sed 's/"//g' > contaminantlist.txt

awk -F ',' '{print $2}' keptContigs.csv | sed 's/"//g' > keptlist.txt

## combines files that match the allow list into a new “keptcontigs” file
cat keptlist.txt  | while read line; do cat ""$line".fasta" ; done > keptContigs.fasta

## combines file that match the deny list into a new “contaminant contigs” file
cat contaminantlist.txt  | while read line; do cat ""$line".fasta" ; done > contaminantContigs2.fasta


## move your kept contigs and contaminant contigs files out of the directory. and delete the split out contigs.
mv contaminantContigs.fasta ./../.
mv keptContigs.fasta ./../.
rm -r test















