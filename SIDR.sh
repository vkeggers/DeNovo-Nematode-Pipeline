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


#########################
#########################
########RUN_BLAST########
#########################
#########################

if [ "$RUN_BLAST" = "No" ]; then
    sbatch sidrblash.sh --genome $GENOME --nt $NT


###################################
###################################
########GENERATE_ALIGNMENTS########
###################################
###################################

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


##############################
##############################
########MAKE_BAM_FILES########
##############################
##############################

#pacbio hifi
samtools view -Sb ./PBaln.sam -o ./PBaln.bam
samtools sort -o ./PBaln_sorted.bam ./PBaln.bam
samtools index ./PBaln_sorted.bam ./PBaln_index.bai

#oxford nanopore
samtools view -Sb ./ONTaln.sam -o ./ONTaln.bam
samtools sort -o ./ONTaln_sorted.bam ./ONTaln.bam
samtools index ./ONTaln_sorted.bam ./ONTaln_index.bai

#RNA
samtools view -Sb ./RNAaln.sam -o ./RNAaln.bam
samtools sort -o ./RNAaln_sorted.bam ./RNAaln.bam
samtools index ./RNAaln_sorted.bam ./RNAaln_index.bai


##############################
##############################
########GENERATE_STATS########
##############################
##############################

cd ..
mkdir stats

#calculate GC content for reads over contig region
bash contig_GC.sh


#GC percent of contigs
bash contig_GC.sh


#plus and minus strand counts
bash plusANDminusCounts.sh


#various forms of coverage
bash coverage.sh


##########################
##########################
########MAKE_TABLE########
##########################
##########################

cd stats
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


########################
########################
########RUN_SIDR########
########################
########################

python xgboost.py


################################
################################
#######MAKE_OUTPUT_FILES########
################################
################################

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















