# DeNovo-Nematode-Pipeline
This pipeline is in reference to and building off of [CRE@UA](https://github.com/BamaComputationalBiology/CRE-UA/blob/main/CRE-Pipeline.md)
Everything was done using the FIU HPC. Press on arrows to expand contents.


<details>
<summary><b>Upload raw data to NCBI SRA</b></summary>

It is a good idea to upload your raw reads to the Sequence Read Archive (SRA) so that it is stored off your system and you can come back and download it if needed. You may also place an embargo on it so that the data will not be public until your paper is published. This is also a good idea because it may take a month to process and you don't want to be worried about this while also trying to publish (most journals require the raw data to be available during the time of review).

</details>



<details>
<summary><b>Data may be obtained through NCBI SRA</b></summary>

**If you have your own data, skip this part!**

<details>
<summary>Nanopore</summary>
```
module load sratoolkit-3.0.0
```

Go to NCBI SRA and search _Oscheius_. use the filters at the side to narrow it down to genome and nanopore reads. Find the sra ID for _Oscheius_ sp.G, the number is **SRR16242712**
```
fasterq-dump SRR16242712
#this will take a while and give you no feedback so just believe it will work.
```

If successful you should have a file named SRR16242712.fastq with 18G of data. Type ls -lh to see this.
</details>



<details>
<summary>Illumina</summary>
```
fastq-dump --split-3 SRR16242711
```

If successful you should have a file named SRR16242711_1.fastq and SRR16242711_2.fastq both with 5.4G of data. Type ls -lh to see this.
</details>



</details>



<details>
<summary><b>Check the Raw Data</b></summary>

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a program to assess the quality of raw reads and give some basic stats. 

```
module load fastqc-0.11.7-gcc-8.2.0-gia624n
fastqc /path/to/reads.fastq
```
The output will be a .html report. To download this from the hpc, exit ssh and log in via sftp:
```
exit
sftp username@hpclogin01.fiu.edu
get *.html
exit
```
You will then need to navigate to your home directory on your computer and open the file in a browser.

**Some questions to pay attention to:**

1. What is the smallest read?
2. What is the largest read?
3. What is the median read length?
4. What is the theoretical coverage of the genome (Size of genome/(median read legth * number of reads)) OR (size of genome/file size)
5. Are there adapter sequences you may need to trim? If there are you should use a software like [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) or [Trimgalore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) to trim the adapter sequences off and run fastqc again. Trimgalore is a module on the HPC but Trimmomatic will have to be installed by the user.

</details>



</details>



<details>
<summary><b>Assembly</b></summary>

It is a good idea to try multiple assembly methods and compare to choose the 'best' one. Best typically means most complete and contiguous. You could try with different softwares, different input data, and different amounts of input data. You can then use that 'best' one for annotation. 

We have tried assembly with:

* [flye](https://github.com/fenderglass/Flye) with [canu correct](https://canu.readthedocs.io/en/latest/quick-start.html#quickstart)

* [nextdenovo](https://github.com/Nextomics/NextDenovo)

* [verkko](https://github.com/marbl/verkko)

* [hifiasm](https://hifiasm.readthedocs.io/en/latest/faq.html)


Flye is an all-rounder, meaning it can assemble with pretty much anything you give it. nextdenovo uses ONT reads primarily. hifiasm uses pacbio hifi reads primarily. Verkko uses pacbio and ONT. 


--------------------------------------------------------------------------------------------------------------------------------------------------------------

<details>
<summary><i>nextDenovo</i></summary>

Between flye and nextdenovo, we find nextDenovo to generally be better and more contiguous.
```
#create the input file
ls SRR16242712.fastq > input.fofn
```

```
#create the configuration file for assembly
vi run.cfg
```

Press[i] for insert and copy and paste the below section (this was obtained by going to nextDenovo documentation and copying the run.cfg file. Then we correct a few lines for our data, like genome size for example. If you don't know the genome size you can estimate it from a related species or use the option auto.

```
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = raw
read_type = ont # clr, ont, hifi
input_fofn = input.fofn
workdir = PB127

[correct_option]
read_cutoff = 1k
genome_size = 120M # estimated genome size, I know because I've already assembled this one
sort_options = -m 20g -t 15
minimap2_options_raw = -t 8
pa_correction = 3
correction_options = -p 15

[assemble_option]
minimap2_options_cns = -t 8
nextgraph_options = -a 1
```

Save by pressing [esc], type ':wq' and press [enter]

```
#create the script to run nextDenovo and create an assembled genome
vi assemble.sh
```

Press [i] for insert mode and copy the below script

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH --output=out_%assemble.log
#SBATCH --mail-user=username@email.com 	#use your own email instead
#SBATCH --mail-type=ALL

module load nextDenovo-2.5.0

nextDenovo run.cfg
```

Save by pressing [esc], type ':wq' and press [enter]


Run the script with: 
```
sbatch assemble.sh
```

To see if your job is running type the following command:
```
squeue --me
```

There is a common issue some face and you may need to load modules before you run the script. In which case use:
```
module load nextDenovo-2.5.0
sbatch assemble.sh
```

The final assembly result is at 03.ctg_graph/nd.asm.fasta

Basic statistics for the assembly are at 03.ctg_graph/nd.asm.fasta.stat
</details>

<details>
<summary><i>Flye</i></summary>

The Canu module is available on HPC but I run into a problem with java when trying to use the module. Additionally, Flye is not available, so let's just conda install them:

Get Canu
```
module load mamba/23.1.0-4
conda create -n canu
source activate canu
conda install -c bioconda canu
```

Create script
```
vi canu_correction.sh
```
Hit [i] for insertion mode and copy/paste the following:

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH --output=out_%canu_correct.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

source activate canu

canu -correct -p PB127_canu -d canu_out genomeSize=120M useGrid=false -nanopore-raw ./SRR16242712.fastq
```

Save by pressing [esc], type ':wq' and press [enter]


Run the script with: 
```
sbatch canu_correction.sh
```

To see if your job is running type the following command:
```
squeue --me
```

This job took 2.5 days to finish, but could be sped up by giving it more resources. Try adding "#SBATCH -n 8" and "#SBATCH --mem=128G" to the script.

The output is in canu_out. The corrected reads are the file: *.correctedReads.fasta.gz

However, an error is thrown because some of the read names match in the first column. To fix this we unzip the file, and replace the spaces with underscores so that the whole column is one long name. 

```
gunzip *.correctedReads.fasta.gz
cat *.correctedReads.fasta | sed 's/ /_/g' > correctedReads2.fasta
```

-----------------------------------------------------------------------------------------------------------------------------------------------------------

Get Flye
```
module load mamba/23.1.0-4
conda create -n flye
source activate flye
conda install -c bioconda flye
```

Create the script
```
vi flye_assemble.sh
```

Hit [i] for insertion mode and copy/paste the following:
```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos node
#SBATCH --partition node
#SBATCH --output=out_%assembly.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

source activate flye    

flye --nano-corr ./canu_out/PB127_canu.correctedReads2.fasta -o flye_assembly -t 8 --genome-size 120M
```
Save and exit by pressing [esc], typing ":wq" and then [enter]

Run the script with: 
```
sbatch flye_assemble.sh
```

To see if your job is running type the following command:
```
squeue --me
```

The final assembly is in ./flye_assembly/assembly.fasta

This took approximately 4hrs to assemble a worm genome ~100Mb

</details>



Before starting Verkko or hifiasm, you may want to select for ultra-long ONT reads (50kb and up). You can do this with awk:
```
awk 'BEGIN {RS = "@"; ORS = ""} NR > 1 {getline seq; getline sep; getline qual; if (length(seq) >= MIN_SIZE) print "@"$0, seq, sep, qual}' MIN_SIZE=50000 ontReads.fastq > filteredONT.fastq
```

<details>
<summary><i>Verkko</i></summary>

**Verkko does not do well with little coverage.**

Install Verkko with Conda:
```
conda create -n verkko -c conda-forge -c bioconda -c defaults verkko
source activate verkko
```

Create the script:
```
vi verkko.sh
```

Press[i] for instertion and copy/paste the following:
```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
# Number of nodes
#SBATCH -N 1

# Number of tasks
#SBATCH -n 16

#SBATCH --output=out_verkko.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

source activate verkko 

export VERKKO=/your/path/to/verkko/bin

verkko -d <work-directory> --hifi <hifi-fastq-files> --nano <ont-fastq-files>
```
Save and exit by pressing [esc], typing ":wq" and then [enter]

Run the script with: 
```
sbatch verkko.sh
```

To see if your job is running type the following command:
```
squeue --me
```

The SnakeMake script that Verkko runs on specifies 4CPUs, thus an error will occur if you are trying to run it on HPC without specifying the cores (#SBATCH -n 16)

The assembly is in the output directory and named assembly.fasta

This takes about 2 hours to complete on a worm genome (~100Mb)
 
</details>



<details>
<summary><i>Hifiasm</i></summary>

Install Hifiasm with conda:
```
module load mamba/23.1.0-4
conda create -n hifiasm
source activate hifiasm
conda install -c bioconda hifiasm
```

Create the script:
```
vi hifiasm_assembly.sh
```

Press[i] for instertion and copy/paste the following:
```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH --output=out_hifi.log
#SBATCH --mail-user=username@email.com  #insert your own email
#SBATCH --mail-type=ALL

source activate hifiasm

#pacbio reads only
hifiasm -o sample.asm -t 32 /path/to/hifi_reads.fastq


#pacbio with nanopore reads over 50kb
hifiasm -o sample.asm -t 32 --ul /path/to/filteredONT.fastq /path/to/hifi_reads.fastq

#if you installed with git then you need to include the full path to hifiasm
#make sure to comment out the option you do not want, remember comment out means to put a # at the beginning of the line
```

Save and exit by pressing [esc], typing ":wq" and then [enter]

Run the script with: 
```
sbatch hifiasm_assembly.sh
```

To see if your job is running type the following command:
```
squeue --me
```
Takes about 1 hour on a 100Mb worm genome. 

When complete, your assembly files are .gfa files, which are information about the overlap graphs. To change them into fasta files you can use awk:
```
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa
```

</details>



</details>



<details>
<summary><b>Assembly Polishing</b></summary>

Illumina has a higher base calling accuracy than nanopore (although nanopore may be catching up soon). Therefore we "polish" the assembly by correcting the long read assembly with Illumina short read data. This applies less with HiFi data since it reached a Q20 (99%) quality score. I'm not sure where the community stands on correcting HiFi reads with Illumina. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------

If you assembled with NextDenovo, proceed with NextPolish. If you assembled with Flye, proceed with Pilon.

<details>
<summary><i>NextPolish</i></summary>

https://github.com/Nextomics/NextPolish
```
#create the input file
ls SRR16242711_1.fastq SRR16242711_2.fastq > sgs.fofn
```

Modify the run.cfg file by typing 'vi run.cfg' and hit [i] for insert. Delete the existing code and copy/paste the following:

```
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = /your/path/to/03.ctg_graph/nd.asm.fasta #genome file
genome_size = 120M
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa
```

Create a script for polishing by typing 'vi polish.sh', hit [i] for insert, and copy/paste the following:

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%polish.log
#SBATCH --mail-user=vegge003@fiu.edu   #use your email
#SBATCH --mail-type=ALL

module load nextPolish-1.4.0   #might need to load before running script

nextPolish run.cfg
```

Run the script. The output will be a file with pid***** and a directory named 01_rundir. The directory contains genome.nextpolish.fasta (the polished genome) and genome.nextpolish.fasta.stat (stats about the corrections made). Please rename the file if working with multiple genomes because all will come out with the same name and it could get confusing. 

</details>

<details>
<summary><i>Pilon</i></summary>

https://github.com/broadinstitute/pilon

Create the script:
```
vi pilon.sh
```

Hit [i] for insertion mode and copy/paste the following:
```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos node_name
#SBATCH --partition node_name
#SBATCH --output=out_%pilon.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

module load pilon-1.22-gcc-8.2.0-33xdiwt

FORWARD=[PATH TO FASTQ_1]
REVERSE=[PATH TO FASTQ_2]
LINE_NAME=PB127 ## YOUR LINE
mkdir ./pilon_out/

## ROUND 1 ##
GENOME=[ path to assembled genome]
#index genome 
bwa index ${GENOME}
#align reads
bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
#sam to bam
samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
##Sort and index the BAM 
samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
samtools index ./pilon_out/bwa.sort

##Pilon it 
java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./pilon_out/${LINE_NAME}_pilon1

## ROUND 2 ##
GENOME=./pilon_out/${LINE_NAME}_pilon1.fasta 
#index genome 
bwa index ${GENOME}
#align reads
bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
#sam to bam
samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
##Sort and index the BAM 
samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
samtools index ./pilon_out/bwa.sort


##Pilon it 
java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./pilon_out/${LINE_NAME}_pilon2


## ROUND 3 ##
GENOME=./pilon_out/${LINE_NAME}_pilon2.fasta 
#index genome 
bwa index ${GENOME}
#align reads
bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
#sam to bam
samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
##Sort and index the BAM 
samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
samtools index ./pilon_out/bwa.sort

##Pilon it 
java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./pilon_out/${LINE_NAME}_pilon3



## ROUND 4 ##
GENOME=./pilon_out/${LINE_NAME}_pilon3.fasta 
#index genome 
bwa index ${GENOME}
#align reads
bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
#sam to bam
samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
##Sort and index the BAM 
samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
samtools index ./pilon_out/bwa.sort

##Pilon it 
java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./:pilon_out/${LINE_NAME}_pilon4
```

</details>

</details>

<details>

<summary><b>Quality Check</b></summary>

Assembly quality has various measures. Things like N50, contig number, assembly size, k-mer counting, and gene presence/absence can all be indications of how good an assembly may be. It is a good idea to try multiple assembly methods and use these metrics to compare them. The "best" assembly is usually the most complete and contiguous. QUAST is particularly nice for comparing multiple assemblies at once.

<details>
<summary><i>BUSCO</i></summary>

[BUSCO](https://busco.ezlab.org/) searches for "near-universal" single-copy orthologs from the taxon given. If we are assembling a nematode genome, we use the database nematoda:

```
wget --no-check-certificate https://busco-data.ezlab.org/v5/data/lineages/nematoda_odb10.2020-08-05.tar.gz
tar -xvzf nematoda_odb10.2020-08-05.tar.gz
```

type 'vi busco.sh' to create a script, hit [i], and copy/paste the lines below:

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos node_name
#SBATCH --partition node_name
#SBATCH --output=out_%busco.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL


module load busco/5.4.7 	#might need to load before running script

export AUGUSTUS_CONFIG_PATH="/your/path/to/Augustus"

busco -c 4 -m genome -i /your/path/to/01_rundir/genome.nextpolish.fasta -o busco_PB127 --offline --lineage_dataset /home/data/jfierst/your_username/nematoda_odb10
```
Notice the AUGUSTUS_CONFIG_PATH. We need to copy the augustus directory, give it write permissions, and tell the program the path to that directory. 

```
cp -R /home/data/jfierst/veggers/programs/Augustus/config /your/path/.
cd Augustus
chmod +777 *  #this is a easy but unsafe way to make sure all directories within the directory Augustus each have all permissions. This will take some time.
```
Edit the script to include your path to Augustus and run the script. BUSCO may take multiple hours to run but should not take longer than a day. Your output will be a short_summary*.txt file.

</details>



<details>
<summary><i>QUAST</i></summary>

[QUAST](https://github.com/ablab/quast) calculates assembly metrics like N50, contig #, and assembly size. If given a reference, it can do more, like tell you about misassemblies, however, if the assembly is _de novo_ then you do not have a reference.

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos node_name
#SBATCH --partition node_name
#SBATCH --output=out_%quast.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

module load quast-5.2.0  #may need to load before running script

quast.py -t 4 --eukaryote --plots-format pdf /your/path/to/01_rundir/genome.nextpolish.fasta -o ./PB127_quast/
```

QUAST only takes a minute or two and the output is in the directory PB127_quast. The file report.txt gives you basic genome assembly stats like GC content, N50, # contigs, etc. If you have multiple assemblies your script might look like:

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos node_name
#SBATCH --partition node_name
#SBATCH --output=out_%quast.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

module load quast-5.2.0  #may need to load before running script

quast.py -t 4 --eukaryote --plots-format pdf /your/path/to/assembly1.fasta /your/path/to/assembly2.fasta /your/path/to/assembly3.fasta -o ./species_quast/
```

The html files are files that display the information in a graphical way using icarus viewer. To view these files, you need to download them to your local machine and then click to open them.

OPTION 1: use HPC GUI

If you are using the web browser HPC access, then navigate to the files tab. Type the path to the directory with the file you want to download, select the file, and click download.

OPTION 2: using HPC from computer terminal

```
#logout of the hpc
exit

#connect to hpc using sftp (secure file transfer)
sftp username@hpclogin01.fiu.edu

get /path/to/file.html

#logout of hpc
exit
```

The .html file should now be in your home directory of your local machine.
</details>



<details>
<summary><i>Merqury</i></summary>
	
[Merqury](https://github.com/marbl/merqury) is a k-mer counting tool with a variety of options. It allows you to see which k-mers appear only in the reads, which occur in the assembly and how many times they occur, potential ploidy of the organsim, haplotype phasing, etc.

Merqury is not on the HPC, but it is available through a conda environment. To download:
```
module load mamba/23.1.0-4
conda create -n merqury
source activate merqury
conda install -c conda-forge -c bioconda merqury
```

To test if the installation works:
```
Rscript $MERQURY/plot/plot_spectra_cn.R --help
```
If the help/options dialogue opens, it worked.

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos node_name
#SBATCH --partition node_name
#SBATCH --output=out_%mercury.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL

#many files will be generated, so lets make a seperate directory to work in so that things stay organized:
mkdir species_merqury
cd species_merqury

#count k-mers in read sets
meryl k=31 count reads_1.fastq reads_2.fastq output DF5033.meryl

#compare the read database with the assembly and generate graphs
merqury.sh DF5033.meryl/ assembly.fasta DF5033_merqury
```

Similar to QUAST, you will need to exit ssh and enter sftp to download the .png files to examine:
```
pwd
exit
sftp username@hpclogin01.fiu.edu
cd <copy and paste the output of pwd above>
get *.png
exit
```



</details>



</details>

<details>

<summary><b>Decontamination</b></summary>

 **SIDR**
 
 https://pypi.org/project/SIDR/
 
**It is still in progress and experimental.**

To decontaminate our genomes we use SIDR, a machine learning program genereated by our lab that takes raw fasta/fastq files, runs blast, creates alignments, and generates various statistics about coverage, gc content, and length. This table is then fed into xgboost to predict contaminants. The link above takes you to the 2018 version in python. The code provided in this project uses linux commands and R; however it is in progress to be coded in a different language to make it faster and more efficient. The Fierst lab is also working on updating the machine learning method used.  

The output of SIDR is two fasta files: keptcontigs.fa and contaminantcontigs.fa

</details>

<details>

<summary>Quality Check</summary>

Using the kept contigs.fa, it is good to repeat the QUAST and BUSCO measures for the assembly to make sure there haven't been crazy changes, or if so, then why.

Modify your busco and quast scripts so that instead of /your/path/to/nextpolish.fa, it is changed to /your/path/to/keptcontigs.fa

</details>


<details>

<summary><b>Masking Repeats</b></summary>

**RepeatMasker/Modeler**

https://github.com/Dfam-consortium/RepeatModeler

If the assembly is good to go, we can begin preliminary repeat annotation. This first strp uses repeatModeler/Masker to mask repeat regions of the genome and make gene annotation easier. This is just a first pass, which creates a library of repeats found in the genome. 

Install repeatModeler/Masker with TE-tools container: https://github.com/Dfam-consortium/TETools
RepeatMasker is installed on the HPC, so alternatively you can try module load RepeatMasker-4.1.0

```
#Get the container
curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh

#Activate the container
./dfam-tetools.sh
```


Run commands one by one in container. If I remember correctly, the BuildDatabase command takes the longest (~12 hours on a 100Mb genome, ~20% repeats)
```
#Build the database
BuildDatabase -name [species_name] [genome.fasta]

#Run RepeatModeler for de novo repeat identification and characterization
RepeatModeler -pa 8 -database [species_name]

#Use the queryRepeatDatabase.pl script inside RepeatMasker/util to extract Rhabditida repeats
queryRepeatDatabase.pl -species rhabditida | grep -v "Species:" > Rhabditida.repeatmasker

#Combine the files to create a library of de novo and known repeats
cat RM*/consensi.fa.classified Rhabditida.repeatmasker > [species_name].repeats

#exit the container
exit
```

Mask the repeats from the library you just generated. 
```
#Module load RepeatMasker
module load RepeatMasker-4.1.0

#Mask the genome of known repeats
RepeatMasker -lib [species_name].repeats -pa 8 -xsmall -nolow [keptcontigs.fasta] 
```
-nolow / -l(ow)

With the option -nolow or -l(ow) only interspersed repeats are masked. By default simple tandem repeats and low complexity (polypurine, AT-rich) regions are masked besides the interspersed repeats. For database searches the default setting is recommended, but sometimes, e.g. when using the masked sequence to predict the presence of exons, it may be better to skip the low complexity masking.


-xsmall 

Returns repetitive regions in lowercase (soft masking) instead of replacing with N's (hard masking). Non-repeat regions remain in uppercase.


-pa

Stands for parallel, for multiprocessing, runs 8 sequences at a time. 


The output of RepeatMasker is *.masked

Also, remember that the output of RepeatModeler (custom library) is in RM*/consensi.fa.classified
Or, if you want a broader library, [species].repeats

</details>


<details>

<summary><b>Gene Annotation</b></summary>

<details>
<summary><i>RNA Alignment</i></summary>

You can do RNA alignment with STAR or hisat2. For nematodes, it seems that STAR performs better and leads to more proteins predicted.

<details>
<summary>STAR</summary>

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH --output=out_star
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

module load star-2.7.9a

# Generate genome index
STAR \
    --runThreadN 12 --runMode genomeGenerate --genomeDir [species]_STAR \
    --genomeSAindexNbases 12 --genomeFastaFiles /path/to/keptcontigs.masked

# Map the reads
STAR \
    --runThreadN 12 --genomeDir [species]_STAR --outSAMtype BAM Unsorted --twopassMode Basic \
    --readFilesIn /path/to/rna_1.fastq /path/to/rna_2.fastq
```

The output is an Aligned.out.bam file which will be used in BRAKER, a.k.a the next step.

Takes about 1 and a half hours on a 100Mb genome.

</details>



<details>
<summary>hisat2</summary>

https://daehwankimlab.github.io/hisat2/manual/

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --output=out_%hisat2.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

module load hisat2-2.1.0-gcc-8.2.0-367kzdd
module load bamtools-2.5.2-gcc-11.2.0
module load samtools-1.9-gcc-8.2.0-o53igvd

hisat2-build ./../masked.fasta species

hisat2 -x species -U ./../reads_1.fastq,./../reads_2.fastq -S species.sam

samtools view -bS species.sam > species.bam

```

The output is an Aligned.out.bam file which will be used in BRAKER, a.k.a the next step.

Takes about 1 hour on a 100Mb genome.

</details>



</details>



<details>
<summary><i>Protein Prediction</i></summary>

**BRAKER3**

https://github.com/Gaius-Augustus/BRAKER

BRAKER3 uses genome, RNA, and protein data to predict proteins. It is highly sensitive to the data you give it. 

This is using BRAKER singularity container published 2/28/2024. 

Install the singularity container:

```
module load singularity-3.8.7
singularity build braker3.sif docker://teambraker/braker3:latest
singularity exec braker3.sif braker.pl
```

If installed you should get the braker help menu. You can also test your installation with their test data:

```
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
export BRAKER_SIF=/your/path/to/braker3.sif # may need to modify
bash test1.sh
bash test2.sh
bash test3.sh
```

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos qos_name
#SBATCH --partition partition_name
#SBATCH --output=out_braker3.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

module load singularity-3.8.7
export BRAKER_SIF=/home/data/jfierst/veggers/braker3.sif

singularity exec braker3.sif braker.pl \
 --workingdir=species_braker3 \
 --species=species_name \
 --useexisting \
 --genome=./masked.fasta \
 --prot_seq=./Metazoa.fa \
 --bam=./species_STAR/Aligned.out.bam \   #if you used hisat2: ./species_hisat2/species.bam
 --softmasking

```


**AGAT**

https://github.com/NBISweden/AGAT

AGAT is used to covert the .gtf output file from BRAKER, to a .gff3 file format. It also calculates some basic statistics such as gene count.

Install with conda

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%braker3_agat.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

#combined fasta and protein
agat_convert_sp_gxf2gxf.pl -g /path/to/[sample]_braker.gtf -o [sample]_braker3.gff3
agat_sp_statistics.pl --gff /path/to/[sample]_braker3.gff3 \
        -f /path/to/keptcontigs.masked \
        -o [sample]_AGATstats

```

</details>



<details>
<summary><i>Functional Annotation</i></summary>
**InterProScan**

https://interproscan-docs.readthedocs.io/en/latest/

InterProScan scans the InterPro database to predict protein function and domains.

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH -n 8
#SBATCH --output=out_InterProScan.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL


module load jdk1.8.0_241
module load interproscan-5.55
module load perl-5.34.0-gcc-8.2.0-b5u622f

/home/applications/interproscan/interproscan-5.33-72.0/interproscan.sh -i /path/to/protein.fa -f tsv -dp -goterms -pa
```

</details>



 
</details>




<details>

<summary><b>Repeat Annotation</b></summary>

**EDTA**

https://github.com/oushujun/EDTA

A de novo repeat annotation program that creates a high-quality non-redundant TE consensus library. While it can run given only a genome fasta file, we also provide the input files:

1) CDS, created from BRAKER output
2) BED of known gene positions, created from BRAKER output
3) Curated TE library, created from RepeatModeler output

I did a manual installation of EDTA, however it is conda-able. As of 11/29/2023 it was not available as a module on the hpc.

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_edta.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

perl /path/to/EDTA.pl --genome /path/to/keptcontigs.fasta --cds /path/to/protein.fasta --curatedlib /path/to/[species_name].repeats --exclude /path/to/genes.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10
```


There are many output files generated by EDTA. 

*.mod.EDTA.TElib.fa is the curated consensus library

*.mod.EDTA.TEanno.gff3 contains all TE annotations, including redundant and fragmented sequences

*.mod.EDTA.TEanno.sum is the summary file for whole genome TE annotation

</details>
