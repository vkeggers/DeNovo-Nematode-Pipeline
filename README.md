# DeNovo-Nematode-Pipeline
This pipeline is in reference to and building off of https://github.com/BamaComputationalBiology/CRE-UA/blob/main/CRE-Pipeline.md
Everything can be done as listed using the FIU HPC. 

<details>
<summary>
	
## Assembly
</summary>

We have assembly methods for nextDenovo and Flye programs. However, we find nextDenovo to generally be better and more contiguous. 

<details>
<summary>nextDenovo</summary>

If you have your own data already, skip down and start at the line of code that says #create the input file.
	
```
module load sratoolkit-3.0.0
```

Go to NCBI SRA and search _Oscheius_. use the filters at the side to narrow it down to genome and nanopore reads. Find the sra ID for _Oscheius_ sp.G, the number is **SRR16242712**

```
fasterq-dump SRR16242712
#this will take a while and give you no feedback so just believe it will work.
```

If successful you should have a file named SRR16242712.fastq with 18G of data. Type ls -lh to see this.


--------------------------------------------------------------------------------------------------------------------------------------------------------------

https://github.com/Nextomics/NextDenovo
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

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%assemble.log
#SBATCH --mail-user=vegge003@fiu.edu 	#use your own email instead
#SBATCH --mail-type=ALL

module load nextDenovo-2.5.0

nextDenovo run.cfg
```

Save by pressing [esc], type ':wq' and press [enter]


Run the script with: 
```
sbatch < assemble.sh
```

To see if your job is running type the following command:
```
squeue --me
```

There is a common issue some face and you may need to load modules before you run the script. In which case use:
```
module load nextDenovo-2.5.0
sbatch < assemble.sh
```

The final assembly result is at 03.ctg_graph/nd.asm.fasta

Basic statistics for the assembly are at 03.ctg_graph/nd.asm.fasta.stat
</details>

<details>
	<summary>Flye</summary>
	
If you have your own data already, skip down and start conda installation.
	
```
module load sratoolkit-3.0.0
```

Go to NCBI SRA and search _Oscheius_. use the filters at the side to narrow it down to genome and nanopore reads. Find the sra ID for _Oscheius_ sp.G, the number is **SRR16242712**

```
fasterq-dump SRR16242712
#this will take a while and give you no feedback so just believe it will work.
```

If successful you should have a file named SRR16242712.fastq with 18G of data. Type ls -lh to see this.


------------------------------------------------------------------------------------------------------------------------------------------------

https://canu.readthedocs.io/en/latest/quick-start.html#quickstart

The Canu module is available on HPC but I run into a problem with java when trying to use the module. Additionally, Flye is not available, and we don't use these programs enough to request their download. Thus, I've just created conda environments for these. You can try using the anaconda module on HPC (module load anaconda2), but I downloaded my own anaconda a long time ago. You can get the linux version of anaconda here: https://www.anaconda.com/download
Miniconda or Mamba probably work too, I just haven't tried.


Get Canu
```
conda create -n canu
conda activate canu
conda install -c bioconda canu
```

Create script
```
vi canu_correction.sh
```
Hit [i] for insertion mode and copy/paste the following:

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%canu_correct.log
#SBATCH --mail-user=vegge003@fiu.edu   #use your own email
#SBATCH --mail-type=ALL

#conda activate canu

canu -correct -p PB127_canu -d canu_out genomeSize=120M useGrid=false -nanopore-raw ./SRR16242712.fastq
```

Save by pressing [esc], type ':wq' and press [enter]


Run the script with: 
```
sbatch < canu_correction.sh
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

https://github.com/fenderglass/Flye

Get Flye
```
conda create -n flye
conda activate flye
conda install -c bioconda flye
```

Create the script
```
vi flye_assemble.sh
```

Hit [i] for insertion mode and copy/paste the following:
```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%assembly.log
#SBATCH --mail-user=vegge003@fiu.edu   #use your own email
#SBATCH --mail-type=ALL

#conda activate flye    

flye --nano-corr ./canu_out/PB127_canu.correctedReads2.fasta -o flye_assembly -t 8 --genome-size 120M
```
Save and exit by pressing [esc], typing ":wq" and then [enter]

Run the script with: 
```
sbatch < flye_assemble.sh
```

To see if your job is running type the following command:
```
squeue --me
```

The final assembly is in ./flye_assembly/assembly.fasta

This took approximately 4hrs to assemble a worm genome ~100Mb

</details>

We have also assembled using Verkko, which is the more current best practice. However this requires PacBio and Nanopore reads while the above only require Nanopore and Illumina. If you can, then Verkko is highly recommended. 

<details>
	<summary>Verkko</summary>

https://github.com/marbl/verkko

Install Verkko with Conda:
```
conda create -n verkko -c conda-forge -c bioconda -c defaults verkko
conda activate verkko
```

Create the script:
```
vi verkko.sh
```

Press[i] for instertion and copy/paste the following:
```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
# Number of nodes
#SBATCH -N 1

# Number of tasks
#SBATCH -n 16

#SBATCH --output=out_verkko.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL

#conda activate verkko 

export VERKKO=/home/data/jfierst/veggers/anaconda3/envs/verkko/lib/verkko/bin

verkko -d <work-directory> --hifi <hifi-read-files> --nano <ont-read-files>
```
Save and exit by pressing [esc], typing ":wq" and then [enter]

Run the script with: 
```
sbatch < verkko.sh
```

To see if your job is running type the following command:
```
squeue --me
```

The SnakeMake script that Verkko runs on specifies 4CPUs, thus an error will occur if you are trying to run it on HPC without specifying the cores (#SBATCH -n 16)

The assembly is in the output directory and named assembly.fasta

This takes about 2 hours to complete on a worm genome (~100Mb)
 
</details>

</details>


<details>
<summary>
	
## Assembly Polishing
</summary>

Illumina has a higher base calling accuracy than nanopore (although nanopore may be catching up soon). Therefore we "polish" the assembly by correcting the long read assembly with Illumina short read data. 

Find the Illumina data associated with _Oscheius_ sp.G on NCBI SRA. You should get the asseccion number: SRR16242711. If you already have your own data then skip down and start at the line of code that says #create the input file.

```
fasterq-dump SRR16242711
```
If successful you should have a file named SRR16242711_1.fastq and SRR16242711_2.fastq both with 5.4G of data. Type ls -lh to see this.

---------------------------------------------------------------------------------------------------------------------------------------------------------------

If you assembled with NextDenovo, proceed with NextPolish. If you assembled with Flye, proceed with Pilon.

<details>
<summary>NextPolish</summary>

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
<summary>Pilon</summary>

https://github.com/broadinstitute/pilon

Create the script:
```
vi pilon.sh
```

Hit [i] for insertion mode and copy/paste the following:
```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%pilon.log
#SBATCH --mail-user=vegge003@fiu.edu   #use your own email
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

<summary>
	
## Quality Check
</summary>

We use QUAST and BUSCO to check the quality of our genome assemblies. There are two ways of doing this: module load from the hpc, or creating a conda environment. 

**BUSCO**

We need to download the nematode dataset so that we can run busco in offline mode. 

```
wget --no-check-certificate https://busco-data.ezlab.org/v5/data/lineages/nematoda_odb10.2020-08-05.tar.gz
tar -xvzf nematoda_odb10.2020-08-05.tar.gz
```

type 'vi busco.sh' to create a script, hit [i], and copy/paste the lines below:

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%busco.log
#SBATCH --mail-user=vegge003@fiu.edu  #use your email
#SBATCH --mail-type=ALL


module load quast-5.2.0 	#might need to load before running script

export AUGUSTUS_CONFIG_PATH="/your/path/to/Augustus"

busco -c 4 -m genome -i /your/path/to/01_rundir/genome.nextpolish.fasta -o busco_PB127 --offline --lineage_dataset /home/data/jfierst/your_username/nematoda_odb10
```
Notice the AUGUSTUS_CONFIG_PATH. We need to copy the augustus directory, give it write permissions, and tell the program the path to that directory. 

```
cp -R /home/data/jfierst/veggers/programs/Augustus/ /your/path/.
cd Augustus
chmod +777 *  #this is a easy but unsafe way to make sure all directories within the directory Augustus each have all permissions. This will take some time.
```
Edit the script to include your path to Augustus and run the script. BUSCO may take multiple hours to run but should not take longer than a day. Your output will be a short_summary*.txt file.

**QUAST**


```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%quast.log
#SBATCH --mail-user=vegge003@fiu.edu    #use your own email
#SBATCH --mail-type=ALL

module load quast-5.2.0  #may need to load before running script

quast.py -t 4 --eukaryote --plots-format pdf /your/path/to/01_rundir/genome.nextpolish.fasta -o ./PB127_quast/
```

QUAST only takes a minute or two and the output is in the directory PB127_quast. The file report.txt gives you basic genome assembly stats like GC content, N50, # contigs, etc. The html files are files that display the information in a graphical way using icarus viewer.

</details>

<details>

 <summary>
	 
## Decontamination
</summary>
 
</details>


<details>

 <summary>
	 
## Quality Check
</summary>
</details>
