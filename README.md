# DeNovo-Nematode-Pipeline
Everything can be done as listed using the FIU HPC

<details>
<summary>Assembly</summary>
	
```
module load sratoolkit-3.0.0
```

Go to NCBI SRA and search _Oscheius_. use the filters at the side to narrow it down to genome and nanopore reads. Find the sra ID for _Oscheius_ sp.G, the number is **SRR16242712**

```
fasterq-dump SRR16242712
#this will take a while and give you no feedback so just believe it will work.
```

If successful you should have a file named SRR16242712.fastq with 18G of data. Type ls -lh to see this.

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
#SBATCH --output=out_%assemble .log
#SBATCH --mail-user=vegge003@fiu.edu 	#use your own email instead
#SBATCH --mail-type=ALL

module load nextDenovo-2.5.0

nextDenovo run.cfg
```

Save by pressing [esc], type ':wq' and press [enter]

To see if your job is running type the following command:

```
squeue --me
```

The final assembly result is at 03.ctg_graph/nd.asm.fasta

Basic statistics for the assembly are at 03.ctg_graph/nd.asm.fasta.stat

</details>

<details>
<summary>Assembly Polishing</summary>

Illumina has a higher base calling accuracy than nanopore (although nanopore may be catching up soon). Therefore we "polish" the assembly by correcting the long read assembly with Illumina short read data. 

Find the Illumina data associated with _Oscheius_ sp.G on NCBI SRA. You should get the asseccion number: SRR16242711.

```
fasterq-dump SRR16242711
```
If successful you should have a file named SRR16242711_1.fastq and SRR16242711_2.fastq both with 5.4G of data. Type ls -lh to see this.

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
genome = /home/data/jfierst/veggers/PB127/03.ctg_graph/nd.asm.fasta #genome file
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
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL

module load nextPolish-1.4.0

nextPolish run.cfg
```

The output will be a file with pid***** and a directory named 01_rundir. The directory contains genome.nextpolish.fasta (the polished genome) and genome.nextpolish.fasta.stat (stats about the corrections made). Please rename the file if working with multiple genomes because all will come out with the same name and it could get confusing. 

</details>

<details>

<summary>Quality Check</summary>

We use QUAST and BUSCO to check the quality of our genome assemblies. There are two ways of doing this: module load from the hpc, or creating a conda environment. As of July 18, 2023, the modules were not on the hpc. To check if they are, type 'module avail busco', or 'module avail quast'. To see if you have conda type 'conda --help'. If it the command isn't recognized then you need to get conda. I have downloaded anaconda from source but there is an anaconda module on the hpc and you can use other programs like mamba or miniconda. Take your pick.

**BUSCO**

``
conda create -n busco #-n is telling it to create an environment named busco
conda activate busco #start up the environment. Your username should be replaced by the name of the environment.
conda install -c bioconda busco #install busco in the environment using the bioconda channel
``


</details>

<details>

 <summary>Decontamination</summary>
 
</details>


<details>

 <summary>Quality Check</summary>
</details>
