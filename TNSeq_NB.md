# TNSeq Barcode-Mapping Pipeline Notebook 
Kenneth Lai 

### 08/19/25 

*Create Scripts: RUN FASTQC = fastqc_pre.sh*
* Ran ```scp -r ktalapas:/projects/sperolab/kenlai/Pipeline/TNSeq_Barcode_Map/results/fastqc .``` to download & examine pre-cleaned fastqc results manually 

*Create Scripts: RUN filter_reads.Pl = filter_reads.py*
* This script will... 
	1. Organize reads into a MATCHED vs. UNMATCHED FASTQ file. 
		MATCHED CRITERIA = 
			a. Transposon sequence was succesfully detected (within hamming distance parameters - only certain # of permitted mismatches to still be considered 'matched') in record. [this fxn outputs the transposon cutoff pos. too, where the detected tranposons seq. ended]
			b. QS threshold was also met (within hamming distance parameters)
	2. TRIM TRANPOSON SEQ. (of matched reads), keeping only the genomic DNA 
		- this will, next, undergo another FASTQC run --> alignment 

* today I wrote functions for 1a + 1b. These functions take in a sequence and other parameters, then output information used to determine if the read will be considered matched or unmatched. If matched, the tranposon sequence function will also give the 0-based postiion where the tranposon-seq-removal needs to happen. 
	* We need to make this script also output things like.... 
		- number of matched vs. unmatched reads based on phred score & mismatch criteria 


### 08/18/25 

* I'm thinking we can make a parameters on the .yml of whether or not fastq files have been cat-ed or not yet. this way, if yes, skip, if not, then can combine all into one. 
* Coded the fastq cat into setup2_organize.sh 

TO DO: 
* create scripts: 
	- RUN FASTQC (pre & post filter).  # really need to think abt creating a hamming distance filter that allows leniancy, look at 08/14/25 Linecount Analysis Notes for reason of concern. 
	- filter_reads.Pl --> PYTHON 

### 08/14/25 

**Uzip & cat FASTQ files**
My primary concern right now, is based off of the input file you work with, will we have to demultiplex at this step/downstream or if this is done upstream already? 

- setup2_organize.sh
	- my test data files look like 'H16_S2_L001_R1_001.fastq.gz' in which case fastq files have already been concatenated per samples into this one file. 
	- in Melanie's NONBARCODED test data, she has FOLDERS FOR EACH SAMPLE. In one of the sample folders, she has files named "21567_TAGCTTAT_L001_R1_001.fastq.gz" "21567_TAGCTTAT_L001_R1_002.fastq.gz" "21567_TAGCTTAT_L002_R1_001.fastq.gz". These are PER SAMPLE FASTQs *split across lanes & file chunks* . 
		- In which case.... you want to implement some command like....
		```
		# Merge all lanes/chunks, PRESERVING CORRECT NUMERIC ORDER (_001,_002, etc.)
			# sort -V ensures that files are processed in numeric order (lane order = L001, L002). 
			# xargs avoids errors such as 'argument list too long' since these files are so large, actually allows cl to apply a command action to files previously islolated in the pipe 
		ls 21567_TAGCTTAT_L*_R1_*.fastq.gz | sort -V | xargs zcat | gzip -c > 21567_TAGCTTAT_R1.fastq.gz
		```
		- File structure of these test data: SampleID | Sample barcode | Lane # | Fwd vs. Rvs read | 'File Chunk #' 
			- File Chunk # = there are multiple fastqs per lane, so to preserve the correct order of reads, they are each labeld with a final # to keep track. 

**Linecount Analysis**
I'm trying to think ahead (filter_reads.pl) where I need to recognize the transposon-flanking regions in each read. I did some exploratory analysis of the test file H16_S2_L001_R1_001.fastq.gz. Here are my findings: 

- These findings use model_pKMW7 to recognize the flanking sequences: 
```
nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG
TCGACGGCTTGGTTTCATCAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGC
```
- Linecount Analysis 
```
# Total # of records in fastq file
zcat H16_S2_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}'
419276

# Linecount for reads that have at least one of the flanking delim. sequences (D1 or D2/TBS)
zcat H16_S2_L001_R1_001.fastq.gz | grep -c -e "CGCCCTGCAGGGATGTCCACGAGGTCTC" -e "CGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG"
79004
--> This outputs to ~ 18.84% of reads in fastq have at least one of the flanking delim. sequences (exactly, w/ out hamming distance thresholds). This % may be a filter you should set? I think making the filter more linient may be necessary since this % is so small. 
```

TO DO: 
* create scripts: 
	- unzip & cat FASTQ files 
	- RUN FASTQC (pre & post filter)
	- filter_reads.Pl --> PYTHON 




### 08/13/25 
- I am troubleshooting my setup1_conda.sh script & environment.yml file to 'find' a combindation of packages that are compatible to each other. Currently my issues is that certain packages depend on newer vs. older versions of python. 
	- MultiQC >= Python 3.7 
	- All other packages required ~ Python 3.6 to be compatible, but this doesn't meet MultiQC's requirements. I am conducting an experiment with Python 3.7+ to see if conda can find a combination that works. In this version, environment.yml looked like this:
	```
	name: TNSeq_env2
	channels:
  - bioconda
  - conda-forge
  - defaults
	dependencies:
  - python>=3.7
  - fastqc
  - multiqc
  - bowtie2
  - samtools
  - yq
	```
	- ! MODULARIZING CONDA ENVs! --> After experimenting with multiple iterations of python versions ( 3.6-3.12 ), I decided to just create a separate conda environment ONLY for multiqc. This way we can have all the packages necessary while avoiding dependncy conflicts. We will need to activate/deactivate between the TNSeq_env and MultiQC_env during the pipeline at some point then. 
	- Turns out conda couldn't solve the dependency issue, neede to use pip to install multiqc into the Mutliqc_env conda env. 

3PM I finally got setup1_conda.sh to work when run manually. Moving onto setup1_organize.sh 

6PM I fixed the yq command in setup1_organize.sh and the input data moved succesfully. master_pipeline.sh ran successfully with no problems as well. 

TO DO: 
* create scripts: 
	- unzip & cat FASTQ files 
	- RUN FASTQC (pre & post filter)
	- filter_reads.Pl --> PYTHON 

	


