# REPROHACK

The goal is to reproduce parts of the analysis described in paper https://pubmed.ncbi.nlm.nih.gov/23861464/ . They performed RNA-Seq in samples from patients with uveal melanoma, some samples are mutated in SF3B1. We want to analyze this data in order to find differentially expressed genes, i.e. genes that are more (or less) expressed in one condition (SF3B1 mutated samples) compared to another (SF3B1 non mutated samples). We need to design Nextflow workflow that must be reproducible. The work will be done in small group.
1. This workflow contains containers (Dockerfiles);
2. Workflow code (Nextflow)
3. README.md + run.sh with all instructions to reproduce the analysis

# Dependencies : 
The pipeline runs on nextflow a domain-specific language created to automate data-analysis pipelines whilst maximising reproducibility. Nextflow enables scientists to focus on their analyses, isolating different parts of the pipeline into processes whose dependencies can be dealt with using containers and virtual environments with technologies such as Docker, Singularity, and Anaconda.

# Hardware requirements :
A machine with at least 32 GB of FREE RAM (to create the index and the mapping on the reference genome) and at least 14 threads. Recommended configuration is 64 GB, by default the mapping process is configured to use 50 GB. Number of threads may be specified in script.nf file. 

# Executing The Workflow :
1 - Clone the repo to your machine or download the run.sh file manually.

git clone https://github.com/braaist/REPROHACK.git
cd REPROHACK 

2 - Run the wokflow (in the tmux or nohup)

bash run.sh

3 – Final results may be found in work/ subdirectory of Stats process.

# Caution 
- A good internet connection is required for the recovery of fastq.

- The workflow will inevitably fail if you attempt to create the genome index and mapping on a machine with less than ~30 GB of available RAM.

- Execution of the workflow on a machine with 14 threads and 64 GB of RAM takes ~6hrs, that's why it's strongly recommended to run run.sh script in tmux or nohup utilities. 


