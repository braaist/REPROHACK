# REPROHACK

The goal is to reproduce parts of the analysis described in paper https://pubmed.ncbi.nlm.nih.gov/23861464/ . They performed RNA-Seq in samples from patients with uveal melanoma, some samples are mutated in SF3B1. We want to analyze this data in order to find differentially expressed genes, i.e. genes that are more (or less) expressed in one condition (SF3B1 mutated samples) compared to another (SF3B1 non mutated samples). We need to design Nextflow workflow that must be reproducible. The work will be done in small group.
1. This workflow contains containers (Dockerfiles);
2. Workflow code (Nextflow)
3. README.md + run.sh with all instructions to reproduce the analysis
