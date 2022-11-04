#!/bin/bash
## groupe 1 Projet hackaton


##installing the images

docker pull delaugustin/rna-star:2.7.10a
docker pull delaugustin/subread:2.0.3
docker pull delaugustin/fastqc:0.11.9
docker pull delaugustin/r_with_desqeq2:4.2.1
ocker pull delaugustin/sra-toolkit:2.11.3

## running command 
mkdir nextflow
nextflow REPROHACK/NEXTFLOW/script.nf -resume
