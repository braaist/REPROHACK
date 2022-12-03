#!/bin/bash

## groupe 1 Projet hackaton
## create working dir
mkdir REPROHACK
cd REPROHACK/


##check that nextflow installed
if ! command -v nextflow &> /dev/null
then
    echo "Nextflow installation"
    sudo apt install default-jre
    curl -fsSL get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin
    exit
fi

if ! command -v docker &>? /dev/null/
then
    echo "Docker installation"
    sudo apt-get install docker-ce docker-ce-cli containerd.io docker-compose-plugin
    exit
fi

##download script.nf and nextflow.config

if test -f script.nf && test -f nextflow.config; then
    echo "script.nf and nextflow.config exists."
else
    wget https://raw.githubusercontent.com/braaist/REPROHACK/main/NEXTFLOW/script.nf
    wget https://raw.githubusercontent.com/braaist/REPROHACK/main/NEXTFLOW/nextflow.config
fi

##installing the images

docker pull delaugustin/rna-star:2.7.10a
docker pull delaugustin/subread:2.0.3
docker pull delaugustin/fastqc:v0.11.9
docker pull delaugustin/r_with_desqeq2:4.2.1
docker pull delaugustin/sra-toolkit:2.11.3
docker pull delaugustin/samtools:v1.16.1

## running command 
nextflow script.nf -resume
