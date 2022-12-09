#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'
## groupe 1 Projet hackaton
## create working dir

echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}|                  Creating working directory                                   |${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
mkdir groupe1_hackaton
cd groupe1_hackaton/
echo -e "${GREEN}successfully done.${NC}\n"

echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}|                  Step 1/ : SOFTWARE INSTALLATION                              |${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
##check that tmux is installed

if ! command -v tmux &> /dev/null
then
    echo "tmux installation"
    sudo apt install tmux
    exit
fi


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
echo -e "${GREEN}successfully done.${NC}\n"



echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}|                  Step 2/ : Downloading script                                 |${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}---------------------------------------------------------------------------------${NC}"

##download script.nf, stat.R and nextflow.config
if test -f script.nf && test -f nextflow.config && test -f stat.R; then
    echo "script.nf and nextflow.config exists."
else
    wget https://raw.githubusercontent.com/braaist/REPROHACK/main/NEXTFLOW/script.nf
    wget https://raw.githubusercontent.com/braaist/REPROHACK/main/NEXTFLOW/nextflow.config
    wget https://raw.githubusercontent.com/braaist/REPROHACK/main/NEXTFLOW/stat.R

fi

echo -e "${GREEN}successfully done.${NC}\n"


echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}|                  Step 3/ : Pulling docker images                              |${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}---------------------------------------------------------------------------------${NC}"

##installing the images

docker pull delaugustin/rna-star:2.7.10a
docker pull delaugustin/subread:v2.0.3
docker pull delaugustin/fastqc:v0.11.9
docker pull delaugustin/r-desqeq2:v4.2
docker pull delaugustin/sra-toolkit:2.11.3
docker pull delaugustin/samtools:v1.16.1
echo -e "${GREEN}successfully done.${NC}\n"



echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}|                  Step 4/ : Running the pipeline                               |${NC}"
echo -e "${RED}}|                                                                               |${NC}"
echo -e "${RED}}---------------------------------------------------------------------------------${NC}"
## running command 
tmux new -s groupe1_hackaton
nextflow script.nf -resume
echo -e "${GREEN}successfully done.${NC}\n"
