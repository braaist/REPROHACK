FROM ubuntu:22.04

#Installation des dépendances
RUN apt-get update --fix-missing
RUN DEBIAN_FRONTEND="noninteractive" TZ="Europe" apt-get install -y tzdata
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev
RUN apt install -y r-base

# On demande à installer les packages BiocManager (ddépendence de DESeq2)
#et DESeq2
RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") '
RUN R -e 'BiocManager::install("DESeq2")'DESeq2")'
