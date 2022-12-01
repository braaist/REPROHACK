#install ubuntu
FROM ubuntu:22.04

#install packages and specify versions for reproducibility
RUN apt-get update
RUN apt-get install wget unzip default-jre libfindbin-libs-perl -y
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && unzip fastqc_v0.11.9.zip -d /home/fastqc && chmod 777 /home/fastqc/FastQC/fastqc
RUN export PATH=$PATH:/home/fastqc/FastQC/
ENV PATH="${PATH}:/home/fastqc/FastQC/"


