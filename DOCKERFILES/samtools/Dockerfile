FROM ubuntu:22.04

#installation of wget and missing commands
RUN apt-get update --fix-missing \
&& apt-get install -y wget \
&& cd /usr/local/ \

#installation of required packages
RUN apt-get update 
RUN apt-get install -y gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev bzip2

#Dowload the compressed file, decompression and make samtools executable
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
&& tar -xjvf samtools-1.16.1.tar.bz2 \
&& rm -rf samtools.tar.bz2 \
&& cd samtools-1.16.1 \
&& ./configure \
&& make \
&& make install