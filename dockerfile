FROM ubuntu:22.04
RUN apt-get update
RUN apt-get install rna-star=2.7.10a+dfsg-1 -y