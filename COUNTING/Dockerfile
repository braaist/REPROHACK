FROM ubuntu:22.04

# installation de wget et des commandes absentes
RUN apt-get update --fix-missing \
&& apt-get install -y wget \
&& cd /usr/local/ \

#telechargement du fichier
&& wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz \

#Decompression de l'executable
&& tar -xf subread-2.0.3-Linux-x86_64.tar.gz \

#Deplacer le fichier decompression dans bin
&& mv subread-2.0.3-Linux-x86_64 /usr/local/bin/ \

#rendre le fichier executable
&& chmod +x /usr/local/bin/subread-2.0.3-Linux-x86_64 \

#supriession du fichier compressé
&& rm subread-2.0.3-Linux-x86_64.tar.gz

ENTRYPOINT ["/usr/local/bin/bin/subread-2.0.3-Linux-x86_64"]
