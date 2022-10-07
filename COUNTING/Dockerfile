FROM ubuntu:22.04
RUN apt-get update --fix-missing \
&& apt-get install -y wget \
&& cd /usr/local/ \
&& wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz \ #telechargement du fichier
&& tar -xf subread-2.0.3-Linux-x86_64.tar.gz \ #Decompression de l'executable
&& sudo mv subread-2.0.3-Linux-x86_64 /usr/local/bin/ \ #Deplacer le fichier decompression dans bin
&& chmod +x /usr/local/bin/subread-2.0.3-Linux-x86_64 \ #rendre le fichier executable
&& rm subread-2.0.3-Linux-x86_64.tar.gz #supriession du fichier compressé

ENTRYPOINT ["/usr/local/bin/bin/subread-2.0.3-Linux-x86_64"]
