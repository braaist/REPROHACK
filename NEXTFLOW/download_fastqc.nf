//specify site for Fasta download as separated prefix
download_prefix="ftp://ftp.sra.ebi.ac.uk/"

//process for downloading reference chromosomes
process DownloadGFF {
	executor = "local"
	publishDir '/home/ubuntu/nextflow/'
	
	input:
	
	output:
	file "*.gtf.gz"
	
	script:
	"""
	wget -O Homo_sapiens.GRCh38.101.chr.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	"""
}

process DownloadRef {

        executor = "local"
        publishDir '/home/ubuntu/nextflow/'
        input:
        tuple val(name)

	output:
	file "ref.fa"
	
	script:
        """
	wget -O chromosome_${name}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${name}.fa.gz
	gunzip -c chromosome_${name}.fa.gz > ref.fa
	rm chromosome_${name}.fa.gz
        """
}

process DownloadFasta {

        executor = "local"
        publishDir '/home/ubuntu/nextflow/'
        output:
	file "*fastq.gz"
	input:
        tuple val(name), val(path)

        script:
        """
        wget -O ${path[0].Name} ${download_prefix}${path[0]}
        wget -O ${path[1].Name} ${download_prefix}${path[1]}
        """
}

process FastQC {

	container "fastqc"
        publishDir '/home/ubuntu/nextflow/'

        input:
	val fastq_file

        output:
        val true
        
	script:
        """
	echo "{fastq_file[0].simpleName}.zip"
        fastqc ${fastq_file[0]}
	fastqc ${fastq_file[1]}
        """

}


workflow {
	//run DownloadGFF
        DownloadGFF()
	//run DownloadRef for 22 chromosomes
        DownloadRef(Channel.from(1..22))
	// run DownloadFasta for SRA
        fasta_files = DownloadFasta(Channel.fromSRA("SRA062359"))
	FastQC(fasta_files)
}
