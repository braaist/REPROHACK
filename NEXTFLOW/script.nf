//specify site for Fasta download as separated prefix
download_prefix="ftp://ftp.sra.ebi.ac.uk/"

//process for downloading reference chromosomes
process DownloadGFF {
	executor = "local"
	
	input:
	
	output:
	
	script:
	"""
	wget -O ${PWD}/Homo_sapiens.GRCh38.101.chr.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	"""
}

process DownloadRef {

        executor = "local"

        input:
        tuple val(name)

	output:
        val true

	script:
        """
	wget -O ${PWD}/chromosome_${name}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${name}.fa.gz
	gunzip -c ${PWD}/chromosome_${name}.fa.gz > ${PWD}/ref.fa
	rm ${PWD}/chromosome_${name}.fa.gz
        """
}

process DownloadFasta {

        executor = "local"

        output:

	input:
        tuple val(name), val(path)

        script:
        """
        wget -O ${PWD}/${path[0].Name} ${download_prefix}${path[0]}
        wget -O ${PWD}/${path[1].Name} ${download_prefix}${path[1]}
        """
}

// Process for creating genome index
process CreatingIndex {

        container 'delaugustin/rna-star'

        input:
        val ready

	output:

        """
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${PWD}/ --genomeFastaFiles ref.fa
        """
}

workflow {
	//run DownloadGFF
        DownloadGFF()

	// Pipeline indexation and mapping
        DownloadRef(Channel.from(1))
        CreatingIndex(DownloadRef.out)

}
