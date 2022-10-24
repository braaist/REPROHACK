//specify site for Fasta download as separated prefix
download_prefix="ftp://ftp.sra.ebi.ac.uk/"

//process for getting gene anotations
process DownloadGFF {
	executor = "local"
	
	input:
	
	output:
        val true
	
	script:
	"""
	wget -O ${PWD}/Homo_sapiens.GRCh38.101.chr.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	"""
}

//process for downloading reference chromosomes
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

// Process for downloading the patient's genomes
process DownloadFastq {

        executor = "local"

        output:
        val true

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
        val true

        """
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${PWD}/ --genomeFastaFiles ${PWD}/ref.fa
        """
}

//Create the mapping for the RNA-seq data
process Mapping {
   
    container 'delaugustin/rna-star'
      
    input:
    val ready
        
    output:
    //val true
           
    """
    gunzip *.gz
    
    STAR  --runThreadN 14 \
    	  --outFilterMultimapNmax 10 \
    	  --genomeDir ${PWD} \
    	  --readFilesIn "*.fastq"  \
    	  --outSAMtype BAM SortedByCoordinate
    """
    
}

workflow {
	
        // ===========Pipeline for getting gene anotation===================
        //run DownloadGFF
        DownloadGFF()
	

        // ===========Pipeline for downloading the patient's genes=================
        fastq_files = DownloadFastq(Channel.fromSRA("SRA062359"))
	

        // ============Pipeline indexation and mapping===========================
        // Run DownloadRef with the channel
        DownloadRef(Channel.from(1)) //1 chromosome for the moment

        // Run CreatingIndex process with DownloadRef's output as input
        CreatingIndex(DownloadRef.out)

        // Run Mapping process with CreatingIndex's output, DownloadFastq's output and DownloadGFF's output as input
        Mapping(CreatingIndex.out,DownloadFastq.out,DownloadGFF.out)
}