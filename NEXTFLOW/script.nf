//specify site for Fasta download as separated prefix
download_prefix = "ftp://ftp.sra.ebi.ac.uk/"

//process for getting gene anotations
process DownloadGFF {

        publishDir "/home/ubuntu/nextflow"
        executor = "local"

        input:

        output:
        val true

        script:
        """
        wget -O Homo_sapiens.GRCh38.101.chr.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
        """
}

//process for downloading reference chromosomes
process DownloadRef {

        publishDir "/home/ubuntu/nextflow"
        executor = "local"

        input:
        tuple val(name)

        output:
        path "ref.fa"

        script:
        """
        wget -O chromosome_${name}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${name}.fa.gz
        gunzip -c chromosome_${name}.fa.gz > ref.fa
        rm chromosome_${name}.fa.gz
        """
}

// Process for downloading the patient's genomes
process DownloadFastq {

	publishDir "/home/ubuntu/nextflow"
        executor = "local"

        output:
        path "*.fastq"

        input:
        tuple val(name), val(path)

        script:
        """
        wget -O ${path[0].Name} ${download_prefix}${path[0]}
        wget -O ${path[1].Name} ${download_prefix}${path[1]}
	gunzip ${path[0].Name}
	gunzip ${path[1].Name}
        """
}

process Mapping {
	container = "delaugustin/rna-star:2.7.10a"
	
	input:
	file file_ref
	tuple file(fastq_file1), file(fastq_file2) 
        
	output:
	val true
    
	script:    
	"""
	STAR --runThreadN 24 --runMode genomeGenerate --genomeDir /genome/ --genomeFastaFiles ${file_ref}
	
	STAR --outSAMstrandField intronMotif \
		--outFilterMismatchNmax 4 \
		--outFilterMultimapNmax 10 \
		--genomeDir /genome/ \
		--readFilesIn ${fastq_file1} ${fastq_file2} \
		--runThreadN 24 \
		--outSAMunmapped None \
		--outSAMtype BAM SortedByCoordinate \
		--outStd BAM_SortedByCoordinate \
		--genomeLoad NoSharedMemory \
		--limitBAMsortRAM 80000000000 \
		> ${fastq_file1.SimpleName}.bam
	"""
}

// Create the process for counting
process Counting {

	container 'delaugustin/subread'

	input:
	val ready

	output:
	val true

	"""
	featureCounts -T 8 -t gene -g gene_id -s 0 -a input.gtf -o output.counts input.bam
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
        file_ref = DownloadRef(Channel.from(1)) //1 chromosome for the moment
        // Run CreatingIndex process with DownloadRef's output as input
	Mapping(file_ref, fastq_files)
}
