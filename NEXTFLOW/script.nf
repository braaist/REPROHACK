//specify site for Fasta download as separated prefix
download_prefix = "ftp://ftp.sra.ebi.ac.uk/"
ref_path = "/genome/"


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
        path "chromosome_*.fa"

        script:
        """
        wget -O chromosome_${name}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${name}.fa.gz
        gunzip chromosome_${name}.fa.gz
        """
}

// Process for downloading the patient's genomes
process DownloadFastq {

	publishDir "/home/ubuntu/nextflow"
        executor = "local"

        output:
        path "*.fastq"

        input:
        tuple val(name), val(pathway)

        script:
        """
        wget -O ${pathway[0].Name} ${download_prefix}${pathway[0]}
        wget -O ${pathway[1].Name} ${download_prefix}${pathway[1]}
	gunzip ${pathway[0].Name}
	gunzip ${pathway[1].Name}
        """
}

process CreatingIndex {
        container = "delaugustin/rna-star:2.7.10a"
	publishDir = "/home/ubuntu/nextflow/"

	input:
        file file_ref
        
        output:
        file "genome/"

        script:
        """
	STAR --runThreadN 14 --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles ${file_ref}
	"""

}

process Mapping {
	container = "delaugustin/rna-star:2.7.10a"
	publishDir = "/home/ubuntu/nextflow/"

	input:
	val index
	file file_ref
	tuple file(fastq_file1), file(fastq_file2) 
        
	output:
	val true
    
	script:    
	"""
	STAR --outSAMstrandField intronMotif \
		--outFilterMismatchNmax 4 \
		--outFilterMultimapNmax 10 \
		--genomeDir /${ref_path}/ \
		--readFilesIn ${fastq_file1} ${fastq_file2} \
		--runThreadN 14 \
		--outSAMunmapped None \
		--outSAMtype BAM SortedByCoordinate \
		--outStd BAM_SortedByCoordinate \
		--genomeLoad NoSharedMemory \
		--limitBAMsortRAM 80000000000 \
		> ${fastq_file1.SimpleName}.bam
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
        file_ref = DownloadRef(Channel.from(1..22)).collectFile(name: 'ref.fa')
        // Run CreatingIndex process with DownloadRef's output as input
	CreatingIndex(file_ref).view()
	//index = CreatingIndex(file_ref)
	//Mapping(index, file_ref, fastq_files)
}
