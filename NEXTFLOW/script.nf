//specify site for Fasta download as separated prefix, NCBI key for requests and SRR list
download_prefix = "ftp://ftp.sra.ebi.ac.uk/"
params.ncbi_api_key = '5aba0d52e8608f675c9fa96c9bd0a5d7ca09'
params.ID_list = ["SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589"]

//process for getting gene anotations
process DownloadGFF {
        executor = "local"

        output:
        path "Homo_sapiens.GRCh38.101.chr.gtf"

        script:
        """
        wget -O Homo_sapiens.GRCh38.101.chr.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
        gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz
        """
}

//process for downloading reference chromosomes
process DownloadRef {
        executor = "local"

        input:
        val(name)

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
process fastqc {
    container = " delaugustin/fastqc:v0.11.9"
    
    input:
    file fastq_files
    
    output :
    path "*html" 
    
    script :
    """
    #!/usr/bin/env bash
 
    fastqc *.fastq 
    echo "Quality control done "
    
    """
}
process CreatingIndex {
        container = "delaugustin/rna-star:2.7.10a"
	cpus 14

	input:
        file file_ref
        
        output:
        path "genome/"

        script:
        """
	mkdir genome
	chmod +x genome
	STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir genome --genomeFastaFiles ${file_ref}
	"""

}

process Mapping {
	container = "delaugustin/rna-star:2.7.10a"
	cpus 14

	input:
	path index
	tuple file(fastq_file1), file(fastq_file2) 
        
	output:
	path "*bam"
    
	script:    
	"""
	STAR --outSAMstrandField intronMotif \
		--outFilterMismatchNmax 4 \
		--outFilterMultimapNmax 10 \
		--genomeDir ${index} \
		--readFilesIn ${fastq_file1} ${fastq_file2} \
		--runThreadN $task.cpus \
		--outSAMunmapped None \
		--outSAMtype BAM SortedByCoordinate \
		--outStd BAM_SortedByCoordinate \
		--genomeLoad NoSharedMemory \
		--limitBAMsortRAM 80000000000 \
		> ${fastq_file1.SimpleName}.bam
	"""
}

process Bam_files_indexation {
	container = "delaugustin/samtools:v1.16.1"
        
	input:
	path bam_files

        output:
        path "*.bam"
    
	script:    
	"""
        samtools index ${bam_files}
	"""
}

process Counting {
	container = "delaugustin/subread:2.0.3"
	cpus 8
        
	input:
	path anotations
        path bam_index_output

        output:
        path "count_tab.tx"
    
	script:    
	"""
        featureCounts -T $task.cpus -t gene -g gene_id -s 0 -a ${anotations} -o count_tab.tx ${bam_index_output}
	"""
}

process Stat_analysis {
	container = "delaugustin/r_with_desqeq2:4.2.1"
        
	input:
	path count_tab

        output:


	script:    
	"""
	"""
}

workflow {
	//run DownloadGFF
        anotations = DownloadGFF()

        fastq_files = DownloadFastq(Channel.fromSRA(params.ID_list, apiKey : params.ncbi_api_key))

        // Getting the human reference chromosome by chromosome ang gather the sequence in a unique file
	file_ref = DownloadRef(Channel.from(1..22)).collectFile(name: 'ref.fa')
        // Quality control 
        html = fastqc(fastq_files)

        // Getting the indices for human genome
	index = CreatingIndex(file_ref)

	// Alignment of the paient genes with on the reference génome
	bam_files = Mapping(index, fastq_files)

        // 
        bam_index_output = Bam_files_indexation(bam_files)

        // Getting the count table
        count_tab = Counting(anotations, bam_index_output)
        
        // Make the statiscal analysis with the results
        //results = Stat_analysis(count_tab)
}
