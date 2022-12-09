//specify site for Fasta download as separated prefix, NCBI key for requests and SRR list
download_prefix = "ftp://ftp.sra.ebi.ac.uk/"
params.ncbi_api_key = '5aba0d52e8608f675c9fa96c9bd0a5d7ca09'
params.ID_list = ["SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589"]
params.outdir = "${PWD}/results"
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
    container = "delaugustin/fastqc:v0.11.9"
    publishDir params.outdir
    input:
    file fastq_files
    
    output :
    path "*html" 
    
    script :
    """
    #!/usr/bin/env bash
 
    fastqc *.fastq 
    echo "Quality control done"
    
    """
}
process CreatingIndex {

	executor = "local"
	
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
	publishDir "/home/ubuntu/REPROHACK/"

	input:
	tuple file(fastq_file1), file(fastq_file2), path(index) 
        
	output:
	path "*.bam"
    
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
	--limitBAMsortRAM 80000000000 > ${fastq_file1.simpleName}.bam 
	"""
}

process Counting {
	container = "delaugustin/subread:2.0.3"
	cpus 14
        
	input:
        file(bam_files)
        path(anotations)

    output:
        path "*count_tab.txt"
    
	script:    
	"""
        featureCounts -T $task.cpus -p -t gene -g gene_id -s 0 -a ${anotations} -o out_count_tab.txt ${bam_files}	
        """
}

process Stat_analysis {
	container = "delaugustin/r-desqeq2:v4.2"
        publishDir params.outdir
	input:
        path script_stats
        path count_tab   


	script:    
	"""
   
        #!/usr/bin/env bash
	Rscript ${PWD}/${script_stats} ${count_tab}
    
	"""
}

workflow {
	//run DownloadGFF
        anotations = DownloadGFF()

        fastq_files = DownloadFastq(Channel.from(params.ID_list, apiKey : params.ncbi_api_key))

        // Getting the human reference chromosome by chromosome ang gather the sequence in a unique file
	file_ref = DownloadRef(Channel.from(1..22)).collectFile(name: 'ref.fa')
        // Quality control 
        html = fastqc(fastq_files)

        // Getting the indices for human genome
	index = CreatingIndex(file_ref)

	// Alignment of the paient genes with on the reference génome
	bam_files = Mapping(fastq_files.combine(index))
        bam_files = bam_files.collect()
	bam_files.view()

        // Getting the count table
        count_tab = Counting(bam_files, anotations)
        count_tab.view()
        // Make the statiscal analysis with the results
        stat = Channel.fromPath("stat.R")
        results = Stat_analysis(stat,count_tab)
}
