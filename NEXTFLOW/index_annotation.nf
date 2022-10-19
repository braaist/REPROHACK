// Process for creating genome index
process CreatingIndex {

        container 'delaugustin/rna-star'

        input:

	output:

        """
        STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${PWD}/ --genomeFastaFiles ref.fa
        """
}

// Process for Getting genome annotations
process GettingAnnotations {
	executor = "local"
	
	input:
        
	
	output:
	
	script:
	"""
	wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	"""
}

workflow {
        // run CreatingIndex
        CreatingIndex()
        // run GettingAnnotations
        GettingAnnotations()
}
