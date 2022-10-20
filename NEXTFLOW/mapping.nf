
//Create the mapping for the RNA-seq data

process Mapping {
   

    container 'delaugustin/rna-star'
      

	input:
        

    output:
        
    
   
    """
    gunzip *.gz
    
    STAR  --runThreadN 14 \
    	  --outFilterMultimapNmax 10 \
    	  --genomeDir ${PWD} \
    	  --readFilesIn "*.fastq"  \
    	  --outSAMtype BAM SortedByCoordinate
    """
    
}


workflow{
    // run Mapping 	
    Mapping() 
}






























