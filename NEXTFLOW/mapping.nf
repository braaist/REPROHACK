
//Create the mapping for the RNA-seq data

process Mapping {
   
    container 'delaugustin/rna-star'
      

    input:
        

    output:
           
    """
    STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir SAindex 
     --readFilesCommand zcat --readFilesIn ${fastq_files[1][0]} ${fastq_files[1][1]} --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate
     --genomeLoad NoSharedMemory --limitBAMsortRAM 50 Aligned.sortedByCoord.out.bam ${fastq_files[0]}.bam
    """
    
    
    
}

workflow{
    // run Mapping 	
    Mapping() 
}
