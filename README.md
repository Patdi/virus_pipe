# virus_discovery
Automated script for virus discovery
Supported inputs: fasta, fastq or fastq.gz. 

Run on command line as: virus_pipe_accurate.sh your_file.fastq genome_to_map_to

Scripts automates Trimmomatic-0.36, fastqc, Trinity, bwa, samtools, fasta_formatter, cap3, blastn, blastx, pulling out virus reads and adding taxonomic information from ICTV and NCBI. Files are gzipped along the way to save physical memory space.

Additional scripts:
build_dictionary.sh: will create a genomic index that can be added to virus_pipe_accurate.sh to automate additional geneomes.
crawler.sh: provides an example of moving within the output file structure to pull out, delete or change files. 
