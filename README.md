# metaamp
The MetaAmp is designed for analyzing universally conserved phylogenetic marker genes (e.g. 16S/18S rRNA) of metagenomic data. It clusters amplicon reads into OTUs and generates diversity indices and taxonomic composition of the input communities . It accepts single-read or paired-end sequencing in fasta or fastq format from various sequencing platforms.

A MetaAmp can be accessed at: http://ebg.ucalgary.ca/metaamp/

A MetaAmp analysis output demo page can be found at: https://xiaoli-dong.github.io/metaamp/

MetaAmp can also be run through command line, here is the example command:

metaamp.pl -map your_mapping_file -an your_analysis_name -seqtype paired -seqformat fastq -s 0.97 -minoverlen 50 -maxdiffs 0  -trunclen 350 -maxee 1 -oligos oligos.txt -pdiffs 0 -email your_email@gmail.com
