# FRAGTE
A completeness-independent method for pre-selection of closely related genomes for species delineation in prokaryotes

# Please Cite
If you use FRAGTE in your publication, please cite: 

# Note
FRAGTE and TETRA programs are in Bin directory; all tested data are in Data directory; other scripts are in Scripts directory.

# Usage
## FRAGTE
1. Fragmenting phase for references/Queries
	
	perl Bin/RepZvalue.pl [GenomeInfo file] [output]
	
	[GenomeInfo file] has 7 fields, including:
	Assembly_assession\tspecies_taxid\torganism_name\tinfraspecific_name\tassembly_level\tGenome size\tfile for Genome

	Note: file for Genome,the file containing sequences in fasta format; one line for one genome.

2. Determining phase
	
	perl Bin/Pairs_byFRAGTE.pl [Ref RepZvalue][Query RepZvalue][output]
	
	[Ref RepZvalue] the output file for reference genome info file from step 1.
	
	[Query RepZvalue] the output file for Query genome info file from step 1.

## TETRA
1. Calculate z-values for references/Queries
	perl Bin/Zvalue.pl [GenomeInfo file] [output]
	
	[GenomeInfo file] has 7 fields, including:
	Assembly_assession\tspecies_taxid\torganism_name\tinfraspecific_name\tassembly_level\tGenome size\tfile for Genome

	Note: file for Genome,the file containing sequences in fasta format; one line for one genome.

2. Calculate PCCDs for queries against references
	perl Bin/Pairs_byTETRA.pl [Ref_Zvalue.xls][Query_Zvalue.xls][output]
	[Ref_Zvalue.xls] Zvalue for references
	[Query_Zvalue.xls] Zvalue for queries
	[output] output file
	
