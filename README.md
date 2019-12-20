# FRAGTE
A completeness-independent method for pre-selection of closely related genomes for species delineation in prokaryotes

# Please Cite
If you use FRAGTE in your publication, please cite: 

# Description
FRAGTE and TETRA programs are in Bin directory; all tested data are in Data directory; other scripts are in Scripts directory.

# Version
1.0

# Developer
Yizhuang Zhou (zhouyizhuang3@163.com)

# Affiliation
Guilin Medical University

# Support platform
FRAGTE was developed and maintained on Linux platform.

# Prerequisite
1. Perl5 with CPAN, FindBin, File::Basename,POSIX and Time::Local modules
(./autobuild_auxiliary will use CPAN to check and install the other four Perl modules)

# Installation
FRAGTE was developed in Perl language. Therefore, it just needs to install the obove four Perl modules. If these Perl modules have been installed, you just need to run FRAGTE. Otherwise, please type:
./autobuild_auxiliary

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
	
# Performance assessment
## Simulated genome
1. FRAGTE Performance
	
	perl Scripts/Simulated_FRAGTE_Performance.pl [Data/Simulated_IntraSpecies.xls][indir][output|Sensitivity][output|NumStat][specificity]
	
	Note: the output file of determining phase should be named as: Ref_P*_Query_P#_Pairs_byFRAGTE.xls (*, % of reference genomes; #, % of query genome). For example, for references extrating 10% of genomes and queries extrating 20% of genomes, the output file of determining phase is Ref_P10_Query_P20_Pairs_byFRAGTE.xls. All Pairs by FRAGTE are within [indir]. 
	
2. TETRA Performance
	
	perl Scripts/Simulated_TETRA_Performance.pl [Data/Simulated_IntraSpecies.xls][TETRA.list][FRAGTE_Sensitivity.xls][output|Sensitivity][output|numstat][output|specificity]
	
	Note: the output file of TETRA should be named as: Ref_P*_Query_P#_TETRA.xls (*, % of reference genomes; #, % of query genome). For example, for references extrating 10% of genomes and queries extrating 20% of genomes, the output file of determining phase is Ref_P10_Query_P20_TETRA.xls.
	
	[TETRA.list] contaning the files for all output files for TETRA; one line for one file.
	
	[FRAGTE_Sensitivity.xls] the output file for FRAGTE sensitivity

## Real genome
1. FRAGTE Performance
	
	perl Scripts/RealAssembly_FRAGTE_Performance.pl [Pairs byFRAGTE][Query GenomeInfo][Ref GenomeInfo][output]
	
	[TETRA.list] containing TETRA output file; one line for one file
	
	[Pairs byFRAGTE] The output of FRAGTE
	
	[Query GenomeInfo] The file containing 7 fields for queries
	
	[Ref GenomeInfo] The file containing 7 fields for references
	
2. TETRA Performance	
	
	perl Scripts/RealAssembly_TETRA_Performance.pl [TETRA][Query_GenomeInfo.xls][Ref_GenomeInfo.xls][Num of intraspecies pairs sieved by FRAGTE][output]
	
	[TETRA] the output file for TETRA
	
	[Query GenomeInfo] The file containing 7 fields for queries
	
	[Ref GenomeInfo] The file containing 7 fields for references
	
	[Num] Num of intraspecies pairs sieved by FRAGTE
	
## MAGs
1. FRAGTE Performance
	
	perl Scripts/MAG_FRAGTE_Performance.pl [Data/MAG_IntraSpecies.xls][MAG_Pairs_byFRAGTE.xls][output]

	[MAG_Pairs_byFRAGTE.xls] the output fie of FRAGTE

2. TETRA Performance
	perl Scripts/MAG_TETRA_Performance.pl [Data/MAG_IntraSpecies.xls][FRAGTE_Performance.xls][MAG_Pairs_byTETRA.xls][output]
	
	[FRAGTE_Performance.xls] the output file of FRAGTE performance on MAGs
	
	[MAG_Pairs_byTETRA.xls] the output fie of TETRA
	
# Support
If you need some other scripts and other materials or encounter bugs, problems, or have any suggestions, please contact me via zhouyizhuang3@163.com
