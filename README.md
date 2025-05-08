# The earliest stages of neoplastic transformation in Familial Adenomatous Polyposis

Scripts used to reproduce the analyses of the manuscript.

Please contact Yichen (yw2@sanger.ac.uk) if you have any questions and enquires.

## Dataset
Raw sequenicng data and aligned BAM (mapped to hg38) are deposited in the European Genome-Phenome Archive (EGA) with accession code EGAD00001015471 and EGAD00001015477. 

driver_mutation: somatic mutations in driver genes for Fig1.

phylogenetic_tree: structures of individual phylogenetic trees.

somatic_mutation: mutation matrices and SNP/indel placed on individual phylogenetic trees.

mutational_signature: exposure of mutational sigantures on individual phylogenetic tree branches, used to generate Fig4.

fap_database_withexposure.csv and ctrl_database_withexposure.csv: metadata, mutations and signatures for each crypts in FAP patients and control, used to generate Fig2&3.

## Fig1_Driver_mutation
Scripts for generating Fig1b and Fig1c.

## Fig2_Mutation_burden
Scripts for generating Fig2 and hypothesis testing for disease effect on SV, CNV and retrotransposition.

## Signature_extraction
Scripts for HDP signature extraction and refitting.

## Fig3_Signature_analysis
Scripts for generating Fig3 a,b,c and hypothesis testing for disease effect on SBS/indel (total burdens and signatures).

## Fig4_Phylogeny
Scripts for visualising signatures on tree and estumating the onset of the somatic APC mutations.