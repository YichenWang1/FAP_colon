# The earliest stages of neoplastic transformation in Familial Adenomatous Polyposis

Code and processed data to reproduce the analyses and figures of the manuscript.

Please contact Yichen (yw2@sanger.ac.uk) if you have any questions and enquires.

## Data availability

| Data | Location |
|---|---|
| Whole-genome sequencing CRAM files (GRCh38 ) | EGA, EGAD00001015471 |
| Whole-genome sequencing BAM files (GRCh38), per sample | EGA, EGAD00001015477 |
| Whole-exome sequencing dataset | EGA, EGAD00001016591 |
| Targeted sequencing dataset| EGA, EGAD00001016590 |
| Somatic mutation VCFs | Zenodo, [10.5281/zenodo.23404884](https://doi.org/10.5281/zenodo.23404884) |

Sequencing data are under controlled access through the WTSI CGP Data Access Committee.

## Repository structure

Each analysis folder has R Markdown (`.Rmd`) files with their knitted `.md` and `.html` reports, so the output can be read without re-running the code.

### `data/`

| Path | Contents |
|---|---|
| `fap_database_withexposure.csv` | Master FAP crypt database with sample metadata, APC hits, burdens, signature exposures, telomere content, CNV |
| `ctrl_database_withexposure.csv` | Wildtype crypt database with sample metadata, burdens, signature exposures, telomere content |
| `driver_mutation/somatic_drivers.csv` | Somatic mutations in cancer genes, with tumour-suppressor, hotspot and putative-driver flags |
| `CNV_SV_RT/` | Copy number changes, structural variants and retrotransposition events |
| `pathology/Pathology.csv` | Histology scores and lesion annotations for each crypt |
| `phylogenetic_tree/tree_structure/` | Phylogenetic trees for each individual (`.tree`), with node tables (`.csv`) |
| `phylogenetic_tree/fap_patients.txt`, `ctrl_patients.txt` | FAP and control individuals |
| `phylogenetic_tree/APC_onset_and_trunk_length.csv` | Trunk burdens and estimated onset of somatic APC mutations |
| `somatic_mutation/snp_on_tree/`, `indel_on_tree/` | Substitutions and indels assigned to tree branches |
| `somatic_mutation/*_fap_control_perbranch.txt` | Branch-level 96-channel SBS and 83-channel indel matrices |
| `somatic_mutation/crypt_matrix_*_rounded.txt` | Crypt-level SBS and indel matrices, sensitivity-corrected and rounded |
| `mutational_signature/sigs_on_tree/` | Signature exposures for each tree branch |

### Analysis folders

| Folder | Contents |
|---|---|
| `Signature_extraction/` | HDP signature extraction (SBS and indel) and refitting of signatures to branches|
| `Fig1_Driver_mutation/` | Driver oncoplot (Fig. 1b) and APC lollipop plots for each individual (Fig. 1c) |
| `Fig2_Mutation_burden/` | Burden against age (Fig. 2a), burdens and structural changes by crypt type (Fig. 2b,c), and hypothesis tests for disease effects on burden, structural variation, retrotransposition and telomere content |
| `Fig3_Signature_analysis/` | Signature exposures for each crypt (Fig. 3a), disease effects on signature burdens (Fig. 3b,c), and the corresponding hypothesis tests |
| `Fig4_Phylogeny/` | Signatures on phylogenetic trees, trunk burdens and timing of somatic APC mutations |
