# DATA folder

## `3_stabsel_cpgs_anno.csv`

Original data from Kristine.

|  column name  |  description  |  notes   |
|---------------|---------------|----------|
| `cpg`         | CpG name      |          |
| `select_prob` | selection probability  | |
| `6_cpg_clock` | is this CpG in the minimal clock? | |
| `AddressA_ID` -- `UCSC_RefGene_Name2` | description of probe, from Illumina Manifest | |

## `all_genes_regs.rds`

R-data file with all the genes that overlap with CpGs. Data fetched from ensembl
using {biomaRt}.

## `all_regul_regs.rds`

R-data file with all the regulatory regions and which genes they control. Data
fetched from ensembl using {biomaRt}.

## `genes_controlled_min_clock.txt` and `genes_controlled_other_clock.txt`

List of unique genes controlled by the regulatory regions overlapping with
CpGs in the minimal clock and overlapping with CpGs not in the minimal clock.

## `genes_min_clock_pos.rds` and `genes_other_pos.rds`

R-data file with genes _controlled by_ regulatory regions overlapping CpGs in
the minimal clock and overlapping with CpGs not in the minimal clock. Data
fetched from ensembl using {biomaRt}.

## `old_new_regul_regs_ids.csv`

Mapping the ensembl IDs of the regulatory regions overlapping with CpGs from
GRCh37 to GRCh38 (data manually obtained from ensembl).

## `regul_regs_gene_targets_min_clock.csv` and `regul_regs_gene_targets_other.csv`

Genes controlled by each of the regulatory regions overlapping with CpGs in the
minimal clock and overlapping with CpGs not in the minimal clock. Data manually
gathered from GeneHancer (via GeneCards).

## `string_db_mapping_genes_control_min_clock.tsv`

Mapping of gene names when querying STRING-db.

## DATA downloaded from databases:

### [EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/downloads)

### `EWAS_trait_trait_logP.txt`

Trait-to-trait relationships based on shared CpGs.

*Format:* Matrix of associations between the traits (there are 372 traits here).

### `EWAS_Atlas_associations.tsv`

Table with all collected reported associations between a CpG probe and a
phenotype/trait.

*Format:* data frame with following columns:    
`Association_ID, probe_ID, trait, case_description, case_beta,
control_description, control_beta, correlation, p_value, rank_in_study,
effect_size, study_ID, PMID`

## [EWAS Catalog](http://ewascatalog.org/download)

These files are tab-delimited text files that are linked by the column "StudyID".
NB: Missing information is left blank or as NA.

### `ewascatalog-studies.txt.gz`

This meta-data file has the following pieces of information:

- Author - the first author of the publication (surname then initials).
- Consortium - the name of the or cohorts used.
- PMID - the PubMed ID of the publication.
- Date - the date the paper was published (YY-MM-DD).
- Trait - the name of the trait.
- EFO - the corresponding ontology term(s) for the trait.
- Analysis - description of the analysis performed.
- Source - the table where the result can be found in the paper.
- Outcome - the outcome of the analysis.
- Exposure - the exposure of the analysis.
- Covariates - the covariates adjusted for in the analysis.
- Outcome_Units - the units of the outcome.
- Exposure_Units - the units of the exposure.
- Methylation_Array - the array used to measure the methylation.
- Tissue - the tissue in which the methylation was measured.
- Further_Details - any further details on the analysis.
- N - the total number of participants used in the analysis.
- N_Cohorts - the total number of cohorts used in the analysis.
- Age - the age group participants belonged to.
- Sex - sex of individuals used in the analysis.
- Ethnicity - ethnicity of the individuals used in the analysis.
- StudyID - ID that can be used to link the meta-data to the results.

### `ewascatalog-results.txt.gz`

- CpG - the CpG site.
- Location - the location of the CpG based on hg19 coordinates.
- Chr - the chromosome where the CpG is located based on hg19 coordinates.
- Pos - the position where the CpG is located based on hg19 coordinates.
- Gene - the gene where the CpG is mapped to based on the Ilumina manifest.
- Type - the type of genomic location the CpG is in (e.g. CpG island).
- Beta - the effect estimate.
- SE - the standard error of beta.
- P - p-value.
- Details - any additional details on the analysis (e.g. sub-trait).
- StudyID - ID that can be used to link the meta-data to the results.

