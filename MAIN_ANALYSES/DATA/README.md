## DATA folder

Here, you can find the supplementary tables from the manuscript.

### Supp_Data_1_stabsel.csv

This table includes results from the stability selection run of the complete dataset, with columns and description as below:

column name  | description
-------------|------------------
CpG_ID       | CpG ID from Illumina   
selection_probability | Selection probability  
chromosome   | Chromosomal location
genomic_coordinates | Genomic location 
relation_to_CpG_island | Relation to CpG island 
present_on_450k | Whether this CpG is present on the 450K array or not
Gene_Illumina | Gene annotation by Illumina manifest file
in_Haftorn_clock | Whether this CpG is present in the Haftorn clock or not
in_Bohlin_clock | Whether this CpG is present in the Bohlin clock or not  
in_Knight_clock | Whether this CpG is present in the Knight clock or not

### Supp_Data_2_stabsel_training.csv

This table includes results from the stability selection run of the training dataset, with columns and description as below:

column name  | description
-------------|------------------
CpG_ID       | CpG ID from Illumina    
selection_probability | Selection probability    
chromosome | Chromosomal location
genomic_coordinates | Genomic location 
relation_to_CpG_island | Relation to CpG island
present_on_450k | Whether this CpG is present on the 450K array or not
Gene_Illumina | Gene annotation by Illumina manifest file
