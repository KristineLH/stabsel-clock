# Introduction
This repository contains scripts and supplementary data belonging to the article "Stability selection enhances feature selection and enables accurate prediction of gestational age using only five DNA methylation sites" by Haftorn et al, 2023. It shows how stability selection (Meinshausen and Bühlmann, 2010) can be used to identify CpGs that are stably predictive of gestational age (or any other trait that can be predicted from DNA methylation data), and how to create 'stable epigenetic clocks' for gestational age using those CpGs and GAM regression (Wood, 2017).

## MAIN_ANALYSES
Contains scripts and supplementary data for the main analyses in the article, including stability selection and creating 'stable clocks'.

## FURTHER_ANALYSES
Contains scripts and supplementary data for the downstream analyses of CpGs that were identified as stably predictive.

## STABLE_CLOCKS
Contains the epigenetic gestational age clocks developed in the article and a brief explanation on how to apply them to DNA methylation data.

# References
Haftorn KL, Romanowska J, Lee Y, et al. Stability selection enhances feature selection and enables accurate prediction of gestational age using only five DNA methylation sites. Clin Epigenetics. 2023;15(1):114.

Meinshausen N, Bühlmann P. Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology). 2010;72(4):417-73.

Wood  SN. Generalized Additive Models: An Introduction with R. 2nd ed: Chapman and Hall/CRC; 2017.
