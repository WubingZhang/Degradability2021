# Model-based Analysis of Protein Degradability (MAPD) 

## Background

Targeted protein degradation (TPD) has rapidly emerged as a therapeutic modality to eliminate previously undruggable proteins by repurposing the cell’s endogenous protein degradation machinery. However, the susceptibility of proteins for targeting by TPD approaches, termed “degradability”, is largely unknown. Recent systematic studies to map the degradable kinome have shown differences in degradation between kinases with similar drug-target engagement, suggesting yet unknown factors influencing degradability. We therefore developed a machine learning model, MAPD (Model-based Analysis of Protein Degradability), to predict degradability from protein features that encompass post-translational modifications, protein stability, protein expression and protein-protein interactions. MAPD shows accurate performance in predicting kinases that are degradable by TPD compounds (auPRC=0.759) and is likely generalizable to independent non-kinase proteins. We found five features with statistical significance to achieve optimal prediction, with ubiquitination potential being the most predictive. By structural modeling, we found that E2-accessible ubiquitination sites, but not lysine residues in general, are particularly associated with kinase degradability. Finally, we extended MAPD predictions to the entire proteome to find 964 disease-causing proteins, including 278 cancer genes, that may be tractable to TPD drug development.


## Notes
This repo provides source codes for reproducing the analysis in our study, titled "Machine learning modeling of protein-intrinsic features predicts tractability of targeted protein degradation".

If you want to train a MAPD model for predicting degradability, please visit https://liulab-dfci.github.io/MAPD/.

If you want to check the results presented in our study, please visit http://mapd.cistrome.org/

## Citation

Wubing Zhang, 
