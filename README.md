GuineaPigVMB
============

This git hub includes six main r-code files. Four of the files are to create figures included in the paper “Chlamydia caviae infection alters abundance but not composition of the guinea pig vaginal microbiota”. 

(1)	Bar_plot.R: Generates a stacked bar plot using the relative abundance of guinea pig vaginal microbiota .  
(2)	Logit_evenness.R: Creates a plot of the evenness data. 
(3)	Abundance_3A_3B: 
a.	3A: Generates line graphs of the mean 16S rRNA gene copy number and ompA gene copy number.
b.	3B: Generates line graphs of the top 15 phylotypes and ompA gene copy number for each experimental type. 
(4)	Heatmap.R : Creates a heatmap of the guinea pig and human vaginal microbiota data.

The remaining two scripts are statistical analysis that looks at the relative abundance data and absolute count data.
(1)	gp_zinbin_relAb.R : Using mixed effects Bayesian Markov Chain Monte Carlo (MCMC) models this script looks at the relative abundance of the vaginal microbiota.
(2)	gp_zinbin_abaAB.R : Using mixed effects Bayesian Markov Chain Monte Carlo (MCMC) models this script looks at the absolute abundance of the vaginal microbiota. Multiplying the relative abundance times the 16S rRNA gene count creates the absolute abundance data.
