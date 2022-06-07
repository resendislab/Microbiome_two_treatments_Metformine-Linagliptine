# Microbiome_two_treatments_Metformine-Linagliptine

Objetive: We have two objectives in this study.
- To determine the effect of linagliptin/metformin compared to metformin on the composition of the gut microbiome in treatment-naive patients with prediabetes and compared to Metformin and linagliptin/metformin.
- To determine the relationship of the composition of the gut microbiome with insulin resistance and insulin secretion in patients with prediabetes treated with linagliptin/metformin compared to metformin and treatment-naive patients.

# Data 

Patients included in this study were part of the project Prevention of diabetes with linagliptin, lifestyle, and metformin (PRELLIM). PRELLIM was a double-blind and randomized clinical trial designed to evaluate the effect of the combination of linagliptin + metformin + lifestyle compared to only metformin + lifestyle on glucose metabolism, insulin resistance, and pancreatic islet function for 12 months (ClinicalTrials.gov number, NCT03004612). The trial conformed to the provisions of the Declaration of Helsinki and with the research and Ethics Committee of the Hospital Regional de Alta Especialidad del Bajío (CI-HRAEB-2017-048 and CEI-22-16). Patients included in this study were 167 eligible randomized participants (65 at basal, 77 at six months follow-up, and 63 at 12 months follow-up). At six months of follow-up, 44 patients belong to the metformin and 33 to the linagliptin/metformin group; at 12 months of follow-up, 28 patients belong to the metformin group and 35 to the linagliptin/metformin group. In addition, I followed up on 24 subjects between 6 and 12 months.

# Methodology

The next figure shows the analysis of the research work.

![image](https://user-images.githubusercontent.com/63600444/172265390-164e5c87-f2d2-451f-ad06-86e61799242f.png)

The sequences obtained in the sequencing were processed using the QIIME2 (Quantitative Insights Into Microbial Ecology) analysis platform. The taxonomic classification for end sequence variants was performed using the SILVA database (version 132). The analysis of the GM was performed using the phyloseq object; through this object, obtained the α and β diversity of the study groups in the different months of follow-up. Subsequently, with the abundance of the different bacterial genus, the genus that classifies for each pharmacological treatment in the different months of follow-up were obtained using the machine learning algorithm (Random forest). Once the bacteria that classified the subjects with metformin and linagliptin/metformin treatment were identified, hierarchical linear regression models were performed to eliminate the effect of confounding variables (obesity, age, and gender). In addition, mixed generalized linear models and finally mediation analysis using structural equation models.

# Analysis scripts and discussion for the Treatment_T2D.

Makes use of helper functions such as QIIME2 and phyloseq. All individual steps are contained in their own sub-folder along with documentation. Input and intermediate data files can be found in the data folder, and some plots can be found in the figures folder.

Getting set up

You will need QIIME2 and R (>=4.1.2). 

Install QIIME2 within a conda environment

```
wget https://data.qiime2.org/distro/core/qiime2-2021.4-py38-linux-conda.yml
conda env create -n qiime2-2021.4 --qiime2-2021.4-py38-linux-conda.yml file
```

Optional cleaning 

```
rm qiime2-2021.4-py38-linux-conda.yml
```

Activate the conda environment.

```
conda activate qiime2-2021.4
```

To install this phyloseq package, start R (version "4.1.2") and enter:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```
