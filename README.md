# Glycemic-Project 🍚
A project on Glycemic Index that I have worked on involving bioinformatics.

## About the Data 💡
The data I have been using had been obtained from blood samples acquired from subjects involved in a randomized controlled trial (RCT) test by Dr Sumanto Haldar and his team. The analysis is a follow-up from the following study:
- [High or low glycemic index (GI) meals at dinner results in greater postprandial glycemia compared with breakfast: a randomized controlled trial](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7202752/)

The blood samples were analyzed by an Agilent 6495C Triple Quadrupole LCMS by Dr Leroy S. Pakkiri in Dr Chester Drum's Lab in the National University of Singapore (NUS). The raw data output (.m) from the LCMS were then converted to .mzML files and processed by the [MRMkit program](https://github.com/cssblab/MRMkit) by the 
Computational & Statistical Systems Biology Laboratory (CSSBLab) for automated peak detection and integration of the MS output. The program also corrects for batch and column variances that may arise from the usage of different 96-well plates and columns in the LCMS.

The resulting output data from MRMkit (Eg. _quant tables_) have been well-checked through for potential machine/procedural biases and are then supplemented with relevant metadata that help to categorize subjects for analysis.

## Contrasts of Interests 🔄
Each subject are measured under 4 circumstances, each across 5 time points. Each subject will have thus have 20 samples in total. The statuses/classifications that each subject adopt are coded as such:
- A: P-E
- B: P-F
- C: Q-E
- D: Q-F

We seek to analyze the metabolite contrasts between subjects of status P and Q, as well as between E and F. Subsequently, we can look at the contrasts between P-E and P-F, Q-E and Q-F, P-E and Q-E, and P-F and Q-F.

## Input Files Required 📁
1. 20210817_quant_table_glycemic.txt
  - This is the quant_table output from MRMkit supplemented with relevant categorical metadata. It will be the main data file that is used throughout the analysis.
2. quant_table_pre_correction.txt
  - This 
