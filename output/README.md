# Output

useful.stuff.aa has convenience functions used by some of the scripts in this workflowr

15C_Data and 15D_Data are files published with the journal article as supplements noting raw data for graphs showing n<20.



To use model training code, one must first generate ARTEMIS and DELFI features using the following code repositories: 
<https://github.com/cancer-genomics/delfi3>
<https://github.com/cancer-genomics/artemis_pipeline>

Instructions to train LCr model:
1. Process outputs from artemis and delfi pipelines using cohort_create.r
2. train model with liver_train.r

Please note that for convenience in using model training code in this folder only, the outputs of preprocess.r designate individuals with any liver diseases as "cancer" and all No Known Liver Disease individuals as "healthy". This is a convenience for training a binary classifier - please refer to ../data/All_LCr_Clinical_Metadata_withscores.csv and ../data/All_Cross_Reactivity_Metadata_withscores.csv for true clinical metadata. 

