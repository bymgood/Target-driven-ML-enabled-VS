# Target-driven ML-enabled VS

echo Enter the uniprot id for the target of interest
read uniprot_ID
echo The uniprot id you entered is $uniprot_ID
echo Enter the working directory of this bash file of the screening package
read directory_1
echo The working directory is $directory_1
echo Let us get it started.

# Module 1
cd $directory_1/1_Target_expansion/
python Target_expansion.py -i $uniprot_ID -f $uniprot_ID -p 0.4
cp $uniprot_ID'.csv' $directory_1/2_Compound_retrieving/

# Module 2
cd $directory_1/2_Compound_retrieving/
python Compound_retrieving.py -i $uniprot_ID'.csv' -f $uniprot_ID'_compounds_collection'
cp *compounds_collection* $directory_1/3_Vectorization/

# Module 3
cd $directory_1/3_Vectorization/
python Vectorization.py -i 'actives_'$uniprot_ID'_compounds_collection.csv' -a active -f 'actives_'$uniprot_ID'_compounds_collection_fp'
python Vectorization.py -i 'inactives_'$uniprot_ID'_compounds_collection.csv' -a inactive -f 'inactives_'$uniprot_ID'_compounds_collection_fp'
cp *fp* $directory_1/4_ML_modeling_training/
cp *fp* $directory_1/6_Post_docking_analysis/

# Module 4
cd $directory_1/4_ML_modeling_training/
python ML_model_training.py -a actives* -i inactives* -f $uniprot_ID
cp *.sav $directory_1/5_Virtural_screening/

# Module 5
cd $directory_1/5_Virtural_screening/
python Virtual_screening.py -m *MLP.sav -t MLP -s Enamine_diversity_50K_morgan_1024_FP.csv -f $uniprot_ID'_50K_VS_MLP'
python Virtual_screening.py -m *random_forest.sav -t RF -s Enamine_diversity_50K_morgan_1024_FP.csv -f $uniprot_ID'_50K_VS_RF'
cp *VS* $directory_1/6_Post_docking_analysis/

# Module 6
cd $directory_1/6_Post_docking_analysis/
python Post_docking_analysis.py -r *MLP.csv -a actives* -i inactives*
python Post_docking_analysis.py -r *RF.csv -a actives* -i inactives*
cp *properties_calculated* $directory_1/7_Data_processing/

# Module 7
cd $directory_1/7_Data_processing/
python Data_processing.py -r *RF* -m *MLP* -f $uniprot_ID'target_driven_ML_enabled_VS'




