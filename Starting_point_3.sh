# Target-driven ML-enabled VS

echo Starting point 3 - specify active and inactive compound list
echo Enter the target name
read uniprot_ID
echo The target name you entered is $uniprot_ID
echo Enter the active compound list
read active_list
echo The active compound list you entered is $active_list
echo Enter the inactive compound list
read inactive_list
echo The inactive compound list you entered is $inactive_list
echo Enter the working directory of this bash file of the screening package
read directory_1
echo The working directory is $directory_1
echo Let us get it started.


# Module 3
cd $directory_1/3_Vectorization/
python Vectorization.py -i $active_list -a active -f 'actives_'$active_list'_compounds_collection_fp'
python Vectorization.py -i $inactive_list -a inactive -f 'inactives_'$inactive_list'_compounds_collection_fp'
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




