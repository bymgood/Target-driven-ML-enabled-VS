# Target-driven-ML-enabled-VS (TAME-VS)
Target-driven machine-learning-enabled virtual screening is a machine learning tool developed to accelerate the early-stage hit identification. In this repository, you can find all you need to launch virtual screening against your target protein with different types of input information:
- Mimimum info: a single uniport ID of your target protein (starting point 1); 
- Intermediate info: Compound datasets containing active and inactive molecules against homologies of your target of interest (starting point 2);
- Advanced info: Compound datasets containing active molecules directly against your target of interest (starting point 3).    
For more information, please refer to our paper (Yuemin Bian): https://www.frontiersin.org/articles/10.3389/fmolb.2023.1163536/full


![alt text](https://www.frontiersin.org/files/Articles/1163536/fmolb-10-1163536-HTML/image_m/fmolb-10-1163536-g001.jpg)
Scheme of the workflow implemented in the TAME-VS platform.

# About this repository
The Starting_point[1-3].sh are the production scripts for launching vritual screening with different types of inputs. Clone the entire repository to your local machine prior to start.

# Prerequisties
To run the virtual screening, you need to install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Then you need to 
- Setup a conda virtual environment with python (3.7) and activate it
```bash
conda create -n TAME_VS python=3.7
conda activate TAME_VS
```
- Install all third-party packages 
```bash
pip install -r requirements.txt
```

# Preparing molecular fingerprints for virtual screening
To use the included Enamine 50k compound library for final ML virtual screening, please run the followig command from **5_Virtural_screening**
```bash
python Library_preparation.py -i Enamine_diversity_50K.csv -s 1 -c 2  -f Enamine_diversity_50K_morgan_1024_FP
```

# Run target-driven machine-learning-enabled VS
The following example uses starting point 1 as an example.
1. Search for the uniport ID of your target proteins (e.g. P09238);
2. Launch the ```Starting_point_1.sh``` script and provide uniport ID and working directory;
```bash
bash Starting_point_1.sh
```
