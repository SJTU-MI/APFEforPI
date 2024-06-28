# APFEforPI
Automated physical feature engineering for polymer informatics (APFEforPI), which has been utilized for the exploitation of high thermal conductivity amorphous polymers.
## Description
Automated physical feature engineering for polymer informatics, including descriptor generation, statistical filtering, and machine learning model selection. The obtained optimized descriptors are available for the description of polymers to enable the construction of machine learning models with desired properties. It has been combined with high-throughput calculations to explore high thermal conductivity polymers and to understand the underlying physical mechanisms.  Please refer to our work "High-Throughput Screening of Amorphous Polymers with High Intrinsic Thermal Conductivity via Automated Physical Feature Engineering" for additional details.
![workflow](https://github.com/SJTU-MI/APFEforPI/blob/main/workflow.jpg)
## Installation
### Files loading and environment setup:

To download, clone this repository:<br>
````
git clone https://github.com/SJTU-MI/APFEforPI.git
````

To run most code in this repository, the relevant anaconda environment can be installed from environment.yml. To build this environment, run：<br>
````
cd ./APFEforPI
conda env create -f environment.yml
conda activate apfeforpi
````
However, some large data files are downloaded from external release repositories. To build this environment, run：<br>
````
wget -i file.txt
chmod 777 file.sh
./file.sh
````
To calculate properties such as thermal conductivity of amorphous polymers, an additional [RadonPy](https://github.com/RadonPy/RadonPy) toolkit is required, run：<br>
````
conda install -c psi4 -c conda-forge rdkit psi4 resp mdtraj matplotlib lammps
git clone -b main https://github.com/RadonPy/RadonPy.git
export PYTHONPATH=<Path-to-RadonPy>:$PYTHONPATH
````

### Try the desired parts of the project:
**01Des_monomer.py**: Polymer monomer descriptors, including van der Waals volume (VDW), moleculamr weight (MW), monomer length (Monomer_length) and ratio of backbone molecular weight to monomer molecular weight (MW_ratio)<br>
**02Des_Mordred.py**: Descriptors of polymer repetitive unit properties calculated by [Mordred](https://github.com/mordred-descriptor/mordred) software<br>
**03Des_MD.py**: Force field descriptors extracted from the MD file after GAFF2 force field assignment<br>
**04Des_merger.py**: Merging polymer descriptors from various sources<br>
**05Statistical_filtering.ipynb**: Statistical filtering of physical descriptors, including removing descriptors with low variance, feature filtering combined with evaluation of Pearson, Spearman, Distance and maximum information coefficients<br>
**06RFEforOpt_descriptors.ipynb**: Physical descriptors optimization realized by recursive feature elimination (RFE) with a hyperparameter-optimized random forest model<br>
**07ML_Opt.ipynb**: Construction of optimized descriptors and machine learning model pairs<br>
**08Virtual_Screening.ipynb**: High throughput virtual screening of high thermal conductivity polymers in the Dataset B and Dataset C<br>
**09SHAP_RF.ipynb**: Feature importance analysis based on SHapley Additive exPlanations([SHAP](https://github.com/slundberg/shap))<br>
**10TC_Cal.py**: A case for modeling and calculating thermal conductivity of the amorphous polymer, taking polyethylene as an example<br>
**Des_genrate.sh**: Automatic generation of physical descriptors for polymers. We have provided a demo for verification.<br>
**Simulation_data.csv**: 104 polymers with simulation data via [RadonPy](https://github.com/RadonPy/RadonPy) software.

## Datasets

The training dataset and the screening datasets in this work are described in the following Table. We utilized 1051 polymer data with calculated TC by nonequilibrium molecular dynamics (NEMD) simulations as the training dataset A, which was sampled by polymer types from PoLyInfo database. Dataset B was collected manually from the PoLyInfo database and retained after de-duplication with dataset A. Dataset C is a specific dataset, and all polymers are polyimides, which were formed by dianhydride and diamine/diisocyanate from PubChem.
| **Name** | **Description** | **Included in Github** | **Source** |
|:--------:|:--------:|:--------:|:--------:|
|**Dataset A**| Benchmark data for ML model training| DatasetA.pkl | [Sampling by polymer types from the PoLyInfo database](https://doi.org/10.1038/s41524-022-00906-4)  | 
|**Dataset B**| PoLyInfo | DatasetB.pkl | [Collected manually from the PoLyInfo database](https://polymer.nims.go.jp/)  | 
|**Dataset C**| Polyimides| DatasetC.pkl | [Hypothetical polyimides formed by dianhydride and diamine/diisocyanate from PubChem](https://data.caltech.edu/records/hend5-jzt61) | 

## Authors

| **AUTHORS** |Xiang Huang, Shengluo Ma, Shenghong Ju            |
|-------------|--------------------------------------------------|
| **VERSION** | V1.0 / April,2023                               |
| **EMAILS**  | shenghong.ju@sjtu.edu.cn                         |

## Publications
1. X. Huang, S. Ma, Y. Wu, C. Wan, C.Y. Zhao, H. Wang, S. Ju, "High-throughput screening of amorphous polymers with high intrinsic thermal conductivity via automated physical feature engineering," Journal of Materials Chemistry A (2023) [[Link](https://pubs.rsc.org/en/content/articlelanding/2023/ta/d3ta03370h)].
2. X. Huang, S. Ju, Tutorial: AI-assisted exploration and active design of polymers with high intrinsic thermal conductivity[J]. Journal of Applied Physics. (2024) [[Link](https://doi.org/10.1063/5.0201522)].


## Attribution
This work is under BSD-2-Clause License. Please, acknowledge use of this work with the appropiate citation to the repository and research article.
