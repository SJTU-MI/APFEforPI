# APFEforPI
Automated physical feature engineering for polymer informatics (APFEforPI), which has been utilized for the exploitation of high thermal conductivity amorphous polymers
## Description
Automated physical feature engineering for polymer informatics, including descriptor generation, statistical filtering, and machine learning model selection. The obtained optimized descriptors are available for the description of polymers to enable the construction of machine learning models with desired properties. It has been combined with high-throughput calculations to explore high thermal conductivity polymers and to understand the underlying physical mechanisms.  Please refer to our work "Informatics Algorithm-Assisted Exploration of High Thermal Conductivity Amorphous Polymers Using Automated Physical Feature Engineering and Symbolic Regression" for additional details.
![workflow](https://github.com/SJTU-MI/APFEforPI/blob/main/workflow.jpg)
## Installation
To download, clone this repository<br>
````
git clone https://github.com/SJTU-MI/APFEforPI.git
````
## Datasets

The training dataset and the screening datasets in this work are described in the following Table. We utilized 1051 polymer data with calculated TC by NEMD simulations as the training dataset A, which was sampled by polymer types from PoLyInfo database. Dataset B was collected manually from the PoLyInfo database and retained after de-duplication with dataset A. Dataset C is a specific dataset, and all polymers are polyimides, which were formed by dianhydride and diamine/diisocyanate from PubChem.
| **Name** | **Description** | **Included in Github** | **Source** |
|:--------:|:--------:|:--------:|:--------:|
|**Dataset A**| Benchmark data for ML model training | DatasetA.pkl | Sampling by polymer types from the PoLyInfo database | 
|**Dataset B**| PoLyInfo | DatasetB.pkl | Collected manually from the PoLyInfo database | 
|**Dataset C**| Polyimides | DatasetC.pkl | Hypothetical polyimides formed by dianhydride and diamine/diisocyanate from PubChem | 


## Authors

| **AUTHORS** |Xiang Huang, Shengluo Ma, Shenghong Ju            |
|-------------|--------------------------------------------------|
| **VERSION** | V1.0 / April,2023                               |
| **EMAILS**  | shenghong.ju@sjtu.edu.cn                         |

## Attribution
This work is under BSD-2-Clause License. Please, acknowledge use of this work with the appropiate citation to the repository and research article.
