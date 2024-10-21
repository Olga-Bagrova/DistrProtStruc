# DistrProtStruc
>Method for research at [Lomonosov Moscow State University](https://msu.ru/en) focused on the distribution of secondary protein structure

<h3> Analysis of Secondary Structure Distributions along Polypeptide Chains of Proteins within Different Functional Classes, Homologous, and Topological Proteins </h3>

**Supervisor:** Ekaterina Belova


## Introduction

Proteins are a highly diverse group of biomolecules that are responsible for a wide range of processes in living organisms. It is thought that the function of proteins is determined by their structure. Moreover, the spatial structure of a protein may be more conserved than its amino acid sequence. In this context, the study of the relationship between protein structure and function can have a significant impact on both fundamental biophysical aspects related to folding and practical applications concerning the prediction of functions and other properties of polypeptides.

The group's research has focused on various groups of proteins and includes several related projects. Overall, this study **aims to identify the structural features of different protein groups**.

A method for calculating the distributions of secondary structures along  polypeptide chains has been developed to identify structural patterns in various protein groups. 

To be able to compare proteins of different lengths, the polypeptide chains are normalized. The results are presented as diagrams showing the distribution of secondary structures at normalized lengths, indicating the frequencies of occurrence of different econdary structure types at corresponding positions.

A large dataset of proteins with available PDB structures was analyzed, dividing it into groups based on molecular function, homology and topology, , in an effort to discover common fundamental patterns in their secondary structures. The results of our research are presented in [publications](#papers).


## How to install

```
git clone https://github.com/Olga-Bagrova/DistrProtStruc.git
```


## Usage

You can use a similar approach to extend your analysis. All you need is a list of PDB IDs and positions that need to be processed. To learn how to perform the analysis, please examine the `Usage_example.ipynb` file.


## Repo content
- `README.md`
- `Usage_example.ipynb` - notebook with usage example
- `distrprotstruc.py` - source code
- `example.xlsx` - input example
- `requirements.txt` - requirements
  

## Papers

1.	**Bagrova O**, Lapshina K, Sidorova A, Shpigun D, Lutsenko A, Belova E. (2024). Secondary structure analysis of proteins within the same topology group. Biochemical and Biophysical Research Communications, 734: 150613. [doi.org/10.1016/j.bbrc.2024.150613](https://doi.org/10.1016/j.bbrc.2024.150613)
2.	Tverdislov VA, Sidorova AE, **Bagrova OE**, Belova EV, Bystrov VS, Levashova NT, Lutsenko AO, Semenova EV, Shpigun DK. (2022). Chirality As a Symmetric Basis of Self-Organization of Biomacromolecules. Biophysics, 67(5): 673-691. [doi.org/10.1134/S0006350922050190](https://doi.org/10.1134/S0006350922050190)
3.	Sidorova AE, Malyshko EV, Lutsenko AO, Shpigun DK, **Bagrova OE**. (2021). Protein Helical Structures: Defining Handedness and Localization Features. Symmetry, 13(5): 879. [doi.org/10.3390/sym13050879](https://doi.org/10.3390/sym13050879)
4.	Malyshko EV, Semenova EV, **Bagrova OE**, Murtazina AR, Tverdislov VA. (2021). Chiral Dualism as a Unifying Principle in Molecular Biophysics. Biophysica, 1(1): 22-37. [doi.org/10.3390/biophysica1010003](https://doi.org/10.3390/biophysica1010003)
5.	Malyshko EV, **Bagrova OE**, Tverdislov VA. (2020). The Relationship between Hierarchical Chiral Structures of Proteins and Their Functions. Biophysics, 65(3): 368-373. [doi.org/10.1134/S0006350920030148](https://doi.org/10.1134/S0006350920030148)
