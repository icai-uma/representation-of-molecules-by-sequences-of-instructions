# IsalChem: Representation of Molecules by Sequences of Instructions

This is a repository for **IsalChem**, a novel formal language intended to represent chemical compounds, described in the paper [Representation of Molecules by Sequences of Instructions](https://doi.org/10.1021/acs.jcim.5c00354)

This repository contains the source code and scripts employed in the study that introduce a novel approach to chemical nomenclature. The article presents a method based on a reduced instruction set that defines a language in which every string represents a valid molecule. Notably, minor alterations in these strings typically correspond to small modifications in the corresponding molecular structures. Computational experiments demonstrate a strong correlation between distances and standard chemical similarity measures in this string space. This innovative approach facilitates integration with state-of-the-art computational intelligence systems, including deep learning models, for advanced chemical information processing.

# Installation Instructions:

To set up the required environment, please execute the following commands:

```
conda install python=3.11

conda install pip numpy networkx python-pptx

pip install rdkit

mkdir sequence
```
