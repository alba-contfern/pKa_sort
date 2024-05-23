# Acidity comparison

Welcome to the repository for the acidity comparison package. The objective of creating this package was to have a program able to sort a list of given molecules from least to most acid. This can be useful for first year organic chemists, when learning about acidity and how to classify acids from weakest to strongest. 

This package mainly uses rdkit and pubchempy packages, pubchem being the source for all the pka values that the package outputs. 

# Installation

To install this package, you need to make sure that the following packages are installed:
- rdkit
- pubchempy
- pandas
- requests

You can install them by entering the following lines into your terminal:

```bash
pip install rdkit==2022.9.5
pip install pubchempy==1.0.4
pip install pandas==2.1.4
pip install requests==2.31.0
```

Note: these are the versions of the packages that were used to build this code, it could potentially work with other versions of these packages.

# Usage

In order to use this package, you will need to import the package into your code, you can do so using the following code:

```bash
import pka_sort as pk
```

You can then call on any of the functions using `pk.function`.