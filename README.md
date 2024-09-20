# RXNRECer
> This repo contains source codes for a Reaction Prediction tool namely RXNRECer.
> RXNRECer is a deep learning based tool for annotaion enzyme function to predict the reaction outcome using protein sequences. The tool is trained on a large dataset of enzyme-catalyzed reactions and can predict the outcome of new reactions with high accuracy.


# Getting Started
## Prerequisites
Conda environment with the following packages:
- PyTorch
- NumPy
- Pandas
- ....
The full list of required packages can be found in rxnrecer.yaml file.

## Installation
``` bash
git clone https://github.com/kingstdio/RXNRECer.git
cd RXNRECer
conda env create -f rxnrecer.yaml
```

## Usage
``` bash
 python production/production.py -i production/sample10-ecoli.fasta -o production/sample10-ecoli-res.json -f json
 ```


## Stargazers over time

[![Stargazers over time](https://starchart.cc/kingstdio/RXNRECer.svg)](https://github.com/kingstdio/RXNRECer/)

