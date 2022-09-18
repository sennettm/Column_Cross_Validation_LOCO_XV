# Leave-One-Column-Out Cross-Validation (LOCO_XV)
**still in beta** 
## Overview

**What is LOCO_XV?**

Leave-One-Column-Out Cross-Validation (LOCO_XV) is a model-selection criteria for evolutionary models in phylogenetics, similar to the Akaike Information Criterion (AIC).
It is a way to quantify how predictive model parameters are. The higher the LOCO_XV, the more predictive a model.

**What is the motivation behind LOCO_XV?**

Common model-selection criteria in phylogenetics, like AIC or BIC, require that you specificy the number of model parameters and/or the number of data points. These
are ambigous in a phylogenetic analysis. How many parameters is it to model evolution as a tree? How many data points are in your dataset? Is it the number of
columns, the number of sequences, the number of amino acids? 

We don't really know the answers to these questions, but by using LOCO_XV we can totally avoid them. 

In addition, phylogenetic analyses often violate assumptions in model-selection criteria like the AIC. For example, the selected model should be "close" to the "true"
model.

Again, we don't know that this is true, but we needn't worry about this when it comes to LOCO_XV.

**How does LOCO_XV work?**

A single column is removed from an alignment. Model parameters and a phylogenetic tree are inferred for that column. The model parameters and tree are used to 
calculate the Ln-likelihood of that removed column. The process is repeated for each column in the alignment.

**What can we do with LOCO_XV?**

We can calculate the predictive Ln-likelihood for different models and determine which is most predictive. 

## Implementing

**Requirements for LOCO_XV**

IQ-Tree v2.2: http://www.iqtree.org/#download
 
***Please note ESR is currently only compatible with IQ-Tree files.***

Dendropy: https://dendropy.org/downloading.html
```
python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git
```

NumPy: https://numpy.org/install/
```
pip install numpy
```

Biopython: https://biopython.org/wiki/Download
```
pip install biopython
```

## Output

siteloglik.txt; a file containing the individual site Ln-likelihoods 
IQTree_Data; a folder containing all the cross-validation alignments, iqtree outputs, etc...

## Running

Copy the scripts and one of the test data sets into some test folder. 

Usage
```
./lco_xv_v0.1.sh -h 
```

```
required:
-a input alignment in fasta format  
-m evolutionary model in IQ-Tree format (e.g. LG+FO+G4) 
-N total number of threads available to work with 
-n number of threads for individual IQ-Tree run

optional: 
-t input starting tree in newick format
```

Example: Perform LOCO_XV on an Apicomplexan L/MDH dataset on a computer with with 8 threads
```
./lco_xv_v0.1.sh -a Apico2020_seqs.fasta -t Apico2020_seqs.fasta.treefile -m LG -N 12 -n 4
```
