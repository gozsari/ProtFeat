<p align="center">
<img width="850" height="100"  src="image/protFeat_banner.png">

![Linux](https://svgshare.com/i/Zhy.svg)
![version](https://img.shields.io/badge/version-1.0-blue)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
![Python](https://img.shields.io/badge/python-v3.9-red)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/Naereen/StrapDown.js/blob/master/LICENSE)
[![Open Source Love svg1](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)

## Description: 

ProtFeat is designed to extract the protein features by employing POSSUM and iFeature python-based tools. ProtFeat includes a total of 39 distinct protein feature extraction methods using 21 PSSM-based protein descriptors from POSSUM and 18 protein descriptors from iFeature.

[POSSUM](https://academic.oup.com/bioinformatics/article/33/17/2756/3813283) (Position-Specific Scoring matrix-based feature generator for machine learning),
a versatile toolkit with an online web server that can generate 21 types of PSSM-based feature descriptors,
thereby addressing a crucial need for bioinformaticians and computational biologists.

 [iFeature](https://academic.oup.com/bioinformatics/article/34/14/2499/4924718), a versatile Python-based toolkit for generating various numerical feature representation schemes for
 both protein and peptide sequences. iFeature is capable of calculating and extracting a comprehensive spectrum
 of 18 major sequence encoding schemes that encompass 53 different types of feature descriptors.


## Installation

ProtFeat is a python package for feature extracting from protein sequences written in Python 3.9. ProtFeat was developed and tested in Ubuntu 20.04 LTS. Please make sure that you have **Anaconda** installed on  your computer and  run the below commands to install requirements. Dependencies are available in requirements.txt file.

```
conda create -n protFeat_env python=3.9
conda activate protFeat_env
```

## How to run ProtFeat to extract the protein features 
Run the following commands in the given order:
### To use ProtFeat as a python package:
```
pip install protFeat
```
Then, you may use protFeat as the following in python:
```
import protFeat
from protFeat.feature_extracter import extract_protein_feature, usage
usage()
extract_protein_feature(protein_feature, place_protein_id, input_folder, fasta_file_name)
```
For example, 
```
extract_protein_feature("AAC", 1, "input_folder", "sample")
```

### To use ProtFeat from terminal:
Clone the Git Repository.
```
git clone https://github.com/gozsari/ProtFeat
```
In terminal or command line navigate into **protFeat** folder.
```
cd ProtFeat
```
Install the requirements by the running the following command.
```
pip install -r requirements.txt
```
Altenatively you may run ProtFeat from the terminal as the following:
```
cd src
python protFeat_command_line.py --pf protein_feature --ppid place_protein_id --inpf input_folder --fname fasta_file_name
```
For example, 
```
python protFeat_command_line.py --pf AAC --ppid 1 --inpf input_folder --fname sample
```
## Explanation of Parameters

**protein_feature: {string}, (default = 'aac_pssm'):** one of the protein descriptors in POSSUM and iFeature.

POSSUM descriptors: 
```
aac_pssm, d_fpssm, smoothed_pssm, ab_pssm, pssm_composition, rpm_pssm,
s_fpssm, dpc_pssm, k_separated_bigrams_pssm, eedp, tpc, edp, rpssm,
pse_pssm, dp_pssm, pssm_ac, pssm_cc, aadp_pssm, aatp, medp , or all_POSSUM
```

Note: all_POSSUM extracts the features of all (21) POSSUM protein descriptors.

iFeature descriptors:
```
AAC, PAAC, APAAC, DPC, GAAC, CKSAAP, CKSAAGP, GDPC, Moran, Geary,
NMBroto, CTDC, CTDD, CTDT, CTriad, KSCTriad, SOCNumber, QSOrder, or all_iFeature
```

Note: all_iFeature extracts the features of all (18) iFeature protein descriptors.

**place_protein_id: {int}, (default = 1):** It indicates the place of protein id in fasta header.
e.g. fasta header: >sp|O27002|....|....|...., seperate the header wrt. '|' then >sp is
in the zeroth position, protein id in the first(1) position.

**input_folder: {string}, (default = 'input_folder'}:** it is the path to the folder that contains the fasta file.

**fasta_file_name: {string}, (default ='sample'):** it is the name of the fasta file exclude the '.fasta' extension.


## Input file 
It must be in fasta format.

## Output file
The extracted feature files will be located under
**feature_extraction_output** 
folder with the name: **fasta_file_name_protein_feature.txt** (e.g. sample_AAC.txt).

The content of the output files: 
  * The output file is *tab-seperated*.
  * Each row corresponds to the extracted features of the protein sequence.
  * The first column of each row is [UniProtKB](https://www.uniprot.org/) id of the proteins, 
      the rest is extracted features of the protein sequence.
      
## Tables of the available protein descriptors

Table 1: Protein descriptors obtained from the POSSUM tool.

| Descriptor group | Protein descriptor                  | Number of dimensions
-------------------|-------------------------------------|---------------------
Row Transformations| AAC-PSSM<br/>D-FPSSM<br/>smoothed-PSMM<br/>AB-PSSM<br/>PSSM-composition<br/>RPM-PSSM<br/>S-FPSSM |20<br/>20<br/>1000<br/>400<br/>400<br/>400<br/>400
Column Transformation| DPC-PSSM <br/> k-seperated-bigrams-PSSM&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br/> tri-gram-PSSM <br/> EEDP <br/> TPC                 | 400 <br/> 400 <br/> 8000 <br/> 4000 <br/>  400
Mixture of row and column transformation | EDP <br/> RPSSM <br/> Pre-PSSM <br/> DP-PSSM <br/> PSSM-AC <br/> PSSM-CC | 20 <br/> 110 <br/> 40 <br/> 240 <br/> 200 <br/> 3800
Combination of above descriptors| AADP-PSSSM <br/> AATP <br/> MEDP | 420 <br/> 420 <br/> 420

<br/>
Table 2: Protein descriptors obtained from the iFeature tool.

| Descriptor group | Protein descriptor| Number of dimensions
-------------------|--------------------|--------------------
Amino acid composition| Amino acid composition (AAC) <br/> Composition of k-spaced amino acid pairs (CKSAAP) <br/>Dipeptide composition (DPC)|20<br/>2400<br/>400
Grouped amino acid composition | Grouped amino acid composition (GAAC) <br/> Composition of k-spaced amino acid group pairs (CKSAAGP) <br/> Grouped dipeptide composition (GDPC)| 5 <br/> 150 <br/> 25
Autocorrelation| Moran (Moran) <br/> Geary (Geary) <br/> Normalized Moreau-Broto (NMBroto) | 240 <br/> 240 <br/> 240 
C/T/D| Composition (CTDC) <br/> Transition (CTDT) <br/> Distribution (CTDD)                                             | 39 <br/> 39 <br/> 195
Conjoint triad| Conjoint triad (CTriad) <br/> Conjoint k-spaced triad (KSCTriad) | 343 <br/> 343*(k+1) 
Quasi-sequence-order| Sequence-order-coupling number (SOCNumber) <br/> Quasi-sequence-order descriptors (QSOrder)     | 60 <br/> 100 
Pseudo-amino acid composition| Pseudo-amino acid composition (PAAC) <br/> Amphiphilic PAAC (APAAC)                                     | 50 <br/> 80 


## License

MIT License

ProtFeat Copyright (C) 2020 CanSyL

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
