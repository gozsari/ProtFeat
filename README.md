## Description: 

The feature_extracter class is designed to extract the protein features by employing POSSUM and iFeature python-based tools.

POSSUM (Position-Specific Scoring matrix-based feature generator for machine learning),
a versatile toolkit with an online web server that can generate 21 types of PSSM-based feature descriptors,
thereby addressing a crucial need for bioinformaticians and computational biologists.

 iFeature, a versatile Python-based toolkit for generating various numerical feature representation schemes for
 both protein and peptide sequences. iFeature is capable of calculating and extracting a comprehensive spectrum
 of 18 major sequence encoding schemes that encompass 53 different types of feature descriptors.

<br/>Table.2: Protein descriptors obtained from the POSSUM tool.

| Descriptor group | Protein descriptor                                                                               | Number of dimensions
----------------------------------------|--------------------------------------------------------------------------------------------------|-----------------------------
Row Transformations| AAC-PSSM<br/>D-FPSSM<br/>smoothed-PSMM<br/>AB-PSSM<br/>PSSM-composition<br/>RPM-PSSM<br/>S-FPSSM |20<br/>20<br/>1000<br/>400<br/>400<br/>400<br/>400
Column Transformation| DPC-PSSM <br/> k-seperated-bigrams-PSSM <br/> tri-gram-PSSM <br/> EEDP <br/> TPC                 | 400 <br/> 400 <br/> 8000 <br/> 4000 <br/>  400
Mixture of row and column transformation | EDP <br/> RPSSM <br/> Pre-PSSM <br/> DP-PSSM <br/> PSSM-AC <br/> PSSM-CC | 20 <br/> 110 <br/> 40 <br/> 240 <br/> 200 <br/> 3800
Combination of above descriptors| AADP-PSSSM <br/> AATP <br/> MEDP | 420 <br/> 420 <br/> 420

<br/>Table.3: Protein descriptors obtained from the İFeature tool.

| Descriptor group | Protein descriptor                                                    | Number of dimensions
----------------------------------------|-----------------------------------------------------------------------|-----------------------------
Amino acid composition| Amino acid composition (AAC) <br/> Composition of k-spaced amino acid pairs (CKSAAP) <br/>Dipeptide composition (DPC)|20<br/>2400<br/>400
Grouped amino acid composition | Grouped amino acid composition (GAAC) <br/> Composition of k-spaced amino acid group pairs (CKSAAGP) <br/> Grouped dipeptide composition (GDPC)| 5 <br/> 150 <br/> 25
Autocorrelation| Moran (Moran) <br/> Geary (Geary) <br/> Normalized Moreau-Broto (NMBroto) | 240 <br/> 240 <br/> 240 
C/T/D| Composition (CTDC) <br/> Transition (CTDT) <br/> Distribution (CTDD)                                             | 39 <br/> 39 <br/> 195
Conjoint triad| Conjoint triad (CTriad) <br/> Conjoint k-spaced triad (KSCTriad) | 343 <br/> 343*(k+1) 
Quasi-sequence-order| Sequence-order-coupling number (SOCNumber) <br/> Quasi-sequence-order descriptors (QSOrder)     | 60 <br/> 100 
Pseudo-amino acid composition| Pseudo-amino acid composition (PAAC) <br/> Amphiphilic PAAC (APAAC)                                     | 50 <br/> 80 


## Explanation of Parameters
Here, we explain about
**protein_feature: {string}, (default = 'aac_pssm'):** one of the 21 PSMM-based protein descriptors in POSSUM.

* **POSSUM descriptors:** aac_pssm, d_fpssm, smoothed_pssm, ab_pssm, pssm_composition, rpm_pssm,
s_fpssm, dpc_pssm, k_separated_bigrams_pssm, eedp, tpc, edp, rpssm,
pse_pssm, dp_pssm, pssm_ac, pssm_cc, aadp_pssm, aatp, medp , or all_POSSUM
<br/><br/>
  * all_POSSUM: it extracts the features of all (21) POSSUM protein descriptors

    
* or one of the 18 protein descriptors in iFeature.

* **iFeature descriptors:** AAC, PAAC, APAAC, DPC, GAAC, CKSAAP, CKSAAGP, GDPC, Moran, Geary,
NMBroto, CTDC, CTDD, CTDT, CTriad, KSCTriad, SOCNumber, QSOrder, or all_iFeature

  * all_iFeature: it extracts the features of all (18) iFeature protein descriptors
<br/><br/>

**place_protein_id: {int}, (default = 1):** It indicates the place of protein id in fasta header.
e.g. fasta header: >sp|O27002|....|....|...., seperate the header wrt. '|' then >sp is
in the zeroth position, protein id in the first(1) position.

**input_folder: {string}, (default = 'input_folder'}:** it is the path to the folder that contains the fasta file.

**fasta_file_name: {string}, (default ='sample'):** it is the name of the fasta file exclude the '.fasta' extension.


## Input file 

* It must be in fasta format
* A sample is also given as **feature_extraction_module/input_files/sample.fasta**


## Downloading Position Specific Scoring Matrices (PSSMs) (optional)
Since extracting PSSMs takes time for POSSUM descriptors, this step is to accelerate 
the feature extraction process. We strongly recommend you to download PSSMs, 
if you have more than 100 proteins in your fasta file(s). 

We have extracted PSSMs for 558,419 proteins available in UniProtKB / SwissProt database.
When you run the module to extract the POSSUM features, an option will be given whether you 
prefer to download PSSM file. If you choose **Y** (Yes), then, the program will automatically 
download and copy the related PSSMs to **feature_extraction_module/pssm_files** folder.

Each PSSM file in the folder is named as **proteinID.pssm**.


## Output file

* The extracted feature files will be located under
**feature_extraction_module/output_folder** 
folder with the name: **fasta_file_name_protein_feature.txt** (e.g. sample_AAC.txt)
* The content of the output files: 
  * Each row corresponds to the extracted features of the protein sequence
  * The first column of each row is protein id (in UniProtKB), 
    the rest is extracted features of the protein sequence.

## License

MIT License

ProFab Copyright (C) 2020 CanSyL

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
