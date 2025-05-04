
A detailed guide to the core bioinformatics and computational biology concepts behind the ProtFeat project. This file is intended to help new users—especially those coming from outside biology—understand the terminology and logic behind the protein feature extraction techniques used in ProtFeat.

---

## 1. Central Dogma Refresher

The central dogma of molecular biology describes the flow of genetic information within biological systems:

**DNA → RNA → Protein**

### Key Steps:
- **Transcription**: DNA is transcribed into messenger RNA (mRNA).
- **Translation**: mRNA is read in triplets (codons) to assemble a chain of amino acids, forming a polypeptide.
- **Folding**: The polypeptide folds into a 3D structure to become a functional protein.

### Why it matters for ProtFeat:
ProtFeat begins with the protein—the final product of this process—and uses its amino acid sequence as the raw input for numerical feature extraction.

---

## 2. Protein Primary Structure & FASTA Format

### Protein Primary Structure

The **primary structure** of a protein is its linear sequence of amino acids, typically represented using one-letter codes:

**Example:**

``` MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQ... ```


Each letter corresponds to one of the 20 standard amino acids (e.g., M = Methionine, A = Alanine, etc.).

### FASTA Format

FASTA is a plain-text format commonly used for storing protein and nucleotide sequences.

**Example:**

``` >sp|Q9H9K5|SYT1_HUMAN Synaptotagmin-1 OS=Homo sapiens... ```
```MSEQNNTEMTFQIQRIY... ```


- The line starting with `>` is the **header**, containing metadata such as UniProt ID and description.
- The second line contains the actual sequence.

### Why it matters for ProtFeat:
- ProtFeat reads protein sequences from FASTA files.
- The `--ppid` parameter allows you to select which word from the header line is used as the unique protein ID.

---

## 3. Amino Acid Physicochemical Properties

Amino acids differ by their chemical and physical characteristics, which influence protein structure and function.

| Property        | Example Residues | Description                            |
|-----------------|------------------|----------------------------------------|
| Hydrophobicity  | L, I, V, F        | Avoid water                            |
| Charge (+)      | K, R              | Positively charged                     |
| Charge (–)      | D, E              | Negatively charged                     |
| Polarity        | S, T, N, Q        | Ability to form hydrogen bonds         |
| Size            | G, A vs. W, Y     | Small vs. bulky side chains            |

### Why it matters for ProtFeat:
Descriptors like CTD or PseAAC use these groupings to extract biologically meaningful statistics about the sequence.

---

## 4. Sequence Alignment Basics

Sequence alignment compares protein sequences to find regions of similarity.

| Concept             | Description                                                |
|---------------------|------------------------------------------------------------|
| Global Alignment    | Aligns full-length sequences (Needleman–Wunsch)            |
| Local Alignment     | Aligns best-matching subsequences (Smith–Waterman)         |
| Substitution Matrix | Scores residue substitutions (e.g., BLOSUM62, PAM)         |

### Why it matters for ProtFeat:
PSSM-based descriptors depend on multiple sequence alignments, typically generated using **PSI-BLAST**, which produces conservation scores used in feature extraction.

---

## 5. Position-Specific Scoring Matrix (PSSM)

A **PSSM** is a matrix of size *L × 20*, where *L* is the sequence length and 20 columns represent the 20 standard amino acids.

Each cell contains a **log-odds score** indicating the likelihood of a particular amino acid occurring at a specific position based on aligned homologous sequences.

**Simplified example:**

| Position | A  | R  | N  | D  | ... | Y  |
|----------|----|----|----|----|-----|----|
| 1        | 4  | -1 | 0  | 0  | ... | -2 |
| 2        | -2 | 6  | 0  | -1 | ... | 0  |

### Why it matters for ProtFeat:
ProtFeat includes multiple **PSSM-based descriptors** (via the POSSUM toolkit) that derive numeric features from evolutionary conservation profiles.

---

## 6. Protein Feature Descriptors

Descriptors are numerical summaries of protein sequences. ProtFeat supports both sequence-based and PSSM-based descriptors.

| Code      | Full Name                        | Description                                      | Typical Vector Length |
|-----------|----------------------------------|--------------------------------------------------|------------------------|
| AAC       | Amino Acid Composition           | Frequency of each amino acid                    | 20                     |
| DPC       | Dipeptide Composition            | Frequency of each 2-letter amino acid pair      | 400                    |
| CTD       | Composition/Transition/Dist.     | Group-wise counts based on properties           | 147                    |
| PseAAC    | Pseudo Amino Acid Composition    | Sequence + physicochemical + order info         | Varies                 |
| PSSM-AC   | PSSM Autocorrelation             | Correlation between positions in the PSSM       | Typically 60           |

### Why it matters for ProtFeat:
These descriptors convert variable-length protein sequences into fixed-length feature vectors for input to machine learning models.

---

## 7. Databases & Identifiers

ProtFeat relies on standard biological identifiers used in well-known databases.

| Database   | Example ID     | Description                                  |
|------------|----------------|----------------------------------------------|
| UniProtKB  | Q9H9K5          | Curated protein knowledgebase                |
| NCBI nr    | NP_000537.3     | Non-redundant protein sequence database       |
| Pfam       | PF00096         | Protein domain/family annotations             |

### Why it matters for ProtFeat:
- FASTA files often come from UniProt or NCBI.
- These IDs are used to track protein features and link sequences to known biological metadata.

---

## 8. Feature Engineering in Proteomics

**Feature engineering** transforms protein sequences into numerical representations suitable for predictive modeling.

### Common ML applications:
- Predicting enzyme/non-enzyme
- Identifying subcellular localization
- Predicting protein–protein interactions
- Classifying disease-causing variants

### Typical workflow:

```
Protein FASTA file
↓
(Optional) PSI-BLAST → PSSM
↓
ProtFeat feature extraction
↓
Tabular feature matrix (.tsv)
↓
Machine learning model
```
### Why it matters for ProtFeat:
ProtFeat provides the core feature extraction step in this pipeline—converting biological sequences into structured data.

---

## Further Reading

**Books:**
- Baxevanis & Ouellette, *Bioinformatics: A Practical Guide to the Analysis of Genes and Proteins*
- Jones & Pevzner, *An Introduction to Bioinformatics Algorithms*

**Online Resources:**
- [UniProt](https://www.uniprot.org)
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov)
- [Pfam](http://pfam.xfam.org)

---

