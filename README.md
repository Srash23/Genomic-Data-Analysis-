# DNA-to-Protein Translation and Amino Acid Frequency Analysis

## Introduction
This project focuses on the computational translation of DNA sequences into proteins, highlighting the fundamental process by which genetic information (genotype) manifests as physical and biochemical traits (phenotype). Using bioinformatics algorithms, this tool transcribes DNA into RNA and subsequently translates RNA into proteins, enabling researchers to explore amino acid frequency differences between gene sequences.

By comparing protein sequences across species, the pipeline aids in disease research, evolutionary studies, and biotechnology applications. The tool generates visual and textual reports, making it a valuable asset for understanding genetic mechanisms and protein function.

## Why Is This Important?
Understanding how genetic information is translated into proteins is crucial in:

**1. Disease Research –** Identifies mutations affecting protein synthesis.

**2. Evolutionary Biology –** Compares protein sequences across species to study evolutionary divergence.

**3. Biotechnology & Genetic Engineering –** Helps design synthetic genes for industrial and medical applications.

## Key Features

**1. DNA-to-RNA Transcription –** Converts DNA sequences into messenger RNA (mRNA).

**2. RNA-to-Protein Translation –** Uses a codon table to generate protein sequences.

**3. Amino Acid Frequency Analysis –** Compares frequency differences across species.

**4. Common Protein Identification –** Finds conserved protein sequences between organisms.

**5. Data Visualization –** Generates bar charts to visualize amino acid differences.

## Installation and Usage

#### Dependencies

**1. Install the required Python libraries using:**
   
pip install matplotlib json

**3. Ensure you have the required input files:**
   
DNA sequences in .txt format.
Codon table in .json format.

## Pipeline Workflow
![Pipeline Workflow - visual selection](https://github.com/user-attachments/assets/7b2ed60c-83aa-4d7b-a086-1009e8f89ca5)

**Step 1: DNA-to-RNA Transcription**


This function reads a DNA sequence from a file and transcribes it into RNA:

RNA = dna[startIndex:].replace('T', 'U')

Input: DNA sequence (e.g., ATGCGTACG...)

Output: RNA sequence (e.g., AUGCGUACG...)


**Step 2: RNA-to-Protein Translation**


The RNA sequence is split into codons (triplets) and mapped to amino acids using a codon table.

Protein.append(codonD[codon]) 

Input: RNA sequence (e.g., AUG-CGU-ACG...)

Output: Protein sequence (e.g., Methionine-Arginine-Threonine...)


**Step 3: Identifying Common Proteins**


This function compares protein sequences across species:

Set_of_Common_Proteins = Unique_Proteins_List_1.intersection(Unique_Proteins_List_2)

Input: Protein lists from two species.

Output: List of shared proteins.


**Step 4: Amino Acid Frequency Analysis**


The relative frequency of each amino acid is calculated:

Frequuency_of_aa_List_1 = aa_count_in_List_1.get(u, 0) / Total_aa_List_1

Input: Protein sequences from different species.

Output: Amino acid frequency table.


**Step 5: Visualization**


The function creates a bar chart comparing amino acid frequencies:

plt.bar(Position_bar_1, freqList1, width = 0.3, label = label1, color = '#0F52BA') 

Output: A bar chart displaying amino acid frequency variations.

## Usage
To run the program on two species' DNA sequences, execute:

runFullProgram()

_**Make sure your files are structured correctly:**_

_**1. DNA Sequences:** data/human_p53.txt, data/elephant_p53.txt_

_**2. Codon Table:** data/codon_table.json_

## Key Insights
**1. Codon-to-Amino Acid Mapping –** Successfully translates RNA into protein sequences.
**2. Conserved Proteins –** Identifies common proteins shared across species.
**3. Amino Acid Variability –** Highlights differences in amino acid composition.
**4. Visual Reports –** Graphically represents mutation effects on protein sequences.

## Biological Significance
1. Differences in amino acid frequencies may indicate functional adaptations in proteins.
2. Conserved sequences across species point to highly essential genes.
3. Variability in key amino acids could be linked to evolutionary divergence or disease mutations.

## Conclusion
This DNA-to-Protein Translation and Amino Acid Frequency Analysis Pipeline provides a bioinformatics framework for studying protein synthesis, genetic variations, and evolutionary relationships. By offering visualization and textual reports, the tool aids in molecular biology, disease research, and biotechnology applications.
