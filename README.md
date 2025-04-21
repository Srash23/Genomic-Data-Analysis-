# Protein Synthesis and Comparative Amino Acid Analysis Pipeline

This repository contains a modular pipeline designed to synthesize proteins from DNA sequences and analyze amino acid composition across different organisms. It enables ORF detection, translation using a codon table, protein synthesis, and amino acid frequency comparison. The pipeline also provides visual comparison via bar plots and identifies significantly differing amino acids.

## Project Overview
This project aims to:
- Translate DNA into proteins using a custom codon table
- Synthesize all possible proteins from a DNA sequence
- Compare synthesized proteins between two species
- Identify shared proteins and quantify amino acid frequency differences
- Highlight amino acids with significant divergence and visualize them

## Pipeline Structure

### 1. Input Parsing
- `readFile()`: Reads DNA from `.txt` file and cleans newlines
- `makeCodonDictionary()`: Converts codon JSON file to lookup dictionary

### 2. Protein Synthesis
- `dnaToRna()`: Converts DNA to RNA and segments into codons
- `generateProtein()`: Translates codons to amino acids using dictionary
- `synthesizeProteins()`: Scans full DNA for multiple ORFs and synthesizes protein chains

### 3. Comparative Analysis
- `commonProteins()`: Finds shared protein sequences
- `findAminoAcidDifferences()`: Detects amino acids with > cutoff frequency difference

### 4. Visualization
- `createChart()`: Plots bar chart of amino acid frequency per species
- `makeEdgeList()`: Highlights bars with the most significant differences

## Sample Output
- Count of total DNA bases and used/unused regions
- Number of synthesized proteins per species
- List of common proteins between two species
- Table of amino acids with largest frequency differences
- Bar plot showing comparative frequency distribution

## Input Files
| Filename | Description |
|----------|-------------|
| `human_p53.txt` | DNA sequence of human p53 gene |
| `elephant_p53.txt` | DNA sequence of elephant p53 gene |
| `codon_table.json` | Codon to amino acid mapping |

## Example Applications
- Studying amino acid evolution across species
- Identifying functionally important amino acid differences
- Teaching foundational bioinformatics and comparative genomics

## Tech Stack
- **Language**: Python
- **Libraries**: Biopython, Matplotlib, JSON

## Acknowledgements
Built as an instructional pipeline for basic bioinformatics and protein comparison workflows.
