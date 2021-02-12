# Motif Mark

Motif Mark is a Python3-based script that allows the user to visualize genes and where motifs occur on them. It creates a color-coded image (.svg) using a fasta file and text file of motif sequences.

## Requirements

Motif Mark requires two files:   
• a fasta file, where every exon is denoted by upper case letters introns by lower case letters  
• a text file of motifs, where every motif sequence is on a new line and follows IUPAC notation.  

An example fasta file and motif file is found in this repo under ex_fasta.fasta and ex_motif.txt, respectively.

## Installation

To run Motif Mark on your own computer, you will need the following Python modules:  
• argparse  
• re  
• cairo  
• math  

