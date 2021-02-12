#!/usr/bin/env python

# Motif Mark Script
# Nikki Szczepanski

# Import necessary modules
import argparse
import re
import cairo
import math
import random

# Argparse function
def get_args():
    """Define the arguments needed for the script."""
    parser = argparse.ArgumentParser(description='Given a fasta file of genes and a text file of motifs, outputs an vector image\
		of a gene and its motifs to scale.')
    parser.add_argument("-f", "--fasta", type=str, required=True, help='Absolute path to fasta file of genes, where introns are\
		in lower case and exons are in upper case.')
    parser.add_argument("-m", "--motifs", type=str, required=True, help='Absolute path to text file containing a list of motifs,\
		where each motif is on a new line. Accepts all IUPAC degenerate base symbols.')
	#parser.add_argument("-c", "--color", type=str, required=False, help='Optional: Absolute path to text file containing a list of\
		#motifs, where the intron color is the first line, the exon color is the second line, and the motifs are ')
    args = parser.parse_args()
    return parser.parse_args()
parse_args = get_args()

# Rename variables from argparse arguments
fasta_file = parse_args.fasta
motifs_file = parse_args.motifs

# Define the ouput image file to be created (with the same path as the input fasta file)
image_file = fasta_file.split(".fa")[0] + ".svg"

def find_motifs():
	"""Read in the motifs file to create a list of motifs named 'motif_list'."""
	motif_list = []
	with open(motifs_file, 'r') as motifs:
		for motif in motifs:
			motif = motif.strip("\n")
			motif_list.append(motif)
		# Sort the motifs from longest to shortest, which will aid in the visualization later on.
		motif_list=sorted(motif_list,key=len,reverse=True)
	return motif_list
motif_list = find_motifs()

# For each IUPAC symbol, define a regex expression.
def degenerate_bases():
    """
	Creates a dictionary in which the keys are degenerate bases and the values are the possible ATCGU base letters.
	This dictionary 'bases' will be used to create regex expressions for motifs including degenerate bases, according to IUPAC.
	"""
    bases = {
		'A': 'A',
		'C': 'C',
		'G': 'G',
		'T': '[TU]',
		'U': '[TU]',
		'W': '[ATU]',
		'S': '[CG]',
		'M': '[AC]',
		'R': '[AG]',
		'Y': '[CTU]',
		'K': '[GTU]',
		'B': '[CGTU]',
		'D': '[AGTU]',
		'H': '[ACTU]',
		'V': '[ACG]',
		'N': '[ACGTU]'}
    return bases
bases = degenerate_bases()

def convert_degenerates():
	"""
	If degenerate bases are present in a motif, create a regex expression that only uses the letters ACTGU, according to IUPAC, using the pre-defined dictionary 'bases'.
	Saves each motif:regex pair in a dictionary 'motif_convmot'. 
	Also creates a dictionary 'converted_motifs' where the key is the regex expression and the value is the length of the original motif.
	Note: Not case sensitive (i.e., can interpret both lower and upper case motifs).
	"""
	motif_convmot = {}
	converted_motifs = {}
	# Iterate through each motif in the list of motifs
	for motif in motif_list:
		regex = ''
		# Iterate through each letter in the motif
		for char in motif:
			if char.upper() in bases.keys():
				# Add each letter's corresponding regex expression to a string 
				regex += bases[char.upper()]
			else:
				print("Error: Unrecognized base in motif. Please include only IUPAC symbols.")
		# Once the string of regex expressions is complete, set this as the value in the motif:regex dictionary
		motif_convmot[motif]=regex
		# Add to the dictionary of regex:original motif lengths
		converted_motifs[regex]=len(motif)
	return motif_convmot,converted_motifs
motif_convmot,converted_motifs = convert_degenerates()

def make_gene_seq_dict():
	"""
	Creates a dictionary 'gene_seq_dict' using the input fasta file in which the keys are the header lines (gene names) and the values are the concatenated strings of the gene's entire sequence.
	"""
	with open(fasta_file, 'r') as fasta:
		seq=''
		gene=''
		gene_seq_dict = {}
		for line in fasta:
			line = line.strip('\n')
			# To ensure that an empty sequence is not set as a value (occurs for the first gene in the file):
			if seq != '':
				gene_seq_dict[gene]=seq
			# Identify headers (gene names):
			if line.startswith('>'):
				# Remove the '>' symbol
				gene=line.split('>')[1]
				# Re-set sequence variable at the start of a new gene
				seq=''
			# Concatenate lines of a sequence
			else:
				seq+=line
		# This last 'if' statement is only for the very last line of the file; 
		# otherwise the dictionary does not incorporate the last line of the last sequence.
		if fasta.readline()=='':
			gene_seq_dict[gene]=seq
	return gene_seq_dict
gene_seq_dict = make_gene_seq_dict()

def find_motifs():
	"""
	For each gene in the fasta file, searches for each motif. For every match of the motif to the gene sequence, saves the starting indices in a list.
	Returns a list 'match_list' of tuples where each tuple includes the gene name, motif regex expression, length of original motif, and list of starting indices. If there are no matches of motif to gene, tuple is not added to the list. 
	Refers to the pre-defined dictionary 'gene_seq_dict'.
	"""
	match_list=[]
	for gene,seq in gene_seq_dict.items():
		for motif,length in converted_motifs.items():
			# Regex expression that allows for overlapping matches. Saves starting index for each match. Case insensitive.
			match = [m.start(0) for m in re.finditer(rf'(?=({motif}))', seq, re.IGNORECASE)]
			# Append gene, motif, motif length and match sites to a list as long as there are motif-gene matches.
			if match != []:
				match_list.append((gene,motif,length,match))
	return match_list
match_list = find_motifs()

def find_exons():
	"""
	Searches for exons in the input fasta file (as denoted by capital letters) and returns a dictionary 'exon_coords_dict' in which the keys are the header lines (gene names) and the values are tuples of the exon's start and end indices (as a tuple) and the length of the gene's entire sequence (introns+exon).
	Refers to the pre-defined dictionary 'gene_seq_dict'.
	"""
	exon_coords_dict={}
	for gene,seq in gene_seq_dict.items():
		# Regex expression that saves the starting and ending indices for capital letter sequences (aka the exon)
		exons = [(e.span(0)) for e in re.finditer('[ATCGUN]+', seq)]
		# Add gene name: exon start and end indices, length of entire gene sequence. Only the first exon is saved.
		exon_coords_dict[gene]=(exons[0],len(seq))
	return exon_coords_dict
exon_coords_dict = find_exons()

def define_boundaries():
	"""
	Determines the dimensions of the output svg based on the number of genes in the fasta file and the length of the longest gene sequence.
	Refers to the pre-defined dictionary 'exon_coords_dict'.
	Outputs the variables 'width' and 'height'.
	"""
	num_genes = 0
	longest_gene = 0
	for gene,values in exon_coords_dict.items():
		num_genes+=1
		if values[1] > longest_gene:
			longest_gene=values[1]
	width = longest_gene + 100
	height = 100*num_genes + 200
	return width,height

width,height = define_boundaries()

def color_scheme():
	"""
	From a list of pre-defined colors, associates each motif to a color.
	Creates a dictionary in which the keys are the motifs as regex expressions and the values are the colors in RGB notation.
	"""
	# Below is the list of default colors that will represent motifs:
	                #blue       magenta     orange      teal        gray        cyan
	color_list=[(0.1,0.4,0.7),(1,0.2,0.4),(1,0.5,0),(0,0.7,0.6),(0.8,0.8,0.8),(0,0.7,1)]
	color_dict={}
	i=0
	for conv_motif in converted_motifs:
		color_dict[conv_motif]=color_list[i]
		i+=1
	return color_dict

color_dict=color_scheme()

def draw_genes():
	"""
	Outputs one svg of genes and motifs using Pycairo, where each gene is shown one after the other in the same order they appear in the fasta file. The motifs are overlaid on each gene where the sequences match. Inludes a legend that displays the motif color key.
	"""
	surface =  cairo.SVGSurface(image_file,width,height)
	context = cairo.Context(surface)
	# Set the width of the intron line
	context.set_line_width(2)
	# Set starting position of file (set margins)
	ypos=50
	indent=50
	# Move along the list of genes, drawing each gene one at a time
	for gene,values in exon_coords_dict.items():
		motif_counter=0
		exon_start=values[0][0]
		exon_length=values[0][1]-values[0][0]
		exon_height=len(motif_list)*5+5
		# Move to x,y position
		context.move_to(indent,ypos-exon_height)
		context.set_source_rgb(0,0,0)
		context.show_text(gene)
		context.move_to(indent,ypos)
		# Draw a line that is to scale to the entire gene sequence (the length of the entire gene is denoted by values[1])
		context.set_line_width(2)
		context.set_source_rgb(0,0,0)
		context.line_to(indent+values[1],ypos)
		context.stroke()
		# Key: context.rectangle(x start,y start,width,height)
		context.rectangle(indent+exon_start,ypos-exon_height/2,exon_length,exon_height)
		context.fill()
		for item in match_list:
			# Refer to one gene at a time
			if gene==item[0]:
				for conv_motif,length in converted_motifs.items():
					# Refer to one motif at a time
					if conv_motif==item[1]:
						motif_counter+=1
						box_height=exon_height-3*motif_counter
						# Set color to associated motif color
						R,G,B=color_dict[conv_motif]
						# Draw each motif match at a time
						for match in item[3]:
							context.rectangle(indent+match,ypos-box_height/2,length,box_height)
							context.set_source_rgb(R,G,B)
							context.fill_preserve()
							context.set_line_width(0.3)
							context.set_source_rgb(0,0,0)
							context.stroke()
		ypos += 100
	# Create box for legend
	box_top = ypos-50
	context.set_source_rgb(0,0,0)
	context.rectangle(indent,box_top,len(max(motif_list,key=len))+100,motif_counter*20+10)
	context.set_line_width(2)
	context.stroke()
	# Add motifs and motif colors to legend
	x=0
	for motif in motif_list:
		# Print motif name in legend
		x+=1
		context.move_to(indent+30,box_top+15*x)
		context.set_source_rgb(0,0,0)
		context.show_text(motif)
		# Associate each motif to its color value
		R,G,B=color_dict[motif_convmot[motif]]
		context.set_source_rgb(R,G,B)
		context.rectangle(indent+10,box_top+15*x-7,10,10)
		context.fill()
	# Finish drawing
	surface.finish()
	return
draw_genes()