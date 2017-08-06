#!/usr/bin/python 
# Looks like this script was written by (maybe?) Alex Williams in 2013 or thereabouts
# This script extracts fields from a cleaned constitutive gtf annotation file and
# writes it to a gene table for use in USEQ 
# gtf format:
# seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
# source - name of the program that generated this feature, or the data source (database or project name)
# feature - feature type name, e.g. Gene, Variation, Similarity
# start - Start position of the feature, with sequence numbering starting at 1.
# end - End position of the feature, with sequence numbering starting at 1.
# score - A floating point value.
# strand - defined as + (forward) or - (reverse).
# frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
# example GTF line:
# chr3	protein_coding	exon	157197323	157197438	.	+	.	 gene_id "ENSMUSG00000028180"; transcript_id "ENSMUST00000106063"; exon_number "E1.2"; gene_name "Zranb2"; gene_biotype "protein_coding"; transcript_name "Zranb2-202";
# example table for USEQ
# (uniqueName1 name2(optional) chrom strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts (commaDelimited)exonEnds).
# ENSG00000183888 C1orf64 chr1 + 16203317 16207889 16203385 16205428 2 16203317,16205000 16203467,16207889
#
# 0       chr3
# 1       protein_coding
# 2       exon
# 3       157199298
# 4       157199350
# 5       .
# 6       +
# 7       .
# 8       gene_id
# 9       "ENSMUSG00000028180";
# 10      transcript_id
# 11      "ENSMUST00000106063";
# 12      exon_number
# 13      "E2.1";
# 14      gene_name
# 15      "Zranb2";
# 16      gene_biotype
# 17      "protein_coding";
# 18      transcript_name
# 19      "Zranb2-202";

import sys

def gtf_to_tab_fst_exn(inl):
	outline = []
	outline.append(inl[9])	# gene_id 								0
	outline.append(inl[15])	# gene name 							1
	outline.append(inl[0])	# chromosome 							2
	outline.append(inl[6])	# strand 								3
	outline.append(inl[3])	# txStart 								4
	outline.append(inl[4])	# txEnd 								5
	outline.append(inl[3])	# cdsstart 								6
	outline.append(inl[4])	# cdsEnd  these lst 2 won't be correct 	7
	outline.append(1)		# exonCount this is an int				8
	outline.append(inl[3])	# exonStart 							9
	outline.append(inl[4])	# exonEnd 								10
	return outline




def gtf_next_exon(inl, outl):
	#if + strand
	if inl[6] == "+":
		outl[5] = inl[4]					# txEnd
		outl[7] = inl[4]					# cdsEnd
		outl[8] = outl[8] + 1 				# add exon
		outl[9] = outl[9] + "," + inl[3] 	# add exon start
		outl[10] = outl[10] + "," + inl[4]	# add exon stop
	else: 
	# - strand
		outl[4] = inl[3]					# txStart
		outl[6] = inl[3]					# cdsStart
		outl[8] = outl[8] + 1 				# add exon	
		outl[9] = inl[3] + "," + outl[9]	# add exon start
		outl[10] = inl[4] + "," + outl[10]
	return outl




def main():
	if len(sys.argv) < 3:
		print "Usage: convert-constitutive-gft-to-gene-table.py <input.gtf> <output file>"
		sys.exit()
	outline = []
	lastline = []
	with open(sys.argv[1], 'r') as gtf:
		with open(sys.argv[2], 'w') as tab:
			for line in gtf:
				gtfsplit = line.rsplit()
				#remove unwanted characters
				for i in range(len(gtfsplit)):
					gtfsplit[i] = gtfsplit[i].translate(None, '";')
				# check for first line
				if len(outline) != 0:
					# check for new gene
					if gtfsplit[9] != outline[0]:
						outline[8] = str(outline[8])
						#write new table line
						tab.write("\t".join(outline))
						tab.write("\n")
						outline = gtf_to_tab_fst_exn(gtfsplit)
					else:
						# next exon
						outline = gtf_next_exon(gtfsplit, outline)
				else: 
					#very first line
					outline = gtf_to_tab_fst_exn(gtfsplit)
			# write last line to table file
			outline[8] = str(outline[8])
			tab.write("\t".join(outline))
			# tab.write("\n")
	print "finished converting file"



if __name__ == '__main__':
	main()
