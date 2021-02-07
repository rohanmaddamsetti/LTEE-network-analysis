#!/usr/bin/python

## fix-Cooper2003-geneIDs.py by Rohan Maddamsetti.
## This script attaches locus_tags to rows in
## ../results/renamed-Cooper2003-expression-array.csv.
## The gene names in these data are not consistent with what I like to use,
## and this script attaches locus_tag data for happy merging.

## This script prints out a list of all genes that map one-to-one between
## REL606.gbk and renamed_arrays.txt (the gene expression data).

## Usage: python expression_label_fixer.py > ../results/thermostability/reformatted-Cooper2003-expression-array.csv

from Bio import SeqIO

## From stackoverflow, for how to find duplicate elements in a list.
## Useful for finding different loci with the same gene name.
def find_duplicate_elements(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set([x for x in seq if x in seen or seen_add(x)])
    # turn the set into a list (as requested)
    return list(seen_twice)

## First, make a list of triples of locus_tags, gene names, and expression names.
## Avoid duplicate entries by hashing by locus_tag.
triples = []

for genome in SeqIO.parse("../data/REL606.8-RM.gbk", "genbank"):
    for feature in genome.features:
        if feature.type != "CDS":
            continue
        else:
            try:
                raw_note = feature.qualifiers["note"][0]
                if raw_note.startswith("b"):
                    note = raw_note[0:5] ## get the blattner ID.
                else:
                    note = ""
            except KeyError:
                note = ""
            try:
                gene = feature.qualifiers["gene"][0]
            except KeyError:
                gene = ""
            locus_tag = feature.qualifiers["locus_tag"][0]
            my_tuple = (locus_tag, gene, note)
            triples.append(my_tuple)


expression_file = open("../results/thermostability/renamed-Cooper2003-expression-array.csv", "r")

for i, line in enumerate(expression_file):
    line = line.strip()
    if i == 0: ## skip the header
        print(line + ",locus_tag")
    else:
        locus_tag = "NA"
        cooper_gene_name = line.split(",")[0]
        for triple in triples:
            if cooper_gene_name in triple:
                locus_tag = triple[0]
                break
        print(line + "," + locus_tag)
