#!/usr/bin/env python

'''
ProteomeVis-to-REL606.py by Rohan Maddamsetti

This script prints out a csv file mapping the genes column in
"../results/thermostability/Ecoli-ProteomeVis-data/Ecoli-proteomevis_inspect.csv"
to the Gene column in
"../results/REL606_IDs.csv".
'''

import os

outfile = "../results/thermostability/REL606-to-ProteomeVis.csv"

REL606_file = "../results/REL606_IDs.csv"
ProteomeVis_file = "../results/thermostability/Ecoli-proteomevis_inspect.csv"

## make a list of Gene,blattner tuples using the REL606 data.
REL606_tuples = []
with open(REL606_file,"r") as REL606_fh:
    for i,line in enumerate(REL606_fh):
        if i == 0: continue ## skip the header
        ldata = line.split(',')
        REL606gene, REL606blattner = ldata[0], ldata[2]
        REL606_tuples.append((REL606gene,REL606blattner))

with open(outfile,"w") as outfh:
    ## write a header.
    ## pdb is the key that will be used for merging REL606 annotation
    ## with ProteomeVis data.
    outfh.write("pdb,Gene,blattner,genes\n")
    with open(ProteomeVis_file, "r") as proteomeVis_fh:
        for i, line in enumerate(proteomeVis_fh):
            if i == 0: continue 
            ldata = line.split(',')
            pdb_field = ldata[1]
            genes_field = ldata[3]
            found_match = False
            ## use REL606 tuples as queries
            for g,b in REL606_tuples:
                if (g in genes_field) or (b in genes_field):
                    row = ','.join([pdb_field, g, b, genes_field])
                    outfh.write(row+'\n')
                    found_match = True
                    continue
            if not found_match:
                print('MISS: ' + genes_field)

