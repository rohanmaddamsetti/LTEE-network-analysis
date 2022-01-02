#!/usr/bin/env python

''' 
filter-ProteomeVis-for-Ecoli.py by Rohan Maddamsetti. 

This script filters the data in ../data/ProteomeVis-data
(This will be just one file: pr)
for the rows that are relevant for E. coli (species == 1),
and then writes to ../results/thermostability.

'''

import os

outputdir = "../results/thermostability"
inputdir = "../data/ProteomeVis-data"
inputfiles = [x for x in os.listdir(inputdir) if x.startswith("proteomevis_inspect")]

for f in inputfiles:
    if ("inspect" in f):
        species_field = -1
    else:
        continue
    print(f)
    out_f = os.path.join(outputdir,"Ecoli-"+f)
    in_f = os.path.join(inputdir,f)
    with open(out_f, "w") as outfh:
        with open(in_f, "r") as in_fh:
            for line in in_fh:
                fields = line.split(',')
                if fields[species_field].startswith('0'): continue
                outfh.write(line)
