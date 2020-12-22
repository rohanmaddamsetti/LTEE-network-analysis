#!/usr/bin/env python

'''
collate-noDeletions-resilience-runs.py by Rohan Maddamsetti

This script makes a big csv file from the output of running
ppi-resilience-analysis2.py on the Duke Compute Cluster
(the files in results/resilience/noDeletions-resilience-analysis-runs)

'''

from os import listdir
from os.path import join


outf = "../results/resilience/collated-noDeletions-resilience-runs.csv"
with open(outf, "w") as outfh:

    header = "dataset,run_type,replicate,strain,resilience"
    outfh.write(header + "\n")
    
    results_dir = "../results/resilience/noDeletions-resilience-analysis-runs"
    for fname in listdir(results_dir):
        ## first, parse the fname to get the dataset, type of run, and replicate.
        fname_metadata = fname.split('.')[0]
        fname_fields = fname_metadata.split('_')
        dataset = fname_fields[0]
        run_type = '_'.join(fname_fields[2:-1])
        replicate = fname_fields[-1].split('Rep')[-1]

        full_fname = join(results_dir,fname)
        with open(full_fname, 'r') as fh:
            for i,line in enumerate(fh):
                if i == 0: continue ## skip header
                line = line.strip()
                fields = line.split(',')
                strain = fields[-2]
                resilience = fields[-1]
                output_line = ','.join([dataset,run_type,replicate,strain,resilience])
                outfh.write(output_line + "\n")

