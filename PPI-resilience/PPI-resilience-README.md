## README by Rohan Maddamsetti

## the genbank reference file REL606.8-RM.gbk was manually edited with the following
changes: 
b0092 is now ddlB rather than ddl.
b0381 is now ddlA rather than ddl.

### to make a table of the mutations in the LTEE metagenomics dataset.
conda activate ltee-metagenomix
cd LTEE-metagenomic-repo
python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv

### To install SNAP.py using the conda package manager, run the following:
conda install -c snap-stanford snap-stanford
