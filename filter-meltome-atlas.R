## filter-meltome-atlas.R by Rohan Maddamsetti.

## This script filters the cross-species dataset for E. coli data,
## and writes to file.

## downloaded from:
## http://meltomeatlas.proteomics.wzw.tum.de:5003/

## Jarzab, A., Kurzawa, N., Hopf, T. et al.
## Meltome atlas—thermal proteome stability across the tree of life.
## Nat Methods 17, 495–503 (2020).
## https://doi.org/10.1038/s41592-020-0801-4

library(tidyverse)

Ecoli.meltome <- read.csv("../data/cross-species-meltome-atlas.csv") %>%
    filter(str_detect(run_name, "Escherichia coli"))
write.csv(Ecoli.meltome,file="../results/thermostability/Ecoli-meltome.csv")
