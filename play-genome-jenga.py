#!/usr/bin/env python

"""
play-genome-jenga.py by Rohan Maddamsetti.

This script generates a minimal metabolic network,
using the method described in Pal et al. (2006) in Nature,
and referenced by Hosseini et al. in PNAS.

"""

import cobra
import pandas as pd
import numpy as np
from scipy import stats
from os import path
import random


def play_genome_jenga1(basic_model, cutoff = 0.01):
    """

    This version completely removes genes from the model.

    While there exist genes that have a beneficial effect when knocked out:
    - randomly order the genes to remove.
    - generate a KO. If deleterious, then reject the change and continue.
      If the change is beneficial, then accept and update the model.
    - once the entire DFE of possible KOs is deleterious,
    return the evolved model.
    """
    random.seed() ## seed the PRNG with the current system time.
    evolved_model = basic_model.copy()
    evolved_model.id = "jenga_evolved"
    
    steps = 0 ## number of genes KO'ed so far
    while(True):
        cur_biomass = evolved_model.slim_optimize()
        KO_results = cobra.flux_analysis.single_gene_deletion(evolved_model)
        num_genes = len(KO_results.ids)
        print("num_genes", num_genes)
        
        """ From the Methods of Pal et al. (2006) in Nature:
        A gene was classified as having no fitness effect 
        if the biomass production rate of the knockout strain 
        was reduced by less than a given cutoff; different cutoffs
        led to very similar results.
        
        filter potential knockouts based on this criterion.
        Pal et al. used cutoff = 0.01 and cutoff = 0.1."""

        filtered_KO_results = KO_results[KO_results.growth >
                                         (cur_biomass - cutoff)]
        num_candidate_genes = len(filtered_KO_results.ids)
        fract_purifying = float(num_genes - num_candidate_genes)/float(num_genes)
        print("num_candidate_genes", num_candidate_genes)
        print("fraction genes under purifying selection: ", fract_purifying)
        if num_candidate_genes: ## there exist genes to delete
            idx_to_delete = random.sample(range(num_candidate_genes),1).pop()
            genes_to_delete = filtered_KO_results.ids.tolist()[idx_to_delete]
            print("deleting: ", genes_to_delete)
            steps += 1
            print("steps <= " + str(steps))
            cobra.manipulation.remove_genes(evolved_model, genes_to_delete)
        else: ## no more genes to delete
            break ## leave the while loop.  
    print("steps <= " + str(steps))
    return evolved_model


def play_genome_jenga2(basic_model, cutoff = 0.01):
    """

    This version inactivates genes in the model, but does not
    turn them off.

    While there exist genes that have a beneficial effect when knocked out:
    - randomly order the genes to remove.
    - generate a KO. If deleterious, then reject the change and continue.
      If the change is beneficial, then accept and update the model.
    - once the entire DFE of possible KOs is deleterious,
    return the evolved model.
    """
    random.seed() ## seed the PRNG with the current system time.
    evolved_model = basic_model.copy()
    evolved_model.id = "jenga_evolved"
    
    steps = 0 ## number of genes KO'ed so far
    while(True):
        cur_biomass = evolved_model.slim_optimize()
        genes_available_to_KO = {x for x in evolved_model.genes if x.functional}
        KO_results = cobra.flux_analysis.single_gene_deletion(evolved_model,
                                                              genes_available_to_KO)
        num_genes = len(genes_available_to_KO)
        print("num_genes", num_genes)
        
        """ From the Methods of Pal et al. (2006) in Nature:
        A gene was classified as having no fitness effect 
        if the biomass production rate of the knockout strain 
        was reduced by less than a given cutoff; different cutoffs
        led to very similar results.
        
        filter potential knockouts based on this criterion.
        Pal et al. used cutoff = 0.01 and cutoff = 0.1."""

        filtered_KO_results = KO_results[KO_results.growth >
                                         (cur_biomass - cutoff)]
        num_candidate_genes = len(filtered_KO_results.ids)
        fract_purifying = float(num_genes - num_candidate_genes)/float(num_genes)
        print("num_candidate_genes", num_candidate_genes)
        print("fraction genes under purifying selection: ", fract_purifying)
        if num_candidate_genes: ## there exist genes to delete
            idx_to_delete = random.sample(range(num_candidate_genes),1).pop()
            genes_to_delete = list(filtered_KO_results.ids.tolist()[idx_to_delete])
            print("deleting: ", genes_to_delete)
            steps += 1
            print("steps <= " + str(steps))
            cobra.manipulation.delete_model_genes(evolved_model, genes_to_delete)
        else: ## no more genes to delete
            break ## leave the while loop.  
    print("steps <= " + str(steps))
    return evolved_model


def main():

    ## use Gurobi. It makes a big difference!
    cobra_config = cobra.Configuration()
    cobra_config.solver = "gurobi"
    print(cobra_config)
    
    BiGG_model_dir = "../data/BiGG-models"

    ## The simplest well-curated model: E. coli core metabolism.
    core_model = cobra.io.load_json_model(path.join(
        BiGG_model_dir, "e_coli_core.json"))
    core_model.id = "Ecoli-core"
    
    # the E. coli K-12 iJO1366 model: this is the best curated complete E. coli model.
    K12_model = cobra.io.load_json_model(path.join(
        BiGG_model_dir, "iJO1366.json"))
    K12_model.id = "K12"

    ## the E. coli REL606 iECB_1328 model: most relevant to LTEE,
    ## but has stochiometric inconsistencies in hydrogen/proton conservation.
    REL606_model = cobra.io.load_json_model(path.join(
        BiGG_model_dir, "iECB_1328.json"))
    REL606_model.id = "REL606"

    evolved_model = play_jenga1(REL606_model)
    return


## play genome Jenga!
main()
