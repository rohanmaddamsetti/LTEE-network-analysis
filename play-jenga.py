#!/usr/bin/env python

"""
play-jenga.py by Rohan Maddamsetti.

This script generates minimal metabolic networks,
using the method described in Pal et al. (2006) in Nature,
and referenced by Hosseini et al. in PNAS.

"""

import cobra
import pandas as pd
import numpy as np
from scipy import stats
from os import path
import random


def generate_LTEE_clone_cobra_models(KOed_genes_df, basic_model,
                                     using_locus_tag=False):
    LTEE_cobra_models = {}
    nonmutator_pops = ["Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5"]
    mutator_pops = ["Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6"]
    LTEE_pop_vec = nonmutator_pops + mutator_pops
    model_genes = [x.id for x in basic_model.genes]
    for LTEE_pop in LTEE_pop_vec:
        cur_KO_model = basic_model.copy()
        cur_KO_model.id = LTEE_pop + "_50K_A_clone"
        is_cur_pop = (KOed_genes_df["Population"] == LTEE_pop)
        
        if using_locus_tag:
            KOed_genes = [x for x in KOed_genes_df[is_cur_pop].locus_tag]
        else:
            KOed_genes = [x for x in KOed_genes_df[is_cur_pop].blattner]
            
        genes_to_remove = [cur_KO_model.genes.get_by_id(x)
                           for x in KOed_genes if x in model_genes]
        if genes_to_remove:
            cobra.manipulation.delete_model_genes(
                cur_KO_model, genes_to_remove)
        
        LTEE_cobra_models[LTEE_pop] = cur_KO_model
    return LTEE_cobra_models


def play_jenga1(basic_model, cutoff = 0.01):
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

def play_jenga2(basic_model, cutoff = 0.01):
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
    # the E. coli K-12 iJO1366 model: this is the best curated complete E. coli model.
    K12_model = cobra.io.load_json_model(path.join(
        BiGG_model_dir, "iJO1366.json"))
    K12_model.id = "K12"

    ## the E. coli REL606 iECB_1328 model: most relevant to LTEE,
    ## but has stochiometric inconsistencies in hydrogen/proton conservation.
    REL606_model = cobra.io.load_json_model(path.join(
        BiGG_model_dir, "iECB_1328.json"))
    REL606_model.id = "REL606"

    ## The simplest well-curated model: E. coli core metabolism.
    core_model = cobra.io.load_json_model(path.join(
        BiGG_model_dir, "e_coli_core.json"))
    core_model.id = "Ecoli-core"

    nonmutator_pops = ["Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5"]
    mutator_pops = ["Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6"]
    LTEE_pop_vec = nonmutator_pops + mutator_pops

    all_KOed_genes_in_50K_A_clones = pd.read_csv(
        "../results/metabolic-enzymes/KOed-genes-in-LTEE-50K-A-clones.csv")

    ## Generate models from K-12.
    iJO1366_LTEE_models = generate_LTEE_clone_cobra_models(
        all_KOed_genes_in_50K_A_clones, K12_model)
    ## Generate models from REL606.
    iECB_1328_LTEE_models = generate_LTEE_clone_cobra_models(
        all_KOed_genes_in_50K_A_clones, REL606_model, using_locus_tag=True)

    ## play_jenga functions work fine with the core_model.
    ##evolved_model = play_jenga1(core_model)
    ##evolved_model = play_jenga2(core_model)


    ## play_jenga1 and play_jenga2 both fail, on K12 and REL606 models, on both
    ## my laptop and the DCC.
    ## then rewrite my code to avoid whatever unknown block
    ## that is failing.
    evolved_model = play_jenga1(K12_model)

    ## play_jenga1 function fails on REL606_model
    ##evolved_model = play_jenga2(REL606_model)
    return

## run the main program
main()
