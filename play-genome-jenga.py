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
import random
import os


## if running as part of a SLURM job array, get the task ID.
## otherwise set it to the default value '9999'.
try:
    TASK_ID = os.environ['SLURM_ARRAY_TASK_ID']
except:
    TASK_ID = "9999"


def outf_path(outdir, outf_name, taskID):
    outf = outf_name + "_Rep" + taskID + ".csv"
    outf_path = os.path.join(outdir,outf)
    return outf_path


def play_genome_jenga1(basic_model, cutoff = 0.00000003):
    """ This version completely removes genes from the model."""
    random.seed() ## seed the PRNG with the current system time.
    evolved_model = basic_model.copy()
    evolved_model.id = "jenga_evolved"
    while(True):
        cur_biomass = evolved_model.slim_optimize()
        KO_results = cobra.flux_analysis.single_gene_deletion(evolved_model)        
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
        if num_candidate_genes: ## there exist genes to delete
            idx_to_delete = random.sample(range(num_candidate_genes),1).pop()
            genes_to_delete = filtered_KO_results.ids.tolist()[idx_to_delete]
            print("deleting: ", genes_to_delete)
            cobra.manipulation.remove_genes(evolved_model, genes_to_delete)
        else: ## no more genes to delete
            break ## leave the while loop.  
    return evolved_model


def play_genome_jenga2(basic_model, cutoff = 0.00000003):
    """ This version inactivates genes in the model, but does not remove them."""
    random.seed() ## seed the PRNG with the current system time.
    evolved_model = basic_model.copy()
    evolved_model.id = "jenga_evolved"
    while(True):
        cur_biomass = evolved_model.slim_optimize()
        genes_available_to_KO = {x for x in evolved_model.genes if x.functional}
        KO_results = cobra.flux_analysis.single_gene_deletion(evolved_model,
                                                              genes_available_to_KO)
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
        if num_candidate_genes: ## there exist genes to delete
            idx_to_delete = random.sample(range(num_candidate_genes),1).pop()
            genes_to_delete = list(filtered_KO_results.ids.tolist()[idx_to_delete])
            print("deleting: ", genes_to_delete)
            cobra.manipulation.delete_model_genes(evolved_model, genes_to_delete)
        else: ## no more genes to delete
            break ## leave the while loop.  
    return evolved_model


def write_locus_tags(evolved_model, gene_outf):
    with open(gene_outf, "w") as fh:
        print("locus_tag", file=fh)
        for locus_tag in evolved_model.genes:
            if locus_tag.functional:
                print(locus_tag, file=fh)
    fh.close()
    return


def write_essential_genes(evolved_model, essential_outf):
    essential_genes = [x for x in cobra.flux_analysis.find_essential_genes(evolved_model)]
    with open(essential_outf, "w") as fh:
        print("locus_tag", file=fh)
        for locus_tag in essential_genes:
            if locus_tag.functional: ## just as an extra check.
                print(locus_tag, file=fh)
    fh.close()
    return


def write_reactions(evolved_model, reaction_outf):
    with open(reaction_outf, "w") as fh:
        print("Reaction", file=fh)
        for rxn in evolved_model.reactions:
            if rxn.functional:
                ## important: print rxn.name, not rxn.
                print(rxn.name, file=fh)
    fh.close()
    return


def main():

    ## Use Gurobi 9.5.
    cobra_config = cobra.Configuration()
    cobra_config.solver = "gurobi"
    
    BiGG_model_dir = "../data/BiGG-models"
    ## The simplest well-curated model: E. coli core metabolism.
    core_model = cobra.io.load_json_model(os.path.join(
        BiGG_model_dir, "e_coli_core.json"))
    core_model.id = "Ecoli-core"
    
    # the E. coli K-12 iJO1366 model: this is the best curated complete E. coli model.
    K12_model = cobra.io.load_json_model(os.path.join(
        BiGG_model_dir, "iJO1366.json"))
    K12_model.id = "K12"

    ## the E. coli REL606 iECB_1328 model: most relevant to LTEE,
    ## but has stochiometric inconsistencies in hydrogen/proton conservation.
    REL606_model = cobra.io.load_json_model(os.path.join(
        BiGG_model_dir, "iECB_1328.json"))
    REL606_model.id = "REL606"

    ## to simulate DM25, we need to add an excess of thiamine to the
    ## default minimal glucose-limited medium
    ## (no thiamine is in the core_model).
    K12_medium = K12_model.medium
    K12_medium['EX_thm_e'] = 1000.0
    K12_model.medium = K12_medium

    REL606_medium = REL606_model.medium
    REL606_medium['EX_thm_e'] = 1000.0
    REL606_model.medium = REL606_medium
    
    ## effective LTEE population size (see Wiser et al. 2013 Supplement)
    Neff = 33000000.0
    neutral_cutoff = 1.0/Neff
    ## play Genome Jenga!
    evolved_model = play_genome_jenga1(REL606_model, neutral_cutoff)
    ## and write the results to file.
    outdir = "../results/metabolic-enzymes/genome-jenga-runs"
    gene_outf = outf_path(outdir, "minimal_genome", TASK_ID)
    write_locus_tags(evolved_model, gene_outf)
    essential_outf = outf_path(outdir, "essential_genes", TASK_ID)
    write_essential_genes(evolved_model, essential_outf)
    rxn_outf = outf_path(outdir, "minimal_reactions", TASK_ID)
    write_reactions(evolved_model, rxn_outf)
    ## signal that the program completed successfully.
    print("SUCCESS!")
    return


## let's play Genome Jenga!
if __name__ == "__main__":
    main()
