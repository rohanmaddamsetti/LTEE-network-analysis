#!/usr/bin/env python

'''
ppi-resilience-analysis.py by Rohan Maddamsetti.

This script does a resilience analysis of PPI networks in the LTEE genomes 
published in Tenaillon et al. (2016).

A critical assumption of this data analysis is that the level of analysis
is focused on sets of genes that are KO'ed-- and not sets of mutations.

I can either choose to keep the number of KO'ed genes fixed, or the number of
mutations fixed. I choose the former.

What are the consequences of this assumption? For one, multi-gene deletions
break consecutive genes. When sampling genes without replacement to be KO'ed
in randomized genomes, it is unlikely to sample a block of KO'ed genes.
This choice for how to randomize genomes may affect the results, 
by breaking up the 'block' structure of multi-gene deletions.
'''

## check if we are on the Duke Compute Cluster.
## if so, then import the path to SWIG so that I can import snap.
from os.path import join, exists
import sys
dcc_swig_path = "/opt/apps/rhel7/swig-4.0.2/swig"
if exists(dcc_swig_path): sys.path.append(dcc_swig_path)

## if running as part of a SLURM job array, get the task ID.
## otherwise set it to the default value '999'.
import os
try:
    taskID = os.environ['SLURM_ARRAY_TASK_ID']
except:
    taskID = '999'

import snap
from math import log
import random
import numpy as np
from scipy.integrate import simps
import pandas as pd
import argparse
from pprint import pprint
from itertools import groupby


def outf_path(outdir, outf_name, taskID):
    outf = outf_name + "_Rep" + taskID + ".csv"
    outf_path = outdir + outf
    return outf_path


def get_REL606_column_set(col):
    REL606_ID_file = "../results/REL606_IDs.csv"
    column_set = set()
    with open(REL606_ID_file, "r") as REL606_fh:
        for i,l in enumerate(REL606_fh):
            if i == 0: continue
            ldata = l.split(',')
            column_set.add(ldata[col])
    return column_set


def get_REL606_gene_set():
    ''' return a set of all gene names in REL606.'''
    return get_REL606_column_set(0)


def get_REL606_column_dict(key_col,val_col):
    REL606_ID_file = "../results/REL606_IDs.csv"
    col_dict = {}
    with open(REL606_ID_file, "r") as REL606_fh:
        for i,l in enumerate(REL606_fh):
            if i == 0: continue
            ldata = l.split(',')
            col_dict[ldata[key_col]] = ldata[val_col]
    return col_dict


def get_REL606_blattner_to_gene_dict():
    return get_REL606_column_dict(2,0)


def get_LTEE_genome_knockout_muts(no_deletions=False):
    '''
    I downloaded this table from Jeff's Shiny web app interface
    to the 264 LTEE genomes dataset. 
    knockout mutations are: nonsense SNPs, small indels, mobile element insertions, and large deletions.
    250 out of 264 genomes have knockout mutations.
    '''
    LTEE_nonsense_indel_MOB_deletions_in_genomes_f = "../data/LTEE-264-genomes-SNP-nonsense-small-indel-MOB-large-deletions.csv"
    LTEE_mut_df = pd.read_csv(LTEE_nonsense_indel_MOB_deletions_in_genomes_f)
    ## filter out all intergenic mutations.
    LTEE_knockout_muts = LTEE_mut_df[~LTEE_mut_df['gene_position'].str.contains("intergenic",na=False)]
    if no_deletions:
        ## then filter out large deletions.
        LTEE_knockout_muts = LTEE_knockout_muts[~LTEE_knockout_muts['mutation_category'].str.contains("large_deletion")]
    return LTEE_knockout_muts


def make_LTEE_strain_to_pop_dict(LTEE_knockout_muts):
    ''' return a dict of LTEE strain to the population it comes from.'''
    LTEE_strain_to_pop = dict()
    LTEE_genome_metadata = LTEE_knockout_muts[['treatment','population','time','strain','clone','mutator_status']].drop_duplicates()
    LTEE_strains = list(LTEE_genome_metadata['strain'])
    matching_pops = list(LTEE_genome_metadata['population'])
    LTEE_strain_to_pop = {x:y for x,y in zip(LTEE_strains, matching_pops)}
    return LTEE_strain_to_pop


def LTEE_strain_to_KO_genes(LTEE_knockout_muts):
    ''' return a dict of LTEE strain to set of knocked out genes in that strain. '''
    LTEE_strain_KO_dict = {}
    ## 250 out of 264 genomes have knockout mutations.
    LTEE_genome_metadata = LTEE_knockout_muts[['treatment','population','time','strain','clone','mutator_status']].drop_duplicates()
    LTEE_strains = list(LTEE_genome_metadata['strain'])
    for clone in LTEE_strains:
        clone_knockout_muts = LTEE_knockout_muts[LTEE_knockout_muts['strain']==clone]
        ## remove square brackets from each gene string in the column, split into a list of entries (which are lists),
        KO_list_of_lists = [x.replace('[', '').replace(']','').split(',') for x in clone_knockout_muts['gene_list']]
        ## and then flatten the list of lists into a set of knocked out genes.
        knocked_out_genes = {item for sublist in KO_list_of_lists for item in sublist}
        LTEE_strain_KO_dict[clone] = knocked_out_genes
    return LTEE_strain_KO_dict


def make_LTEE_pop_to_KO_dict(LTEE_strain_to_KO, LTEE_strain_to_pop):
    ''' 
    return a dict of LTEE pops to the union of KO mutations found in clones
    isolated in that population. The values are a set because of the phylogenetic
    correlations among clones over time (so that identical by descent mutations
    are not double-counted).
    '''
    LTEE_pop_to_KO_dict = {pop:set() for pop in LTEE_strain_to_pop.values()}
    ## add the KO'ed genes in the genomes to the dictionary.
    for clone, pop in LTEE_strain_to_pop.items():
        clone_KOset = LTEE_strain_to_KO.get(clone)
        ## if the strain doesn't have any KO mutations, then it has an empty set.
        if clone_KOset is None: clone_KOset = set()
        ## update the set with the KO'ed genes in the genomes.
        LTEE_pop_to_KO_dict[pop].update(clone_KOset)
    return LTEE_pop_to_KO_dict


def create_graph_and_dicts(edgeFile, col1=0, col2=1, sep=None, nodeSet=None):
    ''' 
    inputs:
    edgeFile: a file containing a list of edges.
    For example, see interactions in Zitnik dataset.

    nodeSet: a set of gene names to filter on,
    OR a dict of blattner IDs to gene names.
    Use the dict for the Zitnik data, and the
    set for the Cong data.

    output: a tuple containing the graph,
    a dictionary of genes to nodes,
    and a dictionary of nodes to genes.
    '''
    G = snap.TUNGraph.New()
    gene_to_node = {}
    node_to_gene = {}
    cur_unused_node_int = 0
    
    with open(edgeFile, "r") as in_fh:
        for l in in_fh:
            l = l.strip()
            if sep is None: ## whitespace
                ldata = l.split()
            else: ## for now, sep must be a comma for csv files.
                assert sep == ','
                ldata = l.split(sep)
            c1 = ldata[col1]
            c2 = ldata[col2]

            ## filter on nodeSet (or keys in nodeSet if a dict)
            if nodeSet is not None:
                if (c1 not in nodeSet) or (c2 not in nodeSet):
                        continue

            if type(nodeSet) == type(dict()):
                ## then we should be dealing with a dict from
                ## blattner to gene name.
                p1 = nodeSet[c1]
                p2 = nodeSet[c2]
            elif type(nodeSet) == type(set()):
                p1 = c1
                p2 = c2
            else:
                raise AssertionError("variable nodeSet must be a set or dict")
            
            if p1 not in gene_to_node:
                p1_int = cur_unused_node_int
                G.AddNode(p1_int)
                gene_to_node[p1] = p1_int
                node_to_gene[p1_int] = p1
                cur_unused_node_int += 1
            else:
                p1_int = gene_to_node[p1]
                assert node_to_gene[p1_int] == p1

            if p2 not in gene_to_node:
                p2_int = cur_unused_node_int
                G.AddNode(p2_int)
                gene_to_node[p2] = p2_int
                node_to_gene[p2_int] = p2
                cur_unused_node_int += 1
            else:
                p2_int = gene_to_node[p2]
                assert node_to_gene[p2_int] == p2
            G.AddEdge(p1_int,p2_int)
    return (G, gene_to_node, node_to_gene)


def GraphComponentDistributionDict(G):
    ''' return a dict of cardinality : number of strongly-connected components 
    with that size in the graph G.'''
    ComponentDist = snap.TIntPrV()
    snap.GetSccSzCnt(G, ComponentDist)
    return {comp.GetVal1():comp.GetVal2() for comp in ComponentDist}


def ComponentDistributionEntropy(component_dict, N):
    ''' Inputs: 
    component_dict: a dictionary of component size to number of components 
    with that size in the graph.
    N: the number of nodes in the original graph.
     Output: (normalized) graph component entropy. See Zitnik et al. (2019) 
    supplement for details.
    '''
    
    ## assert that the num of nodes in component_dict is consistent with N.
    node_num_check = sum([k*v for k,v in component_dict.items()])
    assert node_num_check == N
    
    ''' calculate the p*log(p) entries of the summation, where p = c/N,
    such that c is the number of nodes in the component,
    and p is the fraction of nodes in the component.
    since multiple components have the same size, multiply by that number.
    '''
    H_components = [c_multiplicity * (c/N) * log(c/N) for c, c_multiplicity in component_dict.items()]
    ''' entropy is defined as H = -sum(p*log(p)).
     normalize by log(N) per Zitnik et al. (2019). '''
    H = -sum(H_components)/log(N)
    return H


def GraphResilience(G):
    ''' calculate the resilience of the graph G, 
    using method in Zitnik et al. (2019), described in the supplementary methods.'''

    nodes = G.GetNodes()
    ## get distribution of strongly connected component sizes for the starting graph.
    sscdict = GraphComponentDistributionDict(G)
    ## calculate the entropy of the set of strongly connected components.
    H_0 = ComponentDistributionEntropy(sscdict, nodes)
    ## initialize a dictionary from failure rate to entropy.
    failure_rate_to_entropy = {0 : H_0}

    ## copy G by making a subgraph G2 with the same nodes.
    NIdV = snap.TIntV()
    for n in G.Nodes():
        NId = n.GetId()
        NIdV.Add(NId)
    G2 = snap.GetSubGraph(G, NIdV)
    
    ## iteratively remove edges from random nodes to fragment G2.
    ## calculate H for a range of failure rates.
    ## adaptation of C++ code in Zitnik paper. See:
    ## https://github.com/mims-harvard/life-tree/blob/master/compute-net-stats/analyze.cpp
    ## make sure that we're always working on the G2 object in memory.
    num_deleted = 0
    node_order = [n.GetId() for n in G.Nodes()]
    ## shuffle the order of the nodes. this is done in-place.
    random.shuffle(node_order)
    for fail_rate_p in range(1,100+1): ## failure rate as a percentage from 1 to 100.
        failed_frac = fail_rate_p / 100 ## fraction of failed nodes, ranging from 0 to 1.
        cur_deleted_total = nodes * failed_frac
        while (num_deleted < cur_deleted_total):
            NId = node_order[num_deleted] ## get the next random node
            num_deleted = num_deleted + 1
            G2.DelNode(NId) ## delete it to remove its edges.
            G2.AddNode(NId) ## add it back so that its edges are gone but the node remains.
        ## add entry to failure rate : component entropy dict.
        cur_ssc_dict = GraphComponentDistributionDict(G2)
        failure_rate_to_entropy[fail_rate_p] = ComponentDistributionEntropy(cur_ssc_dict, nodes)
    ## calculate resilience for the graph.
    ## This is 1 - AUC of the interpolated function.
    ## Use Simpson's rule to approximate the integral.
    x = np.array([i/100 for i in range(0,100+1)])
    y = np.array([j for j in failure_rate_to_entropy.values()])
    AUC = simps(y, x)
    resilience = 1 - AUC
    return resilience


def resilience_df(LTEE_strain_to_KO_dict, G, g_to_node, KO_list=None, reps=100):
    '''
    Input:
    a dict of LTEE strains to set of knocked out genes,
    a starting graph G representing nodes and edges in the LTEE ancestor REL606, 
    a dictionary of genes to nodes,
    an optional list of genes, possibly containing repeated entries, to draw KOs from.

    1) Generate a graph for each genome: take the starting graph,
    draw X KO'ed genes from the list,
       and generate a subgraph that contains nodes for all other genes.
    2) calculate network resilience for each genome.
    3) return a DataFrame with the timepoint (Generations), population,
    name, and resilience statistic.
    '''
    REL606_resilience = sum((GraphResilience(G) for x in range(reps)))/float(reps)

    if (KO_list is not None):
        ''' 
        then sort the KO_list, and generate probabilities for each
        item in the list, based on its frequency. See the solution here:
        https://stackoverflow.com/questions/2161752/how-to-count-the-frequency-of-the-elements-in-an-unordered-list
        And the solution here:
        https://medium.com/ibm-watson/incredibly-fast-random-sampling-in-python-baf154bd836a
        '''
        KO_list.sort()
        KO_count_tuples = [(key, len(list(group))) for key, group in groupby(KO_list)]
        KO_elts = [x for (x,_) in KO_count_tuples]
        KO_frequencies = [float(c)/len(KO_list) for (_,c) in KO_count_tuples]    

    clone_col = []
    resilience = []
    for clone, actual_KO_set in LTEE_strain_to_KO_dict.items():
        clone_col.append(clone) ## to ensure that clone and resilience values match up.
        KOsamplesize = len(actual_KO_set)
        if (KO_list is not None):
            ''' then sample KO'ed genes from KO_list to make randomized genomes 
            (preserving the same number of KOs as in the actual genomes).'''
            knocked_out_genes = np.random.choice(KO_elts, KOsamplesize,
                                                 p=KO_frequencies,replace=False)
        elif (KO_list is None):
            ## then calculate the network resilience for the actual LTEE genomes.
            knocked_out_genes = actual_KO_set
        else:
            raise AssertionError("Error: resilience_df() called incorrectly.")
        ''' now get the nodes corresponding to the KO'ed genes to filter them
            from the REL606 graph. '''
        knocked_out_nodes = [g_to_node[x] for x in knocked_out_genes if x in g_to_node]
        ## make a subgraph G2 that omits KO'ed nodes.
        NIdV = snap.TIntV()
        for n in G.Nodes():
            NId = n.GetId()
            if NId not in knocked_out_nodes:
                NIdV.Add(NId)
        G2 = snap.GetSubGraph(G, NIdV)
        
        ## calculate this clone's resilience. default is 100 replicates.
        ## save memory by using a generator comprehension.
        my_resilience = sum((GraphResilience(G2) for x in range(reps)))/float(reps)
        print(my_resilience)
        resilience.append(my_resilience)
    strain_col = ['REL606'] + clone_col
    resilience_col = [REL606_resilience] + resilience
    resilience_results = pd.DataFrame.from_dict({'strain': strain_col, 'resilience': resilience_col})
    return resilience_results


def main():

    parser = argparse.ArgumentParser(description='Provide integer for analysis to run.')
    parser.add_argument('--dataset',type=str)
    parser.add_argument('--analysis',type=int)
    parser.add_argument('--noDeletions', type=bool, default=False)
    args = parser.parse_args()
    assert args.dataset in ["zitnik", "cong"] ## only allowed values.
    
    random.seed() ## seed the random number generator.

    if args.noDeletions == True:
        outdir = "../results/resilience/noDeletions-resilience-analysis-runs/"
    else:
        outdir = "../results/resilience/resilience-analysis-runs/"

    REL606_genes = get_REL606_gene_set()
    LTEE_knockouts_df = get_LTEE_genome_knockout_muts(args.noDeletions)
    LTEE_strain_to_knockouts = LTEE_strain_to_KO_genes(LTEE_knockouts_df)
    LTEE_strain_to_pop =  make_LTEE_strain_to_pop_dict(LTEE_knockouts_df)
    LTEE_pop_to_knockouts = make_LTEE_pop_to_KO_dict(LTEE_strain_to_knockouts, LTEE_strain_to_pop)
    LTEE_genomes_KO_gene_list = [KO for KOset in LTEE_pop_to_knockouts.values() for KO in KOset]    
    all_LTEE_KO_genes = list(set(LTEE_genomes_KO_gene_list))
        
    
    if args.dataset == "zitnik":
        ''' Import protein-protein interaction network from 
        Marinka Zitnik paper in PNAS. '''
        K12_edge_file = "../results/thermostability/Ecoli-Zitnik-data/Ecoli-treeoflife.interactomes/511145.txt"
        ## filter K-12 graph based on blattner IDs found in REL606.
        blattner_to_gene_dict = get_REL606_blattner_to_gene_dict()
        G, g_to_node, node_to_g = create_graph_and_dicts(K12_edge_file,
                                                         nodeSet=blattner_to_gene_dict)
    else:
        ''' Import protein-protein interaction network from Cong et al. (2019).
        Take 2683 interactions from the Cong dataset. these are high-quality 
        coevolution predictions, plus coevolution interactions that are known,
        plus gold standard interactions in PDB and Ecocyc. the information in
        this file is produced by calc-ppi-network-statistics.py. '''
        good_edge_f = "../results/thermostability/Cong-good-interaction-set.tsv"
        ## filter Cong PPI graph based on genes in REL606.
        ## 1191 nodes, 1787 edges in this graph.
        G, g_to_node, node_to_g = create_graph_and_dicts(good_edge_f,
                                                         nodeSet=REL606_genes)
    if args.analysis == 1: 
        ## analyze the resilience of evolved LTEE genomes.
        if args.dataset == "zitnik":
            outf = outf_path(outdir, "Zitnik_PPI_LTEE_genome_resilience", taskID)
        else:
            outf = outf_path(outdir, "Cong_PPI_LTEE_genome_resilience", taskID)
        results = resilience_df(LTEE_strain_to_knockouts, G, g_to_node)
    elif args.analysis == 2:
        ## calculate randomized resilience of genomes, set of all genes in REL606.
        if args.dataset == "zitnik":
            outf = outf_path(outdir, "Zitnik_PPI_all_genes_randomized_resilience", taskID)
        else:
            outf = outf_path(outdir, "Cong_PPI_all_genes_randomized_resilience", taskID)
        results = resilience_df(LTEE_strain_to_knockouts, G, g_to_node,
                                KO_list=list(REL606_genes))
    elif args.analysis == 3:
        ''' calculate randomized resilience of genomes,
        using the list of all genes KO'ed in the LTEE genomics data.
        This preferential samples genes based on the number of times 
        that KO mutations in that gene was observed across populations.
        '''
        if args.dataset == 'zitnik':
            outf = outf_path(outdir, "Zitnik_PPI_across_pops_weighted_randomized_resilience", taskID)
        else:
            outf = outf_path(outdir, "Cong_PPI_across_pops_weighted_randomized_resilience", taskID)
        results = resilience_df(LTEE_strain_to_knockouts, G, g_to_node,
                                KO_list=LTEE_genomes_KO_gene_list)
    elif args.analysis == 4:
        ''' calculate randomized resilience of genomes,
        using the set of all genes KO'ed in the LTEE genomics data. '''
        if args.dataset == "zitnik":
            outf = outf_path(outdir, "Zitnik_PPI_across_pops_randomized_resilience", taskID)
        else:
            outf = outf_path(outdir, "Cong_PPI_across_pops_randomized_resilience", taskID)
        results = resilience_df(LTEE_strain_to_knockouts, G, g_to_node,
                                KO_list=all_LTEE_KO_genes)
    else:
        raise AssertionError("analysis parameter ranges from 1-4.")
    ## write results to file.
    results.to_csv(outf)

if __name__ == "__main__":
    main()
