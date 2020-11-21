#!/usr/bin/env python

'''
snap-ppi-resilience-analysis.py by Rohan Maddamsetti.

This script does a resilience analysis of PPI networks in the LTEE genomes 
published in Tenaillon et al. (2016).
'''

import snap
from math import log
from os.path import join
import random
import numpy as np
from scipy.integrate import simps
import pandas as pd
import gc


def getREL606_column_set(col):
    REL606_ID_file = "../results/REL606_IDs.csv"
    column_set = set()
    with open(REL606_ID_file, "r") as REL606_fh:
        for i,l in enumerate(REL606_fh):
            if i == 0: continue
            ldata = l.split(',')
            column_set.add(ldata[col])
    return column_set


def getREL606_gene_set():
    ''' return a set of all gene names in REL606.'''
    return getREL606_column_set(0)


def getREL606_column_dict(key_col,val_col):
    REL606_ID_file = "../results/REL606_IDs.csv"
    col_dict = {}
    with open(REL606_ID_file, "r") as REL606_fh:
        for i,l in enumerate(REL606_fh):
            if i == 0: continue
            ldata = l.split(',')
            col_dict[ldata[key_col]] = ldata[val_col]
    return col_dict


def getREL606_blattner_to_gene_dict():
    return getREL606_column_dict(2,0)

def CopyGraph(G):
    ## copy graph object into G2 so that we don't mess up G as a side-effect.
    G2 = snap.TUNGraph.New(G.GetNodes(),G.GetEdges())
    ## add nodes.
    for n in G.Nodes():
        G2.AddNode(n.GetId())
    ## add edges.
    for e in G.Edges():
        G2.AddEdge(e.GetSrcNId(),e.GetDstNId())
    return G2

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
    graph_and_dictionaries = (G, gene_to_node, node_to_gene)
    return(graph_and_dictionaries)


def GraphComponentDistributionDict(G):
    ''' return a dict of cardinality : number of strongly-connected components 
    with that size in the graph G.'''
    ComponentDist = snap.TIntPrV()
    snap.GetSccSzCnt(G, ComponentDist)
    component_size_to_count = {comp.GetVal1():comp.GetVal2() for comp in ComponentDist}
    return component_size_to_count


def ComponentDistributionEntropy(component_dict, N):
    ''' Inputs: 
    component_dict: a dictionary of component size to number of components 
    with that size in the graph.
    N: the number of nodes in the original graph.
     Output: (normalized) graph component entropy. See Zitnik et al. (2019) 
    supplement for details.
    '''
    
    ## assert that the number of nodes in component_dict is consistent with N.
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
    return(H)


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
    
    ## copy graph object into G2 so that we don't mess up G as a side-effect.
    G2 = CopyGraph(G)
    
    ## seed random number generator.
    random.seed()
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
            G2.AddNode(NId) ## add it back, so that its edges are gone but the node remains.
        ## add entry to failure rate : component entropy dict.
        cur_ssc_dict = GraphComponentDistributionDict(G2)
        failure_rate_to_entropy[fail_rate_p] = ComponentDistributionEntropy(cur_ssc_dict, nodes)

    del G2 ## explicitly delete G2 from memory now that we don't need it.
        
    ## calculate resilience for the graph.
    ## This is 1 - AUC of the interpolated function.
    ## Use Simpson's rule to approximate the integral.
    x = np.array([i/100 for i in range(0,100+1)])
    y = np.array([j for j in failure_rate_to_entropy.values()])
    AUC = simps(y, x)
    resilience = 1 - AUC
    return resilience

def get_LTEE_genome_knockout_muts():
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
    return LTEE_knockout_muts

def LTEE_strain_to_KO_genes(LTEE_knockout_muts):
    ''' return a dict of LTEE strain to set of knocked out genes in that strain.'''
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
    

def resilience_analysis_of_LTEE_genomes(LTEE_strain_to_KO_dict, G, g_to_node, reps=100):
    '''
    Input: 
    a dict of LTEE strain to set of knocked out genes,
    a starting PPI graph G, 
    and a dictionary of genes to nodes.

    1) Import KO mutations in the genomes in the Tenaillon dataset.
    2) Generate a graph for each genome: take the starting PPI graph,
       and delete nodes/genes hit by nonsense SNPs, indels, mobile elements.
    3) calculate PPI network resilience for each genome.
    4) return a pandas DataFrame of the PPI network resilience of each strain.
    '''

    REL606_resilience = np.mean([GraphResilience(G) for x in range(reps)])

    clone_col = []
    LTEE_strain_resilience = []
    for clone, knocked_out_genes in LTEE_strain_to_KO_dict.items():
        clone_col.append(clone) ## to ensure that clone and resilience match up.
        ## get the nodes to remove from the REL606 graph.
        knocked_out_nodes = []
        for x in knocked_out_genes:
            if x in g_to_node:
                knocked_out_nodes.append(g_to_node[x])
        ## copy starting graph G into G2 so that we don't mess up G as a side-effect.
        ## when we do so, exclude the genes affected by knockout mutations.
        G2 = snap.TUNGraph.New(G.GetNodes(),G.GetEdges())
        ## add nodes.
        for n in G.Nodes():
            nId = n.GetId()
            if nId not in knocked_out_nodes:
                G2.AddNode(nId)
        ## add edges.
        for e in G.Edges():
            srcN = e.GetSrcNId()
            dstN = e.GetDstNId()
            if (srcN not in knocked_out_nodes) and (dstN not in knocked_out_nodes):
                G2.AddEdge(srcN, dstN)
        ## calculate this clone's resilience. default is 100 replicates.
        ## save memory by using a generator comprehension.
        my_resilience = sum((GraphResilience(G2) for x in range(reps)))/float(reps)
        print(my_resilience)
        LTEE_strain_resilience.append(my_resilience)
    strain_col = ['REL606'] + clone_col
    resilience_col = [REL606_resilience] + LTEE_strain_resilience
    resilience_results = pd.DataFrame.from_dict({'strain': strain_col, 'resilience': resilience_col})
    return resilience_results


def resilience_randomized_over_gene_set(LTEE_strain_to_KO_dict, G, g_to_node, KO_set, reps=100):
    '''
    Input:
    a dict of LTEE strains to set of knocked out genes,
    a starting PPI graph G, 
    a dictionary of genes to nodes,
    and a set of genes to draw KOs from.

    This set of genes will be:
    A) all genes in REL606.
    B) KOs in the entire LTEE, from the LTEE metagenomics data.

    1) Generate a graph for each genome: take the starting PPI graph,
    draw X KO'ed genes from the bag,
       and delete nodes/genes from the graph.
    2) calculate PPI network resilience for each genome.
    3) return a DataFrame with the timepoint (Generations), population,
    name, and resilience statistic.
    '''
    REL606_resilience = np.mean([GraphResilience(G) for x in range(reps)])
    ## seed random number generator.
    random.seed()
    clone_col = []
    randomized_resilience = []
    for clone, gset in LTEE_strain_to_KO_dict.items():
        clone_col.append(clone) ## to ensure that clone and randomized resilience match up.
        KOsamplesize = len(gset)
        ## sample the genes to remove from the REL606 graph.
        knocked_out_genes = random.sample(KO_set, KOsamplesize)
        knocked_out_nodes = []
        for x in knocked_out_genes:
            if x in g_to_node:
                knocked_out_nodes.append(g_to_node[x])
        ## copy starting graph G into G2 so that we don't mess up G as a side-effect.
        ## when we do so, exclude the genes affected by knockout mutations.
        G2 = snap.TUNGraph.New(G.GetNodes(),G.GetEdges())
        ## add nodes.
        for n in G.Nodes():
            nId = n.GetId()
            if nId not in knocked_out_nodes:
                G2.AddNode(nId)
        ## add edges.
        for e in G.Edges():
            srcN = e.GetSrcNId()
            dstN = e.GetDstNId()
            if (srcN not in knocked_out_nodes) and (dstN not in knocked_out_nodes):
                G2.AddEdge(srcN, dstN)
        ## calculate this clone's resilience. default is 100 replicates.
        ## save memory by using a generator comprehension.
        my_resilience = sum((GraphResilience(G2) for x in range(reps)))/float(reps)
        print(my_resilience)
        randomized_resilience.append(my_resilience)
    strain_col = ['REL606'] + clone_col
    resilience_col = [REL606_resilience] + randomized_resilience
    resilience_results = pd.DataFrame.from_dict({'strain': strain_col, 'resilience': resilience_col})
    return resilience_results

def get_LTEE_metagenomics_knockouts():
    ''' 
    Import csv of LTEE metagenomics data, and
    return a list of genes affected by knockout mutations.
    knockout mutations are: nonsense SNPs, small indels, mobile element insertions, and large deletions.
'''
    LTEE_metagenomics_df = pd.read_csv("../results/LTEE-metagenome-mutations.csv")
    KO_mut_classes = ["nonsense","indel","sv"]
    KO_df = LTEE_metagenomics_df[LTEE_metagenomics_df.Annotation.isin(KO_mut_classes)]
    KO_genes = {x for x in KO_df['Gene']} - {"intergenic"} ## intergenic is not a Gene.
    return KO_genes

def resilience_randomized_within_LTEE_pops(LTEE_strain_to_KO_dict, LTEE_strain_to_pop_dict, G, g_to_node, KO_dict, reps=100):
    '''
    Input:
    a dict of LTEE strains to set of knocked out genes,
    a dict of LTEE strains to their source populations,
    a starting PPI graph G, 
    a dictionary of genes to nodes,
    and a dict from LTEE population to the set of genes
    with KOs in that population, from the LTEE metagenomics data

    1) Generate a graph for each genome: take the starting PPI graph,
    draw X KO'ed genes from the bag,
       and delete nodes/genes from the graph.
    2) calculate PPI network resilience for each genome.
    3) return a DataFrame with the timepoint (Generations), population,
    name, and resilience statistic.
    '''
    REL606_resilience = np.mean([GraphResilience(G) for x in range(reps)])
    ## seed random number generator.
    random.seed()
    clone_col = []
    randomized_resilience = []
    for clone, gset in LTEE_strain_to_KO_dict.items():
        clone_col.append(clone) ## to ensure that clone and randomized resilience match up.
        KOsamplesize = len(gset)
        ## sample the genes to remove from the REL606 graph.
        pop = LTEE_strain_to_pop_dict[clone]
        KO_set = KO_dict[pop]
        knocked_out_genes = random.sample(KO_set, KOsamplesize)
        knocked_out_nodes = []
        for x in knocked_out_genes:
            if x in g_to_node:
                knocked_out_nodes.append(g_to_node[x])
        ## copy starting graph G into G2 so that we don't mess up G as a side-effect.
        ## when we do so, exclude the genes affected by knockout mutations.
        G2 = snap.TUNGraph.New(G.GetNodes(),G.GetEdges())
        ## add nodes.
        for n in G.Nodes():
            nId = n.GetId()
            if nId not in knocked_out_nodes:
                G2.AddNode(nId)
        ## add edges.
        for e in G.Edges():
            srcN = e.GetSrcNId()
            dstN = e.GetDstNId()
            if (srcN not in knocked_out_nodes) and (dstN not in knocked_out_nodes):
                G2.AddEdge(srcN, dstN)
        ## calculate this clone's resilience. default is 100 replicates.
        ## save memory by using a generator comprehension.
        my_resilience = sum((GraphResilience(G2) for x in range(reps)))/float(reps)
        print(my_resilience)
        randomized_resilience.append(my_resilience)
    strain_col = ['REL606'] + clone_col
    resilience_col = [REL606_resilience] + randomized_resilience
    resilience_results = pd.DataFrame.from_dict({'strain': strain_col, 'resilience': resilience_col})
    return resilience_results

def make_LTEE_strain_to_pop_dict(LTEE_knockout_muts):
    ''' return a dict of LTEE strain to the population it comes from.'''
    LTEE_strain_to_pop = {}
    LTEE_genome_metadata = LTEE_knockout_muts[['treatment','population','time','strain','clone','mutator_status']].drop_duplicates()
    LTEE_strains = list(LTEE_genome_metadata['strain'])
    matching_pops = list(LTEE_genome_metadata['population'])
    LTEE_strain_to_pop = {x:y for x,y in zip(LTEE_strains, matching_pops)}
    return LTEE_strain_to_pop

def make_LTEE_pop_to_KO_dict(LTEE_strain_to_KO, LTEE_strain_to_pop):
    ''' 
    return a dict of LTEE pops to the set of knockout muts in that pop in the 
    metagenomics dataset. Then add the union of KO mutations found in clones
    isolated in that population.
    '''
    LTEE_pop_to_KOset = {}
    LTEE_metagenomics_df = pd.read_csv("../results/LTEE-metagenome-mutations.csv")
    KO_mut_classes = ["nonsense","indel","sv"]
    KO_df = LTEE_metagenomics_df[LTEE_metagenomics_df.Annotation.isin(KO_mut_classes)]
    LTEE_pops = set(KO_df['Population'])
    for p in LTEE_pops:
        ## 'intergenic' is not a gene, so remove.
        KO_set = {x for x in KO_df[KO_df['Population']==p].Gene} - {"intergenic"}
        LTEE_pop_to_KOset[p] = KO_set
    ## now add the KO mutations in the genomics.
    for clone, pop in LTEE_strain_to_pop.items():
        clone_KOset = LTEE_strain_to_KO[clone]
        LTEE_pop_to_KOset[pop].update(clone_KOset)
    return LTEE_pop_to_KOset

def main():

    ''' Import protein-protein interaction networks from:
    1) Marinka Zitnik paper in PNAS.
    2) Qian Cong paper in Science. '''
    
    K12_edge_file = "../results/thermostability/Ecoli-Zitnik-data/Ecoli-treeoflife.interactomes/511145.txt"
    ## filter K-12 graph based on blattner IDs found in REL606.
    blattner_to_gene_dict = getREL606_blattner_to_gene_dict()
    G1, g_to_node1, node_to_g1 = create_graph_and_dicts(K12_edge_file, nodeSet=blattner_to_gene_dict)

    ''' take 2683 interactions from the Cong dataset.
    these are high-quality coevolution predictions, plus
    coevolution interactions that are known, plus
    gold standard interactions in PDB and Ecocyc.
    the information in this file is produced by
    calc-ppi-network-statistics.py.
    '''
    good_edge_f = "../results/thermostability/Cong-good-interaction-set.tsv"
    REL606_genes = getREL606_gene_set()
    ## filter Cong PPI graph based on genes in REL606.
    ## 1191 nodes, 1787 edges in this graph.
    G2, g_to_node2, node_to_g2 = create_graph_and_dicts(good_edge_f, nodeSet=REL606_genes)

    LTEE_knockouts_df = get_LTEE_genome_knockout_muts()
    LTEE_strain_to_knockouts = LTEE_strain_to_KO_genes(LTEE_knockouts_df)
    
    ## Run the resilience analysis, using the graph from Zitnik paper.
    resilience_outf1 = "../results/resilience/Zitnik_PPI_LTEE_genome_resilience.csv"
    resilience_results1 = resilience_analysis_of_LTEE_genomes(LTEE_strain_to_knockouts, G1, g_to_node1, reps=100)

    '''

    resilience_results1.to_csv(resilience_outf1)
    ## Run the resilience analysis, using the graph from Cong paper.
    resilience_outf2 = "../results/resilience/Cong_PPI_LTEE_genome_resilience.csv"
    resilience_results2 = resilience_analysis_of_LTEE_genomes(LTEE_strain_to_knockouts, G2, g_to_node2, reps=100)
    ## write results to file.
    resilience_results2.to_csv(resilience_outf2)
    
    ## calculate randomized resilience of genomes, set of all genes in REL606.
    ## for Zitnik network.
    all_genes_randomized_outf1 = "../results/resilience/Zitnik_PPI_all_genes_randomized_resilience.csv"
    all_genes_randomized_resilience1 = resilience_randomized_over_gene_set(LTEE_strain_to_knockouts, G1, g_to_node1, REL606_genes)
    all_genes_randomized_resilience1.to_csv(all_genes_randomized_outf1)
    ## for Cong network.
    all_genes_randomized_outf2 = "../results/resilience/Cong_PPI_all_genes_randomized_resilience.csv"
    all_genes_randomized_resilience2 = resilience_randomized_over_gene_set(LTEE_strain_to_knockouts, G2, g_to_node2, REL606_genes)
    all_genes_randomized_resilience2.to_csv(all_genes_randomized_outf2)

    ## calculate randomized resilience of genomes,
    ## using the set of all genes KO'ed in BOTH the LTEE metagenomics data
    ## and the LTEE genomics data.
    LTEE_genomics_KO_genes = set()
    for clone_KOset in LTEE_strain_to_knockouts.values():
        LTEE_genomics_KO_genes.update(clone_KOset)
    LTEE_metagenomics_KO_genes = get_LTEE_metagenomics_knockouts()
    all_LTEE_KO_genes = set.union(LTEE_metagenomics_KO_genes,LTEE_genomics_KO_genes)

    ## for Zitnik network.
    acrosspops_randomized_outf1 = "../results/resilience/Zitnik_PPI_across_pops_randomized_resilience.csv"
    acrosspops_randomized_resilience1 = resilience_randomized_over_gene_set(LTEE_strain_to_knockouts, G1, g_to_node1, all_LTEE_KO_genes)
    acrosspops_randomized_resilience1.to_csv(acrosspops_randomized_outf1)
    ## for Cong network.
    acrosspops_randomized_outf2 = "../results/resilience/Cong_PPI_across_pops_randomized_resilience.csv"
    acrosspops_randomized_resilience2 = resilience_randomized_over_gene_set(LTEE_strain_to_knockouts, G2, g_to_node2, all_LTEE_KO_genes)
    acrosspops_randomized_resilience2.to_csv(acrosspops_randomized_outf2)

    ## calculate randomized resilience of genomes,
    ## using the set of within population KO mutations.
    LTEE_strain_to_pop =  make_LTEE_strain_to_pop_dict(LTEE_knockouts_df)
    ## include the set of KOs in the genomics with the set of KOs in the metagenomics
    ## for each population.
    LTEE_pop_to_KO = make_LTEE_pop_to_KO_dict(LTEE_strain_to_knockouts,LTEE_strain_to_pop)
    ## for Zitnik network.
    withinpop_randomized_outf1 = "../results/resilience/Zitnik_PPI_within_pops_randomized_resilience.csv"
    withinpop_randomized_resilience1 = resilience_randomized_within_LTEE_pops(LTEE_strain_to_knockouts, LTEE_strain_to_pop, G1, g_to_node1, LTEE_pop_to_KO)
    withinpop_randomized_resilience1.to_csv(withinpop_randomized_outf1)
    ## for Cong network.
    withinpop_randomized_outf2 = "../results/resilience/Cong_PPI_within_pops_randomized_resilience.csv"
    withinpop_randomized_resilience2 = resilience_randomized_within_LTEE_pops(LTEE_strain_to_knockouts, LTEE_strain_to_pop, G2, g_to_node2, LTEE_pop_to_KO)
    withinpop_randomized_resilience2.to_csv(withinpop_randomized_outf2)

'''

if __name__ == "__main__":
    main()
