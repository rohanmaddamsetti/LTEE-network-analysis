#!/usr/bin/env python

'''
calc-ppi-network-statistics.py by Rohan Maddamsetti.

This script calculates PPI network statistics using the SNAP library and
writes results to file, in order to make comparisons with evolutionary rates
 in the LTEE.

'''

import snap
from math import log
from os.path import join
import random
import numpy as np
from scipy.integrate import simps
import pandas as pd


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


def getCongDataSet(f, col1, col2, col_val_filter=(None,None)):
    '''
    represent PPIs as a set of tuples of gene names.
    col_val_filter is a tuple of the column (0-indexed int) 
    and the value to filter on.
    Use this flag to filter the coevolution findings on what has been
    previously reported, as flagged in Supplementary Table S8 of 
    Cong et al. (2019).
    '''
    pairs = set()
    with open(f,"r") as fh:
        for i,l in enumerate(fh):
            if i == 0: continue ## skip header
            l = l.strip()
            ldata = l.split(',')
            if col_val_filter[0] is not None:
                if ldata[col_val_filter[0]] != col_val_filter[1]:
                    continue
            g1 = ldata[col1]
            g2 = ldata[col2]
            ## add gene pairs in alphabetical order.
            if g1 < g2:
                pairs.add((g1,g2))
            else:
                pairs.add((g2,g1))
    return pairs


def write_Cong_good_interactions(good_interaction_path):
    ## import pairs from Qian Cong dataset.
    cong_data_dir = "../results/thermostability/Ecoli-Cong-data/"
    ## PDB interactions.
    pdb_path = join(cong_data_dir,"S4-Ecoli-PDB-benchmark-set.csv")
    pdb_set = getCongDataSet(pdb_path, col1=1, col2=4)
    ## Ecocyc interactions.
    ecocyc_path = join(cong_data_dir,"S5-Ecocyc-benchmark-set.csv")
    ecocyc_set = getCongDataSet(ecocyc_path, col1=1, col2=4)
    ## yeast two-hybrid interactions.
    Y2H_path = join(cong_data_dir,"S6-Y2H-benchmark-set.csv")
    Y2H_set = getCongDataSet(Y2H_path, col1=1, col2=4)
    ## affinity purification mass-spec interactions.
    APMS_path = join(cong_data_dir,"S7-APMS-benchmark-set.csv")
    APMS_set = getCongDataSet(APMS_path, col1=1, col2=4)
    ## coevolution interactions.
    coevolution_path = join(cong_data_dir,"S8-Ecoli-PPI-by-coevolution.csv")
    ## Note: we are filtering this for the 936 previously known interactions.
    coevolution_set = getCongDataSet(coevolution_path, col1=6 , col2=7, col_val_filter=(5,"yes"))
    ## high-confidence new interactions by coevolution.
    high_conf_new_coevolution_path = join(cong_data_dir,"S10-Confident-novel-interactions-in-coevolution-screen.csv")
    high_conf_new_coevolution_set = getCongDataSet(high_conf_new_coevolution_path, col1=2, col2=3)
    ## 2683 edges in this set.
    good_interaction_set = set.union(pdb_set, ecocyc_set, coevolution_set, high_conf_new_coevolution_set)
    ## write out the good interactions to file.
    with open(good_interaction_path,"w") as outfh:
        good_interaction_list = sorted(list(good_interaction_set))
        for p in good_interaction_list:
            g1, g2 = p
            outfh.write("\t".join([g1,g2])+"\n")
    return


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


def calc_PageRank(Graph, node_to_g):
    prot_to_PageRank = {}
    PRankH = snap.TIntFltH()
    snap.GetPageRank(Graph, PRankH)
    for node in PRankH:
        my_prot = node_to_g[node]
        prot_to_PageRank[my_prot] = PRankH[node]
    return prot_to_PageRank


def calc_HubAndAuthorityScores(Graph, node_to_g):
        ## calculate Hub and Authority scores for nodes in the graph.
    prot_to_hub = {}
    prot_to_authority = {}
    NIdHubH = snap.TIntFltH()
    NIdAuthH = snap.TIntFltH()
    snap.GetHits(Graph, NIdHubH, NIdAuthH)
    for node in NIdHubH:
        my_prot = node_to_g[node]
        prot_to_hub[my_prot] = NIdHubH[node]
    for node in NIdAuthH:
        my_prot = node_to_g[node]
        prot_to_authority[my_prot] = NIdAuthH[node]
    return (prot_to_hub, prot_to_authority)


def calc_ClosenessCentrality(Graph, node_to_g):
    prot_to_closeness_centrality = {}
    for NI in Graph.Nodes():
        my_prot = node_to_g[NI.GetId()]
        CloseCentr = snap.GetClosenessCentr(Graph, NI.GetId())
        prot_to_closeness_centrality[my_prot] = CloseCentr        
    return prot_to_closeness_centrality


def calc_BetweenessCentrality(Graph, node_to_g):
    prot_to_betweeness_centrality = {}
    Nodes = snap.TIntFltH()
    Edges = snap.TIntPrFltH()
    snap.GetBetweennessCentr(Graph, Nodes, Edges, 1.0)
    for node in Nodes:
        my_prot = node_to_g[node]
        prot_to_betweeness_centrality[my_prot] = Nodes[node]
    return prot_to_betweeness_centrality


def calc_EigenvectorCentrality(Graph, node_to_g):
    prot_to_eigenvector_centrality = {}
    NIdEigenH = snap.TIntFltH()
    snap.GetEigenVectorCentr(Graph, NIdEigenH)
    for node in NIdEigenH:
        my_prot = node_to_g[node]
        prot_to_eigenvector_centrality[my_prot] = NIdEigenH[node] 
    return prot_to_eigenvector_centrality


def calc_Degrees(Graph, node_to_g):
    prot_to_degree = {}
    for NI in Graph.Nodes():
        my_prot = node_to_g[NI.GetId()]
        in_degree = NI.GetInDeg()
        out_degree = NI.GetOutDeg()
        assert in_degree == out_degree
        prot_to_degree[my_prot] = in_degree
    return prot_to_degree


def calc_DegreeCentrality(Graph, node_to_g):
    prot_to_degree_centrality = {}
    for NI in Graph.Nodes():
        my_prot = node_to_g[NI.GetId()]
        ## degree centrality of the node
        DegCentr = snap.GetDegreeCentr(Graph, NI.GetId())
        prot_to_degree_centrality[my_prot] = DegCentr
    return prot_to_degree_centrality


def getArticulationPoints(Graph, node_to_g):
    ''' A vertex in an undirected connected graph is an articulation point (or cut vertex)
 iff removing it (and edges through it) disconnects the graph.
 Articulation points represent vulnerabilities in a connected network –
 single points whose failure would split the network into 2 or more disconnected components.
 They are useful for designing reliable networks. '''
    prot_to_articulation_points = {}
    ArtNIdV = snap.TIntV()
    snap.GetArtPoints(Graph, ArtNIdV)
    for node in ArtNIdV:
        my_prot = node_to_g[node]
        prot_to_articulation_points[my_prot] = True
    return prot_to_articulation_points


def getStronglyConnectedComponents(Graph, node_to_g):
    prot_to_SCcomponent = {}
    Components = snap.TCnComV()
    snap.GetSccs(Graph, Components)
    for i, CnCom in enumerate(Components):
        for node in CnCom:
            my_prot = node_to_g[node]
            prot_to_SCcomponent[my_prot] = i+1 ##1-index component membership.
    return prot_to_SCcomponent


def write_NetworkStatistics(Graph, node_to_g, outf):
    ## make dictionaries of genes to network statistics.
    p_to_pagerank = calc_PageRank(Graph, node_to_g)
    p_to_hub, p_to_authority = calc_HubAndAuthorityScores(Graph, node_to_g)
    p_to_closeness_centrality = calc_ClosenessCentrality(Graph, node_to_g)
    p_to_betweeness_centrality = calc_BetweenessCentrality(Graph, node_to_g)
    p_to_eigenvector_centrality = calc_EigenvectorCentrality(Graph, node_to_g)
    p_to_degree = calc_Degrees(Graph, node_to_g)
    p_to_degree_centrality = calc_DegreeCentrality(Graph, node_to_g)
    p_to_articulation_points = getArticulationPoints(Graph, node_to_g)
    p_to_SCcomponent = getStronglyConnectedComponents(Graph, node_to_g)
    header = "Gene,Pagerank,HubScore,AuthorityScore,ClosenessCentrality,BetweenessCentrality,EigenvectorCentrality,Degree,DegreeCentrality,IsArticulationPoint,StronglyConnectedComponent"

    with open(outf,"w") as outfh:
        outfh.write(header+"\n")
        for _, p in node_to_g.items():
            f1 = p
            f2 = p_to_pagerank[p]
            f3 = p_to_hub[p]
            f4 = p_to_authority[p]
            f5 = p_to_closeness_centrality[p]
            f6 = p_to_betweeness_centrality[p]
            f7 = p_to_eigenvector_centrality[p]
            f8 = p_to_degree[p]
            f9 = p_to_degree_centrality[p]
            if p in p_to_articulation_points:
                p_art = 1
            else:
                p_art = 0
            f10 = p_art
            if p in p_to_SCcomponent:
                p_SC = p_to_SCcomponent[p]
            else:
                p_SC = -1
            f11 = p_SC
            row = ','.join([str(x) for x in [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11]])
            outfh.write(row+"\n")
    return


def main():

    ''' Import protein-protein interaction networks from:
    1) Marinka Zitnik paper in PNAS.
    2) Qian Cong paper in Science. '''
    
    K12_edge_file = "../results/thermostability/Ecoli-Zitnik-data/Ecoli-treeoflife.interactomes/511145.txt"
    ## filter K-12 graph based on blattner IDs found in REL606.
    blattner_to_gene_dict = getREL606_blattner_to_gene_dict()
    G1, g_to_node1, node_to_g1 = create_graph_and_dicts(K12_edge_file, nodeSet=blattner_to_gene_dict)

    ## write out 2683 interactions in Cong dataset.
    ## these are high-quality coevolution predictions, plus
    ## coevolution interactions that are known, plus
    ## gold standard interactions in PDB and Ecocyc.
    good_edge_f = "../results/thermostability/Cong-good-interaction-set.tsv"
    write_Cong_good_interactions(good_edge_f)
    REL606_genes = getREL606_gene_set()
    ## filter Cong PPI graph based on genes in REL606.
    ## 1191 nodes, 1787 edges in this graph.
    G2, g_to_node2, node_to_g2 = create_graph_and_dicts(good_edge_f, nodeSet=REL606_genes)

    ## write network statistics to file for analysis in thermostability-analysis.R.
    G1file = "../results/thermostability/Zitnik_network_statistics.csv"
    G2file = "../results/thermostability/Cong_network_statistics.csv"
    write_NetworkStatistics(G1, node_to_g1, G1file)
    write_NetworkStatistics(G2, node_to_g2, G2file)

if __name__ == "__main__":
    main()
