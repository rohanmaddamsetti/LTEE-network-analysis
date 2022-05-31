"""
cluster-jenga-genomes-and-run-STIMS.jl by Rohan Maddamsetti.
Make Supplementary Figure S3 and run STIMS to get statistics quickly.
"""

using DataFrames, DataFramesMeta, CSV, Distances, Clustering, Plots, Plots.PlotMeasures
## NOTA BENE: for this project, I added some code to STIMS.jl that at
## present are NOT in the main branch of STIMS development.
include("STIMS.jl")

################################################################################
## Import results from playing Genome Jenga.

minimal_genome_data = CSV.read(
    "../results/metabolic-enzymes/jenga_minimal_genomes.csv", DataFrame)
essential_genes_data = CSV.read(
    "../results/metabolic-enzymes/jenga_essential_genes.csv", DataFrame)
minimal_rxns_data = CSV.read(
    "../results/metabolic-enzymes/jenga_minimal_reactions.csv", DataFrame)

################################################################################
## represent each minimal genome as a vector, and then cluster by similarity.
## prediction: we see multiple clusters.

function GeneToGenomeMatrix(genome_data)
    genes_for_rows = unique(genome_data.locus_tag)
    gene_to_matrix_idx = Dict(gene => i for (i, gene) in enumerate(genes_for_rows))
    samples_for_cols = unique(genome_data.Replicate)

    gene_genome_matrix = zeros(Bool, length(genes_for_rows),length(samples_for_cols))
    ## fill in the values of the matrix with 1s if the gene is present in the genome.
    for (j, cur_replicate) in enumerate(samples_for_cols)
        replicate_data = @rsubset(genome_data, :Replicate == cur_replicate)
        for x = 1:nrow(replicate_data)
            cur_gene = replicate_data.locus_tag[x]
            i = gene_to_matrix_idx[cur_gene]
            gene_genome_matrix[i,j] = true
        end
    end
    return gene_genome_matrix
end


function ReactionToNetworkMatrix(network_data)
    rxns_for_rows = unique(network_data.Reaction)
    rxn_to_matrix_idx = Dict(rxn => i for (i, rxn) in enumerate(rxns_for_rows))
    samples_for_cols = unique(network_data.Replicate)

    rxn_network_matrix = zeros(Bool, length(rxns_for_rows),length(samples_for_cols))
    ## fill in the values of the matrix with 1s if the rxn is present in the network.
    for (j, cur_replicate) in enumerate(samples_for_cols)
        replicate_data = @rsubset(network_data, :Replicate == cur_replicate)
        for x = 1:nrow(replicate_data)
            cur_rxn = replicate_data.Reaction[x]
            i = rxn_to_matrix_idx[cur_rxn]
            rxn_network_matrix[i,j] = true
        end
    end
    return rxn_network_matrix
end


function JaccardSortMatrix(M)
    ## calculate Jaccard distances between columns (genomes).
    col_dist_matrix = pairwise(Jaccard(), M, dims=2)
    ## and cluster.
    col_hclust = hclust(col_dist_matrix, branchorder = :optimal)

    ## calculate Jaccard distances between rows (genes)
    row_dist_matrix = pairwise(Jaccard(), M, dims=1)
    ## and cluster.
    row_hclust = hclust(row_dist_matrix)#, branchorder = :optimal)
    ## and sort the matrix by the clustering.
    sorted_M = M[row_hclust.order, col_hclust.order]
    return sorted_M
end


function RemoveUninformativeRows(M)
    ## remove rows which are all 1s (present in every single column/sample)
    nrows_of_M, ncols_of_M = size(M)
    rows_to_keep = []
    for i in 1:nrows_of_M
        if sum(M[i,:]) < ncols_of_M
            append!(rows_to_keep, i)
        end
    end
    return M[rows_to_keep, :]
end


sorted_minimal_genome_matrix = minimal_genome_data |>
    GeneToGenomeMatrix |>
    JaccardSortMatrix

sorted_essential_gene_matrix = essential_genes_data |>
    GeneToGenomeMatrix |>
    JaccardSortMatrix

sorted_minimal_rxns_matrix = minimal_rxns_data |>
    ReactionToNetworkMatrix |>
    JaccardSortMatrix

## Supplementary Figure 3A.
S3FigA = heatmap(sorted_minimal_genome_matrix,
                 c = cgrad(:blues, 2, categorical = true),
                 colorbar = false,
                 xlabel = "Genomes", ylabel = "Genes",
                 fontfamily = "Helvetica",
                 xrotation = 45)
## Supplementary Figure 3B.
S3FigB = heatmap(sorted_essential_gene_matrix,
                 c = cgrad(:blues, 2, categorical = true),
                 colorbar = false,
                 xlabel = "Genomes", ylabel = "Essential genes",
                 fontfamily = "Helvetica",
                 xrotation = 45)
## Supplementary Figure 3C.
S3FigC = heatmap(sorted_minimal_rxns_matrix,
                 c = cgrad(:blues, 2, categorical = true),
                 colorbar = false,
                 xlabel = "Genomes", ylabel = "Reactions",
                 fontfamily = "Helvetica",
                 xrotation = 45)

S3Fig = plot(S3FigA, S3FigB, S3FigC, layout=(1,3), legend=false)
savefig(S3Fig, "../results/metabolic-enzymes/S3Fig.pdf")

################################################################################
## Run STIMS statistics, quickly.


## statistics for BiGG core.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/BiGG-core-genes.csv",
               "../results/metabolic-enzymes/BiGG-core-plot.pdf")
## results for BiGG core:
## Ara-5: 0.3154
## Ara-6: 0.2143
## Ara+1: 0.7846
## Ara+2: 0.4091
## Ara+4: 0.6205
## Ara+5: 0.2449
## Ara-1: 0.9493
## Ara-2: 0.6907
## Ara-3: 0.1294
## Ara-4: 0.9339
## Ara+3: 0.9628
## Ara+6: 0.99774

## now repeat, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/BiGG-core-genes.csv",
               "../results/metabolic-enzymes/hypermutator-epoch-BiGG-core-plot.pdf")
## results for BiGG core, just hypermutator epochs:
## Ara-1: 0.9841
## Ara-2: 0.7578
## Ara-3: 0.1061
## Ara-4: 0.9404
## Ara+3: 0.9682
## Ara+6: 0.996


## statistics for superessential metabolic reactions
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/superessential-genes.csv",
               "../results/metabolic-enzymes/superessential-plot.pdf")
## results for superessential genes:
## Ara-5: 0.6621
## Ara-6: 0.9039
## Ara+1: 0.8534
## Ara+2: 0.6734
## Ara+4: 0.8768
## Ara+5: 0.1662
## Ara-1: 1.0
## Ara-2: 0.6516
## Ara-3: 0.6094
## Ara-4: 0.8993
## Ara+3: 0.5625
## Ara+6: 0.9963

## now repeat, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/superessential-genes.csv",
               "../results/metabolic-enzymes/hypermutator-epoch-superessential-plot.pdf")
## results for superessential genes, just hypermutator epochs:
## Ara-1: 1.0
## Ara-2: 0.6625
## Ara-3: 0.5396
## Ara-4: 0.9146
## Ara+3: 0.5475
## Ara+6: 0.9975


STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/specialist-enzymes.csv",
               "../results/metabolic-enzymes/specialist-enzymes-plot.pdf")
## results for specialist enzymes:
## Ara-5: 0.241
## Ara-6: 0.8504
## Ara+1: 0.9901
## Ara+2: 0.9678
## Ara+4: 0.8211
## Ara+5: 0.9354
## Ara-1: 0.5416
## Ara-2: 0.4553
## Ara-3: 0.9180
## Ara-4: 0.9468
## Ara+3: 0.7823
## Ara+6: 0.9998

## now repeat, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/specialist-enzymes.csv",
               "../results/metabolic-enzymes/hypermutator-epoch-specialist-enzymes.pdf")
## results for specialist enzymes, just hypermutator epochs:
## Ara-1: 0.5611
## Ara-2: 0.3723
## Ara-3: 0.8475
## Ara-4: 0.9359
## Ara+3: 0.7829
## Ara+6: 0.9996


STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/generalist-enzymes.csv",
               "../results/metabolic-enzymes/generalist-enzymes-plot.pdf")
## results for generalist enzymes:
## Ara-5: 0.7716
## Ara-6: 0.4526
## Ara+1: 0.6077
## Ara+2: 0.4343
## Ara+4: 0.3933
## Ara+5: 0.0277
## Ara-1: 0.9990
## Ara-2: 0.6964
## Ara-3: 0.3308
## Ara-4: 0.6898
## Ara+3: 0.7102
## Ara+6: 0.8969

## now repeat, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               "../results/metabolic-enzymes/generalist-enzymes.csv",
               "../results/metabolic-enzymes/hypermutator-epoch-generalist-enzymes.pdf")
## results for generalist enzymes, just hypermutator epochs:
## Ara-1: 0.9981
## Ara-2: 0.6889
## Ara-3: 0.2283
## Ara-4: 0.657
## Ara+3: 0.6983
## Ara+6: 0.897


## statistics for Figure 5, and Supplementary Figure 4.
## run STIMS on core and essential genes in the minimal genomes.

jenga_genome_core_csv = "../results/metabolic-enzymes/jenga-genome-core.csv"
jenga_essential_core_csv = "../results/metabolic-enzymes/jenga-essential-core.csv"

## Run STIMS on the core genome of the minimal genomes.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               jenga_genome_core_csv,
               "../results/metabolic-enzymes/jenga-core-genes.pdf")
## Results for core genome of the minimal genomes:
## Ara-5: 0.9293
## Ara-6: 0.91
## Ara+1: 0.9323
## Ara+2: 0.8415
## Ara+4: 0.9989
## Ara+5: 0.7468
## Ara-1: 1.0
## Ara-2: 0.9662
## Ara-3: 0.8234
## Ara-4: 0.9827
## Ara+3: 0.992
## Ara+6: 1.0

## now repeat on the core genome of the minimal genomes, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               jenga_genome_core_csv,
               "../results/metabolic-enzymes/hypermutator-epoch-jenga-core-genes.pdf")
## results for jenga core, just hypermutator epochs:
## Ara-1: 1.0
## Ara-2: 0.9647
## Ara-3: 0.7241
## Ara-4: 0.979
## Ara+3: 0.9946
## Ara+6: 1.0


## Run STIMS on the core essential genes.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
         "../results/REL606_IDs.csv",
         jenga_essential_core_csv,
         "../results/metabolic-enzymes/jenga-essential-core-genes.pdf")
## Results on core essential genes:
## Ara-5: 0.8985
## Ara-6: 0.8473
## Ara+1: 0.9825
## Ara+2: 0.7857
## Ara+4: 0.9949
## Ara+5: 0.6264
## Ara-1: 1.0
## Ara-2: 0.9763
## Ara-3: 0.8782
## Ara-4: 0.9471
## Ara+3: 0.9969
## Ara+6: 1.0

## now repeat on the core essential genes, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               jenga_essential_core_csv,
               "../results/metabolic-enzymes/hypermutator-epoch-jenga-essential-core-genes.pdf")
## results for jenga essential core, just hypermutator epochs:
## Ara-1: 1.0
## Ara-2: 0.9745
## Ara-3: 0.8236
## Ara-4: 0.9532
## Ara+3: 0.9981
## Ara+6: 1.0


## statistics for Figure 6, and Supplementary Figure 5.
## run STIMS on essential genes for growth on
## glucose and citrate in the ancestral REL606 FBA model.

glucose_essential_csv = "../results/metabolic-enzymes/glucose_FBA_essential.csv"
## Run STIMS on essential genes in the REL606 genome FBA model.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               glucose_essential_csv,
               "../results/metabolic-enzymes/glucose-essential-genes.pdf")
## Results for glucose essential genes in the REL606 FBA model:
## Ara-5: 0.9376
## Ara-6: 0.6595
## Ara+1: 0.9533
## Ara+2: 0.8343
## Ara+4: 0.9995
## Ara+5: 0.7281
## Ara-1: 1.0
## Ara-2: 0.9654
## Ara-3: 0.9488
## Ara-4: 0.9467
## Ara+3: 0.9486
## Ara+6: 1.0

## now repeat on glucose essential genes, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               glucose_essential_csv,
               "../results/metabolic-enzymes/hypermutator-epoch-glucose-essential-genes.pdf")
## results for glucose essential genes, just hypermutator epochs:
## Ara-1: 1.0
## Ara-2: 0.9776
## Ara-3: 0.9218
## Ara-4: 0.9331
## Ara+3: 0.9445
## Ara+6: 1.0


citrate_essential_csv = "../results/metabolic-enzymes/citrate_FBA_essential.csv"
## Run STIMS on essential genes in the Citrate REL606 genome FBA model.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               citrate_essential_csv,
               "../results/metabolic-enzymes/citrate-essential-genes.pdf")
## Results for citrate essential genes in the REL606 FBA model:
## Ara-5: 0.9058
## Ara-6: 0.787
## Ara+1: 0.9598
## Ara+2: 0.8405
## Ara+4: 0.9992
## Ara+5: 0.7414
## Ara-1: 1.0
## Ara-2: 0.975
## Ara-3: 0.9861
## Ara-4: 0.9498
## Ara+3: 0.9435
## Ara+6: 1.0

## now repeat on citrate essential genes, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               citrate_essential_csv,
               "../results/metabolic-enzymes/hypermutator-epoch-citrate-essential-genes.pdf")
## results for citrate essential genes, just hypermutator epochs:
## Ara-1: 1.0
## Ara-2: 0.9838
## Ara-3: 0.9693
## Ara-4: 0.9418
## Ara+3: 0.9396
## Ara+6: 1.0


acetate_essential_csv = "../results/metabolic-enzymes/acetate_FBA_essential.csv"
## Run STIMS on essential genes in the Acetate REL606 genome FBA model.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               acetate_essential_csv,
               "../results/metabolic-enzymes/acetate-essential-genes.pdf")
## Results for acetate essential genes in the REL606 FBA model:
## Ara-5: 0.9227
## Ara-6: 0.7025
## Ara+1: 0.9513
## Ara+2: 0.8597
## Ara+4: 0.9989
## Ara+5: 0.7645
## Ara-1: 1.0
## Ara-2: 0.9835
## Ara-3: 0.9566
## Ara-4: 0.9636
## Ara+3: 0.9771
## Ara+6: 1.0

## now repeat on acetate essential genes, just on hypermutator epochs.
STIMS.RunSTIMS("../results/hypermutator-epoch-LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               acetate_essential_csv,
               "../results/metabolic-enzymes/hypermutator-epoch-acetate-essential-genes.pdf")
## results for acetate essential genes, just hypermutator epochs:
## Ara-1: 1.0
## Ara-2: 0.9884
## Ara-3: 0.937
## Ara-4: 0.9534
## Ara+3: 0.9742
## Ara+6: 1.0


## This test is redundant-- citrate essential genes == glucose+acetate essential genes,
## in the FBA analysis!
glucose_acetate_essential_csv = "../results/metabolic-enzymes/glucose_acetate_FBA_essential.csv"
## Run STIMS on essential genes in the Glucose + Acetate REL606 genome FBA model.
STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               glucose_acetate_essential_csv,
               "../results/metabolic-enzymes/glucose-acetate-essential-genes.pdf")
## Results for glucose + acetate essential genes in the REL606 FBA model:
## Ara-5: 0.9158
## Ara-6: 0.7864
## Ara+1: 0.9635
## Ara+2: 0.8465
## Ara+4: 0.9996
## Ara+5: 0.736
## Ara-1: 1.0
## Ara-2: 0.9708
## Ara-3: 0.9864
## Ara-4: 0.9531
## Ara+3: 0.9435
## Ara+6: 1.0


## let's run STIMS on the genes that are essential for citrate (glucose+acetate),
## but not glucose alone.
glucose_essential_df = CSV.read(glucose_essential_csv, DataFrame)
citrate_essential_df = CSV.read(citrate_essential_csv, DataFrame)
acetate_essential_df = CSV.read(acetate_essential_csv, DataFrame)

## This is only citT!
citrate_not_acetate_essential_genes = setdiff(Set(citrate_essential_df.Gene),
                                              Set(acetate_essential_df.Gene))
## This is gapA, citT, pgk, tpiA.
citrate_not_glucose_essential_genes = setdiff(Set(citrate_essential_df.Gene),
                                              Set(glucose_essential_df.Gene))

citrate_not_glucose_essential_df = @rsubset(
    citrate_essential_df,
    :Gene in citrate_not_glucose_essential_genes)

citrate_not_glucose_essential_csv = "../results/metabolic-enzymes/citrate_not_glucose_FBA_essential.csv"
CSV.write(citrate_not_glucose_essential_csv, citrate_not_glucose_essential_df)

STIMS.RunSTIMS("../results/LTEE-metagenome-mutations.csv",
               "../results/REL606_IDs.csv",
               citrate_not_glucose_essential_csv,
               "../results/metabolic-enzymes/citrate-not-glucose-essential-genes.pdf")
## STIMS results for gapA, citT, pgk, tpiA.
## basically, not enough power, but the signal is promising,
## check out Ara-2 and Ara-3.
## Ara-5: 0.1885
## Ara-6: 0.145
## Ara+1: 0.2781
## Ara+2: 0.1245
## Ara+4: 0.144
## Ara+5: 0.1402
## Ara-1: 0.8635
## Ara-2: 0.9023
## Ara-3: 0.9088
## Ara-4: 0.7547
## Ara+3: 0.5001
## Ara+6: 0.2379
