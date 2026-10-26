process cluster_transmission_chains {
    tag { "Computing Single-Linkage Transmission Clusters" }
    publishDir "${params.output}/transmission_analysis", mode: 'copy'
    
    container null 

    input:
    path matrix_tsv

    output:
    path "transmission_cluster_memberships.tsv", emit: clusters_tsv
    path "active_transmission_carriers.tsv",     emit: active_carriers_tsv // <-- CHANGED: New required output track
    path "transmission_clustering_summary.txt",  emit: summary_log

    script:
    def threshold = params.snp_threshold
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(igraph)

    mat_data <- read.table("${matrix_tsv}", header = TRUE, sep = "\\t", row.names = 1, check.names = FALSE)
    mat <- as.matrix(mat_data)

    adjacency_matrix <- ifelse(mat < ${threshold}, 1, 0)
    diag(adjacency_matrix) <- 0

    network <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", diag = FALSE)
    components <- igraph::components(network)
    
    cluster_ids <- components\$membership
    total_clusters <- components\$no
    cluster_sizes <- components\$csize

    output_df <- tibble(
      SampleID = names(cluster_ids),
      TransmissionClusterID = as.integer(cluster_ids)
    ) %>% arrange(TransmissionClusterID, SampleID)

    write_tsv(output_df, "transmission_cluster_memberships.tsv")

    # =========================================================================
    # CHANGED: EXTRACT ONLY ACTIVE MULTI-ISOLATE TRANSMISSION CARRIERS
    # =========================================================================
    active_cluster_numbers <- which(cluster_sizes >= 2)
    active_carriers_df <- output_df %>% filter(TransmissionClusterID %in% active_cluster_numbers)
    write_tsv(active_carriers_df, "active_transmission_carriers.tsv")

    log_file <- "transmission_clustering_summary.txt"
    multi_isolate_clusters <- length(active_cluster_numbers)
    total_in_transmission <- nrow(active_carriers_df)

    write_lines(
      c(
        "==========================================================================================",
        "GENOMIC EPIDEMIOLOGY SUMMARY REPORT: LEGIONELLA PNEUMOPHILA TRANSMISSION CHAINS",
        "==========================================================================================",
        paste("Total Isolates Analyzed:      ", nrow(mat)),
        paste("Configured Cluster Threshold: Less than ${threshold} SNPs (Distance <= ${threshold - 1})"),
        paste("Total Clusters Identified:     ", total_clusters),
        paste("Active Multi-Isolate Chains:  ", multi_isolate_clusters),
        paste("Total Genomes Part of Chain:  ", total_in_transmission),
        "------------------------------------------------------------------------------------------",
        "Cluster ID Size Breakdown:",
        capture.output(table(TransmissionClusterID = cluster_ids)),
        "=========================================================================================="
      ),
      file = log_file
    )
    """
}
