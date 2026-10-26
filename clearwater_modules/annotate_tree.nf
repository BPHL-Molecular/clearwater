process annotate_tree {
    tag { "Annotating Bootstrap Tree Leaves" }
    publishDir "${params.output}/bootstrap_tree", mode: 'copy'
    
    container null 

    input:
    path tree_file
    path scaled_data
    path metadata_file
    path active_carriers_file 

    output:
    path "annotated_phylogeny.pdf",   emit: annotated_pdf
    path "annotated_tree_data.nexus", emit: annotated_nexus

    script:
    def has_metadata = (metadata_file.name != 'EMPTY_METADATA_PLACEHOLDER' && metadata_file.exists()) ? "TRUE" : "FALSE"
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(ggtree)
    library(treeio)

    # CRASH PROTECTION: Manual override for missing ggplot2 internals
    is.waive <- function(x) inherits(x, "waiver")

    # 1. Load the Newick phylogenetic tree structure
    tree <- treeio::read.tree("${tree_file}")

    # 2. Load regression data and standardize subspecies groups
    reg_df <- read_tsv("${scaled_data}", show_col_types = FALSE) %>%
      select(SampleID = sampleID, lp_subspecies) %>%
      mutate(lp_subspecies = case_when(
        str_detect(lp_subspecies, fixed("fraseri"))     ~ "s_fraseri",
        str_detect(lp_subspecies, fixed("pneumophila")) ~ "s_pneumophila",
        str_detect(lp_subspecies, fixed("raphaeli"))    ~ "s_raphaeli",
        TRUE ~ as.character(lp_subspecies)
      ))

    # 3. Incorporate isolate source and county tracking with defensive column matching
    has_meta_flag <- ${has_metadata}
    if (has_meta_flag) {
      raw_meta <- read_csv("${metadata_file}", show_col_types = FALSE)
      
      # SAFEGUARD: Strip invisible byte signatures, spaces, or carriage returns from headers
      names(raw_meta) <- trimws(names(raw_meta))
      
      # FLEXIBLE MATCHING: Dynamically identifies columns even if case variations or BOM tags persist
      sample_col <- grep("SampleID", names(raw_meta), ignore.case = TRUE, value = TRUE)[1]
      src_col    <- grep("Isolate_Sources", names(raw_meta), ignore.case = TRUE, value = TRUE)[1]
      county_col <- grep("County", names(raw_meta), ignore.case = TRUE, value = TRUE)[1]
      
      if (is.na(sample_col) || is.na(src_col) || is.na(county_col)) {
         stop(paste("CRITICAL METADATA ERROR: Could not find required columns in independent CSV file! Found:", 
              paste(names(raw_meta), collapse=", ")))
      }
      
      meta_df <- raw_meta %>%
        select(SampleID = all_of(sample_col), Isolate_Sources = all_of(src_col), County = all_of(county_col)) %>%
        mutate(SampleID = as.character(SampleID))
        
      annotation_df <- left_join(reg_df, meta_df, by = "SampleID")
    } else {
      annotation_df <- reg_df %>% mutate(Isolate_Sources = NA_character_, County = NA_character_)
    }

    # 4. Load the outbreak roster file to flag active transmission carriers
    carriers_df <- read_tsv("${active_carriers_file}", show_col_types = FALSE) %>%
      select(SampleID) %>%
      mutate(Outbreak_Status = "[OUTBREAK]")

    # Join the active transmission tags onto your foundational tracking matrix smoothly
    annotation_df <- left_join(annotation_df, carriers_df, by = "SampleID") %>%
      mutate(Outbreak_Status = replace_na(Outbreak_Status, ""))

    annotation_df <- annotation_df %>% rename(label = SampleID)

    # 5. Bind traits to tree nodes
    annotated_tree <- treeio::full_join(tree, annotation_df, by = "label")

    # 6. Authoritative color palette
    subspecies_colors <- c(
      "s_fraseri"     = "#E41A1C", 
      "s_pneumophila" = "#377EB8", 
      "s_raphaeli"    = "#4DAF4A"  
    )

    # 7. Initialize foundation ggtree canvas with colored labels
    p <- ggtree(annotated_tree, layout = "rectangular") + 
         geom_tiplab(aes(color = lp_subspecies), size = 2.2, hjust = -0.1, align = TRUE) +
         scale_color_manual(values = subspecies_colors, name = "L. pneumophila Subspecies") +
         theme(legend.position = "right")

    # 8. Render Bootstrap Values right above internal branches
    p <- p + geom_nodelab(
      aes(label = label), 
      size = 1.8, 
      color = "#444444", 
      vjust = -0.6, 
      hjust = 1.2
    )

    # DYNAMIC SPACING CALCULATIONS
    max_x <- max(p\$data\$x, na.rm = TRUE)
    max_y <- max(p\$data\$y, na.rm = TRUE)

    # Calculate absolute baseline coordinates for the columns
    sample_id_x  = max_x * 1.02
    outbreak_x   = max_x * 1.18  
    meta_src_x   = max_x * 1.32  
    county_x     = max_x * 1.55  
    header_y     = max_y + 1.2
    underline_y  = max_y + 1.0  

    # =========================================================================
    # 9. ADD HEADERS AND PRECISE TEXT-MATCHED UNDERLINES
    # =========================================================================
    
    # A. Column 1: "SampleID"
    p <- p + annotate(
      "text", x = sample_id_x, y = header_y, label = "SampleID",
      size = 2.4, fontface = "bold", color = "#111111", hjust = 0
    ) +
    annotate(
      "text", x = sample_id_x, y = underline_y, label = "________",
      size = 2.4, fontface = "bold", color = "#111111", hjust = 0
    )

    # B. Column 2: "Status" Header & Row Alert Tags
    p <- p + geom_text(
      data = filter(p\$data, isTip & Outbreak_Status != ""),
      aes(x = outbreak_x, y = y, label = Outbreak_Status), 
      size = 1.9, fontface = "bold", color = "#E41A1C", hjust = 0 
    ) +
    annotate(
      "text", x = outbreak_x, y = header_y, label = "Status",
      size = 2.4, fontface = "bold", color = "#111111", hjust = 0
    ) +
    annotate(
      "text", x = outbreak_x, y = underline_y, label = "______",
      size = 2.4, fontface = "bold", color = "#111111", hjust = 0
    )

    # C. Columns 3 & 4 (Only populated if metadata csv flag is true)
    if (has_meta_flag) {
      p <- p + geom_text(
        data = filter(p\$data, isTip),
        aes(x = meta_src_x, y = y, label = Isolate_Sources), 
        size = 2.0, color = "#555555", hjust = 0 
      )
      
      p <- p + annotate(
        "text", x = meta_src_x, y = header_y, label = "Isolate_Sources",
        size = 2.4, fontface = "bold", color = "#111111", hjust = 0
      ) +
      annotate(
        "text", x = meta_src_x, y = underline_y, label = "_______________",
        size = 2.4, fontface = "bold", color = "#111111", hjust = 0
      )

      p <- p + geom_text(
        data = filter(p\$data, isTip),
        aes(x = county_x, y = y, label = County), 
        size = 2.0, color = "#555555", hjust = 0 
      )
      
      p <- p + annotate(
        "text", x = county_x, y = header_y, label = "County",
        size = 2.4, fontface = "bold", color = "#111111", hjust = 0
      ) +
      annotate(
        "text", x = county_x, y = underline_y, label = "______",
        size = 2.4, fontface = "bold", color = "#111111", hjust = 0
      )
    }

    # Add a clean evolutionary distance scale bar
    p <- p + geom_treescale(x = 0, y = -1, width = max_x * 0.1)

    # 10. Dynamically expand total canvas right boundary to clear all columns nicely
    extended_boundary <- if (has_meta_flag) max_x * 2.0 else max_x * 1.5
    p <- p + xlim(0, extended_boundary) +
         scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) 

    # 11. Save production-ready analytical outputs
    pdf_width <- if (has_meta_flag) 13 else 9.5 
    ggsave("annotated_phylogeny.pdf", plot = p, width = pdf_width, height = 11, units = "in")

    # Save a fully compatible Nexus format file for external viewer compatibility
    treeio::write.beast(annotated_tree, file = "annotated_tree_data.nexus")
    """
}
