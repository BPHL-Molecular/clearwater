process firth_regression {
    container null 
    tag { "Fitting Global Firth Model" }
    publishDir "${params.output}/firth_regression_analysis", mode: 'copy'

    input:
    path imputed_table

    output:
    path "firth_multinomial_results.txt",   emit: stats_table
    path "firth_multinomial_formatted.txt", emit: formatted_table
    path "scaled_subspecies_data.tsv",      emit: scaled_subspecies_data
    path "firth_model.rds",                 emit: model_object

    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(stringr)
    library(brglm2)

    # 1. Load and clean the global data
    data <- read_tsv("${imputed_table}", show_col_types = FALSE)
    data <- data %>% dplyr::rename(lp_subspecies = subspeciesID_kraken)

    data <- data %>%
      mutate(lp_subspecies = case_when(
        str_detect(lp_subspecies, fixed("subsp. pneumophila")) ~ "subsp. pneumophila",
        str_detect(lp_subspecies, fixed("subsp. fraseri"))     ~ "subsp. fraseri",
        str_detect(lp_subspecies, fixed("subsp. raphaeli"))    ~ "subsp. raphaeli",
        str_detect(lp_subspecies, fixed("subsp. pascullei"))   ~ "subsp. pascullei",
        TRUE ~ lp_subspecies
      )) %>%
      filter(!is.na(lp_subspecies))

    # 2. Rescale numeric columns
    target_cols = c('Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh')

    safe_scale <- function(x) {
      if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
        return(rep(0, length(x))) 
      } else {
        return(as.numeric(scale(x)))
      }
    }

    data <- data %>% mutate(across(all_of(target_cols), safe_scale))
    write_tsv(data, "scaled_subspecies_data.tsv")

    # Dynamic Class Prevalence Inspection
    class_counts <- table(data\$lp_subspecies)
    unique_classes <- names(class_counts[class_counts > 0])
    num_classes <- length(unique_classes)

    # Automatically set the absolute dominant category as the reference tier
    reference_class <- names(which.max(class_counts))
    
    # Format reference name strictly to s._[name] syntax
    s_ref_clean <- str_replace(reference_class, "subsp. ", "s._")

    # =================================================================
    # CASE 1: UNIFORM COHORT GUARD (Survive with valid mock structures)
    # =================================================================
    if (num_classes <= 1) {
      message("CRITICAL GUARD: Uniform cohort found (", reference_class, "). Writing fallback files.")
      
      fake_model <- list(status = "skipped_uniform_cohort", cohort = reference_class)
      saveRDS(fake_model, "firth_model.rds")
      
      fallback_results <- tibble(
        response_level = reference_class, term = "Uniform Cohort Grouping",
        coefficient = 0, odds_ratio = 1, or_conf_low = 1, or_conf_high = 1,
        std_dev = 0, wald_statistic = 0, p_value = 1,
        Decision = paste0("Uniform Cohort Detected: Entire dataset is ", reference_class)
      )
      write_tsv(fallback_results, "firth_multinomial_results.txt")
      
      output_file <- "firth_multinomial_formatted.txt"
      write_lines(c(
        "Firth bias-reduced multinomial logistic regression parameters for Legionella pneumophila subspecies classifications.",
        "====================================================================================================",
        "Predictor Gene\tOdds Ratio\tStd. Error\tZ-Statistic\tP-Value\t95% Confidence Interval\tDecision",
        paste0("Uniform Cohort: ", s_ref_clean, "\tN/A\tN/A\tN/A\tN/A\t[N/A, N/A]\tUniform Cohort Detected: No Comparative Variance"),
        "====================================================================================================",
        paste0("*Note: The baseline reference category assigned via max prevalence is ", reference_class, ".")
      ), file = output_file)
      
      quit(status = 0)
    }

    # =================================================================
    # STANDARD CASE: MULTI-CLASS FIRTH ESTIMATION RUN
    # =================================================================
    data <- data %>%
      mutate(lp_subspecies = factor(lp_subspecies, levels = c(reference_class, setdiff(unique_classes, reference_class))))

    model_formula = as.formula(paste("lp_subspecies ~", paste(target_cols, collapse = " + ")))

    model <- brglm2::brmultinom(
      model_formula,
      data = data,
      type = "AS_mean"
    )

    saveRDS(model, "firth_model.rds")

    # 4. Clean matrix statistics extraction
    model_summary <- summary(model)
    coef_matrix <- model_summary\$coefficients
    se_matrix   <- model_summary\$standard.errors
    
    if (is.vector(coef_matrix)) {
      # DYNAMIC MATRIX MAPPING FIX: Safely extracts target level labels without hardcoded index assumptions
      target_label <- setdiff(unique_classes, reference_class)[1]
      coef_matrix <- matrix(coef_matrix, nrow = 1, dimnames = list(target_label, names(coef_matrix)))
      se_matrix   <- matrix(se_matrix, nrow = 1, dimnames = list(target_label, names(se_matrix)))
    }

    wald_matrix <- coef_matrix / se_matrix
    p_matrix    <- 2 * (1 - pnorm(abs(wald_matrix)))
    response_levels <- rownames(coef_matrix)
    terms           <- colnames(coef_matrix)
    
    stats_list <- list()
    for (resp in response_levels) {
      for (trm in terms) {
        est <- coef_matrix[resp, trm]
        se  <- se_matrix[resp, trm]
        w_stat <- wald_matrix[resp, trm]
        p_val  <- p_matrix[resp, trm]
        conf_low  <- est - (1.96 * se)
        conf_high <- est + (1.96 * se)
        
        row_df <- tibble(
          response_level = resp, term = trm, coefficient = est, 
          std_dev = se, wald_statistic = w_stat, p_value = p_val, 
          conf_low = conf_low, conf_high = conf_high
        )
        stats_list[[length(stats_list) + 1]] = row_df
      }
    }
    
    stats_summary <- bind_rows(stats_list)
    
    original_table <- stats_summary %>% 
      mutate(
        odds_ratio = exp(coefficient), 
        or_conf_low = exp(conf_low), 
        or_conf_high = exp(conf_high), 
        Decision = if_else(p_value < 0.05, "statistically significant", "statistically null")
      ) %>% 
      select(response_level, term, coefficient, odds_ratio, or_conf_low, or_conf_high, std_dev, wald_statistic, p_value, Decision)
    
    write_tsv(original_table, "firth_multinomial_results.txt")
    
    # =================================================================
    # 5. DYNAMICALLY FORMATTED PUBLICATION DOCUMENT WRITER
    # =================================================================
    output_file <- "firth_multinomial_formatted.txt"
    
    # Write publication title boundaries
    write_lines(c(
      "Firth bias-reduced multinomial logistic regression parameters for Legionella pneumophila subspecies classifications.", 
      "====================================================================================================",
      "Predictor Gene\tOdds Ratio\tStd. Error\tZ-Statistic\tP-Value\t95% Confidence Interval\tDecision"
    ), file = output_file)

    # Process and append each comparison stack cleanly row-by-row onto disk
    for (resp in response_levels) {
      # Clean response labels to fit the 's._[name]' notation exactly
      s_target_clean <- str_replace(resp, "subsp. ", "s._")
      
      # Double check safety for accidental string duplicates (e.g. s_s_pneumophila)
      s_target_clean <- str_replace(s_target_clean, "s_s_", "s._")
      s_ref_clean    <- str_replace(s_ref_clean, "s_s_", "s._")
      
      comparison_header <- sprintf("Comparison: Subspecies %s vs. %s (Reference Tier)", s_target_clean, s_ref_clean)
      
      # Append stand-alone comparison partition row line
      write_lines(comparison_header, file = output_file, append = TRUE)
      
      # Extract only rows matching this specific target category split
      resp_data <- original_table %>% filter(response_level == resp)
      
      # Order and format rows matching your desired clean gene mapping parameters
      ordered_genes <- c("Gene_asd", "Gene_flaA", "Gene_mip", "Gene_pilE", "Gene_proA", "Gene_neuA_neuAh", "Gene_mompS")
      
      for (g_name in ordered_genes) {
        row_data <- resp_data %>% filter(term == g_name)
        if (nrow(row_data) == 0) next
        
        # Format names from 'Gene_flaA' down to 'flA' / 'flaA' style sheets
        clean_locus <- case_when(
          g_name == "Gene_asd"           ~ "asd",
          g_name == "Gene_flaA"          ~ "flaA",
          g_name == "Gene_mip"           ~ "MIP",
          g_name == "Gene_pilE"          ~ "pilE",
          g_name == "Gene_proA"          ~ "proA",
          g_name == "Gene_neuA_neuAh"    ~ "neuA.neuAh",
          g_name == "Gene_mompS"         ~ "mompS",
          TRUE                           ~ str_replace(g_name, "Gene_", "")
        )
        
        # Extract and protect numeric fields against overflow errors
        or_val   <- formatC(row_data\$odds_ratio, format = "f", digits = 4)
        se_val   <- formatC(row_data\$std_dev, format = "f", digits = 4)
        z_val    <- formatC(row_data\$wald_statistic, format = "f", digits = 3)
        p_val    <- formatC(row_data\$p_value, format = "f", digits = 4)
        dec_val  <- row_data\$Decision
        
        ci_low   <- row_data\$or_conf_low
        ci_high  <- row_data\$or_conf_high

        # SCIENTIFIC CAP ENFORCEMENT: Cleanly cuts off giant intervals for neat publication presentation
        ci_str   <- sprintf("[%s, %s]", 
          if_else(is.finite(ci_low) && ci_low >= 0.001, formatC(ci_low, format = "f", digits = 3), "0.000"),
          if_else(is.finite(ci_high) && ci_high <= 1000.0, formatC(ci_high, format = "f", digits = 3), ">1000")
        )
        
        # Combine parameters into a tab-separated text row
        gene_row_text <- sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s", clean_locus, or_val, se_val, z_val, p_val, ci_str, dec_val)
        write_lines(gene_row_text, file = output_file, append = TRUE)
      }
    }
    
    # Append bottom footnotes
    write_lines(c(
      "====================================================================================================", 
      paste0("*Note: The baseline reference category dynamically assigned for this execution framework is ", reference_class, ".")
    ), file = output_file, append = TRUE)
    """
}