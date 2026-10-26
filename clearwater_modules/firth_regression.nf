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
    path "firth_model.rds",                 emit: model_object // NEW TRACK

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(brglm2)

    # 1. Load and clean the global data
    data <- read_tsv("${imputed_table}", show_col_types = FALSE)
    data <- data %>% dplyr::rename(lp_subspecies = subspeciesID_kraken)

    data <- data %>%
      mutate(lp_subspecies = case_when(
        str_detect(lp_subspecies, fixed("subsp. pneumophila")) ~ "subsp. pneumophila",
        str_detect(lp_subspecies, fixed("subsp. fraseri"))     ~ "subsp. fraseri",
        str_detect(lp_subspecies, fixed("subsp. raphaeli"))    ~ "subsp. raphaeli",
        TRUE ~ lp_subspecies
      ))

    data <- data %>%
      mutate(lp_subspecies = factor(lp_subspecies)) %>%
      mutate(lp_subspecies = relevel(lp_subspecies, ref = "subsp. pneumophila"))

    # 2. Rescale numeric columns
    target_cols = c('Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh')

    safe_scale <- function(x) {
      if (sd(x, na.rm = TRUE) == 0) {
        return(rep(0, length(x))) 
      } else {
        return(as.numeric(scale(x)))
      }
    }

    data <- data %>% mutate(across(all_of(target_cols), safe_scale))
    write_tsv(data, "scaled_subspecies_data.tsv")

    # 3. Fit Firth bias-reduced multinomial logistic regression
    model_formula = as.formula(paste("lp_subspecies ~", paste(target_cols, collapse = " + ")))

    model <- brglm2::brmultinom(
      model_formula,
      data = data,
      type = "AS_mean"
    )

    # EXPLICITLY EXPORT THE TRAINED MODEL OBJECT FOR DOWNSTREAM USE
    saveRDS(model, "firth_model.rds")

    # 4. Clean matrix statistics extraction (rest of your original script remains here...)
    model_summary <- summary(model)
    coef_matrix <- model_summary\$coefficients
    se_matrix   <- model_summary\$standard.errors
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
        row_df <- tibble(response_level = resp, term = trm, coefficient = est, std_dev = se, wald_statistic = w_stat, p_value = p_val, conf_low = conf_low, conf_high = conf_high)
        stats_list[[length(stats_list) + 1]] = row_df
      }
    }
    stats_summary <- bind_rows(stats_list)
    original_table <- stats_summary %>% mutate(odds_ratio = exp(coefficient), or_conf_low = exp(conf_low), or_conf_high = exp(conf_high), Decision = if_else(p_value < 0.05, "statistically significant", "statistically null")) %>% select(response_level, term, coefficient, odds_ratio, or_conf_low, or_conf_high, std_dev, wald_statistic, p_value, Decision)
    write_tsv(original_table, "firth_multinomial_results.txt")
    
    # 5. Format graphic document output...
    formatted_table <- original_table %>% mutate(response_level = str_replace(response_level, "subsp. ", "s_"), term = str_replace(term, "Gene_", ""), term = case_when(term == "neuA_neuAH" ~ "neuA.neuAh", term == "flaA" ~ "flA", term == "mip" ~ "MIP", TRUE ~ term)) %>% filter(term != "(Intercept)") %>% rowwise() %>% mutate(odds_ratio = formatC(odds_ratio, format = "f", digits = 4), confidence_interval = sprintf("[%s, %s]", formatC(or_conf_low, format = "f", digits = 3), formatC(or_conf_high, format = "f", digits = 3)), wald_statistic = formatC(wald_statistic, format = "f", digits = 3), p_value = formatC(p_value, format = "f", digits = 4)) %>% ungroup() %>% select(SubspeciesComparison = response_level, PredictorLocus = term, OddsRatio = odds_ratio, ConfidenceInterval = confidence_interval, WaldZ = wald_statistic, PValue = p_value, Decision) %>% group_by(SubspeciesComparison) %>% mutate(SubspeciesComparison = if_else(row_number() == 1, SubspeciesComparison, "")) %>% ungroup()
    output_file <- "firth_multinomial_formatted.txt"
    write_lines(c("Multinomial logistic regression analysis (Firth's bias-reduction penalized likelihood) evaluating 7-gene allele counts as predictors of Legionella pneumophila subspecies classification.", "===================================================================================================="), file = output_file)
    write_delim(formatted_table, file = output_file, delim = "\t", col_names = TRUE, append = TRUE)
    write_lines(c("====================================================================================================", "*Note: The baseline reference category for the multinomial regression framework is s_pneumophila."), file = output_file, append = TRUE)
    """
}
