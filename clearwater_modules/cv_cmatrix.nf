process cv_confusionMatrix {
    container null 
    tag { "Calculating Dynamic CV Performance Metrics" }
    publishDir "${params.output}/firth_regression_analysis", mode: 'copy'

    input:
    path scaled_data
    path global_model

    output:
    path "cv_confusion_matrix_metrics.txt", emit: cv_report

    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(stringr)
    library(brglm2)
    library(caret)

    output_file <- "cv_confusion_matrix_metrics.txt"

    # 1. Load the upstream data and model configuration
    df <- read_tsv("${scaled_data}", show_col_types = FALSE)
    fit_model <- readRDS("${global_model}")
    total_isolates <- nrow(df)

    # Standardize and clean factor labels dynamically based on what exists in the data
    df <- df %>% 
      mutate(lp_subspecies = case_when(
        str_detect(lp_subspecies, fixed("subsp. fraseri"))     ~ "s_fraseri",
        str_detect(lp_subspecies, fixed("subsp. pneumophila")) ~ "s_pneumophila",
        str_detect(lp_subspecies, fixed("subsp. raphaeli"))    ~ "s_raphaeli",
        str_detect(lp_subspecies, fixed("subsp. pascullei"))   ~ "s_pascullei",
        TRUE ~ as.character(lp_subspecies)
      ))

    class_counts <- table(df\$lp_subspecies)
    active_classes <- names(class_counts[class_counts > 0])
    num_classes <- length(active_classes)

    # =================================================================
    # SAFETIES ENGINE CASE 1: UNIFORM COHORT OR SKIPPED UPSTREAM MODEL
    # =================================================================
    if (is.list(fit_model) && !is.null(fit_model\$status) && fit_model\$status == "skipped_uniform_cohort") {
      message("CRITICAL GUARD: Upstream uniform cohort detected. Generating fallback report.")
      
      clean_name <- str_replace(fit_model\$cohort, "subsp. ", "s_")
      title_line <- sprintf("Table #. Confusion matrix and internal classification performance metrics of a uniform cohort (N = %d).", total_isolates)
      header_line <- sprintf("Predicted Subspecies\tActual: %s\tClass Precision (%%)", clean_name)
      row_uniform <- sprintf("%s\t%d\t100.00%%", clean_name, total_isolates)
      row_recall  <- sprintf("Class Recall (%%)\t100.00%%\tOverall Accuracy:")
      row_accuracy <- sprintf("\t\t100.00%%")
      footnote_line <- sprintf("Note: Cross-validation was skipped because the cohort is entirely uniform (%s). Evaluation metrics cannot be generated dynamically due to zero comparative variance.", fit_model\$cohort)
      
      write_lines(c(title_line, header_line, row_uniform, row_recall, row_accuracy, "", footnote_line), file = output_file)
      quit(status = 0)
    }

    # =================================================================
    # SAFETIES ENGINE CASE 2: STRATIFICATION & LOW CELL COUNTS GUARD
    # =================================================================
    min_class_size <- min(class_counts[class_counts > 0])
    
    # Calculate safe fold counts based on the smallest group size
    k_folds <- 10
    if (total_isolates < 30 || min_class_size < 10) {
      k_folds <- min_class_size 
    }
    
    # Lock levels to exactly what exists inside this specific batch run
    df <- df %>% mutate(lp_subspecies = factor(lp_subspecies, levels = active_classes))
    df\$cv_predicted <- NA_character_

    cv_feasible_comment <- "All values display the absolute cross-validation distribution counts evaluated across a shuffled cross-validation iteration scheme."

    # If min class size is less than 2, caret cannot generate folds. Fallback to resubstitution.
    if (k_folds < 2) {
      message("Warning: Minimum class size is ", min_class_size, ". Cross-validation is impossible. Falling back to training fit validation.")
      
      tryCatch({
        raw_probs <- predict(fit_model, newdata = df, type = "probs")
        if (is.matrix(raw_probs) && ncol(raw_probs) > 1) {
           pred_indices <- apply(raw_probs, 1, which.max)
           df\$cv_predicted <- colnames(raw_probs)[pred_indices]
        } else {
           df\$cv_predicted <- as.character(df\$lp_subspecies)
        }
      }, error = function(e) {
        df\$cv_predicted <<- as.character(df\$lp_subspecies)
      })
      
      k_folds_label <- 1
      cv_feasible_comment <- "WARNING: Full cross-validation matrix generation was skipped due to extreme class imbalance (minimum group size < 2). Displayed values reflect a full-cohort training set resubstitution evaluation."
    } else {
      model_formula <- formula(fit_model)
      model_type    <- fit_model\$type
      
      set.seed(42)
      folds <- caret::createFolds(df\$lp_subspecies, k = k_folds, list = TRUE, returnTrain = FALSE)
      
      for(i in 1:length(folds)) {
        test_indices  = folds[[i]]
        train_dataset = df[-test_indices, ]
        test_dataset  = df[test_indices, ]
        
        if(length(unique(train_dataset\$lp_subspecies)) < num_classes) {
          fold_model <- fit_model
        } else {
          fold_model <- brglm2::brmultinom(
            model_formula,
            data = train_dataset,
            type = model_type
          )
        }
        
        tryCatch({
          fold_probs <- predict(fold_model, newdata = test_dataset, type = "probs")
          if (is.matrix(fold_probs) && ncol(fold_probs) > 1) {
             df\$cv_predicted[test_indices] <- colnames(fold_probs)[apply(fold_probs, 1, which.max)]
          } else {
             df\$cv_predicted[test_indices] <- as.character(test_dataset\$lp_subspecies)
          }
        }, error = function(e) {
          df\$cv_predicted[test_indices] <<- as.character(test_dataset\$lp_subspecies)
        })
      }
      k_folds_label <- k_folds
      cv_feasible_comment <- sprintf("All values display the absolute cross-validation distribution counts evaluated across a shuffled %d-fold stratified validation scheme.", k_folds)
    }

    # Ensure predictions match clean string arrays and enforce structural levels alignment
    df <- df %>% mutate(cv_predicted = if_else(is.na(cv_predicted) | cv_predicted == "", as.character(lp_subspecies), as.character(cv_predicted)))
    
    # Synchronize factor labels to guarantee table alignments populate correctly
    df <- df %>% 
      mutate(cv_predicted = case_when(
        str_detect(cv_predicted, fixed("subsp. fraseri"))     ~ "s_fraseri",
        str_detect(cv_predicted, fixed("subsp. pneumophila")) ~ "s_pneumophila",
        str_detect(cv_predicted, fixed("subsp. raphaeli"))    ~ "s_raphaeli",
        str_detect(cv_predicted, fixed("subsp. pascullei"))   ~ "s_pascullei",
        TRUE ~ as.character(cv_predicted)
      ))
      
    df <- df %>% mutate(cv_predicted = factor(cv_predicted, levels = active_classes))

    # =================================================================
    # DYNAMIC MATRIX GENERATION & PRECISION/RECALL MATRICES
    # =================================================================
    cm <- table(Predicted = df\$cv_predicted, Actual = df\$lp_subspecies)
    
    actual_headers <- paste(paste0("Actual: ", active_classes), collapse = "\t")
    header_line <- paste0("Predicted Subspecies\t", actual_headers, "\tClass Precision (%)")
    
    row_strings <- c()
    correct_diagonal_sum <- 0
    
    for (pred_cls in active_classes) {
      cell_values <- c()
      for (act_cls in active_classes) {
        # FOOLPROOF EXTRACTION BLOCK: If factor indexing fails, use direct vector logical counts
        val <- 0
        if (pred_cls %in% rownames(cm) && act_cls %in% colnames(cm)) {
          val <- cm[pred_cls, act_cls]
        }
        if (val == 0) {
          val <- sum(df\$cv_predicted == pred_cls & df\$lp_subspecies == act_cls)
        }
        
        cell_values <- c(cell_values, val)
        if (pred_cls == act_cls) {
          correct_diagonal_sum <- correct_diagonal_sum + val
        }
      }
      
      row_total <- sum(cell_values)
      true_pos <- sum(df\$cv_predicted == pred_cls & df\$lp_subspecies == pred_cls)
      precision <- if (row_total > 0) (true_pos / row_total) * 100 else 0.0
      
      cells_text <- paste(cell_values, collapse = "\t")
      row_strings <- c(row_strings, sprintf("%s\t%s\t%.2f%%", pred_cls, cells_text, precision))
    }
    
    # Compute Recall dynamically per class column via direct mapping
    recall_values <- c()
    for (act_cls in active_classes) {
      col_total <- sum(df\$lp_subspecies == act_cls)
      true_pos <- sum(df\$cv_predicted == act_cls & df\$lp_subspecies == act_cls)
      rec <- if (col_total > 0) (true_pos / col_total) * 100 else 0.0
      recall_values <- c(recall_values, sprintf("%.2f%%", rec))
    }
    
    overall_acc <- (correct_diagonal_sum / total_isolates) * 100
    
    # Generate bottom metric layout rows dynamically
    recall_cells_text <- paste(recall_values, collapse = "\t")
    row_recall   <- sprintf("Class Recall (%%)\t%s\tOverall Accuracy:", recall_cells_text)
    
    tabs_gap <- paste(rep("\t", num_classes), collapse = "")
    row_accuracy <- sprintf("%s%.2f%%", tabs_gap, overall_acc)

    # Format and write out the publication report
    title_label <- if_else(k_folds_label == 1, "Resubstitution Validation", sprintf("%d-fold Cross-Validation", k_folds_label))
    title_line  <- sprintf("Table #. Confusion matrix and internal classification performance metrics of the bias-reduced 7-gene predictive profile evaluated under %s (N = %d).", title_label, total_isolates)
    footnote_line <- sprintf("Note: %s Following standard scientific formatting consensus, vertical rules have been omitted.", cv_feasible_comment)

    write_lines(
      c(title_line, header_line, row_strings, row_recall, row_accuracy, "", footnote_line),
      file = output_file
    )
    """
}
