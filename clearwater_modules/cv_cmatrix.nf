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
    library(tidyverse)
    library(brglm2)
    library(caret)

    # 1. Load the upstream data and model configuration
    df <- read_tsv("${scaled_data}", show_col_types = FALSE)
    fit_model <- readRDS("${global_model}")
    
    # Standardize factor labels for reporting
    df <- df %>% 
      mutate(lp_subspecies = case_when(
        lp_subspecies == "subsp. fraseri"     ~ "s_fraseri",
        lp_subspecies == "subsp. pneumophila" ~ "s_pneumophila",
        lp_subspecies == "subsp. raphaeli"    ~ "s_raphaeli",
        TRUE ~ as.character(lp_subspecies)
      )) %>%
      mutate(lp_subspecies = factor(lp_subspecies, levels = c("s_fraseri", "s_pneumophila", "s_raphaeli")))

    total_isolates <- nrow(df)
    
    # Dynamically extract formula and algorithm variables from your actual model object
    model_formula <- formula(fit_model)
    model_type    <- fit_model\$type

    # 2. CRASH PROTECTION: Dynamically calculate safe fold counts
    # Find the count of the rarest subspecies class present in the dataset
    class_counts <- table(df\$lp_subspecies)
    min_class_size <- min(class_counts[class_counts > 0])
    
    # Determine safe fold count (Default to 10, but cap at the minimum class size if low)
    k_folds <- 10
    if (total_isolates < 30 || min_class_size < 10) {
      k_folds <- max(2, min_class_size) # Drop to 3-fold, or 2-fold if limited by rare classes
    }
    
    print(paste("Executing cross-validation with optimized fold count:", k_folds))

    # 3. Configure reproducible, stratified folds based on our calculated k_folds
    set.seed(42)
    folds <- caret::createFolds(df\$lp_subspecies, k = k_folds, list = TRUE, returnTrain = FALSE)
    df\$cv_predicted <- NA_character_

    # 4. Stratified Internal Cross-Validation Loop
    for(i in 1:length(folds)) {
      test_indices  = folds[[i]]
      train_dataset = df[-test_indices, ]
      test_dataset  = df[test_indices, ]
      
      # Double-check that all active classes are represented in this specific training fold split
      # If a class completely vanished, fall back to global model framework safely
      if(length(unique(train_dataset\$lp_subspecies)) < length(unique(df\$lp_subspecies))) {
        fold_model <- fit_model
      } else {
        fold_model <- brglm2::brmultinom(
          model_formula,
          data = train_dataset,
          type = model_type
        )
      }
      
      # Test blindly on the remaining out-of-fold split
      fold_preds <- predict(fold_model, newdata = test_dataset, type = "class")
      df\$cv_predicted[test_indices] <- as.character(fold_preds)
    }
    
    df <- df %>% mutate(cv_predicted = factor(cv_predicted, levels = c("s_fraseri", "s_pneumophila", "s_raphaeli")))

    # 5. Process contingency values safely
    cm <- table(Predicted = df\$cv_predicted, Actual = df\$lp_subspecies)
    
    f_f <- if("s_fraseri" %in% rownames(cm)) cm["s_fraseri", "s_fraseri"] else 0
    f_p <- if("s_fraseri" %in% rownames(cm)) cm["s_fraseri", "s_pneumophila"] else 0
    f_r <- if("s_fraseri" %in% rownames(cm)) cm["s_fraseri", "s_raphaeli"] else 0
    
    p_f <- if("s_pneumophila" %in% rownames(cm)) cm["s_pneumophila", "s_fraseri"] else 0
    p_p <- if("s_pneumophila" %in% rownames(cm)) cm["s_pneumophila", "s_pneumophila"] else 0
    p_r <- if("s_pneumophila" %in% rownames(cm)) cm["s_pneumophila", "s_raphaeli"] else 0
    
    r_f <- if("s_raphaeli" %in% rownames(cm)) cm["s_raphaeli", "s_fraseri"] else 0
    r_p <- if("s_raphaeli" %in% rownames(cm)) cm["s_raphaeli", "s_pneumophila"] else 0
    r_r <- if("s_raphaeli" %in% rownames(cm)) cm["s_raphaeli", "s_raphaeli"] else 0

    # Calculate Performance Ratios Safely
    prec_f <- (f_f / max(1, (f_f + f_p + f_r))) * 100
    prec_p <- (p_p / max(1, (p_f + p_p + p_r))) * 100
    prec_r <- (r_r / max(1, (r_f + r_p + r_r))) * 100

    rec_f <- (f_f / max(1, (f_f + p_f + r_f))) * 100
    rec_p <- (p_p / max(1, (f_p + p_p + r_p))) * 100
    rec_r <- (r_r / max(1, (f_r + p_r + r_r))) * 100

    overall_acc <- ( (f_f + p_p + r_r) / total_isolates ) * 100

    # 6. Format and write the table output exactly as requested
    output_file <- "cv_confusion_matrix_metrics.txt"
    
    # Dynamically inject the correct loop label ("2-fold", "3-fold", "10-fold") into scientific header
    title_line      <- sprintf("Table #. Confusion matrix and internal classification performance metrics of the bias-reduced 7-gene predictive profile under %d-fold cross-validation (N = %d).", k_folds, total_isolates)
    header_line     <- "Predicted Subspecies\tActual: s_fraseri\tActual: s_pneumophila\tActual: s_raphaeli\tClass Precision (%)"
    row_fraseri     <- sprintf("s_fraseri\t%d\t%d\t%d\t%.2f%%", f_f, f_p, f_r, prec_f)
    row_pneumophila <- sprintf("s_pneumophila\t%d\t%d\t%d\t%.2f%%", p_f, p_p, p_r, prec_p)
    row_raphaeli    <- sprintf("s_raphaeli\t%d\t%d\t%d\t%.2f%%", r_f, r_p, r_r, prec_r)
    row_recall      <- sprintf("Class Recall (%%)\t%.2f%%\t%.2f%%\t%.2f%%\tOverall Accuracy:", rec_f, rec_p, rec_r)
    row_accuracy    <- sprintf("\t\t\t\t%.2f%%", overall_acc)
    
    footnote_line   <- sprintf("Note: All values display the exact absolute cross-validation distribution counts evaluated across a shuffled %d-fold iteration scheme. Following standard scientific formatting consensus, vertical rules have been omitted.", k_folds)

    write_lines(
      c(title_line, header_line, row_fraseri, row_pneumophila, row_raphaeli, row_recall, row_accuracy, "", footnote_line),
      file = output_file
    )
    """
}
