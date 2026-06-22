#' Stratify cases based on input variables
#'
#' Compute case subgroups by stratification based on multiple input variables.
#' This function creates strata of cases based on continuous and/or categorical
#' variables, with the option of matching cases to these variables.
#'
#' @export
stratify_control <- function(pheno, strat_cont = NULL, strat_cat = NULL, strat_surv = NULL, K = 5,
                             shared_control = FALSE) {

  # Check input data
  #check_stratify_inputs(pheno, strat_cont, strat_cat, K)

  # Re-format the survival matrix (FID, IID, Time, Event)
  if (!is.null(strat_surv)) {
    # Match
    if(!is.null(strat_cont)){
      strat_surv <- strat_surv[match(strat_cont[, 2], strat_surv[, 2]), ]
    } else {
      strat_cont <- strat_surv[, c(1, 2)]
    }

    # Convert
    event_times <- strat_surv[which(strat_surv[, 4] == 1), 3]
    event_times_ctrls <- strat_surv[which(strat_surv[, 4] == 0), 3]

    # Sample censoring time from event times, and code NA if already censored
    controls <- which(strat_surv[, 4] == 0)
    event_times_ctrls_fake <- sample(event_times, length(controls))
    event_times_ctrls_fake[which(event_times_ctrls_fake > event_times_ctrls)] <- NA

    add_cont <- rep(NA, nrow(strat_cont))
    add_cont[which(strat_surv[, 4] == 1)] <- event_times
    add_cont[which(strat_surv[, 4] == 0)] <- event_times_ctrls_fake

    strat_cont <- cbind(strat_cont, add_cont)
  }

  # Store IDs
  ids <- pheno[, 2]
  control_ids <- pheno[which(pheno[, 3] == 0), 2]

  # Convert to data.frame
  pheno <- as.data.frame(pheno)

  if (!is.null(strat_cont)) {
    strat_cont <- as.data.frame(strat_cont)
    strat_cont <- strat_cont[match(ids, strat_cont[, 2]), ]
    strat_cont[, 1] <- strat_cont[, 2] <- ids
  }
  if (!is.null(strat_cat)) {
    strat_cat <- as.data.frame(strat_cat)
    strat_cat <- strat_cat[match(ids, strat_cat[, 2]), ]
    strat_cat[, 1] <- strat_cat[, 2] <- ids
  }

  # Identify cases
  cases <- pheno[which(pheno[, 3] == 1), 2]

  # Initialize for tracking all strata
  all_strat_info <- list()
  strata <- list()
  K_total <- 0
  all_cases_nostrat <- c()
  var_info <- NULL
  ctrl_rows_matched <- NULL

  # Process continuous stratification variables
  if (!is.null(strat_cont)) {
    n_cont <- ncol(strat_cont) - 2
    ctrl_rows_matched <- NULL

    for(i in 1:n_cont) {
      # Extract column
      strat_col <- strat_cont[, 2 + i]

      # Identify cases with missing stratification
      cases_nostrat <- strat_cont[which(strat_cont[, 2] %in% cases & is.na(strat_col)), 2, drop = TRUE]
      all_cases_nostrat <- c(all_cases_nostrat, cases_nostrat)

      # Extract stratification variable for cases only
      strat_cases <- strat_cont[which(pheno[, 3] == 1 & !(pheno[, 2] %in% cases_nostrat)), c(1, 2, 2 + i)]
      colnames(strat_cases) <- c("FID", "IID", "strat_val")
      strat_cases_vals <- strat_cases[, 3]

      # Define groups
      K_use <- K
      K_total <- K_total + K_use
      result <- assign_to_quantiles(strat_cases_vals, K)
      strat_cases$groups <- result$groups
      strat_cases$order <- result$order
      breaks <- result$breaks

      if (shared_control == TRUE) {
        # Also stratify controls
        ctrl_rows <- strat_cont[which(pheno[, 3] == 0), c(1, 2, 2 + i)]
        colnames(ctrl_rows) <- c("FID", "IID", "strat_val")
        ctrl_rows <- ctrl_rows[!is.na(ctrl_rows$strat_val), ]

        ctrl_groups <- cut(ctrl_rows$strat_val,
                   breaks         = breaks,
                   labels         = 1:K,
                   include.lowest = TRUE)
        if (any(table(ctrl_groups) < 100)) 
            stop("Too few controls to match controls on stratification variable.")

        ctrl_rows$groups <- as.integer(ctrl_groups)
        ctrl_rows$order  <- NA
        ctrl_rows_matched <- ctrl_rows[!is.na(ctrl_rows$groups), ]
      }

      # Store information
      var_info <- rbind(var_info, c(colnames(strat_cont)[i+2], "continuous", K))

      # Store stratification info
      all_strat_info[[length(all_strat_info) + 1]] <- list(
        data = strat_cases,
        data_control = ctrl_rows_matched,
        original = strat_cont[, c(1, 2, 2 + i)],
        type = "continuous",
        K = K_use,
        cases_nostrat = cases_nostrat,
        var_index = i,
        group_start = K_total - K_use + 1,
        group_end = K_total
      )
    }
  }
  
  # Process binary/categorical stratification variables
  if(!is.null(strat_cat)) {
    n_bin <- ncol(strat_cat) - 2
    ctrl_rows_matched_cat <- NULL

    for (i in 1:n_bin) {
      # Extract column
      strat_col <- strat_cat[, 2 + i]
      
      # Identify cases with missing stratification
      cases_nostrat <- strat_cat[which(strat_cat[, 2] %in% cases & is.na(strat_col)), 2]
      all_cases_nostrat <- c(all_cases_nostrat, cases_nostrat)
      
      # Extract stratification variable for cases only
      strat_cases <- strat_cat[which(pheno[, 3] == 1 & !(pheno[, 2] %in% cases_nostrat)), c(1, 2, 2 + i)]
      colnames(strat_cases) <- c("FID", "IID", "strat_val")
      strat_cases_vals <- strat_cases[, 3]
      
      # Get unique categories
      categories <- sort(unique(strat_cases_vals[!is.na(strat_cases_vals)]))
      C <- length(categories)

      message(sprintf("Categorical variable %d: %d categories detected", i, C))

      # Assign groups based on categories
      K_total <- K_total + C
      strat_cases$groups <- match(strat_cases_vals, categories)
      strat_cases$order <- NA

      if (shared_control == TRUE) {
        # Also stratify controls
        ctrl_rows_cat <- strat_cat[which(pheno[, 3] == 0), c(1, 2, 2 + i)]
        colnames(ctrl_rows_cat) <- c("FID", "IID", "strat_val")

        ctrl_rows_cat <- ctrl_rows_cat[!is.na(ctrl_rows_cat$strat_val), ]
        ctrl_rows_cat$groups <- match(ctrl_rows_cat$strat_val, categories)

        if (any(table(ctrl_rows_cat$groups) < 100)) 
            stop("Too few controls to match controls on categorical stratification variable.")

        ctrl_rows_cat$order  <- NA
        ctrl_rows_matched_cat <- ctrl_rows_cat[!is.na(ctrl_rows_cat$groups), ]
      }

      # Store information
      var_info <- rbind(var_info, c(colnames(strat_cat)[i + 2], "categorical", C))

      # Store stratification info
      all_strat_info[[length(all_strat_info) + 1]] <- list(
        data = strat_cases,
        data_control = ctrl_rows_matched_cat,
        type = "categorical",
        K = C,
        categories = categories,
        cases_nostrat = cases_nostrat,
        var_index = i,
        group_start = K_total - C + 1,
        group_end = K_total
      )
    }
  }
  
  # Check that we have at least one stratification variable
  if(length(all_strat_info) == 0) {
    stop("No stratification variables provided. Must provide strat_cont or strat_cat")
  }
  
  # Create a combined info dataframe for backward compatibility
  # This combines all stratification info, but individuals can appear multiple times
  all_cases_combined <- do.call(rbind, lapply(all_strat_info, function(x) {
    df <- x$data
    # Adjust group numbers to be global
    df$groups <- df$groups + x$group_start - 1
    df
  }))

  if (shared_control) {
    all_controls_combined <- do.call(rbind, lapply(all_strat_info, function(x) {
        df <- x$data_control
        # Adjust group numbers to be global
        df$groups <- df$groups + x$group_start - 1
        df
    }))
  }
  
  # Create a multivariate dataframe for linear regression + SumHer
  multi <- cbind(ids, ids, pheno[, 3])
  colnames(multi) <- c("FID", "IID", "pheno")

  for (k in 1:length(all_strat_info)) {
    var_add  <- NULL

    if (shared_control) {
        for (j in all_strat_info[[k]]$group_start:all_strat_info[[k]]$group_end) {
            # Find the matched controls
            add <- rep(NA, length(ids))
            add[which(ids %in% all_controls_combined[all_controls_combined$groups == j, 2])] <- 0
            add[which(ids %in% all_cases_combined[all_cases_combined$groups == j, 2])] <- 1
            var_add <- cbind(var_add, add)
        }
    } else {
        for (j in all_strat_info[[k]]$group_start:all_strat_info[[k]]$group_end) {
            add <- rep(NA, length(ids))
            add[which(ids %in% control_ids)] <- 0
            add[which(ids %in% all_cases_combined[all_cases_combined$groups == j, 2])] <- 1
            var_add <- cbind(var_add, add)
        }
    }

    colnames(var_add) <- paste0(all_strat_info[[k]]$type, "_", all_strat_info[[k]]$var_index, "_", 1:all_strat_info[[k]]$K)

    # In the case of a continuous variable, add this as continuous variable
    if(all_strat_info[[k]]$type == "continuous"){
      add_cont <- rep(NA, length(ids))
      add_cont <- all_strat_info[[k]]$original[match(ids, all_strat_info[[k]]$original[, 2]), 3]
      var_add <- cbind(var_add, add_cont)

      colnames(var_add)[ncol(var_add)] <- paste0(all_strat_info[[k]]$type, "_", all_strat_info[[k]]$var_index)
    }

    multi <- cbind(multi, var_add)
  }

  cat("\n")

  # Return list with information (matching original structure)
  strata[["K"]] <- K
  strata[["y"]] <- pheno
  strata[["multi"]] <- multi
  strata[["info"]] <- all_cases_combined
  strata[["ids"]] <- ids
  strata[["strat_miss"]] <- all_cases_nostrat
  strata[["strat_details"]] <- all_strat_info  # Detailed info about each variable
  strata[["n_cont"]] <- if(!is.null(strat_cont)) ncol(strat_cont) - 2 else 0
  strata[["n_bin"]] <- if(!is.null(strat_cat)) ncol(strat_cat) - 2 else 0
  strata[["var_info"]] <- var_info
  strata[["shared_control"]] <- shared_control
  return(strata)
}