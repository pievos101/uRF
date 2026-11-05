#' Label a clustering forest for prediction on new subjects
#'
#' Assign the majority cluster label of the reference subjects in each terminal
#' node. The tree splits are unchanged. The returned copy is a standard ranger
#' classification or probability model and uses the existing predict methods.
#'
#' @param model A saved clustering forest, such as `result$model` from `uRF()`.
#' @param data Reference predictors, with rows in the same order as `clusters`.
#'   Normally the original training predictors. Every leaf must be reached by
#'   at least one reference subject.
#' @param clusters A vector of cluster labels, one per reference subject.
#' @param num.threads Number of threads for routing the reference subjects.
#' @param probability If TRUE, retain cluster proportions in each leaf and
#'   return a ranger probability model. Default FALSE keeps majority-label voting.
#' @return A ranger model. With probability = FALSE, response predictions are
#'   factors. With probability = TRUE, predictions are a subject-by-cluster
#'   matrix averaging leaf cluster proportions across trees. Use
#'   `colnames(p)[max.col(p, ties.method = "first")]` to obtain labels from that
#'   matrix. No trees are fitted and no clustering is rerun.
#' @details Ties within a leaf are resolved by the first factor level (factor
#'   levels are preserved; other labels are converted with `factor()`). Ties
#'   between trees follow ranger's existing prediction behavior and can be made
#'   reproducible with the `seed` argument to `predict()`.
#'
#'   These predictions approximate the original clustering. Evaluate agreement
#'   with the reference clusters; it need not be perfect. Stored OOB predictions
#'   and error metrics are cleared because they do not assess the new labels.
#' @export
assign_clusters <- function(model, data, clusters, num.threads = NULL,
                            probability = FALSE) {
  if (!is.logical(probability) || length(probability) != 1L || is.na(probability)) {
    stop("probability must be TRUE or FALSE.")
  }
  if (!inherits(model, "ranger") || is.null(model$forest) ||
      !identical(model$treetype, "Clustering") ||
      !identical(model$forest$treetype, "Clustering")) {
    stop("model must be a clustering ranger model with a saved forest.")
  }
  if (is.null(nrow(data)) || nrow(data) < 1L) {
    stop("data must contain at least one reference subject.")
  }
  if (!(is.factor(clusters) || is.character(clusters) ||
        is.numeric(clusters) || is.logical(clusters)) ||
      !is.null(dim(clusters)) || length(clusters) != nrow(data) ||
      anyNA(clusters) || (is.numeric(clusters) && any(!is.finite(clusters)))) {
    stop("clusters must be a label vector with one non-missing label per row of data.")
  }
  clusters <- droplevels(factor(clusters, ordered = FALSE))
  labels <- levels(clusters)
  codes <- as.integer(clusters)
  nodes <- predict(model, data, type = "terminalNodes",
                   num.threads = num.threads)$predictions
  nodes <- matrix(nodes, nrow = nrow(data), ncol = model$num.trees)
  if (probability) {
    model$forest$terminal.class.counts <- vector("list", model$num.trees)
  }

  for (tree in seq_len(model$num.trees)) {
    children <- model$forest$child.nodeIDs[[tree]]
    leaves <- which(children[[1]] == 0 & children[[2]] == 0)
    if (probability) {
      # ranger indexes the probability vectors by node ID, including internal nodes.
      model$forest$terminal.class.counts[[tree]] <-
        rep(list(numeric()), length(children[[1]]))
    }
    for (leaf in leaves) {
      # Native node IDs are zero-based; R vector indices are one-based.
      members <- codes[nodes[, tree] == leaf - 1L]
      if (!length(members)) {
        stop("data does not reach every terminal node; supply the original training predictors.")
      }
      counts <- tabulate(members, nbins = length(labels))
      model$forest$split.values[[tree]][leaf] <- which.max(counts)
      if (probability) {
        model$forest$terminal.class.counts[[tree]][[leaf]] <- counts / sum(counts)
      }
    }
  }

  model$forest$class.values <- seq_along(labels)
  model$forest$levels <- labels
  model$treetype <- if (probability) "Probability estimation" else "Classification"
  model$forest$treetype <- model$treetype
  model$predictions <- NULL
  model$prediction.error <- NA_real_
  model$confusion.matrix <- NULL
  model$call <- match.call()
  model
}
