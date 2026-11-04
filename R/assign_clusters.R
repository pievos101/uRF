#' Label a clustering forest for prediction on new subjects
#'
#' Assign the majority cluster label of the reference subjects in each terminal
#' node. The tree splits are unchanged. The returned copy is a standard ranger
#' classification model and uses the existing predict methods.
#'
#' @param model A saved clustering forest, such as `result$model` from `uRF()`.
#' @param data Reference predictors, with rows in the same order as `clusters`.
#'   Normally the original training predictors. Every leaf must be reached by
#'   at least one reference subject.
#' @param clusters A vector of cluster labels, one per reference subject.
#' @param num.threads Number of threads for routing the reference subjects.
#' @return A ranger classification model. Response predictions are factors with
#'   the supplied cluster labels. No trees are fitted and no clustering is rerun.
#' @details Ties within a leaf are resolved by the first factor level (factor
#'   levels are preserved; other labels are converted with `factor()`). Ties
#'   between trees follow ranger's existing prediction behavior and can be made
#'   reproducible with the `seed` argument to `predict()`.
#'
#'   These predictions approximate the original clustering. Evaluate agreement
#'   with the reference clusters; it need not be perfect. Stored OOB predictions
#'   and error metrics are cleared because they do not assess the new labels.
#' @export
assign_clusters <- function(model, data, clusters, num.threads = NULL) {
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

  for (tree in seq_len(model$num.trees)) {
    children <- model$forest$child.nodeIDs[[tree]]
    leaves <- which(children[[1]] == 0 & children[[2]] == 0)
    for (leaf in leaves) {
      # Native node IDs are zero-based; R vector indices are one-based.
      members <- codes[nodes[, tree] == leaf - 1L]
      if (!length(members)) {
        stop("data does not reach every terminal node; supply the original training predictors.")
      }
      counts <- tabulate(members, nbins = length(labels))
      model$forest$split.values[[tree]][leaf] <- which.max(counts)
    }
  }

  model$forest$class.values <- seq_along(labels)
  model$forest$levels <- labels
  model$forest$treetype <- "Classification"
  model$treetype <- "Classification"
  model$predictions <- NULL
  model$prediction.error <- NA_real_
  model$confusion.matrix <- NULL
  model$call <- match.call()
  model
}
