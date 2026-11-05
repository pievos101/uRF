\name{assign_clusters}
\alias{assign_clusters}
\title{Label a clustering forest for prediction on new subjects}
\usage{
assign_clusters(model, data, clusters, num.threads = NULL, probability = FALSE)
}
\arguments{
\item{model}{A saved clustering ranger model, such as \code{result$model} from \code{uRF()}.}
\item{data}{Reference predictors, normally the original training data, with rows in the same order as \code{clusters}. Every terminal node must be reached by at least one reference subject.}
\item{clusters}{A vector of cluster labels, one per reference subject.}
\item{num.threads}{Number of threads for routing reference subjects.}
\item{probability}{If TRUE, store cluster proportions within each leaf and return a ranger probability model. Default FALSE retains majority-label voting.}
}
\description{
Labels each terminal node with the majority cluster among reference subjects
reaching it, or stores leaf cluster proportions when \code{probability = TRUE}.
Returns a copy as a standard ranger classification or probability model.
Tree splits are preserved; no trees are fitted and no clustering is rerun.
Use the existing \code{predict()} methods to assign new subjects.
}
\details{
Ties within a leaf use the first factor level. Supplied factor levels are
preserved after dropping unused levels; other label vectors are converted with
\code{factor()}. Ties between trees use ranger's existing behavior; use
\code{seed} in \code{predict()} for reproducibility.

Predicted assignments approximate the original clustering and may disagree with
reference labels. Stored OOB predictions and error metrics are cleared because
they do not assess the new labels. Impurity importance still describes the
original forest's splits, not a supervised fit to the cluster labels.

With \code{probability = TRUE}, predictions average the cluster proportions in
the reached leaves across trees. They are not majority-vote fractions or
calibrated probabilities. To obtain labels from a probability matrix \code{p},
use \code{colnames(p)[max.col(p, ties.method = "first")]}. This label can differ
from majority voting across leaves; the first column wins a probability tie.

For multiple views, \code{uRF()} currently returns only the last view's forest
in \code{model}. Supply that view's predictors. Predictions then use that view alone.
}
\value{A ranger model. With \code{probability = FALSE}, response predictions are factors containing the cluster labels. With \code{probability = TRUE}, predictions are a subject-by-cluster probability matrix whose rows sum to one. The original model is unchanged.}
\examples{
set.seed(42)
x <- iris[, 1:4]
result <- uRF(list(x), k = 3, ntrees = 10)
assignment_model <- assign_clusters(result$model, x, result$clusters)
predict(assignment_model, x[1, , drop = FALSE], seed = 42)$predictions
probability_model <- assign_clusters(result$model, x, result$clusters,
                                     probability = TRUE)
p <- predict(probability_model, x[1, , drop = FALSE])$predictions
p
colnames(p)[max.col(p, ties.method = "first")]
}
