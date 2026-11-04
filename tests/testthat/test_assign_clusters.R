library(uRF)

test_that("leaf labels use reference majorities through standard prediction", {
  x <- iris[, 1:4]
  original <- ranger(x = x, clustering = TRUE, num.trees = 7,
                     seed = 42, num.threads = 1, min.bucket = 5)
  snapshot <- serialize(original, NULL)
  labels <- factor(rep(c("group B", "group A"), length.out = nrow(x)),
                   levels = c("group B", "group A", "unused"))
  nodes <- predict(original, x, type = "terminalNodes")$predictions
  assigned <- assign_clusters(original, x, labels)
  expect_identical(serialize(original, NULL), snapshot)
  expect_s3_class(assigned, "ranger")
  expect_identical(assigned$treetype, "Classification")
  expect_identical(assigned$forest$levels, c("group B", "group A"))
  expect_null(assigned$predictions)
  expect_true(is.na(assigned$prediction.error))
  expect_identical(assigned$forest$child.nodeIDs, original$forest$child.nodeIDs)
  expect_identical(assigned$forest$split.varIDs, original$forest$split.varIDs)
  expect_equal(predict(assigned, x, type = "terminalNodes")$predictions, nodes)

  expected <- matrix(NA_integer_, nrow(x), original$num.trees)
  for (tree in seq_len(original$num.trees)) {
    for (id in unique(nodes[, tree])) {
      members <- nodes[, tree] == id
      expected[members, tree] <- which.max(table(labels[members])[1:2])
    }
  }
  votes <- predict(assigned, x, predict.all = TRUE)$predictions
  expect_equal(votes, expected)
  majority <- apply(expected, 1, function(v) which.max(tabulate(v, nbins = 2)))
  response <- predict(assigned, x, seed = 42)$predictions
  expect_identical(as.character(response), levels(labels)[majority])
  expect_identical(predict(assigned$forest, x, seed = 42)$predictions, response)
  expect_identical(predict(assigned, x[1, , drop = FALSE], seed = 42)$predictions,
                   response[1])
  expect_equal(predict(assigned, x[, 4:1], seed = 42)$predictions, response)
  restored <- unserialize(serialize(assigned, NULL))
  expect_identical(predict(restored, x, seed = 42)$predictions, response)
})

test_that("README clusters and numeric labels can label a forest", {
  x <- iris[, 1:4]
  result <- uRF(list(x), k = 3, ntrees = 5)
  assigned <- assign_clusters(result$model, x, result$clusters * 10)
  predicted <- predict(assigned, x[1:3, ], seed = 42)$predictions
  expect_s3_class(predicted, "factor")
  expect_length(predicted, 3)
  expect_true(all(as.character(predicted) %in% c("10", "20", "30")))
})

test_that("a single-node forest resolves leaf ties by factor level", {
  x <- data.frame(a = c(1, 1), b = c(2, 2))
  forest <- ranger(x = x, clustering = TRUE, num.trees = 1, min.node.size = 2)
  assigned <- assign_clusters(forest, x, factor(c("b", "a"), levels = c("b", "a")))
  expect_identical(as.character(predict(assigned, x)$predictions), c("b", "b"))
})

test_that("invalid labels and uncovered leaves fail clearly", {
  x <- iris[, 1:4]
  forest <- ranger(x = x, clustering = TRUE, num.trees = 3, seed = 42)
  expect_error(assign_clusters(forest, x, 1:2), "one non-missing label")
  expect_error(assign_clusters(forest, x, rep(NA, 150)), "one non-missing label")
  expect_error(assign_clusters(forest, x, rep(Inf, 150)), "one non-missing label")
  expect_error(assign_clusters(forest, x, matrix(1, 150, 1)), "one non-missing label")
  expect_error(assign_clusters(forest, x[FALSE, ], integer()), "at least one")
  expect_error(assign_clusters(forest, x[1, , drop = FALSE], "a"), "every terminal node")
  expect_error(assign_clusters(list(), x, rep(1, 150)), "clustering ranger model")
  forest$forest <- NULL
  expect_error(assign_clusters(forest, x, rep(1, 150)), "saved forest")
})
