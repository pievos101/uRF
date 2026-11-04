library(uRF)

test_that("response-free forests split and predict with one or multiple threads", {
  x <- iris[, 1:4]
  forests <- lapply(c(1, 2), function(threads) {
    ranger(x = x, clustering = TRUE, num.trees = 10, min.bucket = 5,
           importance = "impurity", seed = 42, num.threads = threads)
  })
  nodes <- lapply(forests, function(forest) {
    expect_identical(forest$treetype, "Clustering")
    expect_true(is.nan(forest$prediction.error))
    expect_true(all(is.finite(forest$variable.importance)))
    prediction <- predict(forest, x, type = "terminalNodes")$predictions
    expect_equal(dim(prediction), c(150, 10))
    expect_true(all(apply(prediction, 2, function(ids) length(unique(ids))) > 1))
    prediction
  })
  expect_identical(nodes[[1]], nodes[[2]])
})

test_that("README clustering, feature importance and federated workflow run", {
  x <- iris[, 1:4]
  result <- uRF(list(x), k = 3)
  expect_length(result$clusters, nrow(x))
  expect_length(unique(result$clusters), 3)
  expect_equal(dim(result$affinity), c(150, 150))
  expect_equal(result$affinity, t(result$affinity))
  expect_true(all(is.finite(result$affinity)))
  expect_gt(length(unique(as.vector(result$affinity))), 1)

  importance <- fimp(x, clusters = result$clusters)
  expect_equal(dim(importance), c(3, 4))
  expect_true(all(is.finite(importance)))

  clients <- split(seq_len(nrow(x)), rep(1:3, length.out = nrow(x)))
  models <- lapply(clients, function(i) uRF(list(x[i, ]), k = 3)$model)
  for (i in clients) {
    federated <- feduRF(list(x[i, ]), models = models, k = 3)
    expect_length(federated$clusters, length(i))
    expect_length(unique(federated$clusters), 3)
    expect_equal(dim(federated$affinity), c(length(i), length(i)))
    expect_true(all(is.finite(federated$affinity)))
  }
})

test_that("missing supervised responses raise an R error", {
  expect_error(ranger(x = iris[, 1:4]), "Either data or x and y")
  expect_error(ranger(clustering = TRUE), "Either data or x and y")
  expect_error(ranger(x = iris[, 1:4], clustering = TRUE,
                      importance = "permutation"), "requires a response")
})

test_that("clustering splits binary features including zero", {
  x <- data.frame(binary = rep(0:1, each = 20))
  forest <- ranger(x = x, clustering = TRUE, num.trees = 5, seed = 42,
                   num.threads = 1, importance = "impurity")
  nodes <- predict(forest, x, type = "terminalNodes")$predictions
  expect_true(all(nodes[1, ] != nodes[40, ]))
  expect_true(all(is.finite(forest$variable.importance)))
})
