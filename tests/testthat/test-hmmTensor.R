library(hmmTensor)
library(rTensor)

context("hmmTensor")

test_that("toyModel generates valid data", {
  toy <- toyModel("simple", T_len = 200, seed = 42)
  expect_identical(toy$K, 2L)
  expect_identical(toy$N, 3L)
  expect_identical(length(toy$Y), 200L)
  expect_identical(length(toy$X), 200L)
  expect_identical(all(toy$Y %in% 1:3), TRUE)
  expect_identical(all(toy$X %in% 1:2), TRUE)
})

test_that("toyModel types work", {
  for (tp in c("simple", "weather", "leftright")) {
    toy <- toyModel(tp, T_len = 100, seed = 1)
    expect_gt(toy$K, 1L)
    expect_gt(toy$N, 1L)
  }
})

test_that("Forward computes valid probabilities", {
  toy <- toyModel("simple", T_len = 100, seed = 42)
  fwd <- Forward(toy$Y, toy$T_mat, toy$O, toy$pi0)
  expect_identical(nrow(fwd$alpha), toy$K)
  expect_identical(ncol(fwd$alpha), length(toy$Y))
  # Each column should sum to ~1 after scaling
  col_sums <- colSums(fwd$alpha)
  expect_lt(max(abs(col_sums - 1)), 1e-10)
  expect_identical(is.finite(fwd$loglik), TRUE)
})

test_that("Backward returns correct dimensions", {
  toy <- toyModel("simple", T_len = 100, seed = 42)
  fwd <- Forward(toy$Y, toy$T_mat, toy$O, toy$pi0)
  bwd <- Backward(toy$Y, toy$T_mat, toy$O, fwd$scale)
  expect_identical(nrow(bwd), toy$K)
  expect_identical(ncol(bwd), length(toy$Y))
})

test_that("Viterbi decoding is accurate", {
  toy <- toyModel("simple", T_len = 500, seed = 42)
  vit <- Viterbi(toy$Y, toy$T_mat, toy$O, toy$pi0)
  acc <- mean(vit$path == toy$X)
  expect_gt(acc, 0.7)
  expect_identical(length(vit$path), length(toy$Y))
})

test_that("BaumWelch recovers parameters", {
  toy <- toyModel("simple", T_len = 1000, seed = 42)
  bw <- BaumWelch(toy$Y, K = 2, N = 3, num.iter = 50)
  # Log-likelihood should increase
  expect_gt(tail(bw$loglik, 1), bw$loglik[1])
  # Rows of T and O should sum to 1
  expect_lt(max(abs(rowSums(bw$T_mat) - 1)), 1e-10)
  expect_lt(max(abs(rowSums(bw$O) - 1)), 1e-10)
})

test_that("Seq2Prob order=2 returns valid matrix", {
  Y <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
  P2 <- Seq2Prob(Y, N = 3, order = 2)
  expect_identical(dim(P2), c(3L, 3L))
  expect_lt(abs(sum(P2) - 1), 1e-10)
})

test_that("Seq2Prob order=3 returns valid tensor", {
  Y <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
  P3 <- Seq2Prob(Y, N = 3, order = 3)
  expect_identical(P3@modes, c(3L, 3L, 3L))
  expect_lt(abs(sum(P3@data) - 1), 1e-10)
})

test_that("HMM symNMF solver works", {
  toy <- toyModel("simple", T_len = 500, seed = 42)
  res <- HMM(toy$Y, K = 2, N = 3, solver = "symNMF", num.iter = 100)
  expect_identical(nrow(res$T_mat), 2L)
  expect_identical(ncol(res$O), 3L)
  expect_lt(max(abs(rowSums(res$T_mat) - 1)), 1e-10)
  expect_lt(max(abs(rowSums(res$O) - 1)), 1e-10)
})

test_that("HMM SVD solver works", {
  toy <- toyModel("simple", T_len = 500, seed = 42)
  res <- HMM(toy$Y, K = 2, N = 3, solver = "SVD")
  expect_identical(nrow(res$T_mat), 2L)
  expect_lt(max(abs(rowSums(res$T_mat) - 1)), 1e-10)
})

test_that("HMM CP solver works", {
  toy <- toyModel("simple", T_len = 500, seed = 42)
  res <- HMM(toy$Y, K = 2, N = 3, solver = "CP", num.iter = 25)
  expect_identical(nrow(res$T_mat), 2L)
  expect_lt(max(abs(rowSums(res$T_mat) - 1)), 1e-10)
})

test_that("HMM TT solver works", {
  toy <- toyModel("simple", T_len = 500, seed = 42)
  res <- HMM(toy$Y, K = 2, N = 3, solver = "TT")
  expect_identical(nrow(res$T_mat), 2L)
  expect_lt(max(abs(rowSums(res$T_mat) - 1)), 1e-10)
})
