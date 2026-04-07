#' Convert Observation Sequences to Co-occurrence Matrix/Tensor
#'
#' Constructs empirical co-occurrence statistics from observation sequences.
#' These form the basis for matrix/tensor decomposition approaches to HMM.
#'
#' @param Y Integer vector of observations (values in 1:N), or a list of
#'   integer vectors for multiple sequences.
#' @param N Number of distinct observation symbols. If NULL, inferred from data.
#' @param order Co-occurrence order: 2 (pairwise, default) or 3 (triple).
#'   Order 2 gives an N x N matrix \eqn{P_{2,1}(i,j) = P(Y_2=i, Y_1=j)}.
#'   Order 3 gives an N x N x N tensor \eqn{P_{3}(i,j,k) = P(Y_3=i, Y_2=j, Y_1=k)}.
#' @param lag Lag for pairwise co-occurrence (default: 1).
#'   \eqn{\Omega^{(\tau)}(i,j) = P(Y_t=i, Y_{t+\tau}=j)}.
#'   Only used when \code{order = 2}.
#' @param smooth Laplace smoothing pseudo-count (default: 0).
#'   Adds \code{smooth} to all counts before normalization.
#' @return For \code{order = 2}: an N x N matrix (normalized).
#'   For \code{order = 3}: an rTensor Tensor object of dimension N x N x N.
#' @export
#' @examples
#' Y <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
#' P2 <- Seq2Prob(Y, N = 3, order = 2)
#' P3 <- Seq2Prob(Y, N = 3, order = 3)
Seq2Prob <- function(Y, N = NULL, order = 2L, lag = 1L, smooth = 0) {
  if (!is.list(Y)) Y <- list(Y)

  if (is.null(N)) {
    N <- max(unlist(Y))
  }

  if (order == 2L) {
    # Pairwise co-occurrence: P(Y_{t+lag} = i, Y_t = j)
    counts <- matrix(smooth, N, N)
    for (y in Y) {
      T_len <- length(y)
      if (T_len <= lag) next
      for (t in 1:(T_len - lag)) {
        counts[y[t + lag], y[t]] <- counts[y[t + lag], y[t]] + 1
      }
    }
    # Normalize to joint probability
    total <- sum(counts)
    if (total > 0) counts <- counts / total
    return(counts)

  } else if (order == 3L) {
    # Triple co-occurrence: P(Y_{t+2} = i, Y_{t+1} = j, Y_t = k)
    counts <- array(smooth, dim = c(N, N, N))
    for (y in Y) {
      T_len <- length(y)
      if (T_len < 3) next
      for (t in 1:(T_len - 2)) {
        counts[y[t + 2], y[t + 1], y[t]] <- counts[y[t + 2], y[t + 1], y[t]] + 1
      }
    }
    total <- sum(counts)
    if (total > 0) counts <- counts / total
    return(rTensor::as.tensor(counts))

  } else {
    stop("order must be 2 or 3")
  }
}
