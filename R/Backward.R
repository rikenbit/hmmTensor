#' Backward Algorithm for HMM
#'
#' Computes backward probabilities \eqn{\beta_t(i) = P(Y_{t+1}, \ldots, Y_T | X_t = i)}
#' using the same scaling factors from the Forward algorithm.
#'
#' @param Y Integer vector of observations (values in 1:N)
#' @param T_mat Transition matrix (K x K), \code{T_mat[i,j] = P(X_t=j | X_{t-1}=i)}
#' @param O Emission matrix (K x N), \code{O[i,j] = P(Y_t=j | X_t=i)}
#' @param scale Scaling factors from \code{\link{Forward}} (length T_len)
#' @return Matrix (K x T_len) of scaled backward probabilities
#' @export
Backward <- function(Y, T_mat, O, scale) {
  K <- nrow(T_mat)
  T_len <- length(Y)

  beta <- matrix(0, K, T_len)

  # Initialization: beta_T(i) = 1 (scaled)
  beta[, T_len] <- 1 / scale[T_len]

  # Recursion
  for (t in (T_len - 1):1) {
    beta[, t] <- T_mat %*% (O[, Y[t + 1]] * beta[, t + 1])
    beta[, t] <- beta[, t] / scale[t]
  }

  beta
}
