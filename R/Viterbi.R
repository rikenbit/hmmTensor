#' Viterbi Algorithm for HMM
#'
#' Finds the most likely state sequence given the observations
#' using the Viterbi algorithm (log-space).
#'
#' @param Y Integer vector of observations (values in 1:N)
#' @param T_mat Transition matrix (K x K), \code{T_mat[i,j] = P(X_t=j | X_{t-1}=i)}
#' @param O Emission matrix (K x N), \code{O[i,j] = P(Y_t=j | X_t=i)}
#' @param pi0 Initial state distribution (length K)
#' @return A list with components:
#'   \describe{
#'     \item{path}{Integer vector of most likely states (length T_len)}
#'     \item{loglik}{Log-likelihood of the best path}
#'   }
#' @export
Viterbi <- function(Y, T_mat, O, pi0) {
  K <- nrow(T_mat)
  T_len <- length(Y)
  eps <- .Machine$double.eps

  log_T <- log(T_mat + eps)
  log_O <- log(O + eps)
  log_pi <- log(pi0 + eps)

  # delta[i, t] = max log-prob of best path ending in state i at time t
  delta <- matrix(-Inf, K, T_len)
  psi <- matrix(0L, K, T_len)  # backpointer

  # Initialization
  delta[, 1] <- log_pi + log_O[, Y[1]]

  # Recursion
  for (t in 2:T_len) {
    for (j in 1:K) {
      candidates <- delta[, t - 1] + log_T[, j]
      psi[j, t] <- which.max(candidates)
      delta[j, t] <- candidates[psi[j, t]] + log_O[j, Y[t]]
    }
  }

  # Backtracking
  path <- integer(T_len)
  path[T_len] <- which.max(delta[, T_len])
  loglik <- delta[path[T_len], T_len]

  for (t in (T_len - 1):1) {
    path[t] <- psi[path[t + 1], t + 1]
  }

  list(path = path, loglik = loglik)
}
