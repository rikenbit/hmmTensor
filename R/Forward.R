#' Forward Algorithm for HMM
#'
#' Computes forward probabilities \eqn{\alpha_t(i) = P(Y_1, \ldots, Y_t, X_t = i)}
#' with log-scaling to avoid underflow.
#'
#' @param Y Integer vector of observations (values in 1:N)
#' @param T_mat Transition matrix (K x K), \code{T_mat[i,j] = P(X_t=j | X_{t-1}=i)}
#' @param O Emission matrix (K x N), \code{O[i,j] = P(Y_t=j | X_t=i)}
#' @param pi0 Initial state distribution (length K)
#' @return A list with components:
#'   \describe{
#'     \item{alpha}{Matrix (K x T_len) of scaled forward probabilities}
#'     \item{loglik}{Log-likelihood of the observation sequence}
#'     \item{scale}{Vector of scaling factors (length T_len)}
#'   }
#' @export
Forward <- function(Y, T_mat, O, pi0) {
  K <- nrow(T_mat)
  T_len <- length(Y)

  alpha <- matrix(0, K, T_len)
  scale <- numeric(T_len)

  # Initialization: alpha_1(i) = pi0(i) * O(i, Y_1)
  alpha[, 1] <- pi0 * O[, Y[1]]
  scale[1] <- sum(alpha[, 1])
  alpha[, 1] <- alpha[, 1] / scale[1]

  # Recursion
  for (t in 2:T_len) {
    alpha[, t] <- (t(T_mat) %*% alpha[, t - 1]) * O[, Y[t]]
    scale[t] <- sum(alpha[, t])
    alpha[, t] <- alpha[, t] / scale[t]
  }

  loglik <- sum(log(scale))

  list(alpha = alpha, loglik = loglik, scale = scale)
}
