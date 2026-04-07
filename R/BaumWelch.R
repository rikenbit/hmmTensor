#' Baum-Welch Algorithm (EM) for HMM
#'
#' Estimates HMM parameters (T, O, pi) from observation sequences
#' using the Expectation-Maximization algorithm.
#'
#' @param Y Integer vector of observations (values in 1:N), or a list of
#'   integer vectors for multiple sequences.
#' @param K Number of hidden states
#' @param N Number of distinct observation symbols
#' @param initT Initial transition matrix (K x K). If NULL, random initialization.
#' @param initO Initial emission matrix (K x N). If NULL, random initialization.
#' @param initPi Initial state distribution (length K). If NULL, uniform.
#' @param num.iter Maximum EM iterations (default: 100)
#' @param thr Convergence threshold on log-likelihood relative change (default: 1e-6)
#' @param verbose Logical (default: FALSE)
#' @return A list with components:
#'   \describe{
#'     \item{T_mat}{Estimated transition matrix (K x K)}
#'     \item{O}{Estimated emission matrix (K x N)}
#'     \item{pi0}{Estimated initial distribution (length K)}
#'     \item{loglik}{Vector of log-likelihoods per iteration}
#'     \item{converged}{Logical}
#'     \item{iter}{Number of iterations}
#'   }
#' @export
BaumWelch <- function(Y, K, N,
                      initT = NULL, initO = NULL, initPi = NULL,
                      num.iter = 100L, thr = 1e-6,
                      verbose = FALSE) {
  eps <- .Machine$double.eps

  # Handle single sequence as list
  if (!is.list(Y)) Y <- list(Y)

  # Initialize parameters
  if (is.null(initT)) {
    T_mat <- matrix(runif(K * K), K, K)
    T_mat <- T_mat / rowSums(T_mat)
  } else {
    T_mat <- initT
  }

  if (is.null(initO)) {
    O <- matrix(runif(K * N), K, N)
    O <- O / rowSums(O)
  } else {
    O <- initO
  }

  if (is.null(initPi)) {
    pi0 <- rep(1 / K, K)
  } else {
    pi0 <- initPi
  }

  loglik_hist <- numeric(num.iter)
  converged <- FALSE

  for (iter in seq_len(num.iter)) {
    # Accumulators
    T_numer <- matrix(eps, K, K)
    T_denom <- rep(eps, K)
    O_numer <- matrix(eps, K, N)
    O_denom <- rep(eps, K)
    pi_accum <- rep(0, K)
    total_loglik <- 0

    for (seq_idx in seq_along(Y)) {
      y <- Y[[seq_idx]]
      T_len <- length(y)

      # E-step
      fwd <- Forward(y, T_mat, O, pi0)
      bwd <- Backward(y, T_mat, O, fwd$scale)
      total_loglik <- total_loglik + fwd$loglik

      # Posterior: gamma_t(i) = P(X_t = i | Y)
      gamma <- fwd$alpha * bwd
      # Correct scaling: alpha is scaled by 1/c_t, beta is scaled by 1/c_t
      # gamma needs normalization
      gamma_sum <- colSums(gamma)
      gamma_sum[gamma_sum == 0] <- eps
      gamma <- sweep(gamma, 2, gamma_sum, "/")

      # xi_t(i,j) = P(X_t=i, X_{t+1}=j | Y)
      for (t in 1:(T_len - 1)) {
        xi <- (fwd$alpha[, t] %o% (O[, y[t + 1]] * bwd[, t + 1])) * T_mat
        xi_sum <- sum(xi)
        if (xi_sum > 0) xi <- xi / xi_sum
        T_numer <- T_numer + xi
        T_denom <- T_denom + gamma[, t]
      }

      # Emission accumulation
      for (t in 1:T_len) {
        O_numer[, y[t]] <- O_numer[, y[t]] + gamma[, t]
        O_denom <- O_denom + gamma[, t]
      }

      pi_accum <- pi_accum + gamma[, 1]
    }

    loglik_hist[iter] <- total_loglik

    # Check convergence
    if (iter > 1) {
      rel_change <- abs(loglik_hist[iter] - loglik_hist[iter - 1]) /
        (abs(loglik_hist[iter]) + eps)
      if (rel_change < thr) {
        converged <- TRUE
        if (verbose) message("Converged at iteration ", iter)
        break
      }
    }

    if (verbose && (iter %% 10 == 0 || iter == 1)) {
      message("Iter ", iter, ": loglik = ", round(total_loglik, 4))
    }

    # M-step
    T_mat <- T_numer / T_denom
    T_mat <- T_mat / rowSums(T_mat)

    O <- O_numer / O_denom
    O <- O / rowSums(O)

    pi0 <- pi_accum / sum(pi_accum)
  }

  list(
    T_mat = T_mat,
    O = O,
    pi0 = pi0,
    loglik = loglik_hist[seq_len(iter)],
    converged = converged,
    iter = iter
  )
}
