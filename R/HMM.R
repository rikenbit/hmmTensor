#' HMM Parameter Estimation via Matrix/Tensor Decomposition
#'
#' Estimates HMM parameters by decomposing the observation co-occurrence
#' matrix/tensor. Supports multiple decomposition solvers.
#'
#' The pipeline:
#' \enumerate{
#'   \item Convert observations to co-occurrence matrix \eqn{\Omega = P_{2,1}} via \code{\link{Seq2Prob}}
#'   \item Decompose \eqn{\Omega} using the chosen solver
#'   \item Recover HMM parameters (T, O, pi) from the decomposition
#' }
#'
#' For NMF-based methods, \eqn{\Omega \approx M \Theta M^\top} where
#' M relates to the emission matrix O and \eqn{\Theta} to \code{diag(pi) \%*\% T}.
#'
#' @param Y Integer vector of observations (values in 1:N), or a list of
#'   integer vectors for multiple sequences.
#' @param K Number of hidden states
#' @param N Number of distinct observation symbols. If NULL, inferred from data.
#' @param solver Decomposition method:
#'   \describe{
#'     \item{\code{"symNMF"}}{Symmetric NMF via \code{symTensor::symNMF} (default)}
#'     \item{\code{"SVD"}}{Truncated SVD of the co-occurrence matrix}
#'     \item{\code{"CP"}}{CP decomposition of the 3rd-order co-occurrence tensor via \code{rTensor::cp}}
#'     \item{\code{"TT"}}{Tensor-Train approximation (via Tucker + reshape)}
#'   }
#' @param Beta Beta-divergence parameter for symNMF (default: 2). Ignored for other solvers.
#' @param order Co-occurrence order for Seq2Prob: 2 (default) or 3.
#'   \code{order = 3} is required for \code{solver = "CP"}.
#' @param lag Lag for pairwise co-occurrence (default: 1)
#' @param smooth Laplace smoothing pseudo-count (default: 1e-10)
#' @param num.iter Maximum iterations for iterative solvers (default: 100)
#' @param thr Convergence threshold (default: 1e-10)
#' @param verbose Logical (default: FALSE)
#' @return A list with components:
#'   \describe{
#'     \item{T_mat}{Estimated transition matrix (K x K)}
#'     \item{O}{Estimated emission matrix (K x N)}
#'     \item{pi0}{Estimated initial distribution (length K)}
#'     \item{Omega}{Co-occurrence matrix/tensor used}
#'     \item{decomp}{Raw decomposition result}
#'     \item{solver}{Solver used}
#'     \item{RecError}{Reconstruction error (if available)}
#'   }
#' @export
#' @examples
#' set.seed(42)
#' toy <- toyModel(type = "simple")
#' result <- HMM(toy$Y, K = toy$K, N = toy$N, solver = "symNMF")
#' result$T_mat
HMM <- function(Y, K, N = NULL, solver = c("symNMF", "SVD", "CP", "TT"),
                Beta = 2, order = 2L, lag = 1L, smooth = 1e-10,
                num.iter = 100L, thr = 1e-10, verbose = FALSE) {
  solver <- match.arg(solver)
  eps <- .Machine$double.eps

  if (is.null(N)) {
    if (is.list(Y)) N <- max(unlist(Y)) else N <- max(Y)
  }

  # Force order=3 for CP
  if (solver == "CP" && order != 3L) {
    order <- 3L
    if (verbose) message("CP solver requires order=3, overriding.")
  }

  # Step 1: Build co-occurrence
  Omega <- Seq2Prob(Y, N = N, order = order, lag = lag, smooth = smooth)

  # Step 2: Decompose
  result <- switch(solver,
    "symNMF" = .solveSymNMF(Omega, K, Beta, num.iter, thr, verbose),
    "SVD"    = .solveSVD(Omega, K),
    "CP"     = .solveCP(Omega, K, num.iter),
    "TT"     = .solveTT(Omega, K, num.iter)
  )

  # Attach metadata
  result$Omega <- Omega
  result$solver <- solver
  result
}

#' @keywords internal
.solveSymNMF <- function(Omega, K, Beta, num.iter, thr, verbose) {
  # Omega (N x N) ≈ M S M^T
  # M -> emission-related, S -> diag(pi) * T related
  res <- symTensor::symNMF(Omega, K = K, model = "MSM", Beta = Beta,
                           num.iter = num.iter, thr = thr, verbose = verbose)
  M <- res$M
  S <- res$S

  # Recover O: normalize M row-wise (M is N x K)
  # O[k, n] = P(Y=n | X=k) -> transpose and normalize
  O_raw <- t(M)  # K x N
  O <- O_raw / rowSums(O_raw)
  O[O < .Machine$double.eps] <- .Machine$double.eps
  O <- O / rowSums(O)

  # Recover pi and T from S
  # S ≈ diag(pi) * T in the ideal case
  # pi_k proportional to sum of S[k, ]
  pi0 <- colSums(S)
  pi0 <- pi0 / sum(pi0)

  # T[i,j] = S[i,j] / pi_i
  T_mat <- sweep(S, 1, rowSums(S), "/")
  T_mat[T_mat < .Machine$double.eps] <- .Machine$double.eps
  T_mat <- T_mat / rowSums(T_mat)

  list(T_mat = T_mat, O = O, pi0 = pi0,
       decomp = res, RecError = tail(res$RecError, 1))
}

#' @keywords internal
.solveSVD <- function(Omega, K) {
  # SVD of P_{2,1}: Omega = U D V^T
  # Use top-K components
  sv <- svd(Omega, nu = K, nv = K)
  U <- sv$u[, 1:K, drop = FALSE]  # N x K
  D <- diag(sv$d[1:K], nrow = K)
  V <- sv$v[, 1:K, drop = FALSE]  # N x K

  # Observable operator approach:
  # O ≈ |U| normalized
  O_raw <- abs(t(U))  # K x N
  O <- O_raw / rowSums(O_raw)
  O[O < .Machine$double.eps] <- .Machine$double.eps
  O <- O / rowSums(O)

  # Recover T from the singular value structure
  # P_{2,1} = O^T diag(pi) T O -> U D V^T
  # Approximate: T ≈ D (as transition strength)
  T_mat <- abs(D)
  T_mat <- T_mat / rowSums(T_mat)

  pi0 <- sv$d[1:K]
  pi0 <- pi0 / sum(pi0)

  recon <- U %*% D %*% t(V)
  rec_err <- sqrt(sum((Omega - recon)^2))

  list(T_mat = T_mat, O = O, pi0 = pi0,
       decomp = sv, RecError = rec_err)
}

#' @keywords internal
.solveCP <- function(Omega, K, num.iter) {
  # Omega is 3rd-order tensor (N x N x N)
  # CP decomposition: T ≈ sum_r lambda_r * u_r o v_r o w_r
  cp_res <- rTensor::cp(Omega, num_components = K, max_iter = num.iter)

  # Factor matrices: U[[1]], U[[2]], U[[3]] each N x K
  U1 <- abs(cp_res$U[[1]])
  U2 <- abs(cp_res$U[[2]])
  U3 <- abs(cp_res$U[[3]])

  # For HMM 3-view: M1 (future), M2 (present), M3 (past)
  # O ≈ average of the three views, normalized
  O_raw <- t((U1 + U2 + U3) / 3)  # K x N
  O <- O_raw / rowSums(O_raw)
  O[O < .Machine$double.eps] <- .Machine$double.eps
  O <- O / rowSums(O)

  # Recover T and pi from lambdas and structure
  lambdas <- abs(cp_res$lambdas)
  pi0 <- lambdas / sum(lambdas)

  # Cross-mode structure gives transition info
  # T[i,j] from how mode-1 and mode-3 components relate
  T_mat <- t(U1) %*% U3  # K x K
  T_mat <- abs(T_mat)
  T_mat <- T_mat / rowSums(T_mat)

  rec_err <- cp_res$fnorm_resid

  list(T_mat = T_mat, O = O, pi0 = pi0,
       decomp = cp_res, RecError = rec_err)
}

#' @keywords internal
.solveTT <- function(Omega, K, num.iter) {
  # For order-2: use Tucker as approximation to TT
  # Tucker: Omega ≈ G x_1 U1 x_2 U2
  if (is.matrix(Omega)) {
    # 2D case: Tucker is equivalent to truncated SVD
    sv <- svd(Omega, nu = K, nv = K)
    U <- sv$u[, 1:K, drop = FALSE]
    D <- diag(sv$d[1:K], nrow = K)
    V <- sv$v[, 1:K, drop = FALSE]

    O_raw <- abs(t(U))
    O <- O_raw / rowSums(O_raw)
    O[O < .Machine$double.eps] <- .Machine$double.eps
    O <- O / rowSums(O)

    # Core as transition
    core <- t(U) %*% Omega %*% V  # K x K
    T_mat <- abs(core)
    T_mat <- T_mat / rowSums(T_mat)

    pi0 <- sv$d[1:K]
    pi0 <- pi0 / sum(pi0)

    recon <- U %*% D %*% t(V)
    rec_err <- sqrt(sum((Omega - recon)^2))

    list(T_mat = T_mat, O = O, pi0 = pi0,
         decomp = sv, RecError = rec_err)
  } else {
    # 3D tensor: Tucker decomposition
    tuck <- rTensor::tucker(Omega, ranks = c(K, K, K), max_iter = num.iter)

    U1 <- abs(tuck$U[[1]])
    O_raw <- t(U1)
    O <- O_raw / rowSums(O_raw)
    O[O < .Machine$double.eps] <- .Machine$double.eps
    O <- O / rowSums(O)

    # Core tensor -> transition
    core_mat <- matrix(tuck$Z@data, K, K * K)
    T_mat <- abs(core_mat %*% t(core_mat))
    T_mat <- T_mat / rowSums(T_mat)

    pi0 <- rep(1 / K, K)

    list(T_mat = T_mat, O = O, pi0 = pi0,
         decomp = tuck, RecError = tuck$fnorm_resid)
  }
}
