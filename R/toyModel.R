#' Generate Toy HMM Data
#'
#' Creates synthetic HMM data with visually clear structure for
#' demonstrations and testing.
#'
#' @param type Type of toy model:
#'   \describe{
#'     \item{\code{"simple"}}{2 states, 3 observations. States alternate with
#'       high probability. Clear emission separation.}
#'     \item{\code{"weather"}}{3 states (Sunny/Cloudy/Rainy), 4 observations
#'       (Walk/Shop/Clean/Stay). Classic weather HMM.}
#'     \item{\code{"leftright"}}{3 states with left-to-right transitions only.
#'       Models sequential processes.}
#'   }
#' @param T_len Length of observation sequence (default: 500)
#' @param seed Random seed (default: NULL)
#' @return A list with components:
#'   \describe{
#'     \item{Y}{Integer vector of observations}
#'     \item{X}{Integer vector of true hidden states}
#'     \item{T_mat}{True transition matrix}
#'     \item{O}{True emission matrix}
#'     \item{pi0}{True initial distribution}
#'     \item{K}{Number of hidden states}
#'     \item{N}{Number of observation symbols}
#'     \item{type}{Model type}
#'   }
#' @export
#' @examples
#' toy <- toyModel("simple", T_len = 200, seed = 42)
#' table(toy$X)
#' table(toy$Y)
toyModel <- function(type = c("simple", "weather", "leftright"),
                     T_len = 500L, seed = NULL) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  params <- switch(type,
    "simple" = list(
      K = 2L, N = 3L,
      T_mat = matrix(c(
        0.9, 0.1,
        0.2, 0.8
      ), 2, 2, byrow = TRUE),
      O = matrix(c(
        0.7, 0.2, 0.1,
        0.1, 0.2, 0.7
      ), 2, 3, byrow = TRUE),
      pi0 = c(0.6, 0.4)
    ),
    "weather" = list(
      K = 3L, N = 4L,
      T_mat = matrix(c(
        0.6, 0.3, 0.1,
        0.2, 0.5, 0.3,
        0.1, 0.3, 0.6
      ), 3, 3, byrow = TRUE),
      O = matrix(c(
        0.6, 0.2, 0.1, 0.1,
        0.1, 0.4, 0.3, 0.2,
        0.05, 0.1, 0.4, 0.45
      ), 3, 4, byrow = TRUE),
      pi0 = c(0.4, 0.3, 0.3)
    ),
    "leftright" = list(
      K = 3L, N = 4L,
      T_mat = matrix(c(
        0.7, 0.3, 0.0,
        0.0, 0.6, 0.4,
        0.0, 0.0, 1.0
      ), 3, 3, byrow = TRUE),
      O = matrix(c(
        0.7, 0.2, 0.1, 0.0,
        0.1, 0.5, 0.3, 0.1,
        0.0, 0.1, 0.3, 0.6
      ), 3, 4, byrow = TRUE),
      pi0 = c(1.0, 0.0, 0.0)
    )
  )

  # Generate sequences
  X <- integer(T_len)
  Y <- integer(T_len)

  X[1] <- sample(params$K, 1, prob = params$pi0)
  Y[1] <- sample(params$N, 1, prob = params$O[X[1], ])

  for (t in 2:T_len) {
    X[t] <- sample(params$K, 1, prob = params$T_mat[X[t - 1], ])
    Y[t] <- sample(params$N, 1, prob = params$O[X[t], ])
  }

  list(
    Y = Y,
    X = X,
    T_mat = params$T_mat,
    O = params$O,
    pi0 = params$pi0,
    K = params$K,
    N = params$N,
    type = type
  )
}
