#' Build ILR Contrast Matrix
#'
#' Constructs the K x (K-1) isometric log-ratio (ILR) contrast matrix **V**
#' satisfying \eqn{V^\top V = I_{K-1}} and \eqn{V^\top \mathbf{1}_K = \mathbf{0}}.
#'
#' The columns of **V** form an orthonormal basis for the subspace of
#' \eqn{\mathbb{R}^K} orthogonal to \eqn{\mathbf{1}_K}, which is the natural
#' tangent space of the probability simplex at its centroid.
#'
#' @param K Integer >= 2. Number of topics (simplex dimension).
#' @return A K x (K-1) numeric matrix.
#' @export
ilr_contrast <- function(K) {
  stopifnot(is.numeric(K), length(K) == 1L, K >= 2L)
  K <- as.integer(K)
  # Helmert-style orthonormal basis for 1_K^\perp
  V <- matrix(0, nrow = K, ncol = K - 1L)
  for (j in seq_len(K - 1L)) {
    # first j entries: 1/sqrt(j*(j+1))
    # (j+1)-th entry:  -j/sqrt(j*(j+1)) = -sqrt(j/(j+1))
    V[seq_len(j), j]  <-  1 / sqrt(j * (j + 1L))
    V[j + 1L,     j]  <- -sqrt(j / (j + 1L))
  }
  V
}
