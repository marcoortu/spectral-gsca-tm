#' Convert ILR scores to topic proportions
#'
#' Applies the ILR inverse (softmax of V z) to a matrix of ILR scores.
#'
#' @param Z Numeric M x (K-1) matrix of ILR scores.
#' @param V Numeric K x (K-1) ILR contrast matrix.
#' @return M x K matrix of topic proportions.
#' @export
ilr_to_proportions <- function(Z, V) {
  scores <- Z %*% t(V)   # M x K
  E <- exp(scores - apply(scores, 1L, max))
  E / rowSums(E)
}

#' Convert topic proportions to ILR scores
#'
#' @param Pi Numeric M x K matrix of topic proportions (rows sum to 1,
#'   all entries > 0).
#' @param V Numeric K x (K-1) ILR contrast matrix.
#' @return M x (K-1) matrix of ILR scores.
#' @export
proportions_to_ilr <- function(Pi, V) {
  log_pi <- log(Pi)
  log_pi %*% V   # M x (K-1)
}

#' Compute Aitchison distance between two sets of compositions
#'
#' @param Pi1,Pi2 Numeric M x K matrices of compositions.
#' @param V Numeric K x (K-1) ILR contrast matrix.
#' @return Numeric vector of length M.
#' @export
aitchison_dist <- function(Pi1, Pi2, V) {
  Z1 <- proportions_to_ilr(Pi1, V)
  Z2 <- proportions_to_ilr(Pi2, V)
  sqrt(rowSums((Z1 - Z2)^2))
}
