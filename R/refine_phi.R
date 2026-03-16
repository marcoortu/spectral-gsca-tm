#' Re-estimate Topic-Word Distributions via One Post-Spectral Refinement Step
#'
#' After fitting sgscatm, the stored topic-term matrix
#' \eqn{\hat{\boldsymbol{\Phi}} = \mathbf{V}\hat{\boldsymbol{\Psi}}
#' + \mathbf{1}_K\bar{\mathbf{w}}^\top} is constrained to a
#' \eqn{(K{-}1)}-dimensional affine ILR subspace. This limits topic
#' exclusivity and diversity because all K distributions lie on the same
#' low-dimensional hyperplane anchored at the corpus mean.
#'
#' This function breaks the subspace constraint by re-estimating
#' \eqn{\hat{\boldsymbol{\Phi}}} as the posterior-weighted average of
#' document term-frequency vectors:
#' \deqn{\hat{\phi}^{(1)}_k \;\propto\;
#'   \sum_i q_{ik}\,\tilde{w}_i + \varepsilon,}
#' where \eqn{q_{ik}} are document-topic assignments and
#' \eqn{\varepsilon} is Laplace smoothing.
#'
#' Two methods are available for computing the assignments \eqn{q_{ik}}:
#' \describe{
#'   \item{`"kmeans"` (default)}{K-means on the ILR scores
#'     \eqn{\sqrt{M}\,\hat{\mathbf{Z}}} (scaled to natural population
#'     variance). Documents are partitioned into K clusters; each cluster
#'     becomes one topic. This directly exploits the spectral structure of
#'     sgscatm and produces **diverse, exclusive** topics.}
#'   \item{`"em"`}{One LDA-style E-step: soft assignments
#'     \eqn{q_{ik} \propto \exp(\sum_n \tilde{w}_{in}\,
#'     \log\hat\phi^{(0)}_{kn})} from the current ILR-based
#'     topic-word distribution \eqn{\hat{\boldsymbol{\phi}}^{(0)} =
#'     \text{softmax}(\mathbf{V}\hat{\boldsymbol{\Psi}})}. Improves
#'     semantic coherence but has limited effect on diversity when
#'     topics are highly similar.}
#' }
#'
#' The document-topic proportions \eqn{\hat{\boldsymbol{\Pi}}}, path
#' coefficients \eqn{\hat{\mathbf{B}}_z}, ILR scores \eqn{\hat{\mathbf{Z}}},
#' and ILR deviations \eqn{\hat{\boldsymbol{\Psi}}} are **unchanged**.
#'
#' An optional temperature \eqn{\tau \in (0,1]} further sharpens the
#' M-step distribution by raising every entry to \eqn{1/\tau}.
#'
#' @param fit An `"sgscatm"` object returned by [sgscatm()].
#' @param W Numeric M x N document-term matrix used to fit `fit`.
#' @param smooth Non-negative numeric. Laplace smoothing added to every entry
#'   of the M-step estimate before row-normalising. Default `1e-4`.
#' @param temp Positive numeric. Temperature for post-M-step sharpening.
#'   `1` (default) leaves the M-step distribution unchanged; values less
#'   than `1` concentrate mass on the top terms.
#' @param method Character. Assignment method: `"kmeans"` (default, better
#'   diversity/exclusivity) or `"em"` (LDA E-step, better coherence).
#' @param nstart Integer. Number of random starts for k-means
#'   (`method = "kmeans"` only). Default `10L`.
#' @param seed Integer or NULL. Random seed for k-means reproducibility.
#'
#' @return The `fit` object with `Phi` replaced by the refined probability
#'   matrix (K x N, rows sum to 1). Fields added: `phi_refined`, `phi_smooth`,
#'   `phi_temp`, `phi_method`.
#'
#' @seealso [topic_word_dist()], [sgscatm()]
#' @export
refine_phi <- function(fit, W, smooth = 1e-4, temp = 1.0,
                       method = c("kmeans", "em"),
                       nstart = 10L, seed = NULL) {
  stopifnot(inherits(fit, "sgscatm"))
  method <- match.arg(method)
  W <- as.matrix(W)
  stopifnot(
    nrow(W) == nrow(fit$Pi),
    is.numeric(smooth), length(smooth) == 1L, smooth >= 0,
    is.numeric(temp),   length(temp)   == 1L, temp   >  0
  )

  # Row-normalise consistently with the fitted model
  if (isTRUE(fit$scale_W)) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W  <- W / rs
  }

  K      <- fit$K
  m_docs <- nrow(fit$Z)

  # ---- compute document-topic assignment weights q (M x K) ------------------
  if (method == "kmeans") {
    # Scale ILR scores to natural population variance: Z_scaled = Z * sqrt(M)
    # (Z_star has unit-norm columns, so Z_scaled has column variance ~1)
    z_scaled <- fit$Z * sqrt(m_docs)            # M x (K-1)

    if (!is.null(seed)) set.seed(seed)
    km       <- kmeans(z_scaled, centers = K,
                       nstart = nstart, iter.max = 50L)
    clusters <- km$cluster                       # M integer vector in 1:K

    # Hard-assignment indicator matrix (M x K)
    q_mat <- matrix(0, m_docs, K)
    for (k in seq_len(K)) {
      members <- which(clusters == k)
      if (length(members) > 0L) q_mat[members, k] <- 1.0
    }

  } else {
    # LDA-style E-step: soft assignments from current ILR topic-word probs
    phi_base <- topic_word_dist(fit, temp = 1.0) # K x N
    log_phi  <- log(pmax(phi_base, 1e-300))       # guard log(0)

    # log_lik[i,k] = sum_n W[i,n]*log(phi_base[k,n])  (cancel normalisation
    # constants via softmax, so only the linear term matters)
    log_lik <- W %*% t(log_phi)                  # M x K
    q_mat   <- .softmax_rows(log_lik)            # M x K
  }

  # ---- M-step: posterior-weighted topic-word distributions -------------------
  phi_new <- crossprod(q_mat, W)                 # K x N

  # Laplace smoothing then row-normalise to proper probabilities
  phi_new <- phi_new + smooth
  phi_new <- phi_new / rowSums(phi_new)

  # Optional temperature sharpening
  if (temp != 1.0) {
    phi_new <- phi_new ^ (1 / temp)
    phi_new <- phi_new / rowSums(phi_new)
  }

  fit$Phi          <- phi_new
  fit$phi_refined  <- TRUE
  fit$phi_smooth   <- smooth
  fit$phi_temp     <- temp
  fit$phi_method   <- method
  fit
}


#' Extract Topic-Word Probability Matrix
#'
#' Returns the K x N topic-word probability matrix for a fitted sgscatm model,
#' with optional temperature sharpening.
#'
#' For models refined via [refine_phi()] the stored `Phi` is returned.
#' For unrefined models the ILR deviation
#' \eqn{\mathbf{V}\hat{\boldsymbol{\Psi}}} is converted to row-wise softmax
#' probabilities. In both cases `temp < 1` can sharpen the distribution.
#'
#' @param fit An `"sgscatm"` object.
#' @param temp Positive numeric. Temperature. `1` (default) returns
#'   the distribution unchanged; values less than `1` sharpen.
#' @return K x N numeric matrix with rows summing to 1.
#' @seealso [refine_phi()]
#' @export
topic_word_dist <- function(fit, temp = 1.0) {
  UseMethod("topic_word_dist")
}

#' @export
topic_word_dist.sgscatm <- function(fit, temp = 1.0) {
  stopifnot(is.numeric(temp), length(temp) == 1L, temp > 0)

  if (isTRUE(fit$phi_refined)) {
    phi_out <- fit$Phi
    if (temp != 1.0) {
      phi_out <- phi_out ^ (1 / temp)
      phi_out <- phi_out / rowSums(phi_out)
    }
  } else {
    phi_dev <- fit$V %*% fit$Psi               # K x N
    phi_out <- .softmax_rows(
      if (temp != 1.0) phi_dev / temp else phi_dev
    )
  }

  phi_out
}
