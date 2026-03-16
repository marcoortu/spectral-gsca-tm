#' @export
print.sgscatm <- function(x, ...) {
  cat("ILR-Spectral-GSCA Structural Topic Model\n")
  cat("  Topics (K)    :", x$K, "\n")
  cat("  Lambda        :", x$lambda, "\n")
  cat("  Documents (M) :", nrow(x$Pi), "\n")
  cat("  Terms (N)     :", ncol(x$Phi), "\n")
  cat("  Covariates (P):", nrow(x$Bz), "\n")
  cat("  Varimax rot.  :", if (isTRUE(x$rotate)) "yes" else "no", "\n")
  ref_str <- if (isTRUE(x$phi_refined))
    sprintf("yes (smooth=%.0e, temp=%.2f)", x$phi_smooth, x$phi_temp)
  else "no"
  cat("  Phi refined   :", ref_str, "\n")
  invisible(x)
}

#' @export
summary.sgscatm <- function(object, n_terms = 10L, ...) {
  cat("ILR-Spectral-GSCA Structural Topic Model\n")
  cat("================================\n")
  rot_str <- if (isTRUE(object$rotate)) "yes" else "no"
  ref_str <- if (isTRUE(object$phi_refined)) "yes" else "no"
  cat("K =", object$K, "  lambda =", object$lambda,
      "  varimax =", rot_str, "  refined =", ref_str, "\n\n")

  cat("Top eigenvalues of S_z:\n")
  print(round(object$eigenvalues, 4))

  cat("\nPath coefficients (Bz):\n")
  print(round(object$Bz, 4))

  invisible(object)
}

#' Extract path coefficients
#' @export
coef.sgscatm <- function(object, ...) object$Bz

#' Extract fitted topic proportions
#' @export
fitted.sgscatm <- function(object, ...) object$Pi

#' Top terms per topic
#'
#' Returns the top `n` terms for each topic, ranked by the corresponding row
#' of the topic-term loading matrix `Phi`.
#'
#' @param x An object of class `"sgscatm"`.
#' @param n Integer. Number of top terms per topic. Default 10.
#' @param vocab Character vector of length N. Term labels. If NULL, column
#'   indices are returned.
#' @return A K x n character (or integer) matrix.
#' @export
top_terms <- function(x, n = 10L, vocab = NULL) {
  UseMethod("top_terms")
}

#' @export
top_terms.sgscatm <- function(x, n = 10L, vocab = NULL) {
  Phi <- x$Phi
  K   <- x$K
  n   <- min(n, ncol(Phi))
  out <- matrix(NA_character_, nrow = K, ncol = n)
  for (k in seq_len(K)) {
    idx <- order(Phi[k, ], decreasing = TRUE)[seq_len(n)]
    out[k, ] <- if (is.null(vocab)) as.character(idx) else vocab[idx]
  }
  rownames(out) <- paste0("Topic", seq_len(K))
  out
}

#' Predict topic proportions for new documents
#'
#' Projects new documents into the fitted ILR topic space using the
#' topic-term loadings.
#'
#' @param object An `"sgscatm"` object.
#' @param newW Numeric matrix M_new x N. New document-term matrix.
#' @param scale_W Logical. Row-normalise newW. Default TRUE.
#' @param ... Ignored.
#' @return M_new x K matrix of predicted topic proportions.
#' @export
predict.sgscatm <- function(object, newW, scale_W = TRUE, ...) {
  newW <- as.matrix(newW)
  if (scale_W) {
    rs <- rowSums(newW); rs[rs == 0] <- 1
    newW <- newW / rs
  }
  # project onto ILR topic axes via Psi
  # W_tilde_new * Psi' * (Psi * Psi')^{-1} gives Z_new
  Psi <- object$Psi                            # (K-1) x N
  w_bar_approx <- colMeans(newW)
  W_tilde_new  <- sweep(newW, 2L, w_bar_approx, "-")
  PsiPsit <- tcrossprod(Psi)                   # (K-1) x (K-1)
  Z_new <- W_tilde_new %*% t(Psi) %*% solve(PsiPsit)  # M_new x (K-1)
  scores <- Z_new %*% t(object$V)
  .softmax_rows(scores)
}
