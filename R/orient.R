#' General-linear orientation of the spectral pilot (internal)
#'
#' Section 3.6(iii): the unit-norm pilot scores are identified only up to an
#' invertible, \eqn{\sqrt M}-scaled linear map, so orientation and scale are
#' supplied from the anchor read-out by a **general-linear** least-squares
#' alignment (never orthogonal Procrustes). Given pilot scores
#' \eqn{Z^\ast} (unit-norm columns) and the anchored read-out scores
#' \eqn{\hat Z_A} on the generative scale,
#' \deqn{\hat A = \arg\min_A \lVert Z^\ast A - H_M\hat Z_A\rVert_F^2,\qquad
#'       \hat Z_0 = Z^\ast\hat A + \mathbf 1\,\bar z_A^\top,}
#' with \eqn{H_M} the centring projector. The \eqn{\sqrt M} rescaling is
#' absorbed into \eqn{\hat A} because \eqn{Z^\ast} is unit-norm while
#' \eqn{\hat Z_A} is on the generative scale.
#' @keywords internal
.sg_gl_align <- function(Zpilot, Ztarget, ridge = 1e-8) {
  a0 <- colMeans(Ztarget)
  Tc <- sweep(Ztarget, 2L, a0)                       # centred target
  A  <- solve(crossprod(Zpilot) + ridge * diag(ncol(Zpilot)),
              crossprod(Zpilot, Tc))
  Z0 <- sweep(Zpilot %*% A, 2L, a0, "+")
  list(Z = Z0, A = A, a0 = a0)
}
