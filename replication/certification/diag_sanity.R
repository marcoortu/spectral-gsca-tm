suppressPackageStartupMessages(devtools::load_all(".", quiet=TRUE))
source("replication/simulation/sim_dgp.R"); sg <- getNamespace("sgscatm")
K<-5L;P<-4L
Bz0 <- matrix(rnorm(P*(K-1),0,0.3),P,K-1)
V <- ilr_contrast(K)
# 1) perm+sign trivial: aligning Bz0 to itself and a relabeled copy -> ~0
cat("self-align mse:", perm_sign_align(Bz0,Bz0,V)$mse, "\n")
Pm <- diag(K)[c(3,1,5,2,4),]; Qp <- crossprod(V, Pm%*%V)
cat("relabeled-align mse:", perm_sign_align(Bz0%*%t(Qp),Bz0,V)$mse, "(should be ~0)\n")
# 2) oracle-frozen@1 at L=1e4, M=2000: bias should be tiny
libfun <- function(m) pmax(rnbinom(m,size=3,mu=1e4),500L)
set.seed(1)
Bs <- vector("list",25); SEs <- vector("list",25)
for(r in 1:25){
  dat <- sim_dgp(M=2000L,N=500L,K=K,P=P,Bz0=Bz0,sigma_eps=0.3,alpha_beta=0.05,
                 doc_length=libfun,seed=99000L+r)
  Wf<-dat$W/rowSums(dat$W); C<-scale(dat$C,TRUE,FALSE)
  fit<-sgscatm(dat$W,dat$C,K=K,lambda=1,rotate=FALSE)
  gl<-sg$.sg_gl_align(fit$Z,scale(dat$Z_true,TRUE,FALSE))
  # frozen@1 with TRUE Phi
  zs<-sg$.sg_z_step(gl$Z,dat$Beta,Wf,V,0,NULL,rep(1e-6,2000L),n_gn=2L,dz_cap=1)
  B<-sg$.sg_b_step(zs$Z,C)
  cv<-perm_sign_coverage(B,Bz0,sg$.sg_sandwich(zs$Z,C,B),V)
  Bs[[r]]<-perm_sign_align(B,Bz0,V)$B; SEs[[r]]<-cv$se
}
arr<-array(unlist(Bs),c(P,K-1,25)); Bbar<-apply(arr,c(1,2),mean)
cat(sprintf("oracle-frozen@1 L=1e4 M=2000: bias_norm=%.4f (||Bz0||=%.3f) SE/SD=%.2f cov=%.3f\n",
  sqrt(sum((Bbar-Bz0)^2)), sqrt(sum(Bz0^2)),
  median(apply((Reduce(`+`,SEs)/25),c(1,2),identity)/apply(arr,c(1,2),sd)),
  mean(mapply(function(b,se) mean(abs(b-Bz0)<=1.96*se), Bs, SEs))))
