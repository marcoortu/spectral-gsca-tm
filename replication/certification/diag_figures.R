suppressPackageStartupMessages(library(ggplot2))
P1 <- readRDS("output/diag_P1_oracle.rds")
dir.create("output/figures", showWarnings=FALSE, recursive=TRUE)
# k-curve: frozen bias vs k, L=1e3 vs 1e4, facet by M
fr <- P1[grepl("^f", P1$tag), ]
fr$k <- as.integer(sub("f","",fr$tag))
fr$L <- factor(fr$L, labels=c("L=1e3","L=1e4"))
p <- ggplot(fr, aes(k, bias_norm, colour=L)) +
  geom_line(linewidth=0.8) + geom_point(size=2) +
  facet_wrap(~M, labeller=label_both, scales="free_y") +
  scale_colour_manual(values=c("L=1e3"="#B7743A","L=1e4"="#534AB7")) +
  labs(x="frozen-Phi sweeps k", y=expression("oracle bias norm ||"*bar(B)-B[z0]*"||"[F]),
       title="A4: frozen bias-vs-k, L-contrast (flat at L=1e4 = depth-cured 1/L bias)") +
  theme_minimal(base_size=11)
ggsave("output/figures/diag_kcurve.pdf", p, width=8, height=3.2)
# A3 lead table
lead <- P1[P1$L==1e4 & P1$tag=="f10", c("M","bias_norm","se_sd","true_cov")]
cat("A3 lead readout (oracle-frozen@10 + plain sandwich, L=1e4):\n"); print(lead, row.names=FALSE, digits=3)
cat("\njoint@conv (L=1e4) for contrast:\n")
print(P1[P1$L==1e4 & P1$tag=="joint", c("M","bias_norm","se_sd","true_cov")], row.names=FALSE, digits=3)
