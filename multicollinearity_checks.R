#check for multicollinearity
rm(list = ls()) #clear memory
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states); track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

Xc <- d[, exo.cols, drop = FALSE]


#Pairwise correlations
cor_mat <- cor(Xc, use = "pairwise.complete.obs", method = "pearson")
print(round(cor_mat, 3))

#list pairs above threshold
thr <- 0.7
cm <- cor_mat
diag(cm) <- NA
idx <- which(abs(cm) >= thr, arr.ind = TRUE)

if (nrow(idx) == 0) {
  cat("\nNo |cor| >=", thr, "pairs found.\n")
} else {
  pairs <- data.frame(
    var1 = rownames(cm)[idx[,1]],
    var2 = colnames(cm)[idx[,2]],
    cor  = cm[idx]
  )
  pairs <- pairs[pairs$var1 < pairs$var2, , drop = FALSE]
  pairs <- pairs[order(-abs(pairs$cor)), , drop = FALSE]
  cat("\nPairs with |cor| >=", thr, ":\n")
  print(pairs, row.names = FALSE)
}

## ---------- Multicollinearity: VIF ----------
vif_manual <- function(X) {
  X <- as.data.frame(X)
  out <- setNames(numeric(ncol(X)), colnames(X))
  for (j in seq_along(out)) {
    y <- X[[j]]
    Z <- X[-j]
    fit <- lm(y ~ ., data = Z)
    r2 <- summary(fit)$r.squared
    out[j] <- 1 / (1 - r2)
  }
  out
}

vif_vals <- vif_manual(Xc)
cat("\nVIF values:\n")
print(round(vif_vals, 3))
cat("\nRules of thumb: VIF > 5 (moderate), VIF > 10 (high).\n")
