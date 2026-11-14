## runall.R -- rerun every experiment in the CLeaR 2027 paper (one self-contained file)
##
##   Table 2  (monotonicity holds)        -> results/table2_rows.tex
##   Table 3  (monotonicity violated)     -> results/table3_rows.tex
##   Section 5 (medical-cost data)        -> results/section5_numbers.tex
##
## Needs only base R (>= 3.6); downloads insurance.csv if it is not in the folder.
##
##   Terminal:  Rscript runall.R          # full run, about 8-10 minutes
##              Rscript runall.R quick    # smoke test, about 1 minute
##   RStudio:   setwd("folder/with/runall.R"); source("runall.R")
##
## Seeds are fixed, so a full run reproduces the numbers in revision_patches.md.

QUICK <- "quick" %in% commandArgs(trailingOnly = TRUE)
OUT   <- "results"
dir.create(OUT, showWarnings = FALSE)
t_start <- Sys.time()

###############################################################################
## 1. Estimators (Appendix D)
##    Local linear CDF estimator (Hall, Wolff & Yao, 1999), product Epanechnikov
##    kernel. Every quantity uses ONE matrix of CDF estimates on ONE shared set of
##    Monte Carlo y-draws, so P_LB - N_LB = CE and P_UB - N_UB = CE hold exactly.
###############################################################################

epan <- function(u) 0.75 * pmax(1 - u^2, 0)

kern_w <- function(Z, z0, h) {
  k <- epan((Z[, 1] - z0[1]) / h[1])
  if (ncol(Z) > 1) for (j in 2:ncol(Z)) k <- k * epan((Z[, j] - z0[j]) / h[j])
  k
}

## LL smoother weights: alpha_hat(y; z0) = sum_i l_i 1{Y_i < y}; NA if singular
ll_weights <- function(Z, z0, h) {
  k <- kern_w(Z, z0, h)
  if (sum(k > 0) <= ncol(Z) + 1) return(rep(NA_real_, nrow(Z)))
  D <- cbind(1, sweep(Z, 2, z0))
  A <- crossprod(D, D * k)
  if (rcond(A) < 1e-12) return(rep(NA_real_, nrow(Z)))
  as.vector(solve(A, t(D * k))[1, ])
}

## Clip to [0, 1] and monotone-rearrange in y (Chernozhukov et al., 2010)
rearrange_rows <- function(Fm, ys) {
  o <- order(ys)
  Fm[, o] <- t(apply(matrix(pmin(pmax(Fm[, o], 0), 1), nrow = nrow(Fm)), 1, sort))
  Fm
}

cdf_matrix <- function(Z, Y, Z0, h, ys) {
  L <- t(apply(Z0, 1, function(z0) ll_weights(Z, z0, h)))
  if (anyNA(L)) stop("singular local design at a target point")
  rearrange_rows(L %*% (outer(Y, ys, "<") * 1), ys)
}

## Eqs. (37)-(42) and (81)-(84); pairs with x0 == x1 contribute 0
pn_from_cdfs <- function(F0, F1, a, b, n_total) {
  I <- function(M) (b - a) * sum(M) / (n_total * ncol(M))
  c(CE   = I(F0 - F1),
    P    = I(pmax(F0 - F1, 0)),
    N    = I(pmax(F1 - F0, 0)),
    P_UB = I(pmin(1 - F1, F0)),
    N_UB = I(pmin(1 - F0, F1)))
}

estimate_pn <- function(X, W, Y, x0, x1, w, h, a, b, ys) {
  Z <- cbind(X, W)
  keep <- x0 != x1
  if (!any(keep)) return(c(CE = 0, P = 0, N = 0, P_UB = 0, N_UB = 0))
  tgt <- function(x) cbind(x, matrix(w, length(x), length(w), byrow = TRUE))
  F0 <- cdf_matrix(Z, Y, tgt(x0[keep]), h, ys)
  F1 <- cdf_matrix(Z, Y, tgt(x1[keep]), h, ys)
  pn_from_cdfs(F0, F1, a, b, n_total = length(x0))
}

## Bandwidth rule h_j = 2 sd_j n^{-1/(p+4)}, p = number of LL regressors
bw_rule <- function(Z, c0 = 2) {
  Z <- as.matrix(Z)
  c0 * apply(Z, 2, sd) * nrow(Z)^(-1 / (ncol(Z) + 4))
}

###############################################################################
## 2. Simulations (Appendix E, Tables 2-3)
##    X := W + 3 U_X - 1, so do(X=0), do(X=2) and the Unif(0, 0.1) / +1.9 shift
##    lie inside the support of X given W = 0.5.
###############################################################################

R_SIM <- if (QUICK) 10 else 100
Ns    <- c(100, 1000, 10000)
N1    <- 10; N2 <- 1000; w0 <- 0.5

gen <- function(n, noise) {
  W  <- runif(n); UX <- runif(n); UY <- runif(n)
  X  <- W + 3 * UX - 1
  Y  <- (0.5 * X + 0.1 * W + 1) * (0.5 - UY) + if (noise) runif(n) else 0
  list(X = X, W = W, Y = Y)
}

one_rep <- function(n, noise) {
  d  <- gen(n, noise)
  hs <- bw_rule(cbind(d$X, d$W))
  a  <- min(d$Y); b <- max(d$Y)
  ys <- runif(N2, a, b)
  x0 <- runif(N1, 0, 0.1); x1 <- x0 + 1.9
  ce <- estimate_pn(d$X, d$W, d$Y, 0, 2, w0, hs, a, b, ys)
  pi <- estimate_pn(d$X, d$W, d$Y, x0, x1, w0, hs, a, b, ys)
  c(CACE = ce[["CE"]], PC = ce[["P"]], NC = ce[["N"]], PC_UB = ce[["P_UB"]], NC_UB = ce[["N_UB"]],
    CPICE = pi[["CE"]], PP = pi[["P"]], NP = pi[["N"]], PP_UB = pi[["P_UB"]], NP_UB = pi[["N_UB"]])
}

run_setting <- function(noise, seed) {
  set.seed(seed)
  res <- lapply(Ns, function(n) {
    t0 <- Sys.time()
    m <- t(replicate(R_SIM, one_rep(n, noise)))
    message(sprintf("  simulation (noise = %s) N = %5d: %.0fs", noise, n,
                    as.numeric(Sys.time() - t0, units = "secs")))
    m
  })
  names(res) <- Ns
  res
}

## LaTeX number: proper minus sign
num <- function(v) { s <- sprintf("%.3f", v); s[s == "-0.000"] <- "0.000"; sub("^-", "$-$", s) }
cell <- function(v) sprintf("%s ([%s, %s])", num(mean(v)), num(quantile(v, .025)), num(quantile(v, .975)))
rows_tex <- function(res, keys, labels, target) {
  vapply(seq_along(keys), function(i)
    sprintf("%s & %s & %s \\\\", labels[i],
            paste(vapply(res, function(m) cell(m[, keys[i]]), ""), collapse = " & "),
            num(target[[keys[i]]])), "")
}

message("Tables 2-3 ...")
res2 <- run_setting(noise = FALSE, seed = 2027)
res3 <- run_setting(noise = TRUE,  seed = 2028)

## True values (Table 2) and true sharp bounds (Table 3) at W = 0.5
target2 <- c(CACE = 0, PC = 0.125, NC = 0.125, CPICE = 0, PP = 0.119, NP = 0.119)
target3 <- c(CACE = 0, PC = 0.106, PC_UB = 0.448, NC = 0.106, NC_UB = 0.448,
             CPICE = 0, PP = 0.101, PP_UB = 0.447, NP = 0.101, NP_UB = 0.447)

tab2 <- rows_tex(res2, c("CACE", "PC", "NC", "CPICE", "PP", "NP"),
                 c("CACE", "P-CACE", "N-CACE", "CPICE", "P-CPICE", "N-CPICE"), target2)
tab3 <- rows_tex(res3, c("CACE", "PC", "PC_UB", "NC", "NC_UB", "CPICE", "PP", "PP_UB", "NP", "NP_UB"),
                 c("CACE", "P-CACE (LB)", "P-CACE (UB)", "N-CACE (LB)", "N-CACE (UB)",
                   "CPICE", "P-CPICE (LB)", "P-CPICE (UB)", "N-CPICE (LB)", "N-CPICE (UB)"), target3)
writeLines(tab2, file.path(OUT, "table2_rows.tex"))
writeLines(tab3, file.path(OUT, "table3_rows.tex"))
saveRDS(list(res2 = res2, res3 = res3), file.path(OUT, "simulations.rds"))

sim_identity <- max(sapply(c(res2, res3), function(m)
  max(abs(m[, "PC_UB"] - m[, "NC_UB"] - m[, "CACE"]), abs(m[, "PP_UB"] - m[, "NP_UB"] - m[, "CPICE"]))))

###############################################################################
## 3. Medical-cost application (Section 5)
##    w = (30 years, non-smoker, male, 1 child); X = BMI; Y = charges ($).
##    Exact matching on sex, smoking, children (n = 130); LL in (BMI, age).
##    CACE: do(BMI = 20) vs do(BMI = 40).
##    CPICE: X^{pi1} = X + 3 if X + 3 <= u(w) (95th pct of BMI | w), else X.
###############################################################################

message("Section 5 ...")
if (!file.exists("insurance.csv"))
  download.file("https://raw.githubusercontent.com/stedy/Machine-Learning-with-R-datasets/master/insurance.csv",
                "insurance.csv", quiet = TRUE)

B     <- if (QUICK) 30 else 1000
N1_RD <- 200; N2_RD <- 4000
shift <- 3; x_lo <- 20; x_hi <- 40; w_age <- 30

raw <- read.csv("insurance.csv")
n_exact <- sum(raw$age == 30 & raw$sex == "male" & raw$smoker == "no" & raw$children == 1)
d  <- raw[raw$smoker == "no" & raw$sex == "male" & raw$children == 1, ]
Xr <- d$bmi; Wr <- d$age; Yr <- d$charges
h  <- bw_rule(cbind(Xr, Wr))
ess <- function(k) sum(k)^2 / sum(k^2)
ess_lo <- ess(kern_w(cbind(Xr, Wr), c(x_lo, w_age), h))
ess_hi <- ess(kern_w(cbind(Xr, Wr), c(x_hi, w_age), h))

wquant <- function(x, wt, p) { o <- order(x); x[o][which(cumsum(wt[o]) / sum(wt) >= p)[1]] }

estimate_once <- function(X, W, Y, seed) {
  set.seed(seed)
  a  <- min(Y); b <- max(Y)
  ys <- runif(N2_RD, a, b)
  kw <- epan((W - w_age) / h[2])
  u  <- wquant(X, kw, 0.95)
  x0 <- sample(X, N1_RD, replace = TRUE, prob = kw)
  x1 <- ifelse(x0 + shift <= u, x0 + shift, x0)
  ce <- estimate_pn(X, W, Y, x_lo, x_hi, w_age, h, a, b, ys)
  pi <- estimate_pn(X, W, Y, x0, x1, w_age, h, a, b, ys)
  c(CACE = ce[["CE"]], P_CACE = ce[["P"]], N_CACE = ce[["N"]], P_CACE_UB = ce[["P_UB"]], N_CACE_UB = ce[["N_UB"]],
    CPICE = pi[["CE"]], P_CPICE = pi[["P"]], N_CPICE = pi[["N"]], P_CPICE_UB = pi[["P_UB"]], N_CPICE_UB = pi[["N_UB"]],
    shifted = mean(x0 != x1), u = u)
}

full <- estimate_once(Xr, Wr, Yr, seed = 1)
set.seed(2027)
boot <- t(sapply(seq_len(B), function(r) {
  i <- sample(nrow(d), replace = TRUE)
  tryCatch(estimate_once(Xr[i], Wr[i], Yr[i], seed = r), error = function(e) rep(NA, length(full)))
}))
nfail <- sum(!complete.cases(boot))
boot  <- boot[complete.cases(boot), , drop = FALSE]
ci    <- apply(boot, 2, quantile, c(.025, .975))
est   <- data.frame(estimate = full, lo = ci[1, ], hi = ci[2, ])[1:10, ]

dollar <- function(v) paste0(ifelse(round(v) < 0, "$-$", ""), "\\$",
                             formatC(abs(round(v)), format = "d", big.mark = "{,}"))
line <- function(k, lab) sprintf("%s: %s (CI: [%s, %s])", lab,
                                 dollar(est[k, 1]), dollar(est[k, 2]), dollar(est[k, 3]))
sec5 <- c(line("CACE", "CACE"), line("P_CACE", "P-CACE"), line("N_CACE", "N-CACE"),
          line("CPICE", "CPICE"), line("P_CPICE", "P-CPICE"), line("N_CPICE", "N-CPICE"),
          line("P_CACE_UB", "Upper bound of P-CACE"), line("N_CACE_UB", "Upper bound of N-CACE"),
          line("P_CPICE_UB", "Upper bound of P-CPICE"), line("N_CPICE_UB", "Upper bound of N-CPICE"),
          sprintf("%% n = %d (exact matches to w: %d); bandwidths BMI %.1f, age %.1f; ESS at BMI 20/40: %.1f/%.1f",
                  nrow(d), n_exact, h[1], h[2], ess_lo, ess_hi),
          sprintf("%% u(w) = %.1f; share shifted = %.2f; bootstrap: %d used, %d singular (discarded)",
                  full[["u"]], full[["shifted"]], nrow(boot), nfail))
writeLines(sec5, file.path(OUT, "section5_numbers.tex"))
saveRDS(list(full = full, boot = boot, h = h, nfail = nfail), file.path(OUT, "section5.rds"))

###############################################################################
## 4. Summary
###############################################################################

cat("\n==== Identity checks (should be ~0) ====\n")
cat(sprintf("simulations: max |P_UB - N_UB - CE| = %.1e\n", sim_identity))
cat(sprintf("Section 5:   P_UB - N_UB - CACE = %.1e ; P_CPICE_UB - N_CPICE_UB - CPICE = %.1e\n",
            full[["P_CACE_UB"]] - full[["N_CACE_UB"]] - full[["CACE"]],
            full[["P_CPICE_UB"]] - full[["N_CPICE_UB"]] - full[["CPICE"]]))
cat("\n==== Table 2 rows ====\n");   cat(tab2, sep = "\n")
cat("\n==== Table 3 rows ====\n");   cat(tab3, sep = "\n")
cat("\n==== Section 5 ====\n");      cat(sec5, sep = "\n")
cat(sprintf("\nfinished in %.1f min; files in ./%s/\n",
            as.numeric(Sys.time() - t_start, units = "mins"), OUT))
if (QUICK) cat("(quick mode: 10 simulation replications, B = 30 -- rerun without 'quick' for the paper)\n")
