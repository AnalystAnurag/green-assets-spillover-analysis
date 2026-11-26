options(width=100, digits=5)
# ---- step01_setup ----
# =====================================================================
# DYNAMIC CONNECTEDNESS AND VOLATILITY SPILLOVERS AMONG CONVENTIONAL AND
# GREEN FINANCIAL MARKETS: A DIEBOLD-YILMAZ ANALYSIS OF THE S&P 500,
# COMMODITIES, BITCOIN AND GREEN ASSETS
# Author: Anurag | Six assets | Daily data 2017-03-07 to 2025-12-31
# (S&P 500 replaces NIFTY 50; USD/INR removed; India-specific events removed)
# Everything below is base R + ggplot2 + (urca, tseries, quantreg, sandwich).
# The VAR, diagnostics, GFEVD, TVP-VAR, quantile and frequency routines are
# written from first principles (adapted unchanged from the earlier seven-asset script (NIFTY and USD/INR version)).
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(urca); library(tseries)
  library(quantreg); library(sandwich); library(scales); library(ggrepel)
})
set.seed(2026)

assets  <- c("SP500", "GOLD", "CRUDE", "BITCOIN", "GREEN_EQ", "GREEN_BOND")
# CRUDE: the file supplied as "BZ=F" contains WTI front-month futures prices (e.g. close 18.27 on 17-Apr-2020 and -37.63 on 20-Apr-2020),
# so the series is treated and reported as WTI crude oil (ticker CL=F).
tickers <- c(SP500 = "^GSPC", GOLD = "GC=F", CRUDE = "CL=F",
             BITCOIN = "BTC-USD", GREEN_EQ = "ICLN", GREEN_BOND = "GRNB")
files   <- c(SP500 = "SP500", GOLD = "GOLD", CRUDE = "CRUDE",
             BITCOIN = "BITCOIN", GREEN_EQ = "GREEN_EQ", GREEN_BOND = "GREEN_BOND")
K <- length(assets)
sample_start <- as.Date("2017-03-07")   # first GRNB price -> binding start date
sample_end   <- as.Date("2025-12-31")
H <- 10; W <- 250; BURN <- 50; BLOCK <- 20
# Bootstrap replications actually used (block bootstrap, resampling rows in blocks):
#   B_STAT = 1000 (static TCI/NET/pairwise, permutation test, TCI difference, incremental value, hedging, MCP);
#   300 (quantile and frequency measures, finite-sample no-spillover benchmark); 500 (200-day event test, block 10).
for (d in c("figures", "tables", "logs")) dir.create(d, showWarnings = FALSE)

# Events: 4 core events (event-window tests) + 3 context events (plot annotation only). No India-specific events.
events <- data.frame(
  id    = 1:7,
  Event = c("COVID-19 pandemic declared", "Russia-Ukraine war", "Fed first rate hike",
            "FTX collapse", "SVB collapse", "US spot Bitcoin ETFs", "US tariff shock"),
  Date  = as.Date(c("2020-03-11", "2022-02-24", "2022-03-16", "2022-11-11",
                    "2023-03-10", "2024-01-10", "2025-04-02")),
  Type  = c("Core", "Core", "Context", "Context", "Core", "Core", "Context"),
  stringsAsFactors = FALSE)
core_ev <- events[events$Type == "Core", ]

pal <- c(SP500 = "#1F77B4", GOLD = "#D4A017", CRUDE = "#333333",
         BITCOIN = "#E66101", GREEN_EQ = "#00A651", GREEN_BOND = "#17BECF")

# ---- step01b_helpers ----
# == plotting helpers (every time-series plot carries the event markers) ==
ev_caption <- function() {
  s <- paste0(events$id, " = ", events$Event, " (", format(events$Date, "%d %b %Y"), ")")
  paste0("Events (red = core, grey = context): ",
         paste(strwrap(paste(s, collapse = "; "), width = 140), collapse = "\n"))
}
add_events <- function() {
  c_ <- events[events$Type == "Core", ]; x_ <- events[events$Type != "Core", ]
  list(
    geom_vline(data = x_, aes(xintercept = Date), colour = "grey45", linetype = "dotted", linewidth = 0.35),
    geom_vline(data = c_, aes(xintercept = Date), colour = "#B2182B", linetype = "dashed", linewidth = 0.4),
    geom_text(data = x_, aes(x = Date, y = Inf, label = id), colour = "grey35", vjust = 1.3, size = 2.7,
              fontface = "bold", inherit.aes = FALSE),
    geom_text(data = c_, aes(x = Date, y = Inf, label = id), colour = "#B2182B", vjust = 1.3, size = 2.9,
              fontface = "bold", inherit.aes = FALSE))
}
theme_dy <- function(bs = 10) theme_minimal(base_size = bs) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 6.5, colour = "grey30"),
        legend.position = "bottom", strip.text = element_text(face = "bold"))
show <- function(p, name, w = 10, h = 5) {
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 150,
         device = ragg::agg_png)
  if (!isTRUE(getOption("knitr.in.progress"))) print(p)   # when knitting, the figure is inserted from the saved PNG
  invisible(p)
}
save_tab <- function(x, name) { write.csv(x, file.path("tables", paste0(name, ".csv")), row.names = TRUE); invisible(x) }

# == time-series helpers ==
skew <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
exkurt <- function(x) { x <- x - mean(x); mean(x^4) / mean(x^2)^2 - 3 }
arch_lm <- function(x, q) {                      # Engle (1982) ARCH-LM on a univariate series
  e2 <- (x - mean(x))^2; X <- embed(e2, q + 1)
  r2 <- summary(lm(X[, 1] ~ X[, -1]))$r.squared; LM <- nrow(X) * r2
  c(LM = LM, p = pchisq(LM, q, lower.tail = FALSE))
}

# == VAR by OLS (equation-by-equation; identical to the vars package) ==
var_fit <- function(Y, p, const = TRUE) {
  Y <- as.matrix(Y); n_all <- nrow(Y); Kv <- ncol(Y)
  Z <- do.call(cbind, lapply(1:p, function(l) Y[(p - l + 1):(n_all - l), , drop = FALSE]))
  X <- if (const) cbind(Z, 1) else Z
  Yt <- Y[(p + 1):n_all, , drop = FALSE]; n <- nrow(Yt)
  B <- solve(crossprod(X), crossprod(X, Yt)); U <- Yt - X %*% B
  A <- lapply(1:p, function(l) t(B[((l - 1) * Kv + 1):(l * Kv), , drop = FALSE]))
  list(A = A, B = B, U = U, X = X, Yt = Yt, n = n, K = Kv, p = p,
       Sigma = crossprod(U) / (n - ncol(X)), SigmaML = crossprod(U) / n,
       cnst = if (const) B[nrow(B), ] else rep(0, Kv))
}
var_ic <- function(Y, pmax) {                    # same effective sample for all p (as vars::VARselect)
  out <- sapply(1:pmax, function(p) {
    f <- var_fit(Y[(pmax - p + 1):nrow(Y), ], p); n <- f$n; Kv <- f$K; d <- log(det(f$SigmaML))
    npar <- p * Kv^2 + Kv
    c(AIC = d + 2 * npar / n, HQ = d + 2 * log(log(n)) * npar / n, SC = d + log(n) * npar / n,
      FPE = ((n + Kv * p + 1) / (n - Kv * p - 1))^Kv * det(f$SigmaML))
  })
  colnames(out) <- 1:pmax; t(out)
}
companion_roots <- function(f) {
  Kv <- f$K; p <- f$p; C <- matrix(0, Kv * p, Kv * p); C[1:Kv, ] <- do.call(cbind, f$A)
  if (p > 1) C[(Kv + 1):(Kv * p), 1:(Kv * (p - 1))] <- diag(Kv * (p - 1))
  Mod(eigen(C, only.values = TRUE)$values)
}
portmanteau <- function(f, h) {                  # adjusted multivariate Portmanteau (Ljung-Box type)
  U <- f$U; n <- nrow(U); C0i <- solve(crossprod(U) / n); Q <- 0
  for (i in 1:h) { Ci <- crossprod(U[(i + 1):n, , drop = FALSE], U[1:(n - i), , drop = FALSE]) / n
    Q <- Q + sum(diag(t(Ci) %*% C0i %*% Ci %*% C0i)) / (n - i) }
  Q <- Q * n^2; df <- f$K^2 * (h - f$p)
  c(stat = Q, df = df, p = pchisq(Q, df, lower.tail = FALSE))
}
arch_mv <- function(f, q) {                      # multivariate ARCH-LM (Luetkepohl 2005, as vars::arch.test)
  U <- f$U; Kv <- ncol(U); ij <- which(lower.tri(diag(Kv), diag = TRUE), arr.ind = TRUE)
  V <- sapply(seq_len(nrow(ij)), function(k) U[, ij[k, 1]] * U[, ij[k, 2]]); n <- nrow(V)
  Yv <- V[(q + 1):n, , drop = FALSE]
  X <- cbind(1, do.call(cbind, lapply(1:q, function(l) V[(q + 1 - l):(n - l), , drop = FALSE])))
  E <- qr.resid(qr(X), Yv); m <- nrow(Yv)
  Om <- crossprod(E) / m; Om0 <- crossprod(scale(Yv, scale = FALSE)) / m
  R2 <- 1 - 2 / (Kv * (Kv + 1)) * sum(diag(Om %*% solve(Om0)))
  LM <- 0.5 * m * Kv * (Kv + 1) * R2; df <- q * (Kv * (Kv + 1) / 2)^2
  c(stat = LM, df = df, p = pchisq(LM, df, lower.tail = FALSE))
}
jb_mv <- function(f) {                           # multivariate Jarque-Bera, Cholesky-standardised (Luetkepohl)
  U <- scale(f$U, scale = FALSE); n <- nrow(U); P <- t(chol(crossprod(U) / n))
  Wm <- t(solve(P, t(U))); b1 <- colMeans(Wm^3); b2 <- colMeans(Wm^4)
  s <- n * sum(b1^2) / 6; k <- n * sum((b2 - 3)^2) / 24; df <- 2 * f$K
  c(skew = s, kurt = k, JB = s + k, df = df, p = pchisq(s + k, df, lower.tail = FALSE))
}
granger_mat <- function(f) {                     # heteroskedasticity-robust (HC3) block Wald tests
  Kv <- f$K; p <- f$p; X <- f$X; XtXi <- solve(crossprod(X)); hat <- rowSums((X %*% XtXi) * X)
  P <- matrix(NA, Kv, Kv, dimnames = list(NULL, NULL))   # P[i, j] = p-value: j does not Granger-cause i
  for (i in 1:Kv) {
    b <- f$B[, i]; u <- f$U[, i]; Xu <- X * (u / (1 - hat)); V <- XtXi %*% crossprod(Xu) %*% XtXi
    for (j in 1:Kv) if (j != i) {
      pos <- j + Kv * (0:(p - 1)); r <- b[pos]; Vr <- V[pos, pos, drop = FALSE]
      P[i, j] <- pchisq(drop(t(r) %*% solve(Vr) %*% r), p, lower.tail = FALSE) }
  }
  P
}

# == Generalised FEVD (Pesaran-Shin 1998), DY (2012) normalisation ==
ma_coefs <- function(A, Hh) {
  Kv <- nrow(A[[1]]); p <- length(A); Phi <- vector("list", Hh); Phi[[1]] <- diag(Kv)
  for (h in 2:Hh) { S <- matrix(0, Kv, Kv); for (l in 1:min(h - 1, p)) S <- S + A[[l]] %*% Phi[[h - l]]; Phi[[h]] <- S }
  Phi
}
gfevd <- function(A, Sigma, Hh = 10) {           # returns row-normalised theta[i, j]: share of i's FEV due to j
  Phi <- ma_coefs(A, Hh); Kv <- nrow(Sigma); num <- matrix(0, Kv, Kv); den <- numeric(Kv)
  for (h in 1:Hh) { PS <- Phi[[h]] %*% Sigma; num <- num + PS^2; den <- den + rowSums(PS * Phi[[h]]) }
  th <- sweep(num, 2, diag(Sigma), "/") / den
  th / rowSums(th)
}
dy_summary <- function(th, nm = assets) {
  Kv <- nrow(th); off <- th; diag(off) <- 0
  FROM <- 100 * rowSums(off); TO <- 100 * colSums(off)
  NPDC <- 100 * (t(th) - th)                     # NPDC[i, j] = net transmission from i to j
  tab <- 100 * th; dimnames(tab) <- list(nm, nm); dimnames(NPDC) <- list(nm, nm)
  names(FROM) <- names(TO) <- nm
  list(theta = th, table = tab, FROM = FROM, TO = TO, NET = TO - FROM, TCI = 100 * sum(off) / Kv, NPDC = NPDC)
}
static_dy <- function(Y, p, Hh = 10) { f <- var_fit(Y, p); dy_summary(gfevd(f$A, f$Sigma, Hh), colnames(Y)) }

# == rolling-window DY ==
summ_array <- function(TH, dates, nm = assets) {
  Kv <- dim(TH)[1]; n <- dim(TH)[3]
  TCI <- numeric(n); TO <- FROM <- matrix(NA, n, Kv, dimnames = list(NULL, nm))
  for (s in 1:n) { o <- TH[, , s]; diag(o) <- 0; TCI[s] <- 100 * sum(o) / Kv; FROM[s, ] <- 100 * rowSums(o); TO[s, ] <- 100 * colSums(o) }
  list(theta = TH, dates = dates, TCI = TCI, TO = TO, FROM = FROM, NET = TO - FROM)
}
rolling_dy <- function(Y, dates, p, Hh, w) {
  n_all <- nrow(Y); idx <- w:n_all; TH <- array(NA, c(ncol(Y), ncol(Y), length(idx)))
  for (s in seq_along(idx)) { t <- idx[s]; f <- var_fit(Y[(t - w + 1):t, ], p); TH[, , s] <- gfevd(f$A, f$Sigma, Hh) }
  summ_array(TH, dates[idx], colnames(Y))
}
pair_series <- function(R, from, to) 100 * (R$theta[match(to, assets), match(from, assets), ] -
                                              R$theta[match(from, assets), match(to, assets), ])
# pair_series(R, i, j) = net pairwise transmission from i to j  (= theta[j, i] - theta[i, j])

# == moving-block bootstrap on the rows of a matrix ==
mbb_idx <- function(n, l = BLOCK) {
  nb <- ceiling(n / l); st <- sample.int(n - l + 1, nb, replace = TRUE)
  head(as.vector(sapply(st, function(s) s:(s + l - 1))), n)
}
q_ci <- function(x, lev = 0.95) quantile(x, c((1 - lev) / 2, 1 - (1 - lev) / 2), na.rm = TRUE)

# ---- step02_load_audit ----
# STEP 2 - Load the OHLC data and audit its quality (before touching it)
raw <- lapply(assets, function(a) {
  d <- read.csv(file.path("data", paste0(files[a], ".csv")), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  # Data repair: a futures price cannot be zero or negative in a log-return model. crude-oil file contains the
  # 20-21 April 2020 bars with negative Open/Low/Close (the WTI May-2020 contract collapse).
  # Bars with any non-positive OHLC value are treated as missing; the common-calendar step then drops those dates for all assets.
  bad <- !(d$Open > 0 & d$High > 0 & d$Low > 0 & d$Close > 0 & d$Adjusted > 0)
  bad[is.na(bad)] <- FALSE
  if (any(bad)) cat("Non-positive price bars set to missing in", a, ":", paste(format(d$Date[bad]), collapse = ", "), "\n")
  d[!bad, ] })
names(raw) <- assets

audit <- do.call(rbind, lapply(assets, function(a) {
  d <- raw[[a]]; ok <- complete.cases(d[, c("Open", "High", "Low", "Close")])
  bad <- ok & (d$High < pmax(d$Open, d$Close) - 1e-9 * d$High | d$Low > pmin(d$Open, d$Close) + 1e-9 * d$Low)
  data.frame(Asset = a, Ticker = tickers[a], First = min(d$Date), Last = max(d$Date), Rows = nrow(d),
             Missing_rows = sum(!ok), Duplicate_dates = sum(duplicated(d$Date)),
             OHLC_inconsistent = sum(bad), Zero_range_days = sum(ok & d$High == d$Low))
}))
rownames(audit) <- NULL
print(audit, row.names = FALSE)
save_tab(audit, "02_data_audit")

# ---- step03_common_sample ----
# STEP 3 - Common calendar: keep only dates on which ALL six assets have a price
px_list <- lapply(assets, function(a) setNames(raw[[a]][, c("Date", "Adjusted")], c("Date", a)))
px_all  <- Reduce(function(x, y) merge(x, y, by = "Date", all = TRUE), px_list)
px_all  <- px_all[px_all$Date >= sample_start & px_all$Date <= sample_end, ]
cat("Union of dates in window:", nrow(px_all), "\n")
miss_tab <- data.frame(Asset = assets, Missing_on_union_calendar = colSums(is.na(px_all[, assets])))
print(miss_tab, row.names = FALSE)
px <- px_all[complete.cases(px_all), ]
rownames(px) <- NULL
dates <- px$Date; Tn <- nrow(px)
P <- as.matrix(px[, assets]); rownames(P) <- as.character(dates)
cat("Common sample:", format(min(dates)), "to", format(max(dates)), "| price obs:", Tn,
    "| return obs:", Tn - 1, "\n")
gap <- as.integer(diff(dates))
cat("Calendar-day gaps between consecutive common dates:\n"); print(table(gap))
stopifnot(all(P > 0), !anyDuplicated(dates))

# Cleaned OHLC aligned to the common calendar (used for range-based volatility)
ohlc <- lapply(assets, function(a) {
  d <- raw[[a]][match(dates, raw[[a]]$Date), c("Open", "High", "Low", "Close")]
  hi <- pmax(d$Open, d$High, d$Low, d$Close); lo <- pmin(d$Open, d$High, d$Low, d$Close)
  d$High <- hi; d$Low <- lo; d }) # repairs any bar with High < max(Open, Close) or Low > min(Open, Close)
names(ohlc) <- assets
stopifnot(!anyNA(unlist(ohlc)))

# Prices rebased to 100 with event markers
pl <- data.frame(Date = rep(dates, K), Asset = factor(rep(assets, each = Tn), levels = assets),
                 Index = as.vector(sweep(P, 2, P[1, ], "/") * 100))
p_px <- ggplot(pl, aes(Date, Index)) + geom_line(linewidth = 0.4, colour = "#1F4E79") +
  add_events() + facet_wrap(~Asset, ncol = 2, scales = "free_y") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "Adjusted prices rebased to 100 (7 March 2017 = 100)", x = NULL, y = "Index level",
       caption = ev_caption()) + theme_dy(9)
show(p_px, "03_prices_events", 10, 8.5)

# ---- step04_returns ----
# STEP 4 - Daily log returns in percent + outlier audit
lP <- log(P); R <- 100 * diff(lP); rdates <- dates[-1]; Tr <- nrow(R)
cat("Return observations:", Tr, "x", ncol(R), "assets\n")
out <- do.call(rbind, lapply(assets, function(a) {
  o <- order(-abs(R[, a]))[1:3]
  data.frame(Asset = a, Date = rdates[o], Return_pct = round(R[o, a], 2),
             Within_60d_of_event = sapply(rdates[o], function(d) { k <- which.min(abs(as.numeric(events$Date - d)))
               if (abs(as.numeric(d - events$Date[k])) <= 60) paste0(events$Event[k], " (", as.integer(d - events$Date[k]), "d)") else "-" })) }))
rownames(out) <- NULL; print(out, row.names = FALSE); save_tab(out, "04_extreme_returns")

rl <- data.frame(Date = rep(rdates, K), Asset = factor(rep(assets, each = Tr), levels = assets), Ret = as.vector(R))
p_ret <- ggplot(rl, aes(Date, Ret)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.25, colour = "#1F4E79") + add_events() +
  facet_wrap(~Asset, ncol = 2, scales = "free_y") + scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "Daily log returns (%)", x = NULL, y = "Return (%)", caption = ev_caption()) + theme_dy(9)
show(p_ret, "04_returns_events", 10, 8.5)

# ---- step05_volatility ----
# STEP 5 - Volatility proxies (daily standard-deviation scale, in %)
#   MAIN      : Garman-Klass (1980) range-based estimator, uses Open-High-Low-Close
#   ROBUSTNESS: Parkinson (1980) high-low estimator;  |r| and r^2 from close-to-close returns
gk_var <- function(d) 0.5 * log(d$High / d$Low)^2 - (2 * log(2) - 1) * log(d$Close / d$Open)^2
pk_var <- function(d) log(d$High / d$Low)^2 / (4 * log(2))
GKv <- sapply(assets, function(a) gk_var(ohlc[[a]]))[-1, ]
PKv <- sapply(assets, function(a) pk_var(ohlc[[a]]))[-1, ]
cat("Negative GK variance days (set to 0):", paste(assets, colSums(GKv < 0), collapse = " | "), "\n")
cat("Zero-range (GK = 0) days           :", paste(assets, colSums(GKv <= 0), collapse = " | "), "\n")
VOL  <- 100 * sqrt(pmax(GKv, 0))           # MAIN volatility proxy
VPK  <- 100 * sqrt(pmax(PKv, 0))           # Parkinson
VAR2 <- R^2                                # squared return (realised-variance proxy)
VABS <- abs(R)                             # absolute return
# Placeholder bars (Open = High = Low = Close) give a zero range although the close-to-close move is non-zero.
# These are data artefacts (stale Yahoo bars), so |return| is used as a floor on those days.
fix <- GKv <= 0 & VABS > 0
cat("Zero-range days repaired with |return| :", paste(assets, colSums(fix), collapse = " | "), "\n")
VOL[fix] <- VABS[fix]; VPK[fix] <- VABS[fix]
rownames(VOL) <- rownames(VPK) <- rownames(VAR2) <- rownames(VABS) <- as.character(rdates)
cat("Correlation of GK-vol with other proxies (average across assets):\n")
print(round(c(PK = mean(sapply(assets, function(a) cor(VOL[, a], VPK[, a]))),
              AbsRet = mean(sapply(assets, function(a) cor(VOL[, a], VABS[, a]))),
              SqRet = mean(sapply(assets, function(a) cor(VOL[, a], VAR2[, a])))), 3))

vl <- data.frame(Date = rep(rdates, K), Asset = factor(rep(assets, each = Tr), levels = assets),
                 Vol = as.vector(VOL * sqrt(252)))
p_vol <- ggplot(vl, aes(Date, Vol)) + geom_line(linewidth = 0.3, colour = "#B35806") + add_events() +
  facet_wrap(~Asset, ncol = 2, scales = "free_y") + scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "Garman-Klass volatility (annualised, % )", x = NULL, y = "Volatility (% p.a.)",
       caption = ev_caption()) + theme_dy(9)
show(p_vol, "05_gk_volatility_events", 10, 8.5)

# ---- step06_descriptives ----
# STEP 6 - Descriptive statistics, normality and serial dependence
desc <- function(M) {
  t(sapply(assets, function(a) { x <- M[, a]
    c(Mean = mean(x), SD = sd(x), Min = min(x), Max = max(x), Skew = skew(x), ExKurt = exkurt(x),
      JB = unname(jarque.bera.test(x)$statistic), JB_p = jarque.bera.test(x)$p.value,
      LB10 = Box.test(x, 10, "Ljung-Box")$statistic, LB10_p = Box.test(x, 10, "Ljung-Box")$p.value,
      LB10_sq = Box.test(x^2, 10, "Ljung-Box")$statistic, LB10_sq_p = Box.test(x^2, 10, "Ljung-Box")$p.value) }))
}
d_ret <- desc(R); d_vol <- desc(VOL)
cat("--- Daily log returns (%) ---\n"); print(round(d_ret, 4))
cat("\n--- Garman-Klass volatility (daily sd, %) ---\n"); print(round(d_vol, 4))
save_tab(d_ret, "06_desc_returns"); save_tab(d_vol, "06_desc_vol")

# Correlation heat-maps (returns and volatility)
cor_long <- function(M, nm) { cm <- cor(M); x <- as.data.frame(as.table(cm)); names(x) <- c("A", "B", "r"); x$Basis <- nm; x }
cl <- rbind(cor_long(R, "Returns"), cor_long(VOL, "GK volatility"))
cl$A <- factor(cl$A, assets); cl$B <- factor(cl$B, rev(assets))
p_cor <- ggplot(cl, aes(A, B, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.8) + facet_wrap(~Basis) + coord_equal() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
  labs(title = "Pearson correlations: returns vs volatility", x = NULL, y = NULL, fill = "r") +
  theme_dy(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
show(p_cor, "06_correlations", 10, 4.8)

# ---- step07_unitroot ----
# STEP 7 - Unit-root, stationarity and cointegration evidence
ur_row <- function(x, type, lab, asset) {
  a <- ur.df(x, type = type, lags = 15, selectlags = "BIC")
  tau <- a@teststat[1, 1]; cv5 <- a@cval[1, "5pct"]
  klag <- sum(grepl("z.diff.lag", rownames(a@testreg$coefficients)))
  k <- ur.kpss(x, type = if (type == "trend") "tau" else "mu", lags = "long")
  kcv5 <- k@cval[1, "5pct"]
  data.frame(Series = lab, Asset = asset, ADF_lags = klag, ADF_stat = tau, ADF_cv5 = cv5,
             ADF_reject_unit_root = tau < cv5, KPSS_stat = k@teststat, KPSS_cv5 = kcv5,
             KPSS_reject_stationarity = k@teststat > kcv5,
             PP_p = suppressWarnings(pp.test(x)$p.value))
}
ur_tab <- rbind(
  do.call(rbind, lapply(assets, function(a) ur_row(lP[, a], "trend", "Log price (trend)", a))),
  do.call(rbind, lapply(assets, function(a) ur_row(R[, a], "drift", "Log return", a))),
  do.call(rbind, lapply(assets, function(a) ur_row(VOL[, a], "drift", "GK volatility", a))))
rownames(ur_tab) <- NULL
print(ur_tab, digits = 3, row.names = FALSE); save_tab(ur_tab, "07_unit_root")

# Johansen trace test on log prices: if prices were cointegrated, a VAR in returns would omit the ECM term
jo <- ca.jo(lP, type = "trace", ecdet = "const", K = 2)
jo_tab <- data.frame(H0 = paste0("r <= ", (K - 1):0), Trace = jo@teststat, CV5 = jo@cval[, "5pct"])
jo_tab$Reject_5pct <- jo_tab$Trace > jo_tab$CV5
print(jo_tab, row.names = FALSE); save_tab(jo_tab, "07_johansen")

# ---- step08_arch ----
# STEP 8 - Volatility clustering: ARCH-LM and Ljung-Box on squared returns
arch_tab <- do.call(rbind, lapply(assets, function(a) {
  a5 <- arch_lm(R[, a], 5); a10 <- arch_lm(R[, a], 10)
  data.frame(Asset = a, ARCH5_LM = a5["LM"], ARCH5_p = a5["p"], ARCH10_LM = a10["LM"], ARCH10_p = a10["p"],
             LB_sq10 = Box.test(R[, a]^2, 10, "Ljung-Box")$statistic, LB_sq10_p = Box.test(R[, a]^2, 10, "Ljung-Box")$p.value) }))
rownames(arch_tab) <- NULL; print(arch_tab, digits = 4, row.names = FALSE); save_tab(arch_tab, "08_arch_lm")

acf_df <- do.call(rbind, lapply(assets, function(a) {
  r1 <- acf(R[, a]^2, lag.max = 25, plot = FALSE); data.frame(Asset = a, Lag = 1:25, ACF = r1$acf[-1]) }))
acf_df$Asset <- factor(acf_df$Asset, assets)
p_acf <- ggplot(acf_df, aes(Lag, ACF)) + geom_col(fill = "#1F4E79", width = 0.6) +
  geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(Tr), linetype = "dashed", colour = "#B2182B") +
  facet_wrap(~Asset, ncol = 4) + labs(title = "ACF of squared daily returns (95% bands dashed)", y = "ACF") + theme_dy(9)
show(p_acf, "08_acf_squared_returns", 10, 5)

# ---- step09_var_select ----
# STEP 9 - VAR lag selection (returns system and GK-volatility system)
sel_r <- var_ic(R, 10); sel_v <- var_ic(VOL, 10)
pick <- function(m) apply(m, 2, which.min)
cat("Returns system - lag chosen by criterion:\n"); print(pick(sel_r))
cat("\nVolatility system - lag chosen by criterion:\n"); print(pick(sel_v))
p_ret <- unname(pick(sel_r)["SC"]); p_vol <- unname(pick(sel_v)["SC"])
cat("\nBaseline lags (Schwarz/BIC, the consistent criterion):  returns p =", p_ret, "| volatility p =", p_vol, "\n")
ic_df <- rbind(data.frame(System = "Returns", Lag = 1:10, sel_r[, c("AIC", "HQ", "SC")]),
               data.frame(System = "GK volatility", Lag = 1:10, sel_v[, c("AIC", "HQ", "SC")]))
ic_long <- reshape(ic_df, direction = "long", varying = c("AIC", "HQ", "SC"), v.names = "Value",
                   timevar = "Criterion", times = c("AIC", "HQ", "SC"), idvar = c("System", "Lag"))
p_ic <- ggplot(ic_long, aes(Lag, Value, colour = Criterion)) + geom_line() + geom_point(size = 1.4) +
  facet_wrap(~System, scales = "free_y") + scale_x_continuous(breaks = 1:10) +
  labs(title = "VAR lag-selection criteria (lower is better)", y = "Criterion value") + theme_dy(9)
show(p_ic, "09_lag_selection", 9, 3.8)

# ---- step10_var_diag ----
# STEP 10 - VAR estimation, stability and residual diagnostics for lags 1..5
var_diag <- function(Y, ps = 1:5) {
  do.call(rbind, lapply(ps, function(p) {
    f <- var_fit(Y, p); pt <- portmanteau(f, 20); ar <- arch_mv(f, 5); jb <- jb_mv(f)
    data.frame(Lag = p, MaxRoot = max(companion_roots(f)), Stable = max(companion_roots(f)) < 1,
               Portmanteau20_p = pt["p"], MV_ARCH5_p = ar["p"], MV_JB_p = jb["p"],
               Resid_corr_max_abs = { cr <- cor(f$U); max(abs(cr[upper.tri(cr)])) }) }))
}
dg_r <- var_diag(R); dg_v <- var_diag(VOL)
cat("--- Returns VAR ---\n"); print(dg_r, digits = 3, row.names = FALSE)
cat("\n--- GK-volatility VAR ---\n"); print(dg_v, digits = 3, row.names = FALSE)
save_tab(dg_r, "10_var_diag_returns"); save_tab(dg_v, "10_var_diag_vol")

# Size (not just significance) of the remaining residual autocorrelation: largest |corr(u_t, u_t-h)| over all variable pairs
resid_acf_max <- function(f, L = 10) sapply(1:L, function(h) { U <- f$U; n <- nrow(U); D <- diag(1 / sqrt(diag(crossprod(U) / n)))
  max(abs(D %*% (crossprod(U[(h + 1):n, , drop = FALSE], U[1:(n - h), , drop = FALSE]) / n) %*% D)) })
fit_r <- var_fit(R, p_ret); fit_v <- var_fit(VOL, p_vol)
cat("\nLargest |residual cross-autocorrelation| at lags 1-10 (returns, baseline):\n"); print(round(resid_acf_max(fit_r), 3))
cat("Largest |residual cross-autocorrelation| at lags 1-10 (volatility, baseline):\n"); print(round(resid_acf_max(fit_v), 3))
cat("\nBaseline returns VAR(", p_ret, ") eigenvalue moduli:", round(sort(companion_roots(fit_r), decreasing = TRUE), 3), "\n")
cat("Baseline volatility VAR(", p_vol, ") largest 5 eigenvalue moduli:",
    round(head(sort(companion_roots(fit_v), decreasing = TRUE), 5), 3), "\n")
cat("Portmanteau h=20 (returns, baseline):", round(portmanteau(fit_r, 20), 4), "\n")
cat("Portmanteau h=20 (volatility, baseline):", round(portmanteau(fit_v, 20), 4), "\n")

# ---- step11_granger ----
# STEP 11 - Predictive causality (HC3-robust block Wald tests): who helps forecast whom?
gc_r <- granger_mat(fit_r); gc_v <- granger_mat(fit_v)
dimnames(gc_r) <- dimnames(gc_v) <- list(Effect_on = assets, Cause = assets)
cat("p-values: row = predicted variable, column = predictor (H0: column does not Granger-cause row)\n")
cat("\n--- Returns ---\n"); print(round(gc_r, 3), na.print = "-")
cat("\n--- GK volatility ---\n"); print(round(gc_v, 3), na.print = "-")
cat("\nSignificant (5%) directional links: returns =", sum(gc_r < 0.05, na.rm = TRUE), "of", K * (K - 1), "| volatility =",
    sum(gc_v < 0.05, na.rm = TRUE), "of", K * (K - 1), "\n")
gcl <- rbind(transform(as.data.frame(as.table(gc_r)), Basis = "Returns"),
             transform(as.data.frame(as.table(gc_v)), Basis = "GK volatility"))
names(gcl)[3] <- "p"; gcl <- gcl[!is.na(gcl$p), ]
gcl$Effect_on <- factor(gcl$Effect_on, rev(assets)); gcl$Cause <- factor(gcl$Cause, assets)
gcl$sig <- cut(gcl$p, c(-1, 0.01, 0.05, 0.1, 2), labels = c("p<0.01", "p<0.05", "p<0.10", "n.s."))
p_gc <- ggplot(gcl, aes(Cause, Effect_on, fill = sig)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3f", p)), size = 2.4) + facet_wrap(~Basis) + coord_equal() +
  scale_fill_manual(values = c("p<0.01" = "#B2182B", "p<0.05" = "#EF8A62", "p<0.10" = "#FDDBC7", "n.s." = "grey95")) +
  labs(title = "Granger-causality p-values (HC3 Wald): column helps predict row", x = "Predictor", y = "Predicted", fill = NULL) +
  theme_dy(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
show(p_gc, "11_granger", 10, 4.8)

# ---- step12_static_dy ----
# STEP 12 - Static Diebold-Yilmaz connectedness (GFEVD, H = 10 days)
dy_print <- function(s) {
  tab <- rbind(s$table, TO = s$TO); tab <- cbind(tab, FROM = c(s$FROM, NA))
  tab <- rbind(tab, NET = c(s$NET, NA)); tab
}
st_r <- static_dy(R, p_ret, H); st_v <- static_dy(VOL, p_vol, H)
cat("=== RETURNS: VAR(", p_ret, "), H =", H, "===\n", sep = ""); print(round(dy_print(st_r), 2), na.print = "")
cat("Total connectedness index (TCI):", round(st_r$TCI, 2), "%\n")
cat("\n=== GK VOLATILITY: VAR(", p_vol, "), H =", H, "===\n", sep = ""); print(round(dy_print(st_v), 2), na.print = "")
cat("Total connectedness index (TCI):", round(st_v$TCI, 2), "%\n")
save_tab(dy_print(st_r), "12_static_returns"); save_tab(dy_print(st_v), "12_static_vol")

# ---- step13_boot_ci ----
# STEP 13 - Moving-block bootstrap confidence intervals (block = 20 days, B = 1000)
B_STAT <- 1000
boot_static <- function(Y, p, Hh = H, B = B_STAT) {
  n <- nrow(Y); Kv <- ncol(Y); TCIb <- numeric(B); NETb <- TOb <- FROMb <- matrix(NA, B, Kv); THb <- array(NA, c(Kv, Kv, B))
  for (b in 1:B) { s <- static_dy(Y[mbb_idx(n), ], p, Hh); TCIb[b] <- s$TCI; NETb[b, ] <- s$NET; TOb[b, ] <- s$TO; FROMb[b, ] <- s$FROM; THb[, , b] <- s$theta }
  list(TCI = TCIb, NET = NETb, TO = TOb, FROM = FROMb, theta = THb)
}
bs_r <- boot_static(R, p_ret); bs_v <- boot_static(VOL, p_vol)
ci_tab <- function(s, bs, nm) {
  data.frame(Basis = nm, Asset = c(assets, "TCI"),
             Estimate = c(s$NET, s$TCI),
             Lower95 = c(apply(bs$NET, 2, quantile, 0.025), quantile(bs$TCI, 0.025)),
             Upper95 = c(apply(bs$NET, 2, quantile, 0.975), quantile(bs$TCI, 0.975)),
             row.names = NULL)
}
ci_all <- rbind(ci_tab(st_r, bs_r, "Returns"), ci_tab(st_v, bs_v, "GK volatility"))
ci_all$Role <- ifelse(ci_all$Asset == "TCI", "-", ifelse(ci_all$Lower95 > 0, "Significant net transmitter",
               ifelse(ci_all$Upper95 < 0, "Significant net receiver", "Not different from zero")))
print(ci_all, digits = 3, row.names = FALSE); save_tab(ci_all, "13_static_ci")

# ---- step14_perm_test ----
# STEP 14 - Is total connectedness larger than chance?  Circular-shift null: every series is shifted by an independent
# random lag, which keeps each series' own dynamics (clustering, fat tails) but destroys all cross-market links.
perm_tci <- function(Y, p, B = B_STAT) {
  n <- nrow(Y); out <- numeric(B)
  for (b in 1:B) { Yb <- sapply(1:ncol(Y), function(j) Y[((seq_len(n) + sample(50:(n - 50), 1) - 1) %% n) + 1, j])
    out[b] <- static_dy(Yb, p, H)$TCI }
  out
}
nl_r <- perm_tci(R, p_ret); nl_v <- perm_tci(VOL, p_vol)
perm_tab <- data.frame(Basis = c("Returns", "GK volatility"), TCI = c(st_r$TCI, st_v$TCI),
                       Null_mean = c(mean(nl_r), mean(nl_v)), Null_95th = c(quantile(nl_r, .95), quantile(nl_v, .95)),
                       Bias_adjusted_TCI = c(st_r$TCI - mean(nl_r), st_v$TCI - mean(nl_v)),
                       p_value = c((1 + sum(nl_r >= st_r$TCI)) / (B_STAT + 1), (1 + sum(nl_v >= st_v$TCI)) / (B_STAT + 1)))
print(perm_tab, digits = 3, row.names = FALSE); save_tab(perm_tab, "14_perm_test")
diff_tci <- st_v$TCI - st_r$TCI
bs_diff <- bs_v$TCI - bs_r$TCI    # paired? (independent resamples) -> use same-index resample below
dd <- numeric(B_STAT); n_ <- nrow(R)
for (b in 1:B_STAT) { ix <- mbb_idx(n_); dd[b] <- static_dy(VOL[ix, ], p_vol, H)$TCI - static_dy(R[ix, ], p_ret, H)$TCI }
cat("\nVolatility TCI minus returns TCI =", round(diff_tci, 2), "pp; paired block-bootstrap 95% CI [",
    round(quantile(dd, .025), 2), ",", round(quantile(dd, .975), 2), "]\n")
tci_diff <- c(diff = diff_tci, lo = unname(quantile(dd, .025)), hi = unname(quantile(dd, .975)))

# ---- step15_static_plots ----
# STEP 15 - Static results in pictures: NET with 95% CI, pairwise shares, and green-asset sources
p_net <- ggplot(ci_all[ci_all$Asset != "TCI", ], aes(reorder(Asset, Estimate), Estimate, fill = Estimate > 0)) +
  geom_col(width = 0.65) + geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.25) + geom_hline(yintercept = 0) +
  coord_flip() + facet_wrap(~Basis, scales = "free_x") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), labels = c("Net receiver", "Net transmitter")) +
  labs(title = "Net connectedness (TO - FROM) with 95% block-bootstrap intervals", x = NULL, y = "NET (percentage points)", fill = NULL) + theme_dy(9)
show(p_net, "15_net_static", 10, 4.2)

pw <- function(s, nm) { x <- as.data.frame(as.table(s$table)); names(x) <- c("Receiver", "Transmitter", "Share"); x$Basis <- nm; x[x$Receiver != x$Transmitter, ] }
hd <- rbind(pw(st_r, "Returns"), pw(st_v, "GK volatility"))
hd$Receiver <- factor(hd$Receiver, rev(assets)); hd$Transmitter <- factor(hd$Transmitter, assets)
p_heat <- ggplot(hd, aes(Transmitter, Receiver, fill = Share)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", Share)), size = 2.7) + facet_wrap(~Basis) + coord_equal() +
  scale_fill_gradient(low = "lightyellow", high = "darkred") +
  labs(title = "Pairwise directional connectedness: share of the receiver's 10-day forecast-error variance due to the transmitter (%)",
       x = "Transmitter", y = "Receiver") + theme_dy(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 9))
show(p_heat, "15_pairwise_heatmap", 11, 4.8)

# RQ3: sources of spillovers to each green asset (with bootstrap CIs)
src <- do.call(rbind, lapply(c("GREEN_EQ", "GREEN_BOND"), function(g) do.call(rbind, lapply(list(list(st_r, bs_r, "Returns"), list(st_v, bs_v, "GK volatility")), function(z) {
  gi <- match(g, assets); others <- setdiff(1:K, gi)
  data.frame(Receiver = g, Basis = z[[3]], Source = assets[others], Share = 100 * z[[1]]$theta[gi, others],
             Lo = 100 * apply(z[[2]]$theta[gi, others, ], 1, quantile, 0.025),
             Hi = 100 * apply(z[[2]]$theta[gi, others, ], 1, quantile, 0.975),
             Own = 100 * z[[1]]$theta[gi, gi]) }))))
rownames(src) <- NULL
print(transform(src, Share = round(Share, 2), Lo = round(Lo, 2), Hi = round(Hi, 2), Own = round(Own, 1)), row.names = FALSE)
save_tab(src, "15_green_sources")
src$Source <- factor(src$Source, assets)
p_src <- ggplot(src, aes(Source, Share, fill = Source)) + geom_col(width = 0.7) + geom_errorbar(aes(ymin = Lo, ymax = Hi), width = 0.25) +
  facet_grid(Basis ~ Receiver, scales = "free_y") + scale_fill_manual(values = pal, guide = "none") +
  labs(title = "RQ3: where do the green assets' forecast-error variances come from? (% share, 95% CI)", x = NULL, y = "Share of FEV (%)") +
  theme_dy(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
show(p_src, "15_green_sources", 9, 5.5)

# ---- step16_rolling ----
# STEP 16 - Rolling-window connectedness (250 trading days, H = 10): returns VAR(p_ret), GK-volatility VAR(p_vol)
rr <- rolling_dy(R, rdates, p_ret, H, W)
rv <- rolling_dy(VOL, rdates, p_vol, H, W)
cat("Rolling windows:", length(rr$TCI), "| first window ends", format(min(rr$dates)), "| last", format(max(rr$dates)), "\n")
cat("Rolling TCI (returns)    : mean", round(mean(rr$TCI), 2), " min", round(min(rr$TCI), 2), " max", round(max(rr$TCI), 2),
    " (peak on", format(rr$dates[which.max(rr$TCI)]), ")\n")
cat("Rolling TCI (volatility) : mean", round(mean(rv$TCI), 2), " min", round(min(rv$TCI), 2), " max", round(max(rv$TCI), 2),
    " (peak on", format(rv$dates[which.max(rv$TCI)]), ")\n")

# Finite-sample bias: a K=7 VAR fitted on only 250 observations produces spurious "spillovers" even for independent series.
# Benchmark = average TCI over random 250-day windows whose columns are independently circular-shifted (no cross-market links).
null_window_tci <- function(Y, p, w = W, B = 300) mean(replicate(B, {
  s <- sample.int(nrow(Y) - w + 1, 1); Yw <- Y[s:(s + w - 1), ]
  static_dy(sapply(1:ncol(Yw), function(j) Yw[((seq_len(w) + sample(10:(w - 10), 1) - 1) %% w) + 1, j]), p, H)$TCI }))
nb_r <- null_window_tci(R, p_ret); nb_v <- null_window_tci(VOL, p_vol)
cat("No-spillover benchmark TCI for a 250-day window: returns", round(nb_r, 2), "| volatility", round(nb_v, 2), "\n")
cat("Bias-adjusted mean rolling TCI: returns", round(mean(rr$TCI) - nb_r, 2), "| volatility", round(mean(rv$TCI) - nb_v, 2), "\n")
tci_df <- rbind(data.frame(Date = rr$dates, TCI = rr$TCI, Basis = "Returns"),
                data.frame(Date = rv$dates, TCI = rv$TCI, Basis = "GK volatility"))
p_tci <- ggplot(tci_df, aes(Date, TCI, colour = Basis)) + geom_line(linewidth = 0.55) + add_events() +
  geom_hline(data = data.frame(Basis = c("Returns", "GK volatility"), y = c(nb_r, nb_v)), aes(yintercept = y, colour = Basis), linetype = "dotted") +
  scale_colour_manual(values = c("Returns" = "#2166AC", "GK volatility" = "#B35806")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Rolling total connectedness (250-day window, H = 10)", subtitle = "Dotted lines: no-spillover benchmark (finite-sample bias of a 6-variable VAR on 250 observations)", x = NULL, y = "TCI (%)", colour = NULL, caption = ev_caption()) + theme_dy(10)
show(p_tci, "16_rolling_tci", 10, 4.8)

net_long <- function(Rl, nm) { x <- as.data.frame(Rl$NET); x$Date <- Rl$dates
  y <- reshape(x, direction = "long", varying = assets, v.names = "NET", timevar = "Asset", times = assets, idvar = "Date")
  y$Basis <- nm; y$Asset <- factor(y$Asset, assets); rownames(y) <- NULL; y }
nl <- rbind(net_long(rr, "Returns"), net_long(rv, "GK volatility"))
p_nets <- lapply(c("GK volatility", "Returns"), function(b) {
  ggplot(nl[nl$Basis == b, ], aes(Date, NET, fill = Asset)) + geom_area(alpha = 0.55) + geom_hline(yintercept = 0, linewidth = 0.3) +
    add_events() + facet_wrap(~Asset, ncol = 2, scales = "free_y") + scale_fill_manual(values = pal, guide = "none") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title = paste0("Rolling net connectedness - ", b, " (positive = net transmitter)"), x = NULL, y = "NET (pp)", caption = ev_caption()) + theme_dy(9) })
show(p_nets[[1]], "16_rolling_net_vol", 10, 8.5)
show(p_nets[[2]], "16_rolling_net_ret", 10, 8.5)

tf_plot <- function(Rl, nm) {
  d <- rbind(
    transform(reshape(data.frame(Rl$TO, Date = Rl$dates), direction = "long", varying = assets, v.names = "Value", timevar = "Asset", times = assets, idvar = "Date"), Dir = "TO others"),
    transform(reshape(data.frame(Rl$FROM, Date = Rl$dates), direction = "long", varying = assets, v.names = "Value", timevar = "Asset", times = assets, idvar = "Date"), Dir = "FROM others"))
  d$Asset <- factor(d$Asset, assets)
  ggplot(d, aes(Date, Value, colour = Dir)) + geom_line(linewidth = 0.4) + add_events() + facet_wrap(~Asset, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = c("TO others" = "#B2182B", "FROM others" = "#2166AC")) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title = paste0("Rolling directional connectedness - ", nm), x = NULL, y = "%", colour = NULL, caption = ev_caption()) + theme_dy(9)
}
show(tf_plot(rv, "GK volatility"), "16_rolling_tofrom_vol", 10, 8.5)
show(tf_plot(rr, "Returns"), "16_rolling_tofrom_ret", 10, 8.5)

# ---- step17_pairs ----
# STEP 17 - Key bilateral channels over time (net pairwise connectedness; positive = first asset transmits to second)
pairs <- list(c("GREEN_EQ", "SP500"), c("GREEN_BOND", "SP500"), c("GREEN_EQ", "GREEN_BOND"),
              c("CRUDE", "GREEN_EQ"), c("BITCOIN", "SP500"), c("GOLD", "GREEN_BOND"))
pair_df <- function(Rl, nm) do.call(rbind, lapply(pairs, function(pp)
  data.frame(Date = Rl$dates, Net = pair_series(Rl, pp[1], pp[2]), Pair = paste(pp[1], "->", pp[2]), Basis = nm)))
pd <- rbind(pair_df(rv, "GK volatility"), pair_df(rr, "Returns"))
pd$Pair <- factor(pd$Pair, unique(pd$Pair))
for (b in c("GK volatility", "Returns")) {
  p <- ggplot(pd[pd$Basis == b, ], aes(Date, Net)) + geom_area(fill = "grey70", alpha = 0.7) + geom_line(linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.3) + add_events() + facet_wrap(~Pair, ncol = 2, scales = "free_y") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title = paste0("Net pairwise connectedness over time - ", b), subtitle = "Positive = the first asset is the net transmitter to the second",
         x = NULL, y = "Net pairwise (pp)", caption = ev_caption()) + theme_dy(9)
  show(p, paste0("17_pairs_", ifelse(b == "Returns", "ret", "vol")), 10, 7.5)
}

# ---- step18_event_windows ----
# STEP 18 - Event-window analysis: 60 trading days before vs 60 after each CORE event (rolling TCI and NET)
evw <- function(vec, dts, d0, w = 60) { pos <- which(dts >= d0)[1]
  if (is.na(pos) || pos - w < 1 || pos + w - 1 > length(vec)) return(c(Pre = NA, Post = NA))
  c(Pre = mean(vec[(pos - w):(pos - 1)]), Post = mean(vec[pos:(pos + w - 1)])) }
ev_tci <- do.call(rbind, lapply(seq_len(nrow(core_ev)), function(i) {
  a <- evw(rr$TCI, rr$dates, core_ev$Date[i]); b <- evw(rv$TCI, rv$dates, core_ev$Date[i])
  data.frame(Event = core_ev$Event[i], Date = core_ev$Date[i], Ret_pre = a[1], Ret_post = a[2], Ret_change = a[2] - a[1],
             Vol_pre = b[1], Vol_post = b[2], Vol_change = b[2] - b[1]) }))
rownames(ev_tci) <- NULL; print(ev_tci, digits = 3, row.names = FALSE); save_tab(ev_tci, "18_event_tci")

ev_net <- do.call(rbind, lapply(seq_len(nrow(core_ev)), function(i) do.call(rbind, lapply(c("GK volatility", "Returns"), function(b) {
  Rl <- if (b == "Returns") rr else rv
  do.call(rbind, lapply(seq_along(assets), function(j) { a <- evw(Rl$NET[, j], Rl$dates, core_ev$Date[i])
    data.frame(Event = core_ev$Event[i], Basis = b, Asset = assets[j], Pre = a[1], Post = a[2], Change = a[2] - a[1]) })) }))))
rownames(ev_net) <- NULL
cat("\nNET change (post - pre), GK volatility:\n")
print(round(reshape(ev_net[ev_net$Basis == "GK volatility", c("Event", "Asset", "Change")], idvar = "Event", timevar = "Asset", direction = "wide")[, -1], 2), row.names = core_ev$Event)
cat("\nNET change (post - pre), Returns:\n")
print(round(reshape(ev_net[ev_net$Basis == "Returns", c("Event", "Asset", "Change")], idvar = "Event", timevar = "Asset", direction = "wide")[, -1], 2), row.names = core_ev$Event)
save_tab(ev_net, "18_event_net")

# Formal test: static TCI on the 200 days before vs the 200 days after each core event (block bootstrap, B = 500).
# VAR(1) is used in both systems because a short window cannot support more parameters without spurious spillovers.
ev_static <- function(Y, dts, d0, p, w = 200, B = 500, bl = 10) {
  pos <- which(dts >= d0)[1]; pre <- Y[(pos - w):(pos - 1), ]; post <- Y[pos:(pos + w - 1), ]
  t0 <- static_dy(pre, p, H)$TCI; t1 <- static_dy(post, p, H)$TCI
  dd <- replicate(B, static_dy(post[mbb_idx(w, bl), ], p, H)$TCI - static_dy(pre[mbb_idx(w, bl), ], p, H)$TCI)
  c(Pre = t0, Post = t1, Diff = t1 - t0, Lo = unname(quantile(dd, .025)), Hi = unname(quantile(dd, .975)))
}
ev_st <- do.call(rbind, lapply(seq_len(nrow(core_ev)), function(i) do.call(rbind, lapply(c("Returns", "GK volatility"), function(b) {
  Y <- if (b == "Returns") R else VOL
  r <- ev_static(Y, rdates, core_ev$Date[i], 1); data.frame(Event = core_ev$Event[i], Basis = b, t(round(r, 2))) }))))
ev_st$Significant <- ifelse(ev_st$Lo > 0 | ev_st$Hi < 0, "Yes", "No"); rownames(ev_st) <- NULL
cat("\nStatic VAR(1) TCI, 200 days before vs 200 days after (95% bootstrap CI of the difference):\n"); print(ev_st, row.names = FALSE)
save_tab(ev_st, "18_event_static_test")

# ---- step19_tvpvar ----
# STEP 19 - TVP-VAR connectedness (Antonakakis, Chatziantoniou & Gabauer 2020): Kalman filter with forgetting factors
tvp_dy <- function(Y, k1 = 0.99, k2 = 0.99, train = 60, Hh = H) {
  Y <- sweep(as.matrix(Y), 2, colMeans(Y)); n <- nrow(Y); Kv <- ncol(Y); nb <- Kv^2
  Z <- Y[1:(train - 1), ]; Y1 <- Y[2:train, ]
  B0 <- t(solve(crossprod(Z) + 0.1 * train * diag(Kv), crossprod(Z, Y1)))     # ridge-OLS prior from the training sample
  beta <- as.vector(t(B0)); P <- 0.01 * diag(nb); Sg <- cov(Y[1:train, ]); TH <- array(NA, c(Kv, Kv, n))
  for (t in 2:n) {
    X <- kronecker(diag(Kv), matrix(Y[t - 1, ], 1)); Rt <- P / k1
    e <- Y[t, ] - drop(X %*% beta); Sg <- k2 * Sg + (1 - k2) * tcrossprod(e)
    G <- Rt %*% t(X) %*% solve(Sg + X %*% Rt %*% t(X)); beta <- beta + drop(G %*% e); P <- (diag(nb) - G %*% X) %*% Rt
    TH[, , t] <- gfevd(list(matrix(beta, Kv, Kv, byrow = TRUE)), Sg, Hh)
  }
  keep <- (train + 1):n; list(S = summ_array(TH[, , keep, drop = FALSE], rdates[keep], colnames(Y)))
}
tv_r <- tvp_dy(R)$S; tv_v <- tvp_dy(VOL)$S
cat("TVP-VAR TCI (returns)    : mean", round(mean(tv_r$TCI), 2), " range [", round(min(tv_r$TCI), 2), ",", round(max(tv_r$TCI), 2), "]\n")
cat("TVP-VAR TCI (volatility) : mean", round(mean(tv_v$TCI), 2), " range [", round(min(tv_v$TCI), 2), ",", round(max(tv_v$TCI), 2), "]\n")
tvdf <- rbind(data.frame(Date = tv_r$dates, TCI = tv_r$TCI, Basis = "Returns"), data.frame(Date = tv_v$dates, TCI = tv_v$TCI, Basis = "GK volatility"))
p_tv <- ggplot(tvdf, aes(Date, TCI, colour = Basis)) + geom_line(linewidth = 0.5) + add_events() +
  scale_colour_manual(values = c("Returns" = "#2166AC", "GK volatility" = "#B35806")) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "TVP-VAR total connectedness (kappa1 = kappa2 = 0.99)", x = NULL, y = "TCI (%)", colour = NULL, caption = ev_caption()) + theme_dy(10)
show(p_tv, "19_tvpvar_tci", 10, 4.8)

tvn <- rbind(transform(net_long(tv_v, "GK volatility")), transform(net_long(tv_r, "Returns")))
p_tvn <- ggplot(tvn[tvn$Basis == "GK volatility", ], aes(Date, NET, fill = Asset)) + geom_area(alpha = 0.55) + geom_hline(yintercept = 0, linewidth = 0.3) +
  add_events() + facet_wrap(~Asset, ncol = 2, scales = "free_y") + scale_fill_manual(values = pal, guide = "none") + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "TVP-VAR net volatility connectedness (positive = net transmitter)", x = NULL, y = "NET (pp)", caption = ev_caption()) + theme_dy(9)
show(p_tvn, "19_tvpvar_net_vol", 10, 8.5)

# robustness to the forgetting factors
kgrid <- list(c(0.96, 0.96), c(0.98, 0.98), c(0.99, 0.99), c(0.99, 0.96))
tv_rob <- lapply(kgrid, function(k) list(r = tvp_dy(R, k[1], k[2])$S, v = tvp_dy(VOL, k[1], k[2])$S))
rob_tab <- do.call(rbind, lapply(seq_along(kgrid), function(i) data.frame(
  kappa1 = kgrid[[i]][1], kappa2 = kgrid[[i]][2],
  Mean_TCI_ret = mean(tv_rob[[i]]$r$TCI), Mean_TCI_vol = mean(tv_rob[[i]]$v$TCI),
  Corr_with_base_ret = cor(tv_rob[[i]]$r$TCI, tv_r$TCI), Corr_with_base_vol = cor(tv_rob[[i]]$v$TCI, tv_v$TCI),
  GREEN_EQ_NET_vol = mean(tv_rob[[i]]$v$NET[, "GREEN_EQ"]), GREEN_BOND_NET_vol = mean(tv_rob[[i]]$v$NET[, "GREEN_BOND"]))))
print(round(rob_tab, 3), row.names = FALSE); save_tab(rob_tab, "19_tvpvar_kappa_robust")
cat("Correlation rolling-window TCI vs TVP-VAR TCI (common dates): returns",
    round(cor(rr$TCI[match(tv_r$dates, rr$dates)], tv_r$TCI, use = "complete.obs"), 3), "| volatility",
    round(cor(rv$TCI[match(tv_v$dates, rv$dates)], tv_v$TCI, use = "complete.obs"), 3), "\n")

# ---- step20_drivers ----
# STEP 20 - What moves connectedness? HAC regressions of the TVP-VAR volatility TCI on market stress and event dummies
Zv <- rowMeans(scale(VOL)); names(Zv) <- as.character(rdates)
stress5 <- stats::filter(Zv, rep(1 / 5, 5), sides = 1)
dum <- function(d0, d1) as.numeric(tv_v$dates >= as.Date(d0) & tv_v$dates <= as.Date(d1))
dd <- data.frame(TCI = tv_v$TCI, Stress = as.numeric(stress5)[match(tv_v$dates, rdates)],
                 COVID = dum("2020-03-11", "2020-06-30"), Ukraine = dum("2022-02-24", "2022-06-30"),
                 Crypto_shock = dum("2022-05-01", "2022-11-30"), Tariff = dum("2025-04-02", "2025-06-30"))
dd <- dd[complete.cases(dd), ]
mod_l <- lm(TCI ~ Stress + COVID + Ukraine + Crypto_shock + Tariff, dd)
mod_d <- lm(diff(dd$TCI) ~ diff(dd$Stress))
ct_l <- coeftest_nw <- function(m, lag) { V <- NeweyWest(m, lag = lag, prewhite = FALSE); b <- coef(m); se <- sqrt(diag(V))
  data.frame(Estimate = b, SE_NW = se, t = b / se, p = 2 * pnorm(-abs(b / se))) }
cat("Levels: TCI(vol, TVP-VAR) on market stress and event dummies; Newey-West (60 lags)\n"); print(round(ct_l(mod_l, 60), 4))
cat("R-squared:", round(summary(mod_l)$r.squared, 3), "\n")
cat("\nFirst differences: change in TCI on change in market stress; Newey-West (10 lags)\n"); print(round(ct_l(mod_d, 10), 4))
save_tab(ct_l(mod_l, 60), "20_drivers_levels")

# ---- step21_quantile ----
# STEP 21 - Connectedness across market states: quantile VAR (Ando, Greenwood-Nimark & Sarisoy 2022)
# Each VAR equation is estimated by quantile regression at tau; GFEVD is built from the tau-quantile coefficient matrix and the
# covariance of the tau-quantile residuals.  Returns: tau = 0.05 bearish, 0.50 normal, 0.95 bullish.
# Volatility: tau = 0.05 calm, 0.50 normal, 0.95 turbulent.
qvar_dy <- function(Y, tau, Hh = H) {
  n <- nrow(Y); Kv <- ncol(Y); Zc <- cbind(1, Y[-n, , drop = FALSE]); Yt <- Y[-1, , drop = FALSE]
  A <- matrix(NA, Kv, Kv); U <- matrix(NA, n - 1, Kv)
  for (i in 1:Kv) { co <- coef(rq.fit.fnb(Zc, Yt[, i], tau = tau)); A[i, ] <- co[-1]; U[, i] <- Yt[, i] - Zc %*% co }
  dy_summary(gfevd(list(A), cov(U), Hh), colnames(Y))
}
taus <- c(0.05, 0.50, 0.95)
qr_r <- lapply(taus, function(t) qvar_dy(R, t)); qr_v <- lapply(taus, function(t) qvar_dy(VOL, t))
B_Q <- 300; nq <- nrow(R)
qb <- lapply(list(R = R, V = VOL), function(Y) { out <- matrix(NA, B_Q, 3); net <- array(NA, c(B_Q, K, 3))
  for (b in 1:B_Q) { ix <- mbb_idx(nq); for (k in 1:3) { s <- qvar_dy(Y[ix, ], taus[k]); out[b, k] <- s$TCI; net[b, , k] <- s$NET } }
  list(TCI = out, NET = net) })
qtab <- do.call(rbind, lapply(c("Returns", "GK volatility"), function(b) {
  res <- if (b == "Returns") qr_r else qr_v; bb <- if (b == "Returns") qb$R else qb$V
  data.frame(Basis = b, Tau = taus, State = if (b == "Returns") c("Bearish", "Normal", "Bullish") else c("Calm", "Normal", "Turbulent"),
             TCI = sapply(res, function(s) s$TCI), Lo = apply(bb$TCI, 2, quantile, .025), Hi = apply(bb$TCI, 2, quantile, .975)) }))
print(qtab, digits = 3, row.names = FALSE); save_tab(qtab, "21_quantile_tci")
cat("\nTCI(tail) - TCI(median), paired block-bootstrap 95% CI:\n")
qd <- do.call(rbind, lapply(list(list("Returns", qr_r, qb$R), list("GK volatility", qr_v, qb$V)), function(z) do.call(rbind, lapply(c(1, 3), function(k) {
  d <- z[[3]]$TCI[, k] - z[[3]]$TCI[, 2]
  data.frame(Basis = z[[1]], Contrast = paste0("tau=", taus[k], " minus tau=0.50"), Diff = z[[2]][[k]]$TCI - z[[2]][[2]]$TCI,
             Lo = unname(quantile(d, .025)), Hi = unname(quantile(d, .975))) }))))
print(qd, digits = 3, row.names = FALSE); save_tab(qd, "21_quantile_diff")
qnet <- do.call(rbind, lapply(c("Returns", "GK volatility"), function(b) { res <- if (b == "Returns") qr_r else qr_v
  do.call(rbind, lapply(1:3, function(k) data.frame(Basis = b, Tau = taus[k], Asset = assets, NET = res[[k]]$NET))) }))
cat("\nNET connectedness by quantile:\n")
for (b in c("Returns", "GK volatility")) { x <- qnet[qnet$Basis == b, ]; cat("\n", b, "\n"); print(round(reshape(x[, c("Tau", "Asset", "NET")], idvar = "Asset", timevar = "Tau", direction = "wide")[, -1], 2), row.names = assets) }
save_tab(qnet, "21_quantile_net")
qtab$Label <- paste0(qtab$State, " (tau=", qtab$Tau, ")")
p_q1 <- ggplot(qtab, aes(factor(Tau), TCI, fill = State)) + geom_col(width = 0.6) + geom_errorbar(aes(ymin = Lo, ymax = Hi), width = 0.2) +
  facet_wrap(~Basis) + geom_text(aes(label = State, y = 1), vjust = 0, size = 3, colour = "white") + scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "Total connectedness across quantiles (95% block-bootstrap interval)", x = "Quantile (tau)", y = "TCI (%)") + theme_dy(10)
qnet$Asset <- factor(qnet$Asset, assets)
p_q2 <- ggplot(qnet, aes(Asset, NET, fill = factor(Tau))) + geom_col(position = "dodge", width = 0.75) + geom_hline(yintercept = 0) + facet_wrap(~Basis, scales = "free_y") +
  scale_fill_manual(values = c("0.05" = "#2166AC", "0.5" = "grey60", "0.95" = "#B2182B")) + labs(title = "Net connectedness by quantile", x = NULL, y = "NET (pp)", fill = "tau") +
  theme_dy(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
show(p_q1 / p_q2, "21_quantile", 10, 8)

# ---- step22_frequency ----
# STEP 22 - Connectedness across investment horizons: frequency decomposition (Barunik & Krehlik 2018)
freq_dy <- function(f, nf = 500, cuts = c(20, 5)) {          # bands: long (>20d), medium (5-20d), short (1-5d)
  Kv <- f$K; Sg <- f$Sigma; om <- pi * (1:nf) / nf; sdg <- diag(Sg)
  Th <- array(NA, c(Kv, Kv, nf)); Gm <- matrix(NA, nf, Kv)
  for (k in 1:nf) {
    Pm <- diag(Kv) + 0i; for (l in seq_along(f$A)) Pm <- Pm - f$A[[l]] * exp(-1i * om[k] * l)
    Ps <- solve(Pm); PS <- Ps %*% Sg; Spec <- Re(diag(PS %*% Conj(t(Ps))))
    Th[, , k] <- sweep(Mod(PS)^2, 2, sdg, "/") / Spec; Gm[k, ] <- Spec }
  Gw <- sweep(Gm, 2, colMeans(Gm), "/")                       # frequency weights Gamma_i(omega)
  for (k in 1:nf) Th[, , k] <- Th[, , k] * Gw[k, ]
  inf <- apply(Th, c(1, 2), mean); rs <- rowSums(inf)          # Theta(infinity) and its row sums
  band <- ifelse(om > pi / cuts[2], "Short (1-5d)", ifelse(om > pi / cuts[1], "Medium (5-20d)", "Long (>20d)"))
  out <- lapply(c("Short (1-5d)", "Medium (5-20d)", "Long (>20d)"), function(bn) {
    Tb <- apply(Th[, , band == bn, drop = FALSE], c(1, 2), sum) / nf / rs; off <- Tb; diag(off) <- 0
    list(Band = bn, Abs = 100 * sum(off) / Kv, Within = 100 * sum(off) / sum(Tb), TO = 100 * colSums(off), FROM = 100 * rowSums(off), Tb = Tb) })
  names(out) <- sapply(out, `[[`, "Band"); out
}
fr_r <- freq_dy(fit_r); fr_v <- freq_dy(fit_v)
ftab <- do.call(rbind, lapply(list(list("Returns", fr_r), list("GK volatility", fr_v)), function(z)
  data.frame(Basis = z[[1]], Band = names(z[[2]]), Absolute_TCI = sapply(z[[2]], `[[`, "Abs"), Within_TCI = sapply(z[[2]], `[[`, "Within"),
             Share_of_total = 100 * sapply(z[[2]], `[[`, "Abs") / sum(sapply(z[[2]], `[[`, "Abs")))))
rownames(ftab) <- NULL
cat("Absolute TCI by band sums to the infinite-horizon TCI:  returns", round(sum(ftab$Absolute_TCI[ftab$Basis == "Returns"]), 2),
    "| volatility", round(sum(ftab$Absolute_TCI[ftab$Basis == "GK volatility"]), 2), "\n")
# bootstrap intervals for the band shares
fb <- lapply(list(R = list(R, p_ret), V = list(VOL, p_vol)), function(z) t(replicate(300, { ix <- mbb_idx(nrow(z[[1]]))
  fr <- freq_dy(var_fit(z[[1]][ix, ], z[[2]]), nf = 200); a <- sapply(fr, `[[`, "Abs"); 100 * a / sum(a) })))
ftab$Share_lo <- c(apply(fb$R, 2, quantile, .025), apply(fb$V, 2, quantile, .025)); ftab$Share_hi <- c(apply(fb$R, 2, quantile, .975), apply(fb$V, 2, quantile, .975))
print(ftab, digits = 3, row.names = FALSE); save_tab(ftab, "22_frequency_tci")
fnet <- do.call(rbind, lapply(list(list("Returns", fr_r), list("GK volatility", fr_v)), function(z)
  do.call(rbind, lapply(z[[2]], function(b) data.frame(Basis = z[[1]], Band = b$Band, Asset = assets, TO = b$TO, FROM = b$FROM, NET = b$TO - b$FROM)))))
rownames(fnet) <- NULL
for (b in c("Returns", "GK volatility")) { x <- fnet[fnet$Basis == b, ]; cat("\nNET by band -", b, "\n")
  print(round(reshape(x[, c("Band", "Asset", "NET")], idvar = "Asset", timevar = "Band", direction = "wide")[, -1], 2), row.names = assets) }
save_tab(fnet, "22_frequency_net")
ftab$Band <- factor(ftab$Band, c("Short (1-5d)", "Medium (5-20d)", "Long (>20d)")); fnet$Band <- factor(fnet$Band, levels(ftab$Band)); fnet$Asset <- factor(fnet$Asset, assets)
p_f1 <- ggplot(ftab, aes(Band, Absolute_TCI, fill = Band)) + geom_col(width = 0.6) + facet_wrap(~Basis) + scale_fill_manual(values = c("Short (1-5d)" = "#D95F02", "Medium (5-20d)" = "#7570B3", "Long (>20d)" = "#1B9E77"), guide = "none") +
  labs(title = "Total connectedness by horizon band (absolute: bands sum to the total)", x = NULL, y = "TCI contribution (pp)") + theme_dy(10)
p_f2 <- ggplot(fnet, aes(Asset, NET, fill = Band)) + geom_col(position = "dodge", width = 0.75) + geom_hline(yintercept = 0) + facet_wrap(~Basis, scales = "free_y") +
  scale_fill_manual(values = c("Short (1-5d)" = "#D95F02", "Medium (5-20d)" = "#7570B3", "Long (>20d)" = "#1B9E77")) + labs(title = "Net connectedness by horizon band", x = NULL, y = "NET (pp)", fill = NULL) +
  theme_dy(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
show(p_f1 / p_f2, "22_frequency", 10, 8)

# ---- step23_incremental ----
# STEP 23 - RQ6: do the green assets change system connectedness and the roles of the conventional assets?
conv <- c("SP500", "GOLD", "CRUDE", "BITCOIN"); gr <- c("GREEN_EQ", "GREEN_BOND"); ci5 <- match(conv, assets); gi <- match(gr, assets)
inc_stat <- function(Y7, p) {
  s7 <- static_dy(Y7, p, H); s5 <- static_dy(Y7[, conv], p, H); th <- s7$theta
  within5 <- 100 * (sum(th[ci5, ci5]) - sum(diag(th[ci5, ci5]))) / length(conv)           # conventional FEV explained by OTHER conventional assets (in the full system)
  from_green <- 100 * mean(rowSums(th[ci5, gi])); to_green_from_conv <- 100 * mean(rowSums(th[gi, ci5]))
  c(TCI_all = s7$TCI, TCI_conv = s5$TCI, dTCI = s7$TCI - s5$TCI, Conv_within = within5, Conv_from_green = from_green, Green_from_conv = to_green_from_conv,
    setNames(s7$NET[ci5], paste0("NETall_", conv)), setNames(s5$NET, paste0("NETconv_", conv)), setNames(s7$NET[ci5] - s5$NET, paste0("dNET_", conv)))
}
inc_obs <- list(Returns = inc_stat(R, p_ret), Vol = inc_stat(VOL, p_vol))
inc_bs <- lapply(list(Returns = list(R, p_ret), Vol = list(VOL, p_vol)), function(z) t(replicate(B_STAT, inc_stat(z[[1]][mbb_idx(nrow(z[[1]])), ], z[[2]]))))
inc_tab <- do.call(rbind, lapply(c("Returns", "Vol"), function(b) data.frame(Basis = ifelse(b == "Vol", "GK volatility", b), Measure = names(inc_obs[[b]]),
  Estimate = inc_obs[[b]], Lo95 = apply(inc_bs[[b]], 2, quantile, .025), Hi95 = apply(inc_bs[[b]], 2, quantile, .975), row.names = NULL)))
inc_tab$Excludes_zero <- ifelse(inc_tab$Lo95 > 0 | inc_tab$Hi95 < 0, "Yes", "No")
print(transform(inc_tab, Estimate = round(Estimate, 2), Lo95 = round(Lo95, 2), Hi95 = round(Hi95, 2)), row.names = FALSE)
save_tab(inc_tab, "23_incremental_value")

# Same comparison through time: rolling TCI of the 4-asset conventional system vs the 6-asset system
rr5 <- rolling_dy(R[, conv], rdates, p_ret, H, W); rv5 <- rolling_dy(VOL[, conv], rdates, p_vol, H, W)
d5 <- rbind(data.frame(Date = rr$dates, TCI = rr$TCI, System = "6 assets (with green)", Basis = "Returns"),
            data.frame(Date = rr5$dates, TCI = rr5$TCI, System = "4 conventional only", Basis = "Returns"),
            data.frame(Date = rv$dates, TCI = rv$TCI, System = "6 assets (with green)", Basis = "GK volatility"),
            data.frame(Date = rv5$dates, TCI = rv5$TCI, System = "4 conventional only", Basis = "GK volatility"))
cat("\nRolling TCI, mean: returns 6-asset", round(mean(rr$TCI), 2), "vs 4-asset", round(mean(rr5$TCI), 2),
    "| volatility 6-asset", round(mean(rv$TCI), 2), "vs 4-asset", round(mean(rv5$TCI), 2), "\n")
p_d5 <- ggplot(d5, aes(Date, TCI, colour = System)) + geom_line(linewidth = 0.5) + add_events() + facet_wrap(~Basis, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("6 assets (with green)" = "#00A651", "4 conventional only" = "#444444")) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "RQ6: rolling total connectedness with and without the green assets", x = NULL, y = "TCI (%)", colour = NULL, caption = ev_caption()) + theme_dy(10)
show(p_d5, "23_rolling_5_vs_7", 10, 6.5)

# ---- step24_hedging ----
# STEP 24 - RQ7: do the green assets help reduce SP500 risk?  (hedging and diversification)
# (a) Minimum-variance hedge ratio h = Cov(SP500, X)/Var(X): hedged return = r_SP500 - h * r_X.
#     Hedging effectiveness HE = 1 - Var(hedged)/Var(SP500).  Four precisely defined versions:
#       HE_in      : same-day, one hedge ratio estimated on the whole sample (in-sample)
#       HE_oos     : same-day, ratio from the previous 250 days, evaluated on the later days (out-of-sample ratio)
#       HE_nextday : SP500(t+1) hedged with X(t), one ratio on the whole sample (lagged hedge, in-sample ratio)
#       HE_next_oos: SP500(t+1) hedged with X(t), ratio from the previous 250 pairs (lagged hedge AND out-of-sample ratio)
# (b) Two-asset minimum-variance portfolio of SP500 and X (Kroner & Ng 1998) with rolling 250-day covariances.
# Timing: the S&P 500 and the US-listed ETFs (ICLN, GRNB) close together, so a same-day hedge is tradeable for them
# (gold and WTI futures settle at slightly different times). The "next-day" rows are a lagged-hedge check:
# S&P 500(t+1) hedged with X(t).
rN <- R[, "SP500"]; others <- setdiff(assets, "SP500")
HE <- function(y, x, h) 1 - var(y - h * x) / var(y)
rollh <- function(y, x, w = W) { n <- length(y); h <- rep(NA, n)
  for (t in (w + 1):n) { ii <- (t - w):(t - 1); h[t] <- cov(y[ii], x[ii]) / var(x[ii]) }; h }
hedge_tab <- do.call(rbind, lapply(others, function(a) {
  x <- R[, a]; h_in <- cov(rN, x) / var(x); hr <- rollh(rN, x); ok <- !is.na(hr)
  he_oos <- 1 - var((rN - hr * x)[ok]) / var(rN[ok])
  hn <- cov(rN[-1], x[-length(x)]) / var(x[-length(x)])                       # next-day (implementable)
  he_next <- HE(rN[-1], x[-length(x)], hn)
  y1 <- rN[-1]; x1 <- x[-length(x)]; hr1 <- rollh(y1, x1); ok1 <- !is.na(hr1)
  he_next_oos <- 1 - var((y1 - hr1 * x1)[ok1]) / var(y1[ok1])
  covid <- rdates >= as.Date("2020-02-20") & rdates <= as.Date("2020-06-30")
  data.frame(Hedge_asset = a, h_in = h_in, HE_in = HE(rN, x, h_in), HE_oos = he_oos, HE_nextday = he_next, HE_next_oos = he_next_oos,
             HE_covid = HE(rN[covid], x[covid], cov(rN[covid], x[covid]) / var(x[covid])), HE_calm = HE(rN[!covid], x[!covid], cov(rN[!covid], x[!covid]) / var(x[!covid])),
             Corr = cor(rN, x)) }))
rownames(hedge_tab) <- NULL
hedge_tab[, -1] <- lapply(hedge_tab[, -1], function(v) round(v, 4)); print(hedge_tab, row.names = FALSE)

# block-bootstrap interval for in-sample HE of each green asset
he_ci <- do.call(rbind, lapply(c("GREEN_EQ", "GREEN_BOND", "GOLD", "CRUDE"), function(a) {
  b <- replicate(B_STAT, { ix <- mbb_idx(Tr); y <- rN[ix]; x <- R[ix, a]; HE(y, x, cov(y, x) / var(x)) })
  data.frame(Hedge_asset = a, HE_in = HE(rN, R[, a], cov(rN, R[, a]) / var(R[, a])), Lo95 = unname(quantile(b, .025)), Hi95 = unname(quantile(b, .975))) }))
cat("\nIn-sample hedging effectiveness with 95% block-bootstrap interval:\n"); print(round(he_ci[, -1], 4), row.names = he_ci$Hedge_asset)
save_tab(hedge_tab, "24_hedge_ratios"); save_tab(he_ci, "24_hedge_ci")

# (b) Diversification: hold (1 - a) SP500 + a * X for a = 10% and 20% (fixed weights, no optimiser). A two-asset minimum-variance
#     optimiser is NOT used because it simply parks most of the money in the lowest-volatility asset (e.g. a bond ETF).
cvar5 <- function(x) mean(x[x <= quantile(x, .05)])
mvp <- do.call(rbind, lapply(others, function(a) do.call(rbind, lapply(c(0.10, 0.20), function(w) {
  port <- (1 - w) * rN + w * R[, a]
  bs <- replicate(B_STAT, { ix <- mbb_idx(Tr); 100 * (1 - var(((1 - w) * rN + w * R[, a])[ix]) / var(rN[ix])) })
  data.frame(Pair = paste0("SP500 + ", a), Allocation_to_X = w, Var_reduction_pct = 100 * (1 - var(port) / var(rN)), Lo95 = unname(quantile(bs, .025)), Hi95 = unname(quantile(bs, .975)),
             CVaR5_SP500 = cvar5(rN), CVaR5_portfolio = cvar5(port)) }))))
rownames(mvp) <- NULL
cat("\nDiversification with a fixed allocation to X (variance reduction vs holding SP500 alone, 95% block-bootstrap interval):\n")
print(transform(mvp, Var_reduction_pct = round(Var_reduction_pct, 2), Lo95 = round(Lo95, 2), Hi95 = round(Hi95, 2), CVaR5_SP500 = round(CVaR5_SP500, 3), CVaR5_portfolio = round(CVaR5_portfolio, 3)), row.names = FALSE)
save_tab(mvp, "24_fixed_allocation")

# dynamic hedge ratios (green assets) with events
hr_df <- do.call(rbind, lapply(c("GREEN_EQ", "GREEN_BOND", "GOLD"), function(a) { h <- rollh(rN, R[, a]); data.frame(Date = rdates, h = h, Asset = a) }))
p_h <- ggplot(na.omit(hr_df), aes(Date, h, colour = Asset)) + geom_line(linewidth = 0.5) + geom_hline(yintercept = 0, linewidth = 0.3) + add_events() +
  scale_colour_manual(values = pal) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Dynamic minimum-variance hedge ratio for SP500 (rolling 250 days)", subtitle = "Units of the hedge asset to short per unit of SP500", x = NULL, y = "Hedge ratio", colour = NULL, caption = ev_caption()) + theme_dy(10)
show(p_h, "24_dynamic_hedge_ratio", 10, 4.8)

# Minimum Connectedness Portfolio (Broadstock et al. 2022): w = C^-1 1 / (1' C^-1 1), C = pairwise connectedness index (PCI) matrix of the rolling VAR.
# (The row-normalised GFEVD matrix cannot be used directly: its rows sum to one, so C^-1 1 = 1 and every weight is 1/K.)
# Weights from window ending t are applied to day t+1 (long-only truncation). Benchmarks: equal weight and minimum variance.
mcp_eval <- function(Rl, Y) {
  idx <- match(Rl$dates, rdates); n <- length(idx) - 1; Wm <- Wv <- matrix(NA, n, K); pr <- matrix(NA, n, 3, dimnames = list(NULL, c("MinConnectedness", "MinVariance", "EqualWeight")))
  for (s in 1:n) {
    th <- Rl$theta[, , s]; C <- matrix(1, K, K)                       # PCI_ij = 2(theta_ij + theta_ji)/(theta_ii + theta_ij + theta_ji + theta_jj), diag = 1
    for (i in 1:K) for (j in 1:K) if (i != j) C[i, j] <- 2 * (th[i, j] + th[j, i]) / (th[i, i] + th[i, j] + th[j, i] + th[j, j])
    w <- solve(C, rep(1, K)); w <- pmax(w, 0); w <- w / sum(w)
    S <- cov(Y[(idx[s] - W + 1):idx[s], ]); v <- solve(S, rep(1, K)); v <- pmax(v, 0); v <- v / sum(v)
    ynext <- Y[idx[s] + 1, ]; pr[s, ] <- c(sum(w * ynext), sum(v * ynext), mean(ynext)); Wm[s, ] <- w }
  list(pr = pr, W = Wm, dates = Rl$dates[-1][seq_len(n)])
}
mc <- mcp_eval(rr, R)
perf <- function(x) c(AnnMean = mean(x) * 252, AnnSD = sd(x) * sqrt(252), Sharpe = mean(x) / sd(x) * sqrt(252), VaR5 = unname(quantile(x, .05)), CVaR5 = mean(x[x <= quantile(x, .05)]))
perf_tab <- t(apply(mc$pr, 2, perf)); cat("\nPortfolio performance (returns in %, out of sample, long-only):\n"); print(round(perf_tab, 3))
sd_diff <- replicate(B_STAT, { ix <- mbb_idx(nrow(mc$pr)); sd(mc$pr[ix, "MinConnectedness"]) - sd(mc$pr[ix, "EqualWeight"]) })
cat("SD(MCP) - SD(EqualWeight) per day:", round(sd(mc$pr[, 1]) - sd(mc$pr[, 3]), 4), " 95% CI [", round(quantile(sd_diff, .025), 4), ",", round(quantile(sd_diff, .975), 4), "]\n")
mw <- colMeans(mc$W); names(mw) <- assets; cat("Average MCP weights:\n"); print(round(mw, 3))
save_tab(perf_tab, "24_portfolio_perf"); save_tab(data.frame(Asset = assets, Weight = mw), "24_mcp_weights")
wl <- data.frame(Date = rep(mc$dates, K), Asset = factor(rep(assets, each = nrow(mc$W)), assets), w = as.vector(mc$W))
p_w <- ggplot(wl, aes(Date, w, fill = Asset)) + geom_area(alpha = 0.85) + add_events() + scale_fill_manual(values = pal) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Minimum-connectedness portfolio weights over time (returns system)", x = NULL, y = "Weight", fill = NULL, caption = ev_caption()) + theme_dy(10)
show(p_w, "24_mcp_weights", 10, 4.8)

# ---- step25_robustness ----
# STEP 25 - Robustness grid: do the headline conclusions survive other specifications?
#   TCI level and the NET role of GREEN_EQ / GREEN_BOND under: horizon, lag, volatility proxy, transformation, sample
LGK <- log(pmax(VOL, 1e-3))                                           # log-volatility reduces the weight of extreme days
exc <- !(rdates >= as.Date("2020-02-15") & rdates <= as.Date("2020-06-30"))
exc2 <- !(rdates >= as.Date("2022-05-01") & rdates <= as.Date("2022-11-30"))   # 2022 crypto crash (Terra/Luna, FTX)
specs <- list(
  "Baseline"                       = list(r = list(R, p_ret, H),   v = list(VOL, p_vol, H)),
  "H = 5 days"                     = list(r = list(R, p_ret, 5),   v = list(VOL, p_vol, 5)),
  "H = 20 days"                    = list(r = list(R, p_ret, 20),  v = list(VOL, p_vol, 20)),
  "VAR lag = 1"                    = list(r = list(R, 1, H),       v = list(VOL, 1, H)),
  "VAR lag = 3"                    = list(r = list(R, 3, H),       v = list(VOL, 3, H)),
  "VAR lag = 5"                    = list(r = list(R, 5, H),       v = list(VOL, 5, H)),
  "Parkinson volatility"           = list(r = NULL,                v = list(VPK, p_vol, H)),
  "Absolute returns"               = list(r = NULL,                v = list(VABS, p_vol, H)),
  "Squared returns"                = list(r = NULL,                v = list(VAR2, p_vol, H)),
  "Log GK volatility"              = list(r = NULL,                v = list(LGK, p_vol, H)),
  "Excl. COVID (15 Feb-30 Jun 20)" = list(r = list(R[exc, ], p_ret, H), v = list(VOL[exc, ], p_vol, H)),
  "Excl. 2022 crypto crash"        = list(r = list(R[exc2, ], p_ret, H), v = list(VOL[exc2, ], p_vol, H)))
rob <- do.call(rbind, lapply(names(specs), function(nm) do.call(rbind, lapply(c("r", "v"), function(b) {
  z <- specs[[nm]][[b]]; if (is.null(z)) return(NULL); s <- static_dy(z[[1]], z[[2]], z[[3]])
  data.frame(Spec = nm, Basis = ifelse(b == "r", "Returns", "Volatility"), TCI = s$TCI, NET_GREEN_EQ = s$NET["GREEN_EQ"], NET_GREEN_BOND = s$NET["GREEN_BOND"],
             Top_transmitter = names(which.max(s$NET)), Top_receiver = names(which.min(s$NET)), row.names = NULL) }))))
rownames(rob) <- NULL; print(transform(rob, TCI = round(TCI, 2), NET_GREEN_EQ = round(NET_GREEN_EQ, 2), NET_GREEN_BOND = round(NET_GREEN_BOND, 2)), row.names = FALSE)
save_tab(rob, "25_robustness_grid")
cat("\nGREEN_EQ is the top net transmitter or net-positive in", sum(rob$NET_GREEN_EQ > 0), "of", nrow(rob), "specifications;",
    "GREEN_EQ NET>0 in returns specs:", sum(rob$NET_GREEN_EQ[rob$Basis == "Returns"] > 0), "of", sum(rob$Basis == "Returns"),
    "| volatility specs:", sum(rob$NET_GREEN_EQ[rob$Basis == "Volatility"] > 0), "of", sum(rob$Basis == "Volatility"), "\n")
cat("GREEN_BOND NET sign: positive in", sum(rob$NET_GREEN_BOND > 0), "of", nrow(rob), "specs\n")

# Rolling-window robustness: 200 / 250 / 300 days
rw <- list(`200` = list(r = rolling_dy(R, rdates, p_ret, H, 200), v = rolling_dy(VOL, rdates, p_vol, H, 200)),
           `250` = list(r = rr, v = rv),
           `300` = list(r = rolling_dy(R, rdates, p_ret, H, 300), v = rolling_dy(VOL, rdates, p_vol, H, 300)))
cd <- Reduce(intersect, lapply(rw, function(z) as.character(z$r$dates)))
get_tci <- function(Rl) Rl$TCI[match(cd, as.character(Rl$dates))]
wtab <- do.call(rbind, lapply(c("r", "v"), function(b) do.call(rbind, lapply(c("200", "300"), function(w) {
  a <- get_tci(rw[["250"]][[b]]); o <- get_tci(rw[[w]][[b]])
  data.frame(Basis = ifelse(b == "r", "Returns", "Volatility"), Window = w, Mean_TCI_250 = mean(a), Mean_TCI_alt = mean(o), Corr = cor(a, o), MAD = mean(abs(a - o)), RMSE = sqrt(mean((a - o)^2)),
             Sign_agree_GREEN_EQ = mean(sign(rw[["250"]][[b]]$NET[match(cd, as.character(rw[["250"]][[b]]$dates)), "GREEN_EQ"]) == sign(rw[[w]][[b]]$NET[match(cd, as.character(rw[[w]][[b]]$dates)), "GREEN_EQ"])),
             Sign_agree_GREEN_BOND = mean(sign(rw[["250"]][[b]]$NET[match(cd, as.character(rw[["250"]][[b]]$dates)), "GREEN_BOND"]) == sign(rw[[w]][[b]]$NET[match(cd, as.character(rw[[w]][[b]]$dates)), "GREEN_BOND"]))) }))))
rownames(wtab) <- NULL; cat("\nRolling-window robustness (common dates):\n"); print(round(wtab[, -c(1, 2)], 3)); print(wtab[, 1:2])
save_tab(wtab, "25_rolling_window_robust")
wd <- do.call(rbind, lapply(c("200", "250", "300"), function(w) rbind(
  data.frame(Date = rw[[w]]$r$dates, TCI = rw[[w]]$r$TCI, Window = w, Basis = "Returns"), data.frame(Date = rw[[w]]$v$dates, TCI = rw[[w]]$v$TCI, Window = w, Basis = "GK volatility"))))
p_w3 <- ggplot(wd, aes(Date, TCI, colour = Window)) + geom_line(linewidth = 0.45) + add_events() + facet_wrap(~Basis, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("200" = "#1B9E77", "250" = "#D95F02", "300" = "#7570B3")) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Rolling TCI under three window lengths", x = NULL, y = "TCI (%)", caption = ev_caption()) + theme_dy(10)
show(p_w3, "25_window_robust", 10, 6.5)

# role stability through time: share of days GREEN_EQ / GREEN_BOND are net transmitters (volatility and returns)
role <- data.frame(Asset = rep(c("GREEN_EQ", "GREEN_BOND"), each = 2), Basis = rep(c("Returns", "Volatility"), 2),
  Share_days_net_transmitter = c(mean(rr$NET[, "GREEN_EQ"] > 0), mean(rv$NET[, "GREEN_EQ"] > 0), mean(rr$NET[, "GREEN_BOND"] > 0), mean(rv$NET[, "GREEN_BOND"] > 0)),
  Mean_NET = c(mean(rr$NET[, "GREEN_EQ"]), mean(rv$NET[, "GREEN_EQ"]), mean(rr$NET[, "GREEN_BOND"]), mean(rv$NET[, "GREEN_BOND"])))
cat("\nRole persistence in rolling windows:\n"); print(role, digits = 3, row.names = FALSE); save_tab(role, "25_role_persistence")

# export data for the independent Python cross-check
write.csv(data.frame(Date = rdates, R), "xval_returns.csv", row.names = FALSE); write.csv(data.frame(Date = rdates, VOL), "xval_vol.csv", row.names = FALSE)
xv <- list(TCI_r = st_r$TCI, TCI_v = st_v$TCI, NET_r = st_r$NET, NET_v = st_v$NET, p_ret = p_ret, p_vol = p_vol, sel_r = sel_r, sel_v = sel_v,
           roll_tci_r_last = tail(rr$TCI, 3), roll_tci_v_last = tail(rv$TCI, 3), pt_r = portmanteau(fit_r, 20), jb_r = jb_mv(fit_r), arch_r = arch_mv(fit_r, 5))
saveRDS(xv, "xval_R.rds")

# ---- step26_scorecard ----
# STEP 26 - Hypothesis scorecard (hypotheses are derived from RQ1-RQ7 of the literature review; each row cites its own test)
ie <- match("GREEN_EQ", assets); ib <- match("GREEN_BOND", assets)
nd_eq_bond <- bs_v$NET[, ie] - bs_v$NET[, ib]; nd_eq_bond_r <- bs_r$NET[, ie] - bs_r$NET[, ib]
ci_row <- function(basis, a) ci_all[ci_all$Basis == basis & ci_all$Asset == a, ]
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
sg <- function(b, g) { s <- src[src$Basis == b & src$Receiver == g, ]; s <- s[order(-s$Share), ]; paste0(s$Source[1], " ", fmt(s$Share[1], 1), "% [", fmt(s$Lo[1], 1), ", ", fmt(s$Hi[1], 1), "]") }
qdv <- qd[qd$Basis == "GK volatility" & grepl("0.95", qd$Contrast), ]; qdr <- qd[qd$Basis == "Returns" & grepl("0.05", qd$Contrast), ]
fr_v_long <- ftab[ftab$Basis == "GK volatility" & ftab$Band == "Long (>20d)", ]; fr_r_short <- ftab[ftab$Basis == "Returns" & ftab$Band == "Short (1-5d)", ]
inc_g <- function(b, m) inc_tab[inc_tab$Basis == b & inc_tab$Measure == m, ]
scorecard <- data.frame(
  Hypothesis = c(
    "H1a (RQ1): system connectedness is significantly above the no-spillover benchmark",
    "H1b (RQ1): volatility connectedness exceeds return connectedness",
    "H2a (RQ2): GREEN_EQ is a net transmitter of volatility",
    "H2b (RQ2): GREEN_BOND is a net receiver of volatility and differs from GREEN_EQ",
    "H3 (RQ3): conventional markets are the main sources of the green assets' variance",
    "H4a (RQ4): connectedness rises around COVID-19 and Russia-Ukraine",
    "H4b (RQ4): the Bitcoin-ETF approval does not raise system connectedness",
    "H5a (RQ5): connectedness is higher in extreme (tail) states than in the median state",
    "H5b (RQ5): the horizon structure differs between returns and volatility",
    "H6 (RQ6): adding the green assets raises system connectedness and changes conventional roles",
    "H7 (RQ7): green assets reduce SP500 risk (hedging / diversification)"),
  Evidence = c(
    paste0("TCI returns ", fmt(st_r$TCI, 1), "%, volatility ", fmt(st_v$TCI, 1), "%; circular-shift null 95th pct ", fmt(perm_tab$Null_95th[1], 1), "% / ", fmt(perm_tab$Null_95th[2], 1), "%; p = ", fmt(perm_tab$p_value[1], 3)),
    paste0("Static difference ", fmt(tci_diff["diff"], 1), " pp, paired 95% CI [", fmt(tci_diff["lo"], 1), ", ", fmt(tci_diff["hi"], 1), "]; rolling bias-adjusted mean TCI returns ", fmt(mean(rr$TCI) - nb_r, 1), " vs volatility ", fmt(mean(rv$TCI) - nb_v, 1)),
    paste0("NET volatility ", fmt(ci_row("GK volatility", "GREEN_EQ")$Estimate, 1), " pp, 95% CI [", fmt(ci_row("GK volatility", "GREEN_EQ")$Lower95, 1), ", ", fmt(ci_row("GK volatility", "GREEN_EQ")$Upper95, 1), "]; positive in ", sum(rob$NET_GREEN_EQ[rob$Basis == "Volatility"] > 0), " of ", sum(rob$Basis == "Volatility"), " volatility specifications"),
    paste0("GREEN_BOND NET volatility ", fmt(ci_row("GK volatility", "GREEN_BOND")$Estimate, 1), " pp, CI [", fmt(ci_row("GK volatility", "GREEN_BOND")$Lower95, 1), ", ", fmt(ci_row("GK volatility", "GREEN_BOND")$Upper95, 1), "]; EQ minus BOND = ", fmt(mean(ci_row("GK volatility", "GREEN_EQ")$Estimate - ci_row("GK volatility", "GREEN_BOND")$Estimate), 1), " pp, CI [", fmt(quantile(nd_eq_bond, .025), 1), ", ", fmt(quantile(nd_eq_bond, .975), 1), "]"),
    paste0("GREEN_EQ volatility variance: own shocks ", fmt(src$Own[src$Receiver == "GREEN_EQ" & src$Basis == "GK volatility"][1], 0), "%, conventional assets ", fmt(sum(src$Share[src$Receiver == "GREEN_EQ" & src$Basis == "GK volatility" & src$Source %in% conv]), 0), "%, GREEN_BOND ", fmt(src$Share[src$Receiver == "GREEN_EQ" & src$Basis == "GK volatility" & src$Source == "GREEN_BOND"], 0), "% (largest single source: ", sg("GK volatility", "GREEN_EQ"), "). GREEN_BOND volatility variance: own ", fmt(src$Own[src$Receiver == "GREEN_BOND" & src$Basis == "GK volatility"][1], 0), "%, conventional ", fmt(sum(src$Share[src$Receiver == "GREEN_BOND" & src$Basis == "GK volatility" & src$Source %in% conv]), 0), "%, GREEN_EQ ", fmt(src$Share[src$Receiver == "GREEN_BOND" & src$Basis == "GK volatility" & src$Source == "GREEN_EQ"], 0), "%; returns: largest source of GREEN_BOND is ", sg("Returns", "GREEN_BOND")),
    paste0("Rolling TCI change (60d post - pre): COVID returns ", fmt(ev_tci$Ret_change[1], 1), ", volatility ", fmt(ev_tci$Vol_change[1], 1), "; Russia-Ukraine ", fmt(ev_tci$Ret_change[2], 1), " / ", fmt(ev_tci$Vol_change[2], 1), "; 200d static test significant for: ", paste(unique(ev_st$Event[ev_st$Significant == "Yes"]), collapse = ", ")),
    paste0("Bitcoin ETF: rolling TCI change returns ", fmt(ev_tci$Ret_change[4], 1), ", volatility ", fmt(ev_tci$Vol_change[4], 1), "; systems with a significant 200d static change: ", paste(ev_st$Basis[ev_st$Event == "US spot Bitcoin ETFs" & ev_st$Significant == "Yes"], collapse = ", ")),
    paste0("Volatility TCI at tau=0.95 minus median: ", fmt(qdv$Diff, 1), " pp, CI [", fmt(qdv$Lo, 1), ", ", fmt(qdv$Hi, 1), "]; returns bearish minus median: ", fmt(qdr$Diff, 1), " pp, CI [", fmt(qdr$Lo, 1), ", ", fmt(qdr$Hi, 1), "]"),
    paste0("Returns: ", fmt(fr_r_short$Share_of_total, 0), "% of connectedness is short-horizon (1-5d); volatility: ", fmt(fr_v_long$Share_of_total, 0), "% is long-horizon (>20d)"),
    paste0("dTCI (6 minus 4 assets): returns ", fmt(inc_g("Returns", "dTCI")$Estimate, 1), " [", fmt(inc_g("Returns", "dTCI")$Lo95, 1), ", ", fmt(inc_g("Returns", "dTCI")$Hi95, 1), "], volatility ", fmt(inc_g("GK volatility", "dTCI")$Estimate, 1), " [", fmt(inc_g("GK volatility", "dTCI")$Lo95, 1), ", ", fmt(inc_g("GK volatility", "dTCI")$Hi95, 1), "]; S&P 500 NET shifts by ", fmt(inc_g("Returns", "dNET_SP500")$Estimate, 1), " pp [", fmt(inc_g("Returns", "dNET_SP500")$Lo95, 1), ", ", fmt(inc_g("Returns", "dNET_SP500")$Hi95, 1), "]"),
    paste0("HE of SP500 hedged with GREEN_EQ: in-sample same-day ", fmt(100 * he_ci$HE_in[1], 1), "% [", fmt(100 * he_ci$Lo95[1], 1), ", ", fmt(100 * he_ci$Hi95[1], 1), "]; out-of-sample ratio, same-day ", fmt(100 * hedge_tab$HE_oos[hedge_tab$Hedge_asset == "GREEN_EQ"], 1), "%; next-day (lagged hedge, in-sample ratio) ", fmt(100 * hedge_tab$HE_nextday[hedge_tab$Hedge_asset == "GREEN_EQ"], 1), "%; next-day with out-of-sample ratio ", fmt(100 * hedge_tab$HE_next_oos[hedge_tab$Hedge_asset == "GREEN_EQ"], 1), "%; GREEN_BOND in-sample ", fmt(100 * he_ci$HE_in[2], 1), "%")),
  # Decisions below were assigned after reviewing the Evidence column of the run on the supplied data (2017-03-07 to 2025-12-31);
  # re-check them if the data or settings change.
  Decision = c("Supported",
               "Not supported",
               "Not supported (sign positive in 11 of 12 specifications, but the CI includes zero; S&P 500 is the significant volatility transmitter)",
               "Supported",
               "Partly supported (green equity: yes; green bond: mostly own shocks)",
               "Supported (rolling TCI); formal 200-day test significant for Russia-Ukraine only",
               "Supported",
               "Supported (volatility); not supported for returns",
               "Supported",
               "Supported",
               "Partly supported (large same-day hedging effectiveness, none from a lagged hedge)"),
  stringsAsFactors = FALSE)
print(scorecard[, c("Hypothesis", "Decision")], right = FALSE, row.names = FALSE)
for (i in seq_len(nrow(scorecard))) cat("\n[", scorecard$Hypothesis[i], "]\n  Evidence: ", scorecard$Evidence[i], "\n", sep = "")
save_tab(scorecard, "26_hypothesis_scorecard")

# ---- step27_export ----
# STEP 27 - Export every result table to one Excel workbook for the appendix
library(openxlsx)
wb <- createWorkbook()
for (f in sort(list.files("tables", pattern = "\\.csv$", full.names = TRUE))) {
  nm <- substr(sub("\\.csv$", "", basename(f)), 1, 31); addWorksheet(wb, nm); writeData(wb, nm, read.csv(f, check.names = FALSE))
}
saveWorkbook(wb, "DY_Six_Asset_Results.xlsx", overwrite = TRUE)
cat("Workbook written with", length(names(wb)), "sheets.\n")

