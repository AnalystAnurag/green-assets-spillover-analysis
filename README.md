# Do Green Assets Really Protect a Stock Portfolio?

**Return and volatility spillovers among the S&P 500, gold, crude oil, Bitcoin and green assets (2017–2025)**

*Diebold–Yilmaz connectedness analysis in R · Term paper, Time Series Analysis, M.Sc. Business Analytics, Manipal Academy of Higher Education*

---

## In one minute

Investors often treat green investments (clean-energy stocks and green bonds) as a diversifier, something that holds up when the stock market falls. We tested that claim on about 2,200 trading days of data across six asset classes.

**Short answer: they don't protect you when it matters.**

- **The S&P 500 drives the system.** It is the only significant net *transmitter* of volatility. Shocks start in US equities and spread to the other assets.
- **Linkages jump in a crisis.** In normal markets, about 20% of each asset's volatility comes from the others. In turbulent markets, that rises to **71%**, so diversification fails exactly when you need it.
- **Clean-energy stocks don't really hedge.** They cut S&P 500 risk by 42% only if you hedge on the same day using information you wouldn't have had in advance. With a realistic hedge set a day ahead, the reduction falls to **0.4%**. Clean-energy equity behaves like more equity exposure, not like protection.
- **Green bonds are different.** They mostly move on their own shocks and *absorb* volatility rather than spreading it. But they are too loosely linked to the S&P 500 to hedge it in a meaningful way (about 4% risk reduction).

---

## The question and the hypotheses

We framed the work as **7 research questions and 11 testable hypotheses**, each with a clear accept or reject rule. The script outputs a final scorecard (`26_hypothesis_scorecard`):

| # | Hypothesis | Result |
|---|---|---|
| H1a | The assets are significantly connected (more than chance would produce) | ✅ Supported |
| H1b | Volatility is more connected than returns | ❌ Not supported |
| H2a | Green equity is a net transmitter of volatility | ❌ Not supported (positive sign, but not significant) |
| H2b | Green bonds are a net receiver and differ from green equity | ✅ Supported |
| H3 | Conventional markets are the main source of green assets' risk | 🟡 Partly (yes for green equity; green bonds mostly move on their own shocks) |
| H4a | Connectedness rises around COVID-19 and the Russia–Ukraine war | ✅ Supported |
| H4b | Approval of the US spot Bitcoin ETFs did not raise connectedness | ✅ Supported |
| H5a | Connectedness is higher in extreme market states | ✅ Volatility / ❌ Returns |
| H5b | Return spillovers are short-lived; volatility spillovers persist | ✅ Supported |
| H6 | Adding green assets raises overall system connectedness | ✅ Supported |
| H7 | Green assets reduce S&P 500 risk | 🟡 Partly (same-day only; no benefit from a realistic lagged hedge) |

---

## Data

| Asset | Proxy | Ticker |
|---|---|---|
| US equities | S&P 500 index | `^GSPC` |
| Gold | Gold futures | `GC=F` |
| Crude oil | WTI crude futures | `CL=F` |
| Bitcoin | Bitcoin / USD | `BTC-USD` |
| Green equity | iShares Global Clean Energy ETF | `ICLN` |
| Green bonds | VanEck Green Bond ETF | `GRNB` |

- **Period:** 7 Mar 2017 – 31 Dec 2025 (the start date is set by GRNB's first price)
- **Frequency:** daily, restricted to dates on which all six assets traded
- **Measures:** log returns, and Garman–Klass volatility built from daily open/high/low/close prices
- **Cleaning:** the April 2020 negative WTI prices are treated as missing, since a log return can't be taken on a negative price

---

## Method

The analysis runs as a 27-step pipeline in one R script. The core VAR, generalized variance decomposition, TVP-VAR, quantile and frequency routines are written from first principles rather than taken from a connectedness package.

1. **Data audit and preparation:** quality checks, a common trading calendar, returns and volatility
2. **Diagnostics:** descriptive statistics, unit-root tests (ADF, KPSS, Johansen), ARCH effects, VAR lag selection and residual tests, robust Granger causality
3. **Static connectedness:** Diebold–Yilmaz (2012) spillover tables from a generalized forecast-error variance decomposition, with **moving-block bootstrap confidence intervals** (1,000 replications) and a **permutation test** against a no-spillover benchmark
4. **Dynamic connectedness:** 250-day rolling windows and a TVP-VAR, annotated with 7 market events (COVID-19, Russia–Ukraine, Fed hikes, FTX, SVB, Bitcoin ETFs, US tariffs)
5. **Event studies:** pre/post-event comparisons around 4 core events
6. **Market states:** quantile connectedness (calm, normal and turbulent markets)
7. **Time horizons:** frequency connectedness (short, medium and long-term spillovers)
8. **Incremental value:** how the system changes when the green assets are added to the four conventional assets
9. **Portfolio application:** hedge ratios and hedging effectiveness (in-sample, out-of-sample, same-day vs next-day, COVID vs calm periods), plus minimum-variance, minimum-connectedness and equal-weight portfolios
10. **Robustness:** 11 alternative specifications, each run on both returns and volatility (forecast horizon, VAR lag order, alternative volatility measures, excluding COVID and the 2022 crypto crash)

---

## Key numbers

| Measure | Returns | Volatility |
|---|---|---|
| Total connectedness (static) | 21.8% | 22.2% |
| S&P 500 net spillover | +4.8 pp (significant) | +14.4 pp (significant) |
| Green equity net spillover | +4.5 pp (significant) | +2.9 pp (not significant) |
| Green bond net spillover | −2.0 pp (significant) | −4.4 pp (significant) |
| Connectedness, median → extreme (95th percentile) state | 21.7% → 22.7% | 19.5% → 71.1% |

| Hedging the S&P 500 with… | Same-day effectiveness | Next-day effectiveness |
|---|---|---|
| Green equity (ICLN) | 41.9% | 0.4% |
| Green bonds (GRNB) | 4.1% | 0.1% |
| Bitcoin | 7.4% | 0.1% |
| Gold | 0.4% | 0.1% |

---

## Repository structure

```
├── Connectedness_Analysis_SP500_Full.R   # full analysis pipeline (27 steps)
├── DY_Six_Asset_Results.xlsx             # all 34 result tables, one sheet per table
├── Term_paper_TSA_Group_34.pdf           # full write-up
├── data/                                 # input price files (see below)
├── figures/                              # created by the script
└── tables/                               # created by the script
```

Each sheet in `DY_Six_Asset_Results.xlsx` is prefixed with the pipeline step that produced it. For example, `12_static_vol` is the static volatility spillover table and `24_hedge_ratios` holds the hedging results.

---

## How to run

**1. Install the R packages**

```r
install.packages(c("ggplot2", "patchwork", "urca", "tseries", "quantreg",
                   "sandwich", "scales", "ggrepel", "ragg", "openxlsx"))
```

**2. Add the data.** Place six daily price files in `data/`, named `SP500.csv`, `GOLD.csv`, `CRUDE.csv`, `BITCOIN.csv`, `GREEN_EQ.csv` and `GREEN_BOND.csv`. Each file needs the columns `Date, Open, High, Low, Close, Adjusted`, which is the standard Yahoo Finance export format for the tickers above.

**3. Run the script** from the repository folder:

```r
source("Connectedness_Analysis_SP500_Full.R")
```

The script creates `figures/` (24 charts) and `tables/` (CSV files), then collects every table into `DY_Six_Asset_Results.xlsx`. The bootstrap steps are compute-heavy, so expect a long run. The random seed is fixed (`set.seed(2026)`), so the results are reproducible.

---

## Limitations

- Only one proxy is used for each green asset class (ICLN and GRNB), so the results may not generalize to other green indices.
- The analysis uses daily data only, so intraday spillovers are not captured.
- The sample starts in 2017 (when GRNB data begins), which limits how many crisis episodes can be studied.

---

## Authors

**Anurag Dash** and **Kratika**, M.Sc. Business Analytics, Manipal School of Commerce & Economics
Guided by **Dr. Satyaban Sahoo**, Assistant Professor

## Key reference

Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. *International Journal of Forecasting*, 28(1), 57–66.
