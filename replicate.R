# Run from the project root with Rscript --vanilla replicate.R.
# foreign is a recommended package distributed with standard R installations.
d <- foreign::read.dta(file.path("data", "FAid_Final.dta"),
                       convert.factors = FALSE)
stopifnot(!anyDuplicated(d[c("obs", "year")]))
previous <- match(paste(d$obs, d$year - 1), paste(d$obs, d$year))
d$lag_wheat <- d$US_wheat_production[previous] / 1000
d$instrument <- d$lag_wheat * d$fadum_avg
d$wheat_aid <- d$wheat_aid / 1000
d <- subset(d, year >= 1971 & year <= 2006)
d$gdp_avg <- ave(d$ln_rgdpch, d$risocode, FUN = function(x) mean(x, na.rm = TRUE))

weather <- paste0(rep(c("all_Precip_", "all_Temp_"), each = 12),
                  rep(tolower(month.abb), 2))
weather_interactions <- paste0(weather, "_faavg")
for (i in seq_along(weather)) {
  d[[weather_interactions[i]]] <- d[[weather[i]]] * d$fadum_avg
}
d$oil_exposure <- d$oil_price_2011_USD * d$fadum_avg
d$income_exposure <- log(d$USA_rgdpch) * d$fadum_avg
d$pres_exposure <- d$US_president_democ * d$fadum_avg
country_averages <- c("gdp_avg", "real_usmilaid_avg", "real_us_nonfoodaid_ecaid_avg",
                      "recipient_pc_cereals_prod_avg", "cereal_pc_import_quantity_avg")
controls <- c(weather, weather_interactions, "oil_exposure", "income_exposure",
              "pres_exposure", country_averages)
required <- c(controls, "intra_state", "wheat_aid", "instrument", "risocode", "wb_region")
s <- d[complete.cases(d[required]) & d$risocode != "" & d$wb_region != "", ]
s$country <- factor(s$risocode)
s$region_year <- interaction(s$wb_region, s$year, drop = TRUE)
s$year_fe <- factor(s$year)
s$time <- s$year - 1971
n <- nrow(s)
g <- nlevels(s$country)
stopifnot(n == 4089, !anyDuplicated(s[c("country", "year")]))

# Partial out the same nuisance regressors from Y, X and Z.
baseline <- reformulate(c("country", "region_year", setdiff(controls, country_averages),
                           paste0("year_fe:", country_averages)), intercept = FALSE)
estimate <- function(specification, label) {
  design <- model.matrix(specification, s)
  design <- sweep(design, 2, pmax(sqrt(colSums(design^2)), 1e-12), "/")
  nuisance <- qr(design, tol = 1e-8)
  residuals <- qr.resid(nuisance, as.matrix(s[c("intra_state", "wheat_aid", "instrument")]))
  y <- residuals[, 1]
  x <- residuals[, 2]
  z <- residuals[, 3]
  df <- n - nuisance$rank - 1
  correction <- (g / (g - 1)) * ((n - 1) / df)
  fit <- function(outcome, regressor, instrument = regressor, adjustment = correction) {
    denominator <- sum(instrument * regressor)
    beta <- sum(instrument * outcome) / denominator
    score <- instrument * (outcome - regressor * beta)
    se <- sqrt(adjustment * sum(rowsum(score, s$country)^2)) / abs(denominator)
    c(beta = beta, se = se)
  }
  fits <- rbind(fit(y, x), fit(y, z), fit(x, z), fit(y, x, z, 1))
  # The published reduced form expresses the production shock in million tons.
  fits[2, ] <- fits[2, ] * 1000
  t_ratio <- fits[, "beta"] / fits[, "se"]
  p <- 2 * pt(-abs(t_ratio), df = g - 1)
  p[4] <- 2 * pnorm(-abs(t_ratio[4]))
  data.frame(specification = label, model = c("OLS", "Reduced form", "First stage", "2SLS"),
             fits, p = p, n = n, clusters = g, residual_df = df,
             first_stage_F = t_ratio[3]^2, row.names = NULL)
}
base <- estimate(baseline, "Baseline")
trends <- estimate(update(baseline, ~ . + country:time), "Country trends")
results <- rbind(base, trends)

# Published Table 2, column 6, panels A, B, D and C, in that order.
published <- rbind(
  c(-0.00005, 0.00017),
  c( 0.00909, 0.00322),
  c( 0.00358, 0.00103),
  c( 0.00254, 0.00088)
)
stopifnot(all(round(as.matrix(base[c("beta", "se")]), 5) == published))

dir.create("results", showWarnings = FALSE)
write.csv(results, "results/estimates.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/session.txt")
table_tex <- function(path, header, rows, column_format) {
  line <- function(x) paste(paste(x, collapse = " & "), "\\\\")
  writeLines(c(paste0("\\begin{tabular}{", column_format, "}"), "\\toprule",
               line(header), "\\midrule", apply(rows, 1, line),
               "\\bottomrule", "\\end{tabular}"), path)
}
table_tex("results/main.tex", c("Table 2, column 6", "Published", "Replication"),
          cbind(base$model, sprintf("%.5f (%.5f)", published[, 1], published[, 2]),
                sprintf("%.5f (%.5f)", base$beta, base$se)), "lrr")
extension_rows <- function(x) {
  iv <- x[x$model == "2SLS", ]
  fs <- x[x$model == "First stage", ]
  c(sprintf("%.5f (%.5f)", iv$beta, iv$se), sprintf("%.4f", iv$p),
    sprintf("%.5f (%.5f)", fs$beta, fs$se), sprintf("%.2f", fs$first_stage_F),
    iv$n, iv$clusters)
}
table_tex("results/extension.tex", c("Statistic", "Baseline", "Country trends"),
          cbind(c("2SLS coefficient (SE)", "2SLS $p$-value", "First stage (SE)",
                  "Clustered first-stage Wald $F$", "Country-years", "Countries"),
                extension_rows(base), extension_rows(trends)), "lrr")

# Figures 3 and 4 in one two-panel chart, using the common regression sample.
s$group <- ifelse(s$fadum_avg < 0.30, "Irregular recipients", "Regular recipients")
annual <- aggregate(cbind(intra_state, lag_wheat) ~ year + group, s, mean)
write.csv(annual, "results/figure_data.csv", row.names = FALSE)
pdf("results/conflict_and_wheat.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.4, 3, 0.7))
for (group_name in unique(sort(annual$group))) {
  a <- subset(annual, group == group_name)
  a$wheat_million <- a$lag_wheat / 1000
  plot(a$wheat_million, a$intra_state, pch = 19, col = "#236578",
       xlim = range(annual$lag_wheat / 1000), ylim = c(0, max(0.4, annual$intra_state)),
       xlab = "Lagged US wheat production (million metric tons)",
       ylab = "Share experiencing civil conflict", main = group_name, cex.lab = 0.85)
  abline(lm(intra_state ~ wheat_million, a), col = "#702C3D", lwd = 2)
}
invisible(dev.off())
print(results, row.names = FALSE)
