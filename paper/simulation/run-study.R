required_wdm_version <- "0.3.0"
if (!requireNamespace("wdm", quietly = TRUE) ||
    packageVersion("wdm") < required_wdm_version) {
    stop("Install wdm >= 0.3.0 before running the simulation.")
}

methods <- c(
    "pearson", "spearman", "kendall", "blomqvist", "hoeffding",
    "chatterjee"
)
significance_level <- 0.05
null_repetitions <- as.integer(Sys.getenv("WDM_NULL_REPS", "2000"))
power_repetitions <- as.integer(Sys.getenv("WDM_POWER_REPS", "1000"))
sample_sizes <- c(50L, 100L, 200L, 500L)

if (is.na(null_repetitions) || null_repetitions < 1L ||
    is.na(power_repetitions) || power_repetitions < 1L) {
    stop("WDM_NULL_REPS and WDM_POWER_REPS must be positive integers.")
}

make_weights <- function(sample_size, weight_profile) {
    if (weight_profile == "equal")
        return(NULL)
    qlnorm(
        (seq_len(sample_size) - 0.5) / sample_size,
        meanlog = 0,
        sdlog = 0.8
    )
}

make_response <- function(predictor, scenario) {
    noise <- rnorm(length(predictor))
    if (scenario == "null")
        return(noise)
    if (scenario == "linear") {
        return(
            0.25 * predictor + sqrt(1 - 0.25^2) * noise
        )
    }
    0.5 * (predictor^2 - 1) / sqrt(2) + sqrt(1 - 0.5^2) * noise
}

simulate_rejections <- function(scenario, sample_size, weight_profile,
                                repetitions, seed) {
    set.seed(seed)
    weights <- make_weights(sample_size, weight_profile)
    rejection_counts <- setNames(numeric(length(methods)), methods)

    for (repetition in seq_len(repetitions)) {
        predictor <- rnorm(sample_size)
        response <- make_response(predictor, scenario)
        p_values <- vapply(methods, function(method) {
            wdm::indep_test(
                predictor,
                response,
                method = method,
                weights = weights,
                alternative = if (method == "chatterjee")
                    "greater"
                else
                    "two-sided"
            )$p_value
        }, numeric(1))
        rejection_counts <- rejection_counts +
            as.numeric(p_values < significance_level)
    }

    effective_sample_size <- if (is.null(weights)) {
        sample_size
    } else {
        sum(weights)^2 / sum(weights^2)
    }
    data.frame(
        scenario = scenario,
        sample_size = sample_size,
        weight_profile = weight_profile,
        effective_sample_size = effective_sample_size,
        method = methods,
        repetitions = repetitions,
        significance_level = significance_level,
        rejection_rate = as.numeric(rejection_counts / repetitions),
        row.names = NULL
    )
}

null_design <- expand.grid(
    weight_profile = c("equal", "lognormal"),
    sample_size = sample_sizes,
    scenario = "null",
    stringsAsFactors = FALSE
)
null_design$seed <- 20260822L + seq_len(nrow(null_design))
power_design <- expand.grid(
    weight_profile = c("equal", "lognormal"),
    sample_size = sample_sizes,
    scenario = c("linear", "quadratic"),
    stringsAsFactors = FALSE
)
power_design$seed <- 20260842L + seq_len(nrow(power_design))

start_time <- proc.time()[["elapsed"]]
null_results <- do.call(rbind, lapply(seq_len(nrow(null_design)), function(i) {
    simulate_rejections(
        null_design$scenario[i],
        null_design$sample_size[i],
        null_design$weight_profile[i],
        null_repetitions,
        null_design$seed[i]
    )
}))
power_results <- do.call(rbind, lapply(seq_len(nrow(power_design)), function(i) {
    simulate_rejections(
        power_design$scenario[i],
        power_design$sample_size[i],
        power_design$weight_profile[i],
        power_repetitions,
        power_design$seed[i]
    )
}))
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
results <- rbind(null_results, power_results)

results_directory <- file.path("paper", "simulation", "results")
if (!dir.exists(results_directory))
    dir.create(results_directory, recursive = TRUE)
write.csv(
    results,
    file.path(results_directory, "rejection-rates.csv"),
    row.names = FALSE
)

rejection_rate <- function(scenario, sample_size, weight_profile, method) {
    sprintf(
        "%.3f",
        results$rejection_rate[
            results$scenario == scenario &
                results$sample_size == sample_size &
                results$weight_profile == weight_profile &
                results$method == method
        ]
    )
}

table_rows <- function(scenario, weight_profile) {
    vapply(sample_sizes, function(sample_size) {
        effective_sample_size <- unique(
            results$effective_sample_size[
                results$scenario == scenario &
                    results$sample_size == sample_size &
                    results$weight_profile == weight_profile
            ]
        )
        paste0(
            sample_size, " & ", sprintf("%.1f", effective_sample_size), " & ",
            paste(
                vapply(methods, function(method) {
                    rejection_rate(
                        scenario, sample_size, weight_profile, method
                    )
                }, character(1)),
                collapse = " & "
            ),
            " \\\\"
        )
    }, character(1))
}
writeLines(
    c(
        "\\begin{table}[t]",
        "  \\centering",
        "  \\small",
        "  \\setlength{\\tabcolsep}{3.5pt}",
        "  \\caption{Empirical rejection probabilities under independence at nominal level $0.05$. The Monte Carlo standard error is at most $0.0049$.}",
        "  \\label{tab:simulation-size}",
        "  \\begin{tabular}{rrcccccc}",
        "    \\toprule",
        "    $n$ & $n_{\\mathrm{eff}}$ & Pear. & Spear. & Kend. & Blom. & Hoef. & Chatt. \\\\",
        "    \\midrule",
        "    \\multicolumn{8}{l}{Equal weights} \\\\",
        paste0("    ", table_rows("null", "equal")),
        "    \\addlinespace",
        "    \\multicolumn{8}{l}{Unequal weights} \\\\",
        paste0("    ", table_rows("null", "lognormal")),
        "    \\bottomrule",
        "  \\end{tabular}",
        "\\end{table}"
    ),
    file.path(results_directory, "null-size.tex")
)

writeLines(
    c(
        "\\begin{table}[p]",
        "  \\centering",
        "  \\small",
        "  \\setlength{\\tabcolsep}{3.5pt}",
        "  \\caption{Empirical power at nominal level $0.05$. The Monte Carlo standard error is at most $0.0158$.}",
        "  \\label{tab:simulation-power}",
        "  \\begin{tabular}{rrcccccc}",
        "    \\toprule",
        "    $n$ & $n_{\\mathrm{eff}}$ & Pear. & Spear. & Kend. & Blom. & Hoef. & Chatt. \\\\",
        "    \\midrule",
        "    \\multicolumn{8}{l}{Linear alternative, equal weights} \\\\",
        paste0("    ", table_rows("linear", "equal")),
        "    \\addlinespace",
        "    \\multicolumn{8}{l}{Linear alternative, unequal weights} \\\\",
        paste0("    ", table_rows("linear", "lognormal")),
        "    \\addlinespace",
        "    \\multicolumn{8}{l}{Quadratic alternative, equal weights} \\\\",
        paste0("    ", table_rows("quadratic", "equal")),
        "    \\addlinespace",
        "    \\multicolumn{8}{l}{Quadratic alternative, unequal weights} \\\\",
        paste0("    ", table_rows("quadratic", "lognormal")),
        "    \\bottomrule",
        "  \\end{tabular}",
        "\\end{table}"
    ),
    file.path(results_directory, "power.tex")
)

writeLines(
    c(
        paste("wdm version:", packageVersion("wdm")),
        paste("R version:", R.version.string),
        paste("null repetitions per design point:", null_repetitions),
        paste("power repetitions per design point:", power_repetitions),
        sprintf("elapsed seconds: %.2f", elapsed_seconds)
    ),
    file.path(results_directory, "run-info.txt")
)

message(sprintf("Simulation completed in %.2f seconds.", elapsed_seconds))
