#' Internal function of did_multiplegt_dyn - bootstrap se
#'
#' Bootstrap standard errors for did_multiplegt_dyn(). The implementation is
#' weight-based: instead of physically duplicating rows for sampled units, we
#' draw a multiplicity per unit (or cluster) and merge that count back into
#' the data as a weight column. Inside the slim engine, the bootstrap weights
#' are multiplied with the user weight (if any) so that controls regressions,
#' residualizing, parametric trends and the final estimation all see a
#' coherent weight.
#'
#' Each call to the slim engine runs in its own callr::r() subprocess. We do
#' this because polars allocates its working buffers in Rust's allocator, and
#' the OS only reclaims those pages when the process exits -- gc()/rm() in R
#' do nothing for them. After the subprocess returns the three numeric
#' vectors we need (effects/ATE/placebos), it dies and its memory pages go
#' back to the OS.
#'
#' For verification (see DID_BOOTSTRAP_SAMPLE_DIR option below) the function
#' can also write/read a CSV of the per-replicate unit selections, so that
#' the legacy row-replication path and the new weight-based path can be
#' driven by the exact same bootstrap draws and their outputs compared.
#'
#' @param df df
#' @param outcome outcome
#' @param group group
#' @param time time
#' @param treatment treatment
#' @param effects effects
#' @param placebo placebo
#' @param ci_level ci_level
#' @param switchers switchers
#' @param trends_nonparam trends_nonparam
#' @param weight weight
#' @param controls controls
#' @param dont_drop_larger_lower dont_drop_larger_lower
#' @param drop_if_d_miss_before_first_switch drop_if_d_miss_before_first_switch
#' @param cluster cluster
#' @param same_switchers same_switchers
#' @param same_switchers_pl same_switchers_pl
#' @param only_never_switchers only_never_switchers
#' @param effects_equal effects_equal
#' @param save_results save_results
#' @param normalized normalized
#' @param predict_het predict_het
#' @param trends_lin trends_lin
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @param bootstrap bootstrap
#' @param bootstrap_seed bootstrap_seed
#' @param base base
#' @returns A list of the final results updated with the bootstrap standard errors
#' @noRd
did_multiplegt_bootstrap <- function(
  df,
  outcome,
  group,
  time,
  treatment,
  effects,
  placebo,
  ci_level,
  switchers,
  trends_nonparam,
  weight,
  controls,
  dont_drop_larger_lower,
  drop_if_d_miss_before_first_switch,
  cluster,
  same_switchers,
  same_switchers_pl,
  only_never_switchers,
  effects_equal,
  save_results,
  normalized,
  predict_het,
  trends_lin,
  less_conservative_se,
  continuous,
  bootstrap,
  bootstrap_seed = NULL,
  base
){

    # ------------------------------------------------------------------
    # Resolve which path to take and where (if anywhere) to mirror the
    # bootstrap draws on disk. Three options drive the behaviour, all
    # OFF by default so production users see no change.
    # ------------------------------------------------------------------
    method <- tolower(getOption("DID_BOOTSTRAP_METHOD", "subprocess"))
    if (!method %in% c("subprocess", "row_replication")) {
        warning("Unknown DID_BOOTSTRAP_METHOD '", method,
                "', falling back to 'subprocess'.")
        method <- "subprocess"
    }
    sample_dir   <- getOption("DID_BOOTSTRAP_SAMPLE_DIR", NULL)
    sample_load  <- isTRUE(getOption("DID_BOOTSTRAP_LOAD_SAMPLES", FALSE))
    if (!is.null(sample_dir)) {
        dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    }

    ## Set seed if provided
    base_seed <- if (!is.null(bootstrap_seed)) bootstrap_seed else as.integer(Sys.time())
    if (!sample_load) set.seed(base_seed)

    bresults_effects <- NULL
    bresults_ATE     <- NULL
    bresults_placebo <- NULL

    n_effects <- nrow(base$Effects)
    bresults_effects <- matrix(NA, nrow = bootstrap, ncol = n_effects)
    if (isFALSE(trends_lin)) {
        bresults_ATE <- matrix(NA, nrow = bootstrap, ncol = 1)
    }
    if (placebo > 0) {
        n_placebo <- nrow(base$Placebos)
        bresults_placebo <- matrix(NA, nrow = bootstrap, ncol = n_placebo)
    }

    bs_group <- if (!is.null(cluster)) cluster else group
    df       <- as.data.frame(df)

    # ------------------------------------------------------------------
    # Pre-compute the unique resampling units. Sampling will draw N_units
    # units with replacement on each replicate, where N_units is the
    # number of distinct values of bs_group in df.
    # ------------------------------------------------------------------
    bs_vals  <- df[[bs_group]]
    keep_bs  <- !is.na(bs_vals)
    unit_ids <- sort(unique(bs_vals[keep_bs]))
    n_units  <- length(unit_ids)

    # Index from unit value -> row positions, used by the legacy
    # row-replication path. Built once.
    if (method == "row_replication") {
        df_pos <- split(seq_len(nrow(df))[keep_bs], match(bs_vals[keep_bs], unit_ids))
    }

    # ------------------------------------------------------------------
    # Identify the slim set of columns the engine actually needs. Anything
    # that is not used by did_multiplegt_main is dropped before we hand
    # the data over to the subprocess: smaller payloads, fewer pages
    # touched in Rust, and a much smaller serialized blob.
    # ------------------------------------------------------------------
    needed_cols <- unique(c(
        outcome, group, time, treatment,
        trends_nonparam, weight, controls, cluster,
        unlist(predict_het[1])
    ))
    needed_cols <- needed_cols[!is.null(needed_cols) &
                               !is.na(needed_cols) &
                               nzchar(needed_cols)]
    needed_cols <- intersect(needed_cols, names(df))
    df_slim     <- df[, needed_cols, drop = FALSE]

    # ------------------------------------------------------------------
    # Helper: draw a vector of unit multiplicities of length n_units that
    # sum to n_units (cluster bootstrap). When DID_BOOTSTRAP_LOAD_SAMPLES
    # is set, the draw is read back from the CSV that the previous run
    # produced -- this is what lets us prove the new path matches the old
    # path on the same draws.
    #
    # NOTE on by-stratification: did_multiplegt_dyn calls this helper once
    # per by-stratum. Each call has its OWN unit_ids (different counties
    # belong to each stratum). We therefore stamp the CSV filenames with a
    # tag derived from unit_ids so two strata never collide. Both the
    # writer and the reader derive the same tag from the same unit_ids, so
    # a load-samples run finds the same files the write run produced.
    # ------------------------------------------------------------------
    boot_weight_col <- "boot_weight_XX"
    stratum_tag <- local({
        # Use rlang::hash if available, else a simple, deterministic
        # fingerprint of (length, sum, min, max). Truncate to 10 chars.
        if (requireNamespace("rlang", quietly = TRUE)) {
            substr(rlang::hash(unit_ids), 1, 10)
        } else {
            sprintf("n%d_%s_%s_%s",
                    length(unit_ids),
                    format(sum(as.numeric(unit_ids)), scientific = FALSE),
                    format(min(unit_ids)),
                    format(max(unit_ids)))
        }
    })
    sample_csv_path <- function(j) {
        file.path(sample_dir, sprintf("rep_%05d_%s.csv", j, stratum_tag))
    }

    draw_multiplicities <- function(j) {
        if (sample_load && !is.null(sample_dir)) {
            f <- sample_csv_path(j)
            if (!file.exists(f)) {
                stop("DID_BOOTSTRAP_LOAD_SAMPLES is on but ", f, " is missing.")
            }
            tab <- utils::read.csv(f, stringsAsFactors = FALSE)
            counts <- integer(n_units)
            m <- match(tab$unit_id, unit_ids)
            counts[m[!is.na(m)]] <- as.integer(tab$count[!is.na(m)])
            return(counts)
        }
        sampled <- sample.int(n_units, n_units, replace = TRUE)
        tabulate(sampled, nbins = n_units)
    }

    write_sample_csv <- function(j, counts) {
        if (is.null(sample_dir) || sample_load) return(invisible(NULL))
        nz <- counts > 0
        utils::write.csv(
            data.frame(unit_id = unit_ids[nz], count = counts[nz]),
            sample_csv_path(j),
            row.names = FALSE
        )
    }

    # ------------------------------------------------------------------
    # The slim engine runs inside callr::r() in subprocess mode. We pull
    # the necessary R packages explicitly so that the worker has
    # everything it needs and nothing else.
    # ------------------------------------------------------------------
    use_callr <- method == "subprocess" && requireNamespace("callr", quietly = TRUE)

    # Build the constant kwargs for the slim engine once.
    smaller_kwargs <- list(
        boot_weight_col = boot_weight_col,
        outcome = outcome, group = group, time = time, treatment = treatment,
        effects = effects, placebo = placebo, ci_level = ci_level,
        switchers = switchers, trends_nonparam = trends_nonparam,
        weight = weight, controls = controls,
        dont_drop_larger_lower = dont_drop_larger_lower,
        drop_if_d_miss_before_first_switch = drop_if_d_miss_before_first_switch,
        cluster = cluster, same_switchers = same_switchers,
        same_switchers_pl = same_switchers_pl,
        only_never_switchers = only_never_switchers,
        effects_equal = effects_equal, save_results = save_results,
        normalized = normalized, predict_het = predict_het,
        trends_lin = trends_lin,
        less_conservative_se = less_conservative_se,
        continuous = continuous
    )

    # Worker for the new (subprocess + weights) path.
    #
    # callr launches a fresh R via fork()+exec(). Under memory pressure or
    # transient scheduler hiccups this can fail with errors like "could not
    # start R, exited with non-zero status, has crashed or was killed",
    # especially after many successful iterations have left the parent's
    # address space large. We retry up to MAX_TRIES, gc()'ing the parent
    # between attempts to release committed pages.
    .MAX_SUBPROC_TRIES <- 3L
    run_subprocess <- function(df_boot) {
        if (!use_callr) {
            # Fallback: same engine, same process. Useful when callr is
            # not available, but does NOT free the polars buffers between
            # iterations.
            return(do.call(did_multiplegt_main_smaller,
                           c(list(df = df_boot), smaller_kwargs)))
        }
        last_err <- NULL
        for (attempt in seq_len(.MAX_SUBPROC_TRIES)) {
            res <- tryCatch(
                callr::r(
                    function(df_boot, kwargs) {
                        # polars must be ATTACHED (not just required), because
                        # internal helpers like pl_over_cols look up `pl` via
                        # lexical scoping which falls through to the search path.
                        suppressPackageStartupMessages({
                            for (pkg in c("polars", "DIDmultiplegtDYN")) {
                                if (!paste0("package:", pkg) %in% search()) {
                                    loadNamespace(pkg)
                                    attachNamespace(pkg)
                                }
                            }
                        })
                        do.call(
                            getFromNamespace("did_multiplegt_main_smaller",
                                             "DIDmultiplegtDYN"),
                            c(list(df = df_boot), kwargs)
                        )
                    },
                    args = list(df_boot = df_boot, kwargs = smaller_kwargs),
                    show = FALSE,
                    spinner = FALSE
                ),
                error = function(e) e
            )
            if (!inherits(res, "error")) return(res)
            last_err <- res
            # Free what we can in the parent before retrying so fork() has
            # room. Small back-off so a transient resource crunch can clear.
            gc(verbose = FALSE)
            Sys.sleep(0.5 * attempt)
        }
        stop("did_multiplegt_bootstrap: subprocess failed after ",
             .MAX_SUBPROC_TRIES, " attempts: ",
             conditionMessage(last_err))
    }

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------
    for (j in seq_len(bootstrap)) {

        counts <- draw_multiplicities(j)
        write_sample_csv(j, counts)

        if (method == "subprocess") {
            # ---- Weight-based path ----
            sel <- which(counts > 0)
            count_lookup <- data.frame(
                .__unit__ = unit_ids[sel],
                .__count__ = counts[sel],
                stringsAsFactors = FALSE
            )
            names(count_lookup)[1] <- bs_group
            names(count_lookup)[2] <- boot_weight_col

            df_boot <- merge(df_slim, count_lookup, by = bs_group, all = FALSE)

            res <- run_subprocess(df_boot)

            n_res_effects <- length(res$Effects)
            if (n_res_effects > 0) {
                n_copy <- min(ncol(bresults_effects), n_res_effects)
                bresults_effects[j, seq_len(n_copy)] <- res$Effects[seq_len(n_copy)]
            }
            if (!is.null(bresults_ATE) && length(res$ATE) >= 1 &&
                !is.na(res$ATE[1])) {
                bresults_ATE[j, 1] <- res$ATE[1]
            }
            if (!is.null(bresults_placebo) && length(res$Placebos) > 0) {
                n_copy <- min(ncol(bresults_placebo), length(res$Placebos))
                bresults_placebo[j, seq_len(n_copy)] <- res$Placebos[seq_len(n_copy)]
            }

            rm(df_boot, count_lookup, res)
            # Release R-side allocations of the merged df_boot etc. so the
            # parent doesn't keep growing across iterations. The polars
            # buffers held by parent are unrelated (they live in Rust); the
            # gc() here is about keeping fork() room available for the
            # next iteration's subprocess.
            invisible(gc(verbose = FALSE))

        } else {
            # ---- Legacy row-replication path ----
            sel <- which(counts > 0)
            idx <- unlist(lapply(sel, function(k) {
                rep(df_pos[[k]], counts[k])
            }))
            df_boot <- df[idx, , drop = FALSE]
            df_boot <- df_boot[order(df_boot[[group]], df_boot[[time]]), ]

            suppressMessages({
            df_est <- did_multiplegt_main(df = df_boot, outcome = outcome, group = group, time = time, treatment = treatment, effects = effects, placebo = placebo, ci_level = ci_level, switchers = switchers, trends_nonparam = trends_nonparam, weight = weight, controls = controls, dont_drop_larger_lower = dont_drop_larger_lower, drop_if_d_miss_before_first_switch = drop_if_d_miss_before_first_switch, cluster = cluster, same_switchers = same_switchers, same_switchers_pl = same_switchers_pl, only_never_switchers = only_never_switchers, effects_equal = effects_equal, save_results = save_results, normalized = normalized, predict_het = predict_het, trends_lin = trends_lin, less_conservative_se = less_conservative_se, continuous = continuous)})

            res <- df_est$did_multiplegt_dyn

            n_res_effects <- nrow(res$Effects)
            if (n_res_effects > 0) {
                n_copy <- min(ncol(bresults_effects), n_res_effects)
                bresults_effects[j, seq_len(n_copy)] <- res$Effects[seq_len(n_copy), 1]
            }
            if (!is.null(bresults_ATE) && !is.null(res$ATE[1])) {
                bresults_ATE[j, 1] <- res$ATE[1]
            }
            if (!is.null(bresults_placebo) && !is.null(res$Placebos)) {
                n_res_placebo <- nrow(res$Placebos)
                if (n_res_placebo > 0) {
                    n_copy <- min(ncol(bresults_placebo), n_res_placebo)
                    bresults_placebo[j, seq_len(n_copy)] <- res$Placebos[seq_len(n_copy), 1]
                }
            }

            rm(res, df_est, df_boot, idx)
        }

        progressBar(j, bootstrap)
    }

    ci_level <- ci_level / 100

    # Fast C++ SD computation for effects
    effect_sds <- bootstrap_compute_sd_cpp(bresults_effects)
    n_eff <- nrow(base$Effects)
    base$Effects[1:n_eff, 2] <- effect_sds[1:n_eff]

    # Fast C++ CI computation for effects
    ci_effects <- bootstrap_compute_ci_cpp(base$Effects[1:n_eff, 1], effect_sds[1:n_eff], ci_level)
    base$Effects[1:n_eff, 3] <- ci_effects$lb
    base$Effects[1:n_eff, 4] <- ci_effects$ub

    if (nrow(base$Effects) == 1) {
        class(base$Effects) <- "numeric"
    }

    # ATE SE computation using C++
    if (!is.null(bresults_ATE) && !is.null(base$ATE[1])) {
        ate_sd <- bootstrap_compute_sd_cpp(bresults_ATE)
        base$ATE[2] <- ate_sd[1]
        ci_ate <- bootstrap_compute_ci_cpp(base$ATE[1], ate_sd[1], ci_level)
        base$ATE[3] <- ci_ate$lb[1]
        base$ATE[4] <- ci_ate$ub[1]
    }

    # Fast C++ SD computation for placebos
    if (!is.null(bresults_placebo)) {
        placebo_sds <- bootstrap_compute_sd_cpp(bresults_placebo)
        n_pl <- nrow(base$Placebos)
        base$Placebos[1:n_pl, 2] <- placebo_sds[1:n_pl]

        # Fast C++ CI computation for placebos
        ci_placebo <- bootstrap_compute_ci_cpp(base$Placebos[1:n_pl, 1], placebo_sds[1:n_pl], ci_level)
        base$Placebos[1:n_pl, 3] <- ci_placebo$lb
        base$Placebos[1:n_pl, 4] <- ci_placebo$ub

        if (nrow(base$Placebos) == 1) {
            class(base$Placebos) <- "numeric"
        }
    }
    return(base)
}

#' Internal function to convert lists to vectors (optimized)
#' @param lis A list
#' @returns A vector
#' @noRd
list_to_vec <- function(lis) {
    unlist(lis, use.names = FALSE)
}
