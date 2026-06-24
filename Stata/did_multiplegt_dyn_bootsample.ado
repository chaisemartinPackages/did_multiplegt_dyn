*! ===========================================================================
*! did_multiplegt_dyn_bootsample.ado
*! ---------------------------------------------------------------------------
*! Standalone cluster sampler used by did_multiplegt_dyn's bsmanual() flow.
*! Draws cluster resamples using `bsample, cluster()` -- the same primitive
*! Stata's built-in `bootstrap, cluster(X) seed(N)` calls internally. Under
*! the same seed, the per-rep cluster-id sequence here is bit-identical to
*! Stata's bootstrap, so a downstream reweighting of N_gt by these counts
*! reproduces the built-in bootstrap distribution.
*!
*! Output CSV (wide):  cluster_id, rep_1, rep_2, ..., rep_REPS
*! ===========================================================================

capture program drop did_multiplegt_dyn_bootsample
program did_multiplegt_dyn_bootsample
    version 12.0
    syntax , Cluster(varname) Reps(integer) Seed(integer) Saving(string)

    quietly {
        * Snapshot the caller's data so we can restore at the end (no
        * `preserve` -- the caller may already have an active preserve and
        * nested preserves are illegal).
        tempfile _bs_caller
        save "`_bs_caller'", replace

        * Build the wide base table: one row per unique cluster, one column
        * per rep initialised to 0 (clusters not drawn in rep b stay 0).
        keep `cluster'
        duplicates drop `cluster', force
        sort `cluster'
        rename `cluster' cluster_id
        forvalues b = 1/`reps' {
            gen long rep_`b' = 0
        }
        tempfile base
        save "`base'", replace

        * Loop reps drawing the cluster resample on the FULL dataset via the
        * same `bsample, cluster()` call Stata's bootstrap uses internally.
        * `set seed` once before the loop -> deterministic identical sequence
        * to what `bootstrap, cluster(`cluster') seed(`seed')` produces.
        set seed `seed'

        forvalues b = 1/`reps' {
            * Load full original data and resample clusters with replacement.
            use "`_bs_caller'", clear
            bsample, cluster(`cluster')

            * Convert the resample into per-cluster counts.
            contract `cluster'
            rename _freq _count_XX
            rename `cluster' cluster_id
            keep cluster_id _count_XX
            tempfile rep_b_data
            save "`rep_b_data'", replace

            * Merge counts onto the base, fill in rep_b, save back.
            use "`base'", clear
            merge 1:1 cluster_id using "`rep_b_data'", nogen
            replace rep_`b' = _count_XX if !missing(_count_XX)
            drop _count_XX
            save "`base'", replace
        }

        * Write wide CSV.
        use "`base'", clear
        order cluster_id rep_1-rep_`reps'
        export delimited using "`saving'", replace

        * Restore caller's original data.
        use "`_bs_caller'", clear
    }
end
