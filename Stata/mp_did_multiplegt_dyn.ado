cap program drop mp_did_multiplegt_dyn
program define mp_did_multiplegt_dyn, eclass
version 12.0
syntax varlist(min=4 max=4 numeric) [if] [in] [, effects(integer 1) placebo(integer 1) switchers(string) only_never_switchers controls(varlist numeric) trends_nonparam(varlist numeric) weight(varlist numeric max=1) dont_drop_larger_lower NORMALIZED cluster(varlist numeric max=1) same_switchers same_switchers_pl effects_equal drop_if_d_miss_before_first_switch trends_lin CONTinuous(integer 0) less_conservative_se bootstrap(string)]

    di ""
    di as text "Step 1: Retrieving dynamic effects"
    qui {
        qui did_multiplegt_dyn `varlist', effects(`effects') switchers(`switchers') `only_never_switchers' `effects_equal' controls(`controls') trends_nonparam(`trends_nonparam') weight(`weight') `dont_drop_larger_lower' `normalized' cluster(`cluster') `same_switchers' `same_switchers_pl' `drop_if_d_miss_before_first_switch' `trends_lin' continuous(`continuous') `less_conservative_se' bootstrap(`bootstrap') graph_off

        mat b_dyn = e(b)[1, 1..l_XX]
        mat b_ATE = e(b)[1, l_XX + 1]
        mat V_dyn = e(V)[1..l_XX, 1..l_XX]
        mat V_ATE = e(V)[1+l_XX, 1+l_XX]
        scalar l_XX_og = l_XX
    }

    di "{hline 80}"
    di _skip(13) "{bf:Estimation of treatment effects: Event-study effects}"
    di "{hline 80}"
    matlist mat_res_XX[1..l_XX, 1..6]
    di "{hline 80}"
    cap qui di scalar(p_equality_effects)
    if _rc == 0 {
        di as text "{it:Test of equality of the effects : p-value =} " scalar(p_equality_effects)
    }

    matrix mat_res_avg_XX=mat_res_XX[l_XX+1, 1..6]
    matrix mat_res_avg_XX=(mat_res_avg_XX, .z )
    matrix colnames mat_res_avg_XX= "Estimate" "SE" "LB CI" "UB CI" "N" "Switch" "x Periods"
    display _newline
    di "{hline 80}"
    di _skip(15) "{bf:Average cumulative (total) effect per treatment unit}"
    di "{hline 80}"
    matlist mat_res_avg_XX, nodotz
    di "{hline 80}"

    di ""
    di as text "Step 2: Running again with `placebo' out of `=max_pl_XX' feasible placebo(s)"
    qui {

        preserve
        add_periods `varlist', add(`=max_pl_gap_XX')
        if `effects' < `placebo' {
            local effects = `placebo'
        }

        qui did_multiplegt_dyn `varlist', effects(`effects') placebo(`placebo') switchers(`switchers') `only_never_switchers' `effects_equal' controls(`controls') trends_nonparam(`trends_nonparam') weight(`weight') `dont_drop_larger_lower' `normalized' cluster(`cluster') `same_switchers' `same_switchers_pl' `drop_if_d_miss_before_first_switch' `trends_lin' continuous(`continuous') `less_conservative_se' bootstrap(`bootstrap') graph_off
        mat b = b_dyn, e(b)[1, l_XX+1..l_XX+l_placebo_XX], b_ATE
        mat V = (V_dyn, J(l_XX_og, l_placebo_XX, 0)) \ (J(l_placebo_XX, l_XX_og, 0), e(V)[l_XX+1..l_XX+l_placebo_XX, l_XX+1..l_XX+l_placebo_XX])
        mat V = (V, J(l_XX_og + l_placebo_XX, 1, 0)) \ (J(1, l_XX_og + l_placebo_XX, 0), V_ATE)
        restore

    }

    di "{hline 80}"
    di _skip(10) "{bf:Testing the parallel trends and no anticipation assumptions}"
    di "{hline 80}"
    matlist mat_res_XX[l_XX+2...,1..6]
    di "{hline 80}"
    cap qui di  scalar(p_jointplacebo)
    if _rc == 0 {
        di as text "{it:Test of joint nullity of the placebos : p-value =} " scalar(p_jointplacebo)
    }

    local rown ""
    forv j = 1/`=scalar(l_XX_og)' {
        local rown  "`rown' Effect_`j'"
    }
    forv j = 1/`=scalar(l_placebo_XX)' {
        local rown "`rown' Placebo_`j'" 
    }
    local rown "`rown' Av_tot_eff"
    mat rownames V = `rown'
    mat colnames V = `rown'
    mat colnames b = `rown'
    ereturn post b V
end

cap program drop add_periods
program define add_periods, rclass
syntax varlist, add(string)
qui {
    if `add' != 0 {
        tokenize `varlist'
        sum `3'
        local top = r(max)
        unique `2'
        forv j = 1/`=r(unique)' {
            local N = _N
            insobs `add'    
            replace `1' = 0 if missing(`2')
            replace `4' = 0 if missing(`2')
            replace `2' = `j' if missing(`2')
            forv i = 1/`add' {
                replace `3' = `top' + `i' in `=`N' + `i''
            }
        }
        sort `2' `3'
    }
}
end