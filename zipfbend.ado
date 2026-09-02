*! version 1.0.2  13july2026 R Bergs
*! Calculates Zipf bend point (Rank - 0.5) and KS-test for city sizes

program define zipfbend, rclass
    version 12.0
    syntax varname(numeric) [if] in [, Generate(string)]

    preserve
    
    * Restrict sample if if/in is specified
    keep `if' `in'
    
    * Exclude missing and non-positive values, then sort city sizes descending
    qui drop if missing(`varlist')
    qui drop if `varlist' <= 0
    qui sort `varlist'
    qui gen long _orig_obs = _n
    qui gsort -`varlist'
    
    * 1. Generate variables (y = ln(size), x = ln(rank - 0.5))
    qui gen long _rank = _n
    qui gen double y = ln(`varlist')
    qui gen double x_obs = ln(_rank - 0.5)
    
    * 2. Linear OLS regression to estimate alpha and beta
    qui regress y x_obs
    scalar s_alpha = _b[_cons]
    scalar s_beta  = _b[x_obs]
    
    * 3. Calculate estimated x-values: x_est = (y - alpha) / beta
    qui gen double x_est = (y - s_alpha) / s_beta
    
    * 4. Calculate difference
    qui gen double diff_x = x_obs - x_est
    
    * Find the maximum positive difference
    qui summarize diff_x
    local max_diff = r(max)
    
    * Find the rank, size, and original ID corresponding to this maximum difference
    qui sort diff_x
    local k_rank = _rank[_N]
    local k_size = `varlist'[_N]
    local k_id   = _orig_obs[_N]
    
    * 5. Prepare Kolmogorov-Smirnov Test 
    qui sort `varlist'
    qui gen double cdf_obs = _n / _N
    
    * Theoretical Pareto CDF based on OLS estimators
    local p_exp = -1 / s_beta
    qui summarize `varlist'
    local s_min = r(min)
    qui gen double cdf_theo = 1 - (`s_min' / `varlist')^(`p_exp')
    
    * Calculate KS test statistic manually (maximum absolute distance)
    qui gen double ks_dist = abs(cdf_obs - cdf_theo)
    qui summarize ks_dist
    local ks_stat = r(max)
    
    * Critical value for KS (Approximation for Alpha=0.05 with large N: 1.36 / sqrt(N))
    local ks_crit = 1.36 / sqrt(_N)
    
    * Display results in the Stata results window
    di as text "{hline 65}"
    di as res " ZIPF DISTRIBUTION & BEND POINT ANALYSIS (Rank - 0.5)"
    di as text "{hline 65}"
    di as text "Number of obs (N):       " as result _N
    di as text "Estimated Zipf coeff (b):   " as result %9.4f s_beta
    di as text "Implied Pareto exponent:    " as result %9.4f `p_exp'
    di as text "{hline 65}"
    di as text "Maximum positive difference:" as result %9.4f `max_diff'
    di as text "-> Occurs at rank:          " as result `k_rank'
    di as text "-> Item size at this obs:   " as result %12.0fc `k_size'
    di as text "{hline 65}"
    di as text "Kolmogorov-Smirnov Test (Goodness-of-Fit for Power Law):"
    di as text "-> KS test statistic (D):   " as result %9.4f `ks_stat'
    di as text "-> Critical value (p=0.05): " as result %9.4f `ks_crit'
    if `ks_stat' > `ks_crit' {
        di as error " H0 (Data follows Pareto) rejected. Non-linear behavior detected."
    }
    else {
        di as result " H0 cannot be rejected. Pure Pareto distribution is plausible."
    }
    di as text "{hline 65}"
    
    * If requested, save the calculated difference variable in the original dataset
    if "`generate'" != "" {
        restore
        qui gen _rank = _n
        qui gen x_obs = ln(_rank - 0.5)
        qui gen y_ln = ln(`varlist')
        qui regress y_ln x_obs
        qui gen x_est = (y_ln - _b[_cons]) / _b[x_obs]
        qui gen `generate' = x_obs - x_est
        di as text "Difference variable '" `generate' "' has been saved to your dataset."
    }
    else {
        restore
    }

end
