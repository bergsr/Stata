program define zipfbend, rclass
    version 12.0
    syntax varname(numeric) [if] in [, Generate(string)]
    
    preserve 
    
    * Exclude missing/non-positive values and sort descending globally
    qui drop if missing(`varlist')
    qui drop if `varlist' <= 0
    qui sort `varlist'
    qui gen long _orig_obs = _n
    qui gsort -`varlist'
    
    * Generate global rank and log variables for the ENTIRE distribution
    qui gen long _rank = _n
    qui gen double y = ln(`varlist')
    qui gen double x_obs = ln(_rank - 0.5)

    * 1. Regression based ONLY on the absolute first and last obs of TOTAL distribution
    qui count
    local global_last = r(N)
    tempvar global_base
    qui gen byte `global_base' = (_rank == 1 | _rank == `global_last')
    qui regress y x_obs if `global_base' == 1
    scalar b_alpha = _b[_cons]
    scalar b_beta = _b[x_obs]  // Stable chord slope preserved here!
    
         * NOW apply user's if/in restrictions or sample truncation safely
    keep `if' `in'

    * 2. Calculate estimated x-values based on the STABLE regression parameters
    qui gen double x_est = (y - b_alpha) / b_beta 
    
    * 3. Calculate difference 
    qui gen double diff_x = x_obs - x_est 
    
    * Find maximum positive difference and corresponding row identifiers
    qui summarize diff_x
    local max_diff = r(max)
    
    * Find the rank, size, and original ID corresponding to this maximum difference
    qui sort diff_x
    local k_rank = _rank[_N]
    local k_size = `varlist'[_N]
    local k_id   = _orig_obs[_N]
	
	* 4. Full Sample Regression (Global)
    qui regress y x_obs
    scalar s_alpha = _b[_cons]
    scalar s_beta = _b[x_obs]
    
    * 5. Prepare Kolmogorov-Smirnov Test (remains based on Full Sample OLS)
    qui sort `varlist'
    qui gen double cdf_obs = _n / _N
    
    * Theoretical Pareto CDF based on full sample OLS estimators
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
    di as text "Maximum positive difference from baseline:" as result %9.4f `max_diff'
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
    
    * If requested, handle the generation of variables before restoring
    if "`generate'" != "" {
        * Change the name of your calculated difference to the user-specified name
        capture drop `generate'
        qui gen double `generate' = diff_x
        
        * Change preserve behavior to change the data tracking on exit
        restore, not
        di as text "Residuals variable '" as result "`generate'" as text "' has been saved to your dataset."
    }
    else {
        * If generate wasn't used, pull the original dataset back unmodified
        restore
    }
    
    * 6. Save in r()-Vektor
    return scalar alpha_full = s_alpha
    return scalar beta_full  = s_beta
    return scalar pareto_exp = `p_exp'
    
    return scalar alpha_base = b_alpha
    return scalar beta_base  = b_beta
    
    return scalar bend_rank  = `k_rank'
    return scalar bend_size  = `k_size'
    return scalar bend_id    = `k_id'
    return scalar max_diff   = `max_diff'
    
    return scalar ks_stat    = `ks_stat'
    return scalar N          = `global_last'

end
