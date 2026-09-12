{smcl}
{* *! version 1.1.0  14July2026}{...}
{cmd:help zipfbend}
{hline}

.-
help for zipfbend                                                (R Bergs, 14July2026)
.-

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{hi:zipfbend} {hline 2}}Calculates the transition point between an upper Pareto and a lower log-normal tail in rank-size distributions (Rank - 0.5) and Kolmogorov-Smirnov test for city sizes{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmd:zipfbend} {varname} {if} in [{cmd:,} {cmdab:g:enerate(}{it:newvar}{cmd:)}]


{title:Description}

{p 4 4 2}
{cmd:zipfbend} estimates the parameters of a Zipf (Power Law) distribution using 
the modified log rank-size OLS regression: 
ln(Size) = alpha + beta * ln(Rank - 0.5).

{p 4 4 2}
The command automatically drops missing and non-positive values before estimation. 
It identifies a structural "bend point" where the empirical log-rank deviates 
maximally (positive distance) from the OLS-predicted stable log-rank baseline 
(chord between obs 1 and obs n). It also performs a manual Kolmogorov-Smirnov (KS) 
goodness-of-fit test against a theoretical Pareto distribution based on the 
estimated alpha and beta coefficients.


{title:Options}

{p 4 8 2}
{cmd:generate(}{it:newvar}{cmd:)} creates a new variable in the original dataset 
containing the calculated difference between the observed log-rank and the 
OLS-estimated log-rank (x_obs - x_est). 


{title:Results}

{p 4 4 2}
{cmd:zipfbend} displays the following:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:}number of observations{p_end}
{synopt:}estimated Zipf coefficient (b){p_end}
{synopt:}implied Pareto exponent (-1/b){p_end}
{synopt:}maximum positive difference from baseline{p_end}
{synopt:}rank corresponding to the bend point{p_end}
{synopt:}item size at the bend point{p_end}
{synopt:}calculated Kolmogorov-Smirnov test statistic (D){p_end}
{synopt:}approximate critical value for KS test (alpha = 0.05){p_end}


{title:Examples}

{p 4 8 2}. zipfbend city_pop in 1/100{p_end}
{p 4 8 2}. zipfbend city_pop if city_pop>10000 in 10/50, generate(zipf_deviation){p_end}


{title:Author}

{p 4 4 2}R Bergs{p_end}
