set more off
clear
log using ../output/analysis.log, replace

import delimited using ../output/small.csv, clear
rename v1 t
rename v2 m
rename v3 n
rename v4 w
rename v5 r
rename v6 x
rename v7 l
rename v8 k
rename v9 pq
rename v10 tfp
rename v11 ns

replace tfp = tfp * 10
save ../output/sample.dta, replace

drop if n == 0
collapse (mean) x w r tfp_mean = tfp (sd) tfp_sd = tfp, by(m t ns)

forv s = 1/100 {
    reg tfp_mean x if ns == `s'
    matrix betas = nullmat(betas) \ e(b)
    reg tfp_sd x if ns == `s'
    matrix alphas = nullmat(alphas) \ e(b)
}
clear
svmat betas
sum
clear
svmat alphas
sum



