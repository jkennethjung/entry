set more off
clear

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

use ../temp/small.dta, clear
reg l n w r x
sum l
reg n w r x

use ../output/sample.dta, clear
* Aux regressions
forv s = 1/100 {
    preserve
    keep if ns == `s'
    reg l n w r x
    matrix b = e(b)
    sum l
    matrix b = b, r(mean), r(sd)^2 
    reg n w r x
    matrix b = b, e(b)
    restore
    matrix b_aux = nullmat(b_aux) \ b 
}
preserve
svmat b_aux
sum
restore

* Simulated productivities
drop if n == 0
collapse (mean) x w r tfp_mean = tfp (sd) tfp_sd = tfp, by(m t ns)

forv s = 1/100 {
    reg tfp_mean x if ns == `s'
    matrix b = e(b)
    matrix beta = b[1,1]

    reg tfp_mean x w r if ns == `s'
    matrix b = e(b)
    matrix beta = beta, b[1,1] 

    reg tfp_sd x if ns == `s'
    matrix b = e(b)
    matrix beta = beta, b[1,1] 

    reg tfp_sd x w r if ns == `s'
    matrix b = e(b)
    matrix beta = beta, b[1,1] 

    matrix betas = nullmat(betas) \ beta
}
clear
svmat betas
sum
