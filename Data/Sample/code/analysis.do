clear
set more off
log using ../output/analysis.log, replace

use ../temp/data_firms.dta, clear
merge m:1 t_state using ../temp/data_state.dta, assert(3) keep(3) ///
    keepusing(N_t M_t) nogen
keep if N_t <= 20 & M_t <= 10
unique t_state
save ../output/data_N20M10.dta, replace

keep if N_t <= 10
unique t_state
save ../output/data_N10M10.dta, replace
sum

log close

