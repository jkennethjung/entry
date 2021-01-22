clear
set more off
log using ../output/analysis.log, replace

use ../temp/data_firms.dta, clear
merge m:1 t_state using ../temp/data_state.dta, assert(3) keep(3) ///
    keepusing(N_t M_t) nogen
drop if labor == 0 & N_m != 0 
keep if N_t <= 20 & M_t <= 10
unique t_state
rename N_m n
rename wage w
rename rent r
rename labor l 
keep t_state cea n w r x l N_t M_t
replace x = log(x)
sum x
replace x = 0.1 * x / r(mean)
sum w
replace w = w / r(mean)
sum r
replace r = r / r(mean)
sum l
replace l = l / r(mean)
save ../temp/raw.dta, replace

collapse (first) n, by(t_state)
gen t = _n
keep t t_state
save ../temp/raw_xw.dta, replace 

use ../temp/raw.dta, clear
merge m:1 t_state using ../temp/raw_xw.dta, assert(3) keep(3) nogen
collapse (first) n t_state, by(t cea)
bysort t: egen n_rank = rank(-n), unique
sort t n_rank
rename n_rank m
list if _n < 20
save ../temp/mt.dta, replace

use ../temp/raw.dta, clear
merge m:1 cea t using ../temp/mt.dta, assert(3) nogen
ds
order t t_state cea m n w r x l
sort t m
sum
save ../output/data.dta, replace
keep t m n w r x l
export delimited using ../output/data.csv, replace novar delim(",")

* Small sample
use ../temp/raw.dta, clear
keep if N_t <= 10
save ../temp/small.dta, replace

collapse (first) n, by(t_state)
gen t = _n
keep t t_state
save ../temp/small_xw.dta, replace 

use ../temp/small.dta, clear
merge m:1 t_state using ../temp/small_xw.dta, assert(3) keep(3) nogen
collapse (first) n t_state, by(t cea)
bysort t: egen n_rank = rank(-n), unique
sort t n_rank
rename n_rank m
save ../temp/small_mt.dta, replace

use ../temp/small.dta, clear
merge m:1 cea t using ../temp/small_mt.dta, assert(3) nogen
ds
order t t_state cea m n w r x l
sort t m
sum
save ../output/small.dta, replace
keep t m n w r x l
export delimited using ../output/small.csv, replace novar delim(",")

log close
