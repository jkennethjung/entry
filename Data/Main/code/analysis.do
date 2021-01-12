clear
set more off
log using ../output/analysis.log, replace

import delimited using ../temp/cea_fips_xw.csv, clear varnames(1)
ds
keep fips name state v26 
rename v26 cea
unique cea
gen ch_fip = regexm(fips, "[a-zA-Z]")
drop if ch_fip == 1
drop ch_fip
destring fips, replace
unique fips
duplicates tag fips, gen(dup_fips)
tab dup_fips
list if dup_fips > 0
drop if dup_fips > 0
save ../temp/xw.dta, replace

import delimited using ../temp/efsy_cbp_2016.csv, clear varnames(1)
keep if regexm(naics, "^23")
keep if naics == "236220"
gen fips = fipstate*1000 + fipscty
drop fipscty
unique fips
rename emp x_cty

merge 1:1 fips using ../temp/xw.dta, keep(3) nogen
bysort cea fipstate: egen x_state_cea = sum(x_cty)
bysort cea: egen x_max_state = max(x_state_cea)
gen is_max_state = x_state_cea == x_max_state
gen cea_state = fipstate if is_max_state
bysort cea: egen t_state = mean(cea_state) 
assert t_state == cea_state if is_max_state
bysort cea: egen x = sum(x_cty)
drop cea_state x_max_state x_state_cea is_max_state
sum 
desc
save ../temp/demand.dta, replace

use ../temp/infogroup.dta, clear
keep primary_naics_code sic_code* fips_code employee_size_location archive 
rename fips_code fips
rename employ labor 
rename archive year
keep if year == 2016
drop year
destring fips, replace

gen rmc = regexm(primary_naics_code, "^32732")
replace rmc = regexm(sic_code, "^3273") if rmc == 0
forv i = 1/4 {
    replace rmc = regexm(sic_code_`i', "^3273") if rmc == 0
}
count if rmc
keep if rmc

merge m:1 fips using ../temp/demand.dta, keep(1 3)
drop if _merge == 1
drop _merge

gen j = 1
bysort cea: egen N_m = sum(j)
drop j
sum N_m
tab N_m
plot N_m x 
save ../output/data_firms.dta, replace

collapse (mean) N_m (sum) labor, by(cea t_state)
save ../output/data_cea.dta, replace

sum
gen M_t = 1 
collapse (sum) N_t = N_m M_t, by(t_state)
save ../output/data_state.dta, replace

sum

