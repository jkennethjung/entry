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
desc
save ../temp/xw.dta, replace

import delimited using ../temp/efsy_cbp_2016.csv, clear varnames(1)
desc
keep if regexm(naics, "^23")
keep if naics == "236220"
gen fips = fipstate*1000 + fipscty
drop fipstate fipscty
unique fips
sum
save ../temp/demand.dta, replace

use ../temp/infogroup.dta, clear
desc
keep primary_naics_code sic_code* fips_code employee_size_location archive 
rename fips_code fips
rename employ labor 
rename archive year
sum
desc
keep if year == 2016
drop year
destring fips, replace

gen rmc = regexm(primary_naics_code, "^3273")
replace rmc = regexm(sic_code, "^3273") if rmc == 0
forv i = 1/4 {
    replace rmc = regexm(sic_code_`i', "^3273") if rmc == 0
}
count if rmc
keep if rmc

merge m:1 fips using ../temp/xw.dta, assert(2 3) keep(3) nogen
merge m:1 fips using ../temp/demand.dta, keep(1 3)
drop if _merge == 1
drop _merge
save ../output/data.dta, replace
sum 
desc

