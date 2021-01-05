clear
set more off
log using ../output/analysis.log, replace

import delimited using ../temp/efsy_cbp_2016.csv, clear varnames(1)
desc
keep if regexm(naics, "^23")
sum

use ../temp/infogroup.dta, clear
sum
desc

gen rmc = regexm(primary_naics_code, "^3273")
replace rmc = regexm(sic_code, "^3273") if rmc == 0
forv i = 1/4 {
    replace rmc = regexm(sic_code_`i', "^3273") if rmc == 0
}

count if rmc
