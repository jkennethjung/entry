clear
set more off
log using ../output/analysis.log, replace

import delimited using ../temp/allhlcn16.csv, clear varnames(1)
desc, fullnames
keep if areatype == "County" 
keep if naics == 1013 | naics == 1012 
destring areacode, replace
rename areacode fips
rename annualaverageweeklywage wage
keep fips wage 
destring wage, replace ignore(",")
collapse (mean) wage, by(fips)
sum
count if wage == 0
save ../temp/wages.dta, replace

import delimited using ../temp/FY2016F_50_RevFinal.csv, clear varnames(1)
desc, fullnames
gen fips = state*1000 + county
rename rent50_0 rent 
rename pop2010 pop
keep fips rent pop
sum
collapse (mean) rent, by(fips)

merge 1:1 fips using ../temp/wages.dta
keep if _merge == 3
drop _merge
save ../output/costs.dta, replace

log close

