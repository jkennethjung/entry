set more off
clear
log using ../output/analysis.log, replace

import delimited using ../temp/data.csv, clear
rename v1 t
rename v2 m 
rename v3 n
rename v4 w
rename v5 r
rename v6 x
rename v7 l
rename v8 k
rename v9 pq 
save ../temp/raw.dta, replace

collapse (first) n, by(t m)
bysort t: egen n_rank = rank(-n), unique 
sort t n_rank
list
save ../temp/mt.dta, replace

use ../temp/raw.dta, clear
merge m:1 m t using ../temp/mt.dta, assert(3) nogen
drop m 
order t n_rank n w r x l k pq
sort t n_rank
export delimited using ../output/data.csv, replace novar delim(",")

log close
