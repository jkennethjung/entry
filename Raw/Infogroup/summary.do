clear
set more off
log using summary.log, replace

use infogroup.dta, clear
sum
desc, fullnames

log close
