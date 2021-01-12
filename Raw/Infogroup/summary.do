clear
set more off
log using summary.log, replace

use infogroup.dta, clear
sum
desc, fullnames

tab company_holding_status
local vars parent_number company_holding_status subsidiary_number abi site_number
foreach v in `vars' {
    destring `v', replace 
    count if missing(`v')
}
list `vars' if _n < 20
unique abi

log close
