/*******************************************************************************
ASSIGNMENT 2: PANEL DATA ANALYSIS OF MICROECONOMIC DECISIONS
DANIEL REDEL S.
DATE: Fall 2022
*******************************************************************************/
clear all

*ssc install blindschemes, replace all 
*ssc install asdoc
*ssc install xtabond2
set scheme tab2, permanently

/*******************************************************************************
1. BINARY CHOICE MODELS
*******************************************************************************/

*-------------*
*--Load Data--*
*-------------*
clear
use "C:\Users\danny\OneDrive\Análisis Cuantitativo y Econometría\Panel Data\Panel Data Analysis - Tilburg\2022 Class\Assignments\Assignment 2\data_ass2_2022.dta"

keep if WAVE <= 8 // My SNR = 2102630 ---> Waves 1 to 8
xtset ID WAVE

*----------------------------------------------------*
*--Question 1.1. Static Random Effects Logit Model --*
*----------------------------------------------------*
gen AGE2 = AGE^2
global covariates AGE AGE2 MARRIED YOUNG_CH NUM_CH EDU_13_15 EDU_15 WHITE

xtlogit EMPL $covariates, re
est store re_logit

*--------------------------------------*
*--Question 1.2. Pooled Logit Model  --*
*--------------------------------------*
logit EMPL $covariates
est store pooled_logit

*----------------------------------------------*
*--Question 1.2.a. Pooled vs Random Effects  --*
*----------------------------------------------*
* See Document

*-----------------------------------*
*--Question 1.3. Marginal Effects --*
*-----------------------------------*
dis (-0.131+2*(0.001)*40)*(0.5)*[1-0.5]

*--------------------------------------------------*
*--Question 1.4. Static Fixed Effects Logit Model--*
*--------------------------------------------------*
xtlogit EMPL $covariates, fe
est store fe_logit
ereturn list

*----------------------------------*
*--Question 1.4.a. HAUSMAN Test  --*
*----------------------------------*
xtlogit EMPL AGE AGE2 YOUNG_CH NUM_CH EDU_13_15 EDU_15, re
est store re_haus
xtlogit EMPL AGE AGE2 YOUNG_CH NUM_CH EDU_13_15 EDU_15, fe
est store fe_haus

hausman fe_haus re_haus

*-------------------------------------------------*
*--Question 1.5. Quasi-Fixed Effects Logit Model--*
*-------------------------------------------------*
** Mundlak Version
egen MEAN_EDU_13_15 = mean(EDU_13_15), by(ID)
egen MEAN_EDU_15 = mean(EDU_15), by(ID)

xtlogit EMPL $covariates MEAN_EDU_13_15 MEAN_EDU_15, re
est store qfe_logit

*------------------------------------------------------------------*
*--Question 1.6. Dynamic Random Effects Logit Model (Wooldridge) --*
*------------------------------------------------------------------*
* Lagged Value *
gen lagged_EMPL=.
replace lagged_EMPL=l.EMPL if WAVE>1

* Initial Value *
gen EMPL_0 = EMPL
by ID, sort: replace EMPL_0 = EMPL[1]

* Model *
xtlogit EMPL $covariates lagged_EMPL EMPL_0 if WAVE > 1, re
* xtlogit EMPL $covariates lagged_EMPL EMPL_0, re // should be the same
est store dynamic_re_logit

*-------------------------------------*
*--Question 1.6.a. State Dependence --*
*-------------------------------------*
* See Document

*---------------------------------------------*
*--Question 1.6.b. Unobserved Heterogeneity --*
*---------------------------------------------*
sum EMPL if WAVE==1
sum EMPL_0 // V(y) == 0.493

dis (2.534)^2*(0.493)*(1-0.493)+(1.629)^2 // New Sigma2_alpha
dis (4.2586154/(4.2586154+(3.14159^2/3))) // New Rho

*--------------------------------------------------*
*--Question 1. FINAL TABLE: BINARY CHOICE MODELS --*
*--------------------------------------------------*	
esttab pooled_logit re_logit fe_logit qfe_logit dynamic_re_logit ///
	using table1.tex, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f)))  ///
	mtitles("Pooled Logit" "RE Logit" "FE Logit" "QFE Logit" "Dynamic") ///
	nonumbers starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) ///
	stats(N ll sigma_u rho, fmt(a4 4 a4) ///
	labels("Observations" "Log-Likelihood" "sigma_u" "rho")) replace

/*******************************************************************************
2. TOBIT MODEL
*******************************************************************************/
est clear

**See Censored
gen censored = 0 
replace censored = 1 if INCOME==0
replace censored = 0 if INCOME == .
tab censored // 46.77 are censored
tab EMPL // 45.39 unemployed

*----------------------------------------------------*
*--Question 2.1. Static Random Effects Tobit Model --*
*----------------------------------------------------*
global covariates_tobit SEX AGE AGE2 WHITE MARRIED EDU_13_15 EDU_15
xttobit INCOME $covariates_tobit, ll(0)
est store tobit1

*-------------------------------------------*
*--Question 2.2. Unobserved Heterogeneity --*
*-------------------------------------------*
* See Document 

*------------------------------------------------------------------*
*--Question 2.3. Static Random Effects Tobit Model: Year Dummies --*
*------------------------------------------------------------------*
xttobit INCOME $covariates_tobit i.YEAR, ll(0)
est store tobit2

*Dummy for each Year
forvalues i = 1(1)8{
 generate year`i' = 0
 replace year`i' = 1 if YEAR == 1992+2*`i'
 }

* RE Tobit
xttobit INCOME $covariates_tobit ///
year2 year3 year4 year5 year6 year7 year8, ll(0)
est store tobit2

*Test whether they are jointly significant
test year2 year3 year4 year5 year6 year7 year8
/* Time dummies are jointly significantly different from zero (at the 5%)*/

*--------------------------------------------------------------------*
*--Question 2.4. Static Random Effects Tobit Model: Marital Status --*
*--------------------------------------------------------------------*
gen MARRIED_SEX = MARRIED*SEX
xttobit INCOME $covariates_tobit MARRIED_SEX, ll(0)
est store tobit3

*-----------------------------------------------------*
*--Question 2.5. Dynamic Random Effects Tobit Model --*
*-----------------------------------------------------*
* Lagged Value *
gen lagged_INCOME=.
replace lagged_INCOME=l.INCOME if WAVE>1

* Initival Value *
gen INCOME_0 = INCOME
by ID, sort: replace INCOME_0 = INCOME[1]

* Model *
xttobit INCOME $covariates_tobit lagged_INCOME INCOME_0 if WAVE > 1, ll(0)
est store tobit4
return list
matrix list r(table)
matrix list e(b)

*---------------------------------------------*
*--Question 2.5.1. Unobserved Heterogeneity --*
*---------------------------------------------*
sum INCOME_0 // V(y_0) = 18522.89

* New Sigma Alpha
dis _b[INCOME_0]^2*(18522.89)^2 + _b[/: sigma_u]^2
* New RHO:
dis "RHO: "(_b[INCOME_0]^2*(18522.89)^2 /// 
	+ _b[/: sigma_u]^2)/(_b[INCOME_0]^2*(18522.89)^2 ///
	+ _b[/: sigma_u]^2 + _b[/:sigma_e]^2) 

*------------------------------------------*
*--Question 2. FINAL TABLE: TOBIT MODELS --*
*------------------------------------------*
esttab tobit1 tobit2 tobit3 tobit4 using table2.tex ///
	, cells(b(star) se(par fmt(%9.3f)))  ///
	mtitles("RE Tobit 1" "RE Tobit 2" "RE Logit 3" "Dynamic T") nonumbers ///
	starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) ///
	stats(N ll sigma_u rho, fmt(a4 4 a4) ///
	labels("Observations" "Log-Likelihood" "sigma_u" "rho")) replace

*-----------------------------------------------------------*
*--Question 2.6. Individual Effects against Initial Value --*
*-----------------------------------------------------------*
test INCOME_0
/* Time dummies are jointly significantly different from zero (at the 5%)*/

/*******************************************************************************
2. ORDERED RESPONSE MODELS
*******************************************************************************/

*-------------------------------------------------------*
*--Question 3.1. Static Random Effects Ordered Logit  --*
*-------------------------------------------------------*
xtologit JOB_SAT AGE INCOME EDU_13_15 EDU_15 if EMPL==1
est store re_ologit

*--------------------------------*
*--Question 3.2. Constant Term --*
*--------------------------------*
* See Document

*---------------------------------------------------*
*--Question 3.3.  Likelihood Ratio Test in Output --*
*---------------------------------------------------*
* See Document

*-----------------------------------------------------*
*--Question 3.4. Static Fixed Effects Ordered Logit --*
*-----------------------------------------------------*
feologit JOB_SAT AGE INCOME EDU_13_15 EDU_15 if EMPL==1
est store fe_ologit

*----------------------------------------------------*
*--Question 3.5. Hausmann Test RE vs FE: Education --*
*----------------------------------------------------*
dis (-0.272--0.345)^2/(0.170^2-0.049^2) // 0.201

*------------------------------------------------*
*--Question 3.6.1. RE ORDERED RESPONSE MODEL 2 --*
*------------------------------------------------*
xtologit JOB_SAT AGE AGE2 INCOME SEX WHITE MARRIED EDU_13_15 EDU_15 if EMPL==1
est store re_ologit1

*-----------------------------------------------*
*--Question 3.6.2. QFE ORDERED RESPONSE MODEL --*
*-----------------------------------------------*
*** Quasi-Fixed Effects: Mundlak Version
egen MEAN_INCOME = mean(INCOME), by(ID)

xtologit JOB_SAT AGE AGE2 INCOME SEX WHITE MARRIED ///
	EDU_13_15 EDU_15 MEAN_INCOME if EMPL==1
est store qfe_ologit

*-------------------------------------*
*--Question 3.6.3. MODEL ASSESSMENT --*
*-------------------------------------*
 ** RE1 Model vs RE2 Model, 4 degrees-of-freedom
 dis 2*(-41024.1569--41047.4076)
 
** RE2 Model vs QFE Model, 1 degree-of-freedom
 dis 2*(-41013.9830--41047.4076)
 
** FE Model vs RE Model
dis (-0.013--0.012)^2/(0.003^2-0.002^2)
dis (-0.021--0.071)^2/(0.120^2-0.043^2)
dis (-0.272--0.345)^2/(0.170^2-0.049^2)
** Income:
dis (-0.00000189--0.00000315)^2/(0.00000062^2-0.00000043^2)

** FE Model vs QFE Model
** Income :
dis (-0.00000189--0.00000198)^2/(0.00000062^2-0.00000051^2) // == 0.065

**format fmt
esttab fe_ologit re_ologit qfe_ologit ///
	, cells(b(star fmt(%9.8f)) se(par fmt(%9.8f))) ///
	mtitles("FE Ologit" "RE Ologit") nonumbers ///
	starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) ///
	modelwidth(12) stats(N ll sigma_u, fmt(a4 4 a4) ///
	labels("Observations" "Log-Likelihood" "sigma_u"))

*-------------------------------------------------*
*--Question 3.7.  Marginal Effects of Education --*
*-------------------------------------------------*
*Using QFE Ologit
dis -(0.5)*(1-0.5)*(-0.240)
dis  (0.5)*(1-0.5)*(-0.240)

*-----------------------------------------------------*
*--Question 3. FINAL TABLE: ORDERED RESPONSE MODELS --*
*-----------------------------------------------------*
esttab re_ologit fe_ologit re_ologit1 qfe_ologit using table3.tex///
	, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f)))  ///
	mtitles("RE Ologit 1" "FE Ologit 1" "RE Ologit 2" "QFE Ologit") nonumbers ///
	starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) ///
	stats(N ll sigma_u, fmt(a4 4 a4) ///
	labels("Observations" "Log-Likelihood" "sigma_u")) replace



