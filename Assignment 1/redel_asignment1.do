/*******************************************************************************
ASSIGNMENT 1: PANEL DATA ANALYSIS OF MICROECONOMIC DECISIONS
DANIEL REDEL S.
DATE: Fall 2022
*******************************************************************************/
clear all

*ssc install blindschemes, replace all 
*ssc install asdoc
*ssc install xtabond2
set scheme tab2, permanently

/*******************************************************************************
1. THE EFFECTS OF AGE, COHORT AND TIME
*******************************************************************************/

*-------------*
*--Load Data--*
*-------------*
clear
use "C:\Users\danny\OneDrive\Análisis Cuantitativo y Econometría\Panel Data\Panel Data Analysis - Tilburg\2022 Class\Assignments\Assignment 1\soep.dta"
xtset persnr year
gen cohort = year - age // We generate the cohort variable

*--------------------------------------*
*--Question 1.a. Collinearity Problem--*
*--------------------------------------*
reg s_life year age cohort
reg s_life year cohort
reg s_life age cohort
 
*------------------------------------------------------*
*--Question 1.b.1. Life Satisfaction Plot, by Cohort --*
*------------------------------------------------------*
preserve
collapse (mean) s_life, by(age cohort)

*----- syntax for graph and graph -----
levelsof(cohort), local(gc)

local i = 1
foreach g of local gc {
    local call `call' || line s_life age if cohort == `g', sort // or "connected"
    local leg `leg' label(`i++' "`g'") // syntax for legend
}

twoway `call' legend(`leg') /// graph
    , xtitle("Age") ytitle("Life Satisfaction") 
 
gr export figure1.png, as(png) replace

restore

*-----------------------------------------------------------*
*--Question 1.b.2. Change in Life Satisfaction, by Cohort --*
*-----------------------------------------------------------*
gen d_s_life = d.s_life

preserve
collapse (mean) d_s_life, by(age cohort)

*----- syntax for graph and graph -----
levelsof(cohort), local(gcx)

local i = 1
foreach g of local gcx {
    local callx `callx' || line d_s_life age if cohort == `g', sort
    local legx `legx' label(`i++' "`g'") 
}

twoway `callx' legend(`legx') /// graph
    , xtitle("Age") ytitle("Change in Life Satisfaction") 
 
gr export figure2.png, as(png) replace

restore

*-----------------------------------------------------------*
*--Question 1.b.3. Implied Relationship --*
*-----------------------------------------------------------*

preserve
collapse (mean) s_life year, by(age cohort)
xtset cohort year

**---POOLED OLS---**
reg s_life age, r
est store pols1
**---CONTROLLING FOR COHORT---***
reg s_life age cohort
est store pols2
**---FE ESTIMATOR---***
tab cohort, gen(byD)
reg s_life age byD*
est store d_fe
xtreg s_life age,  fe
est store fe

* Results Table				
esttab pols1 pols2 d_fe fe, keep(age cohort) cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
mtitles("POLS 1" "POLS 2" "DUMMY FE" "FE") nonumbers ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N)

restore

*------------------------------------*
*--Question 1.c. Life Satisfaction --*
*------------------------------------*
preserve

*db collapse
collapse (mean) s_life, by(age)

twoway line s_life age, sort /// graph
    , xtitle("Age") ytitle("Life Satisfaction")

gr export figure3.png, as(png) replace

restore

/*******************************************************************************
2. ONE DRAW OF SIMULATED DATA
*******************************************************************************/
clear all
*------------------------------------*
*--Question 2.a. Explaining the DGP--*
*------------------------------------*
set seed 345398 
drawnorm alpha_i, n(200) 
expand 5 
drawnorm nu_it e_it, n(1000) 
g x_it=nu_it+alpha_i 
drop nu_it 
g y_it=3+alpha_i+2*x_it+e_it 

*-----------------------------------------------------*
*--Question 2.b. Correlation with Unobserved Effects--*
*-----------------------------------------------------*

asdoc pwcorr, sig

reg y_it x_it
est store pols1

*------------------------------------------------------*
*--Question 2.c. OLS with Unobserved Effects Included--*
*------------------------------------------------------*

reg y_it x_it alpha_i
est store pols2

* Results Table				
esttab pols1 pols2 using table2.tex, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N) replace

/*******************************************************************************
3. MANY DRAWS OF SIMULATED DATA
*******************************************************************************/

*-----------------------------------------------------------*
*--Question 3.a. Estimated vs. Unbiased Standard Deviation--*
*-----------------------------------------------------------*
clear all

set seed 345398
capture program drop mcprog
program mcprog
clear
drawnorm alpha_i, n(200)
gen id = _n
expand 5
drawnorm nu_it e_it, n(1000)
g x_it=nu_it+alpha_i
drop nu_it
g y_it=3+alpha_i+2*x_it+e_it
regress y_it x_it
end

simulate _b _se, reps(100): mcprog
sum

estpost sum
esttab, cells("mean sd min max") nomtitle nonumber replace
esttab using table3.tex, cells("mean sd min max") nomtitle nonumber replace

*----------------------------------*
*--Question 3.b. Biased Estimates--*
*----------------------------------*
set scheme tab2, permanently

kdensity _b_x_it, ///
xline(2.498296, lstyle(foreground) lpattern(dash)) xline(2.0, lcolor(red)) ///
xlabel(1.8(0.2)2.6) ///
title("") ///
xtitle("") ///
ytitle("")

gr export figure4.png, as(png) replace
 
*--------------------------------------------------*
*--Question 3.c.1. Clustered Standard Errors in OLS--*
*--------------------------------------------------*
clear all

set seed 345398		
cap prog drop prog_mc								
program prog_mc, rclass
	drawnorm alpha_i, n(200)
	gen id = _n
	expand 5
	drawnorm nu_it e_it, n(1000)
	g x_it=nu_it+alpha_i
	drop nu_it
	g y_it=3+alpha_i+2*x_it+e_it
	
	bysort id: gen time = _n
	xtset id time

	/* POLS regression */
	reg y_it x_it			
	scalar b1_pols = _b[x_it] 
	scalar seb1_pols = _se[x_it]
									
	/* POLS Clustered regression */
	reg y_it x_it, cluster(id)			
	scalar b1_pols_cl = _b[x_it] 
	scalar seb1_pols_cl = _se[x_it]

	/* RE regression */
	xtreg y_it x_it, re // cluster(id)
	scalar b1_re = _b[x_it]	
	scalar seb1_re = _se[x_it]
	
	drop e_it x_it y_it alpha_i time id 
	
end

/* the following line performs 100 simulations */ 
simulate b1_pols b1_pols_cl b1_re seb1_pols seb1_pols_cl seb1_re, reps(100): prog_mc
renvarlab _sim_1-_sim_6, label /* relabels the results stored by stata */ 
sum b1_pols b1_pols_cl b1_re seb1_pols seb1_pols_cl seb1_re /* we take the average over all 100 estimates*/

estpost summarize
esttab, cells("mean sd min max") nomtitle nonumber replace
esttab using table4.tex, cells("mean sd min max") nomtitle nonumber replace


/*******************************************************************************
4. FIXED EFFECTS AND FIRST DIFFERENCES ESTIMATION
*******************************************************************************/
clear all
set seed 345398
drawnorm alpha_i, n(200)
gen id = _n
expand 5

drawnorm nu_it e_it, n(1000)
g x_it=nu_it+alpha_i
drop nu_it
g y_it=3+alpha_i+2*x_it+e_it

bysort id: gen time = _n
xtset id time

*------------------------------------------------------*
*--Question 4.a.1. Fixed Effects vs First Differences--*
*------------------------------------------------------*

* FE regression 
xtreg y_it x_it, fe
est store fe

* FD regression 
reg D.(y_it x_it), nocons
est store fd

* Results Table				
estout fd fe, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N)
	
*----------------------------------------------*
*--Question 4.a.2. Standard Errors: FE vs. FD--*
*----------------------------------------------*

*The FD does not need adjustment
gen x_s = x_it - l.x_it							// differencing
gen y_s = y_it - l.y_it							// differencing

* FD without SE adjustment
reg y_s x_s
est store fd0

* Results Table				
estout fd0 fd, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N)


*The FE needs adjustment, but xtreg already does that for us
bysort id: egen mn_x = mean(x_it)			// mean per individual
gen x_w = x_it - mn_x							// demeaning
bysort id: egen mn_y = mean(y_it)			// mean per individual
gen y_w = y_it - mn_y							// demeaning

* FE without SE adjustment
reg y_w x_w
est store fe0

* Results Table				
estout fe0 fe, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N) 

*--------------------------------*
*--Question 4.b. Random Effects--*
*--------------------------------*

* RE regression 
xtreg y_it x_it, re
est store re

* Results Table	
esttab fd fe re, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N) replace
			
esttab fd fe re using table5.tex, cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(16) modelwidth(12) stats(r2 N) replace


/*******************************************************************************
5. DYNAMIC MODEL
*******************************************************************************/

*-----------------------------------*
*--Question 5.1. Monte Carlo T = 5--*
*-----------------------------------*
clear all
set seed 345398

capture program drop mcprog
program mcprog
	clear
	drawnorm alpha_i, n(200)
	gen id = _n
	expand 5
	
	drawnorm nu_it e_it, n(1000) // n(T*N)
	g x_it=nu_it+alpha_i
	drop nu_it
	gen y_it = 3 + alpha_i + 2*x_it + e_it 
	
	bysort id: gen time = _n
	xtset id time
	
	/* Make the first year == 0 */
	replace y_it = 3 + alpha_i + 0.5*l.y_it + 2*x_it + e_it if time > 1
	
	/* RE regression */
	xtreg y_it l.y_it x_it, fe
end

simulate _b _se, reps(100): mcprog
sum

*------------------------------------*
*--Question 5.1. Monte Carlo T = 10--*
*------------------------------------*
clear all
set seed 345398

capture program drop mcprog
program mcprog
	clear
	drawnorm alpha_i, n(200)
	gen id = _n
	expand 10
	
	drawnorm nu_it e_it, n(2000) // n(T*N)
	g x_it=nu_it+alpha_i
	drop nu_it
	gen y_it = 3 + alpha_i + 2*x_it + e_it 
	
	bysort id: gen time = _n
	xtset id time
	
	/* Make the first year == 0 */
	replace y_it = 3 + alpha_i + 0.5*l.y_it + 2*x_it + e_it if time > 1
	
	/* RE regression */
	xtreg y_it l.y_it x_it, fe
end

simulate _b _se, reps(100): mcprog
sum

*------------------------------------*
*--Question 5.1. Monte Carlo T = 20--*
*------------------------------------*
clear all
set seed 345398

capture program drop mcprog
program mcprog
	clear
	drawnorm alpha_i, n(200)
	gen id = _n
	expand 20
	
	drawnorm nu_it e_it, n(4000) // n(T*N)
	g x_it=nu_it+alpha_i
	drop nu_it
	gen y_it = 3 + alpha_i + 2*x_it + e_it 
	
	bysort id: gen time = _n
	xtset id time
	
	/* Make the first year == 0 */
	replace y_it = 3 + alpha_i + 0.5*l.y_it + 2*x_it + e_it if time > 1
	
	/* RE regression */
	xtreg y_it l.y_it x_it, fe
end

simulate _b _se, reps(100): mcprog
sum

*------------------------------------*
*--Question 5.1. Monte Carlo T = 50--*
*------------------------------------*
clear all
set seed 345398

capture program drop mcprog
program mcprog
	clear
	drawnorm alpha_i, n(200)
	gen id = _n
	expand 50
	
	drawnorm nu_it e_it, n(10000) // n(T*N)
	g x_it=nu_it+alpha_i
	drop nu_it
	gen y_it = 3 + alpha_i + 2*x_it + e_it 
	
	bysort id: gen time = _n
	xtset id time
	
	/* Make the first year == 0 */
	replace y_it = 3 + alpha_i + 0.5*l.y_it + 2*x_it + e_it if time > 1
	
	/* RE regression */
	xtreg y_it l.y_it x_it, fe
end

simulate _b _se, reps(100): mcprog
sum

/*******************************************************************************
6. INSTRUMENTAL VARIABLES ESTIMATION
*******************************************************************************/
*ssc install xtabond2

clear all
set seed 345398

capture program drop mcprog
program mcprog
	clear
	drawnorm alpha_i, n(200)
	gen id = _n
	expand 5
	
	drawnorm nu_it e_it, n(1000) // n(T*N)
	g x_it=nu_it+alpha_i
	drop nu_it
	gen y_it = 3 + alpha_i + 2*x_it + e_it 
	
	bysort id: gen time = _n
	xtset id time
	
	replace y_it = 3 + alpha_i + 0.5*l.y_it + 2*x_it + e_it if time > 1
	
	/* Arellano-Bond */
	xtabond2 y_it l.y_it x_it, gmm(l.y_it) nolevel
end

simulate _b _se, reps(100): mcprog
sum

estpost summarize
esttab using table7.tex, cells("mean sd min max") nomtitle nonumber replace






