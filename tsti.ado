cap program drop tsti
program define tsti, rclass
	*Syntax: estimate se rope_lb rope_ub, options. df has a second companion 'ghost' option designed to assess whether it has been provied
	
	*Get tokens
	gettoken estimate 0 : 0, parse(" ,")
	gettoken se 0 : 0, parse(" ,")
	gettoken rope_lb 0 : 0, parse(" ,")
	gettoken rope_ub 0 : 0, parse(" ,")
	
	*Confirm that all entries are numbers
	confirm number `estimate'
	confirm number `se'
	confirm number `rope_lb'
	confirm number `rope_ub'
	
	*If se <= 0...
	if (`se' <= 0) {
		
		*... then stop the function
		display "'se' must be strictly greater than zero"
		exit
		
	}
	
	*If rope_lb >= rope_ub...
	if (`rope_lb' >= `rope_ub') {
		
		*... then stop the function
		display "'rope_lb' must be strictly less than 'rope_ub'"
		exit
		
	}
	
	*Check if moremata is installed
	capture which mf_mm_root.hlp
	
	*If moremata is not installed...
	if (_rc == 111) {
		
		*... then stop the function
		display "Command 'tst' requires command mm_root. To install, type 'ssc install moremata'"
		exit
		
	}
	
	syntax [, df(real 10000000000000000000) alpha(real 0.05) power(real 0.8) df2(real 10000000000000000000)]
	
	*************************
	***** ERRORS & PREP *****
	*************************
	
	*If alpha is not between 0 and 0.5...
	if (`alpha' <= 0 | `alpha' >= 0.5) {
		
		*... then stop the function
		display "'alpha' must be between 0 and 0.5"
		exit
		
	}
	
	*If power is not between 0.5 and 1...
	if (`power' <= 0.5 | `power' >= 1) {
		
		*... then stop the function
		display "'power' must be between 0.5 and 1"
		exit
		
	}
	
	*If the estimate is exactly midway between the lower and upper bounds of the ROPE...
	if (`estimate' == (`rope_lb' + `rope_ub')/2) {
		
		*Then select the upper bound as the relevant TOST bound
		local bound = `rope_ub'
		
	}
	
	*Otherwise...
	if (`estimate' != (`rope_lb' + `rope_ub')/2) {
		
		*If the upper bound is closer to the estimate than the lower bound...
		if (abs(`estimate' - `rope_ub') < abs(`estimate' - `rope_lb')) {
			
			*Then select the upper bound as the relevant TOST bound
			local bound = `rope_ub'
			
		}
		
		*Otherwise...
		if (abs(`estimate' - `rope_ub') > abs(`estimate' - `rope_lb')) {
			
			*Select the lower bound as the relevant TOST bound
			local bound = `rope_lb'
			
		}
		
	}
	
	*Generate confidence percentage
	local confidence_pct = round(1 - `alpha', .01)*100
	
	*Generate power percentage
	local power_pct = round(`power', .01)*100
	
	*Generate a test matrix
	matrix test_mat = J(3, 3, .)
	
	*If the estimate is located above the ROPE...
	if (`estimate' > `rope_ub') {
		
		*Then designate the test for bounding above the ROPE as relevant and the other two tests as irrelevant
		mat test_mat[1, 3] = 1
		mat test_mat[2, 3] = 0
		mat test_mat[3, 3] = 0
		
	}
	
	*If the estimate is located below the ROPE...
	if (`estimate' < `rope_lb') {
		
		*Then designate the test for bounding below the ROPE as relevant and the other two tests as irrelevant
		mat test_mat[1, 3] = 0
		mat test_mat[2, 3] = 0
		mat test_mat[3, 3] = 1
		
	}
	
	*If the estimate is located inside the ROPE...
	if (`estimate' <= `rope_ub' & `estimate' >= `rope_lb') {
		
		*Then designate the TOST p-value as relevant and the other two tests as irrelevant
		mat test_mat[1, 3] = 0
		mat test_mat[2, 3] = 1
		mat test_mat[3, 3] = 0
		
	}
	
	*Clear mata
	clear mata
	
	************************
	***** SUB-ROUTINES *****
	************************
	
	*If df is not provided...
	if (`df' == `df2') {
		
		*Generate the bounds of the ECI
		local ECI_LB = `estimate' - invnormal(1 - `alpha')*`se'
		local ECI_UB = `estimate' + invnormal(1 - `alpha')*`se'
		
		*Enter mata to obtain number of standard errors necessary to produce ROSE
		mata: function power_func(epsilon) return(normal(epsilon - invnormal(1 - `alpha')) + normal(-epsilon - invnormal(1 - `alpha')))
		mata: function power_diff(epsilon) return(power_func(epsilon) - `power')
		mata: mm_root(epsilon = ., &power_diff(), 0, invnormal(`power') + invnormal(1 - `alpha'))
		mata: st_numscalar("ROSE_SE", epsilon)
		mata: st_local("ROSE_SE", strofreal(epsilon))
		
		*Generate the bounds of the ROSE
		local ROSE_LB = `estimate' - `ROSE_SE'*`se'
		local ROSE_UB = `estimate' + `ROSE_SE'*`se'
		
		*Store the z-statistic and p-value of the one-sided test for bounding above the ROPE
		mat test_mat[1, 1] = (`estimate' - `rope_ub')/`se'
		mat test_mat[1, 2] = 1 - normal(test_mat[1, 1])
		
		*If the lower bound of the ROPE is the relevant TOST bound...
		if (`bound' == `rope_lb') {
			
			*Then store the z-statistic as estimate - min(ROPE) in standard error units...
			mat test_mat[2, 1] = (`estimate' - `rope_lb')/`se'
			*... and store the p-value of the one-sided test in the upper tail
			mat test_mat[2, 2] = 1 - normal(test_mat[2, 1])
			
		}
		
		*If the upper bound of the ROPE is the relevant TOST bound...
		if (`bound' == `rope_ub') {
			
			*Then store the z-statistic as estimate - max(ROPE) in standard error units...
			mat test_mat[2, 1] = (`estimate' - `rope_ub')/`se'
			*... and store the p-value of the one-sided test in the lower tail
			mat test_mat[2, 2] = normal(test_mat[2, 1])
			
		}
		
		*Store the z-statistic and p-value of the one-sided test for bounding below the ROPE
		mat test_mat[3, 1] = (`estimate' - `rope_lb')/`se'
		mat test_mat[3, 2] = normal(test_mat[3, 1])
		
		*If no p-value is < alpha...
		if (test_mat[1, 2] >= `alpha' & test_mat[2, 2] >= `alpha' & test_mat[3, 2] >= `alpha') {
			
			*Then store the conclusion
			local conclusion "The significance of the estimate is inconclusive."
			
		}
		
		*If the p-value of the one-sided test for bounding above the ROPE < alpha...
		if (test_mat[1, 2] < `alpha') {
			
			*Then store the conclusion
			local conclusion "The estimate is significantly bounded above the ROPE."
			
		}
		
		*If the p-value of the TOST procedure < alpha...
		if (test_mat[2, 2] < `alpha') {
			
			*Then store the conclusion
			local conclusion "The estimate is significantly bounded within the ROPE."
			
		}
		
		*If the p-value of the one-sided test for bounding below the ROPE < alpha...
		if (test_mat[3, 2] < `alpha') {
			
			*Then store the conclusion
			local conclusion "The estimate is significantly bounded below the ROPE."
			
		}
		
		********************
		*** BOUNDS TABLE ***
		********************
		
		disp ""
		disp in smcl in gr "{ralign 59: Approximate bounds}" 																_col(59) " {c |}" 	_col(71) in gr "Lower bound"  		 _col(94) in gr "Upper bound"
		disp in smcl in gr "{hline 60}{c +}{hline 52}"
		disp in smcl in gr "{ralign 59:Region of practical equivalence (ROPE)}"        										_col(59) " {c |} " 	_col(71) as result %9.3f `rope_lb'   _col(94) %9.3f  `rope_ub'
		disp in smcl in gr "{ralign 59:`confidence_pct'% equivalence confidence interval (ECI)}"    						_col(59) " {c |} " 	_col(71) as result %9.3f `ECI_LB'    _col(94) %9.3f  `ECI_UB'      
		disp in smcl in gr "{ralign 59:`confidence_pct'% region of statistical equivalence (ROSE) with `power_pct'% power}" _col(59) " {c |} " 	_col(71) as result %9.3f `ROSE_LB'   _col(94) %9.3f  `ROSE_UB'      
		
		*********************
		*** RESULTS TABLE ***
		*********************
		
		disp ""
		disp in smcl in gr "{ralign 46: Testing results}" 							   _col(47) " {c |} " _col(52) in gr "z-statistic"			  _col(67) in gr "p-value"	    _col(80) in gr "Relevant"
		disp in smcl in gr "{hline 47}{c +}{hline 40}"
		disp in smcl in gr "{ralign 46:Test: Estimate bounded above ROPE (one-sided)}" _col(47) " {c |} " _col(52) as result %9.3f test_mat[1, 1] _col(64) %9.3f test_mat[1, 2]	_col(76) %9.0f  test_mat[1, 3]
		disp in smcl in gr "{ralign 46:Test: Estimate bounded within ROPE (TOST)}"     _col(47) " {c |} " _col(52) as result %9.3f test_mat[2, 1] _col(64) %9.3f test_mat[2, 2]	_col(76) %9.0f  test_mat[2, 3]
		disp in smcl in gr "{ralign 46:Test: Estimate bounded below ROPE (one-sided)}" _col(47) " {c |} " _col(52) as result %9.3f test_mat[3, 1] _col(64) %9.3f test_mat[3, 2]	_col(76) %9.0f  test_mat[3, 3]
		
		*************************
		*** PRINT DISCLAIMERS ***
		*************************
		disp ""
		disp "`conclusion'"
		disp ""
		disp "Asymptotically approximate equivalence confidence intervals (ECIs) and three-sided testing (TST) results reported"
		disp "If using for academic/research purposes, please cite the paper underlying this program:"
		disp "Fitzgerald, Jack (2024). The Need for Equivalence Testing in Economics. Institute for Replication Discussion Paper Series No. 125. https://www.econstor.eu/handle/10419/296190."
		
	}
	
	*If df is provided...
	if (`df' != `df2') {
		
		*If df is not greater than zero...
		if (`df' <= 0) {
			
			*... then stop the function
			display "If 'df' is specified, then it must be greater than zero"
			exit
			
		}
		
		*Generate the bounds of the ECI
		local ECI_LB = `estimate' - invt(`df', 1 - `alpha')*`se'
		local ECI_UB = `estimate' + invt(`df', 1 - `alpha')*`se'
		
		*Enter mata to obtain number of standard errors necessary to produce ROSE
		mata: function power_func(epsilon) return(nt(`df', invt(`df', 1 - `alpha'), epsilon) + nt(`df', invt(`df', 1 - `alpha'), -epsilon))
		mata: function power_diff(epsilon) return(power_func(epsilon) - `power')
		mata: mm_root(epsilon = ., &power_diff(), 0, invnt(`df', invt(`df', 1 - `alpha'), `power'))
		mata: st_numscalar("ROSE_SE", epsilon)
		mata: st_local("ROSE_SE", strofreal(epsilon))
		
		*Generate the bounds of the ROSE
		local ROSE_LB = `estimate' - `ROSE_SE'*`se'
		local ROSE_UB = `estimate' + `ROSE_SE'*`se'
		
		*Store the t-statistic and p-value of the one-sided test for bounding above the ROPE
		mat test_mat[1, 1] = (`estimate' - `rope_ub')/`se'
		mat test_mat[1, 2] = 1 - t(`df', test_mat[1, 1])
		
		*If the lower bound of the ROPE is the relevant TOST bound...
		if (`bound' == `rope_lb') {
			
			*Then store the t-statistic as estimate - min(ROPE) in standard error units...
			mat test_mat[2, 1] = (`estimate' - `rope_lb')/`se'
			*... and store the p-value of the one-sided test in the upper tail
			mat test_mat[2, 2] = 1 - t(`df', test_mat[2, 1])
			
		}
		
		*If the upper bound of the ROPE is the relevant TOST bound...
		if (`bound' == `rope_ub') {
			
			*Then store the t-statistic as estimate - max(ROPE) in standard error units...
			mat test_mat[2, 1] = (`estimate' - `rope_ub')/`se'
			*... and store the p-value of the one-sided test in the upper tail
			mat test_mat[2, 2] = t(`df', test_mat[2, 1])
			
		}
		
		*Store the t-statistic and p-value of the one-sided test for bounding below the ROPE
		mat test_mat[3, 1] = (`estimate' - `rope_lb')/`se'
		mat test_mat[3, 2] = t(`df', test_mat[3, 1])
		
		*If no p-value is < alpha...
		if (test_mat[1, 2] >= `alpha' & test_mat[2, 2] >= `alpha' & test_mat[3, 2] >= `alpha') {
			
			*Then store the conclusion
			local conclusion "The significance of the estimate is inconclusive."
			
		}
		
		*If the p-value of the one-sided test for bounding above the ROPE < alpha...
		if (test_mat[1, 2] < `alpha') {
			
			*Then store the conclusion
			local conclusion "The estimate is significantly bounded above the ROPE."
			
		}
		
		*If the p-value of the TOST procedure < alpha...
		if (test_mat[2, 2] < `alpha') {
			
			*Then store the conclusion
			local conclusion "The estimate is significantly bounded within the ROPE."
			
		}
		
		*If the p-value of the one-sided test for bounding below the ROPE < alpha...
		if (test_mat[3, 2] < `alpha') {
			
			*Then store the conclusion
			local conclusion "The estimate is significantly bounded below the ROPE."
			
		}
		
		********************
		*** BOUNDS TABLE ***
		********************
		
		disp ""
		disp in smcl in gr "{ralign 59: Exact bounds}" 																		_col(59) " {c |}" 	_col(71) in gr "Lower bound"  		 _col(94) in gr "Upper bound"
		disp in smcl in gr "{hline 60}{c +}{hline 52}"
		disp in smcl in gr "{ralign 59:Region of practical equivalence (ROPE)}"        										_col(59) " {c |} " 	_col(71) as result %9.3f `rope_lb'   _col(94) %9.3f  `rope_ub'
		disp in smcl in gr "{ralign 59:`confidence_pct'% equivalence confidence interval (ECI)}"    						_col(59) " {c |} " 	_col(71) as result %9.3f `ECI_LB'    _col(94) %9.3f  `ECI_UB'      
		disp in smcl in gr "{ralign 59:`confidence_pct'% region of statistical equivalence (ROSE) with `power_pct'% power}" _col(59) " {c |} " 	_col(71) as result %9.3f `ROSE_LB'   _col(94) %9.3f  `ROSE_UB'     
		
		*********************
		*** RESULTS TABLE ***
		*********************
		
		disp ""
		disp in smcl in gr "{ralign 46: Testing results}" 							   _col(47) " {c |} " _col(52) in gr "z-statistic"			  _col(67) in gr "p-value"	    _col(80) in gr "Relevant"
		disp in smcl in gr "{hline 47}{c +}{hline 40}"
		disp in smcl in gr "{ralign 46:Test: Estimate bounded above ROPE (one-sided)}" _col(47) " {c |} " _col(52) as result %9.3f test_mat[1, 1] _col(64) %9.3f test_mat[1, 2]	_col(76) %9.0f  test_mat[1, 3]
		disp in smcl in gr "{ralign 46:Test: Estimate bounded within ROPE (TOST)}"     _col(47) " {c |} " _col(52) as result %9.3f test_mat[2, 1] _col(64) %9.3f test_mat[2, 2]	_col(76) %9.0f  test_mat[2, 3]
		disp in smcl in gr "{ralign 46:Test: Estimate bounded below ROPE (one-sided)}" _col(47) " {c |} " _col(52) as result %9.3f test_mat[3, 1] _col(64) %9.3f test_mat[3, 2]	_col(76) %9.0f  test_mat[3, 3]
		
		*************************
		*** PRINT DISCLAIMERS ***
		*************************
		disp ""
		disp "`conclusion'"
		disp ""
		disp "Exact equivalence confidence intervals (ECIs) and three-sided testing (TST) results reported"
		disp "If using for academic/research purposes, please cite the paper underlying this program:"
		disp "Fitzgerald, Jack (2024). The Need for Equivalence Testing in Economics. Institute for Replication Discussion Paper Series No. 125. https://www.econstor.eu/handle/10419/296190."
		
	}
	
	*Return results
	return local estimate = `estimate'
	return local se = `se'
	return local ROPE_LB = `rope_lb'
	return local ROPE_UB = `rope_ub'
	return local alpha = `alpha'
	return local power = `power'
	return local ECI_LB = `ECI_LB'
	return local ECI_UB = `ECI_UB'
	return local ROSE_LB = `ROSE_LB'
	return local ROSE_UB = `ROSE_UB'
	return local ts_above = test_mat[1, 1]
	return local p_above = test_mat[1, 2]
	return local relevant_above = test_mat[1, 3]
	return local ts_TOST = test_mat[2, 1]
	return local p_TOST = test_mat[2, 2]
	return local relevant_TOST = test_mat[2, 3]
	return local ts_below = test_mat[3, 1]
	return local p_below = test_mat[3, 2]
	return local relevant_below = test_mat[3, 3]
	return local conclusion `conclusion'
	
	end
	
