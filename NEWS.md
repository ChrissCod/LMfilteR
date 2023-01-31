# News

## 0.1.3.1
* Minor typos and corrections
* Updated maintainer's email, adding author.
## 0.1.3
* Separation of dependent and independent variables into two separate arguments in PF_lm.R and PF_lm_ss.R
* Inclusion of utils.R auxiliary functions with resampling methods following Li, T., Bolic, M., & Djuric, P. M. (2015). Resampling methods for particle filtering: classification, implementation, and strategies. IEEE Signal processing magazine, 32(3), 70-86.
* Separate arguments for dependent variable (Y) and independent variables (Data1) in PF_lm and PF_lm_ss
* Included parameter lbd for the initial priors when initDisPar is not provided, apply for PF_lm and PF_lm_ss
* The returned list of PF_lm and PF_lm_ss includes a summary of the estimated parameters
* Argument initDisPar in PF_lm and PF_lm_ss now only includes the parameters that are going to be estimated. See Details section.
* Included two algorithms that include evolutionary algorithms-based parameters inside the particle filters version for both linear and non-linear (logistic) models: EPF_L_compl.R and EPF_logist_compl.R

* Updated author's email

## 0.1.2
* Function PF_lm.R has changed the output, now it returns a list of three elements
* Changes in the examples section of PF_lm.R documentation
* Better performance achieved with the new function PF_lm_ss.R which uses the method of simple sampling as resampling method 

## 0.1.1
* Changes in the documentation

## 0.1.0 Initial Release
