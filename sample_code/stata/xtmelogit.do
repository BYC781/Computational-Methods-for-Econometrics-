*Grab the most recent version from the internet
insheet using "https://stats.idre.ucla.edu/stat/data/hdp.csv", comma clear


foreach i of varlist familyhx smokinghx sex cancerstage school {
encode `i', gen(`i'2)
drop `i'
rename `i'2 `i'
}

summarize


tab1 remission cancerstage lengthofstay

cor il6 crp lengthofstay experience

xtmelogit remission il6 crp i.cancerstage lengthofstay experience || did: , intpoints(10)
xtmelogit remission il6 crp i.cancerstage lengthofstay experience || did: lengthofstay, intpoints(10)

/*The variable lists that make up each equation describe how the doctor random effects enter into the model, 
either as random intercepts (constant term) or as random coefficients on regressors in the data. */
