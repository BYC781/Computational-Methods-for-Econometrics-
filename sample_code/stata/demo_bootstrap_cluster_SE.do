 use "C:\Users\ssunr\Dropbox\teaching_NTU\Econ7218\STATA_codes\DDK2011.dta"
reg totalscore tracking, cluster(schoolid)
reg totalscore tracking, cluster(schoolid) vce(bootstrap, reps(1000) bca)
estat bootstrap, all
