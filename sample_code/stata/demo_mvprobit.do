ssc install mvprobit
clear
set seed 12309
set obs 3000
matrix R=(1, .25, .5, .75 \ .25, 1, .75, .5 \ .5, .75, 1, .75 \ .75, .5, .75, 1)
drawnorm u1 u2 u3 u4, corr(R)
correlate u*
generate x1 = uniform()-.5
generate x2 = uniform() + 1/3
generate x3 = 2*uniform() + .5
generate x4 = .5*uniform() - 1/3


ge y1s = .5 + 4*x1 + u1
ge y2s = 3 + .5*x1 - 3*x2 + u2
ge y3s = 1 - 2*x1 + .4*x2 - .75*x3 + u3
ge y4s = -6 + 1*x1 - .3*x2 + 3*x3 - .4*x4 + u4
ge y1 = y1s>0
ge y2 = y2s>0
ge y3 = y3s>0
ge y4 = y4s>0

mvprobit (y1=x1) (y2=x1 x2) (y3 = x1 x2 x3) (y4=x1 x2 x3 x4), dr(75)
