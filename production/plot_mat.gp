set terminal postscript landscape enhanced color lw 1 font 16
set encoding iso_8859_1
set grid lt 0 lw .3
set samples 10000
set fit results
set pointsize 1
set border lw 0.5

set linetype cycle 16
set linetype 1 pt 4
set linetype 2 pt 6
set linetype 3 pt 8
set linetype 4 pt 10
set linetype 5 pt 12 lc 7
set linetype 6 pt 1
set linetype 7 pt 2 lc 8
set linetype 8 pt 3
set linetype 10 pt 7

!rm -f *.eps

betas = "20 30 40 50 60"
mat = "rubrene pentacene C10-DNTT C10-DNBDT DNTT C8-DNTT-C8"
experiment = "8.6 1.45 8.5 12.1"
max_err = 0.10

bin(x,width)=width * (floor(x/width) + .5)

set ylabel "frequency"

binwidth=.05
set boxwidth binwidth

set title "Integrated autocorrelation time. The smaller the better (but always t_{int} >= 0.5). t_{int}=0.5 means no correlation."
set xlabel "t_{int}"
set output "tau_int.eps"
plot "summary_mat.csv" u (bin($9, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "<x>", "" u (bin($12, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "<n>", "" u (bin($15, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "conductivity", "" u (bin($23, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "mobility"

set title "Acceptance should be non-zero. Optimal performance around acc=70%. Any acc > 30% is ok." 
set xlabel "acc"
set output "acc.eps"
plot "summary_mat.csv" u (bin($16, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes not

binwidth=.002
set boxwidth binwidth

set title "Sanity check: expectation value of Boltzmann weight should be 1 for infinite statistics."
set xlabel "<exp(-dH)>"
set output "boltzmann.eps"
plot "summary_mat.csv" u (bin($17, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes not

set xlabel "{/Symbol m}_0 [eV]"
set logscale y

b = 40

set title "electron number; {/Symbol b} = 40 eV^{-1}"
set ylabel "<n>"
set output "n_avg.eps"
plot for [i=1:4] "summary_mat.csv" u ($3/($4==b)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):10:11 w yerrorbars t word(mat, i)

set title "conductivity"
set ylabel "{/Symbol s}"
set output "cond.eps"
plot for [i=1:4] "summary_mat.csv" u ($3/($4==b)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):($13):($14) w yerrorbars t word(mat, i)

unset logscale y
mfac(x) = 0.1519*x
#mfac(x) = 3.59
set ylabel "{/Symbol m} [cm^2/(Vs)]"
set title "mobility; {/Symbol b} = 40 eV^{-1}"
set output "mobility.eps"
plot for [i=1:4] "summary_mat.csv" u ($3/($4==b)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):(mfac($6)*$21):(mfac($6)*$22) w yerrorbars t word(mat, i), for [i=1:4] word(experiment, i)+0 lt i t word(mat, i).", exp"

set title "mobility; {/Symbol m}_0 <= -0.15 (<n> -> 0)"
set xlabel "{/Symbol b} [eV^{-1}]"
set output "mobility_of_beta.eps"
plot for [i=1:4] "summary_mat.csv" u ($4/($3<-0.13)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):(mfac($6)*$21):(mfac($6)*$22) w yerrorbars not, for [i=1:6] "summary_mat_z.csv" u ($4/($3<-0.13)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):(mfac($6)*$21):(mfac($6)*$22) lt i w yerrorbars t word(mat, i)

kB = 8.617333262e-5
set xlabel "T [K]"
set output "mobility_of_temp.eps"
plot for [i=1:4] "summary_mat.csv" u (1./kB/$4/($3<-0.13)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):(mfac($6)*$21):(mfac($6)*$22) w yerrorbars not, for [i=1:6] "summary_mat_z.csv" u (1./kB/$4/($3<-0.13)/(abs($17-1)<max_err)/(stringcolumn(5) eq word(mat, i))):(mfac($6)*$21):(mfac($6)*$22) lt i w yerrorbars t word(mat, i)

##mfac(x) = pi/(.08*x)**2 # this goes to 1 at zero hopping
#mfac(x) = 395./x # this appears to be consistent with Alessandro's paper's formulas
#set ylabel "L^2"
#
#set title "squared transient localization length; {/Symbol m}_0 <= -0.1 (<n> -> 0); {/Symbol b} = 40eV^{-1}"
#set output "lsq_of_theta.eps"
#plot for [fac in "0.005 0.02 0.04 0.06 0.08 0.1"] "summary_mat.csv" u ($6/($3<-0.07)/(abs($17-1)<max_err)/($4==40)/($5==fac)):(mfac($4)*$21):(mfac($4)*$22) w yerrorbars t "J = ".fac."eV"
