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
facs = "0.005 0.02 0.04 0.06 0.08 0.1"
max_err = 1.03

bin(x,width)=width * (floor(x/width) + .5)

set ylabel "frequency"

binwidth=.05
set boxwidth binwidth

set title "Integrated autocorrelation time. The smaller the better (but always t_{int} >= 0.5). t_{int}=0.5 means no correlation."
set xlabel "t_{int}"
set output "tau_int.eps"
plot "summary_angle.csv" u (bin($9, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "<x>", "" u (bin($12, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "<n>", "" u (bin($15, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "conductivity", "" u (bin($23, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes t "mobility"

set title "Acceptance should be non-zero. Optimal performance around acc=70%. Any acc > 30% is ok." 
set xlabel "acc"
set output "acc.eps"
plot "summary_angle.csv" u (bin($16, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes not

binwidth=.002
set boxwidth binwidth

set title "Sanity check: expectation value of Boltzmann weight should be 1 for infinite statistics."
set xlabel "<exp(-dH)>"
set output "boltzmann.eps"
plot "summary_angle.csv" u (bin($17, binwidth)/(abs($17-1)<max_err)):(1.0) smooth freq with boxes not

set xlabel "{/Symbol m}_0 [eV]"
set logscale y

phonon_num(x) = .5 / tanh(x*0.006/2) - .5
lfac(x) = .5 / sqrt(phonon_num(b) + .5)

set key bottom
set title "hopping rescaling, symmetric; f=1 (maximal)"
set ylabel "-{/Symbol l}<x>"
set output "x_avg.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):(-lfac(b)*$7):8 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}"

set title "hopping rescaling, 1D; f=1 (maximal)"
set ylabel "-{/Symbol l}<x>"
set output "x_avg_1d.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0)/($5==0.1)):(-lfac(b)*$7):8 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}"

set title "phonon number deviation, symmetric; f=1 (maximal)"
set ylabel "<m{/Symbol w}_0x^2>-1/2 - <n_p>_{non-interacting}"
set output "quantum_avg.eps"
plot for [i=1:5] "summary_angle.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):($18-phonon_num(word(betas,i))):19 w yerrorbars t "{/Symbol b} = ".word(betas,i)."eV^{-1}"

set title "phonon number deviation, 1D; f=1 (maximal)"
set ylabel "<m{/Symbol w}_0x^2>-1/2 - <n_p>_{non-interacting}"
set output "quantum_avg_1d.eps"
plot for [i=1:5] "summary_angle.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0)/($5==0.1)):($18-phonon_num(word(betas,i))):19 w yerrorbars t "{/Symbol b} = ".word(betas,i)."eV^{-1}"

unset logscale y
set key top
set title "phonon number, symmetric; f=1 (maximal)"
set ylabel "<m{/Symbol w}_0x^2>-1/2"
set output "x2_avg.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):18:19 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] phonon_num(word(betas,i)) lc i t "non-interacting, {/Symbol b} = ".word(betas,i)."eV^{-1}",\
for [i=1:5] 1./(0.006*word(betas,i)) lc i dt 3 t "classical, {/Symbol b} = ".word(betas,i)."eV^{-1}"

set title "phonon number, 1D; f=1 (maximal)"
set ylabel "<m{/Symbol w}_0x^2>-1/2"
set output "x2_avg_1d.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0)/($5==0.1)):18:19 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] phonon_num(word(betas,i)) lc i t "non-interacting, {/Symbol b} = ".word(betas,i)."eV^{-1}",\
for [i=1:5] 1./(0.006*word(betas,i)) lc i dt 3 t "classical, {/Symbol b} = ".word(betas,i)."eV^{-1}"

set logscale y
set key bottom
set title "electron number, symmetric; f=1 (minimal)"
set ylabel "<n>"
set output "n_avg.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):10:11 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):($10*exp($4*$3)) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

set title "electron number, 1D; f=1 (minimal)"
set ylabel "<n>"
set output "n_avg_1d.eps"
plot [][:1] for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0)/($5==0.1)):10:11 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0)/($5==0.1)):($10*exp($4*$3)) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

set title "electron number, 1D, {/Symbol q} = {/Symbol p}; f=1 (minimal)"
set ylabel "<n>"
set output "n_avg_1d_inv.eps"
plot [][:1] for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==3.14159)/($5==0.1)):10:11 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==3.14159)/($5==0.1)):($10*exp($4*$3)) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

set title "conductivity, symmetric; f=1 (maximal)"
set ylabel "{/Symbol s}"
set output "cond.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):($13):($14) w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):($13*exp($4*$3)) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

set title "conductivity, 1D; f=1 (maximal)"
set ylabel "{/Symbol s}"
set output "cond_1d.eps"
plot [][:1] for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0)/($5==0.1)):($13):($14) w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0)/($5==0.1)):($13*exp($4*$3)) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

unset logscale y
set key top
mfac(x) = 3.59
set ylabel "{/Symbol m} [cm^2/(Vs)]"
set title "mobility, symmetric; f=1 (maximal)"
set output "mobility.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):(mfac(b)*$21):(mfac(b)*$22) w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u ($3/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0.942478)/($5==0.1)):(mfac($4)*$21) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

set title "mobility, 1D; f=1 (maximal)"
set output "mobility_1d.eps"
plot for [b in betas] "summary_angle.csv" u ($3/($4==b)/(abs($17-1)<max_err)/($6==0)/($5==0.1)):(mfac(b)*$21):(mfac(b)*$22) w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] "summary_zero.csv" u 3:(mfac($4)*$21/($4==word(betas,i))/(abs($17-1)<max_err)/($6==0)/($5==0.1)) w l lc i t "{/Symbol b} = ".word(betas,i)."eV^{-1}, asymptotic"

set key box opaque
set key width -4
set title "mobility, 1D; {/Symbol b} = 40eV^{-1}"
set output "mobility_1d_inv.eps"
plot for [i=1:5] "summary_angle.csv" u ($3/($4==40)/(abs($17-1)<max_err)/($6==0)/($5==word(facs, i))):(mfac($4)*$21):(mfac($4)*$22) w yerrorbars t "{/Symbol q} = 0, J = ".word(facs, i)."eV",\
for [i=1:5] "summary_angle.csv" u ($3/($4==40)/(abs($17-1)<max_err)/($6==3.14159)/($5==word(facs, i))):(mfac($4)*$21):(mfac($4)*$22) w yerrorbars lc i t "{/Symbol q} = {/Symbol p}, J = ".word(facs, i)."eV"

set key nobox
set key width 0
set title "mobility; {/Symbol m}_0 <= -0.1 (<n> -> 0); f=1 (maximal)"
set xlabel "{/Symbol b} [eV^{-1}]"
set output "mobility_of_beta.eps"
plot for [th=0:10:2] "summary_angle.csv" u ($4/($3<-0.07)/(abs($17-1)<max_err)/(abs($6-th*0.314159)<0.01)/($5==0.1)):(mfac($4)*$21):(mfac($4)*$22) w yerrorbars t "10{/Symbol q}/{/Symbol p} = ".th

#set xtics add ("Rubrene" 0.21) rotate by 40 right
#set xtics add ("Pentacene" 2.54) rotate by 40 right
#set xtics add ("C_{10}-DNTT" 0.678) rotate by 40 right
set xtics add ("C_{10}-DNBDT" 0.83) rotate by 40 right
set title "mobility; {/Symbol m}_0 <= -0.15 (<n> -> 0); {/Symbol b} = 40eV^{-1}"
set xlabel "{/Symbol q}"
set output "mobility_of_theta.eps"
plot for [fac in facs] "summary_angle.csv" u ($6/($3<-0.13)/(abs($17-1)<max_err)/($4==40)/($5==fac)):(mfac($4)*$21):(mfac($4)*$22) w yerrorbars t "J = ".fac."eV"

#mfac(x) = pi/(.08*x)**2 # this goes to 1 at zero hopping
mfac(x) = 395./x # this appears to be consistent with Alessandro's paper's formulas
set ylabel "L^2"

set title "squared transient localization length; {/Symbol m}_0 <= -0.15 (<n> -> 0); {/Symbol b} = 40eV^{-1}"
set output "lsq_of_theta.eps"
plot for [fac in facs] "summary_angle.csv" u ($6/($3<-0.13)/(abs($17-1)<max_err)/($4==40)/($5==fac)):(mfac($4)*$21):(mfac($4)*$22) w yerrorbars t "J = ".fac."eV"
