set terminal postscript landscape enhanced color lw 1 font 16
set encoding iso_8859_1
set grid lt 0 lw .3
set samples 10000
set fit results
set pointsize 1
set border lw 0.5

set linetype 1 pt 4
set linetype 2 pt 6
set linetype 3 pt 8
set linetype 4 pt 10
set linetype 5 pt 1
set linetype 6 pt 2
set linetype 7 pt 2

betas = "20 30 40 50 60"

bin(x,width)=width * (floor(x/width) + .5)

set ylabel "frequency"

binwidth=.05
set boxwidth binwidth

set xlabel "t_{int}"
set output "tau_int.eps"
plot "summary.csv" u (bin($7, binwidth)/($2>30)):(1.0) smooth freq with boxes t "<x>", "" u (bin($10, binwidth)/($2>30)):(1.0) smooth freq with boxes t "<n>", "" u (bin($13, binwidth)/($2>30)):(1.0) smooth freq with boxes t "conductivity"

set xlabel "acceptance"
set output "acc.eps"
plot "summary.csv" u (bin($14, binwidth)/($2>30)):(1.0) smooth freq with boxes not

binwidth=.002
set boxwidth binwidth

set xlabel "<exp(-dH)>"
set output "boltzmann.eps"
plot "summary.csv" u (bin($15, binwidth)/($2>30)):(1.0) smooth freq with boxes not

set xlabel "{/Symbol m}_0 [eV]"
set logscale y

phonon_num(x) = .5 / tanh(x*0.006/2) - .5
lfac(x) = .5 / sqrt(phonon_num(b) + .5)
set title "hopping rescaling; L=10,12,15,18; Nt=32"
set ylabel "-{/Symbol l}<x>"
set output "x_avg.eps"
plot for [b in betas] "summary.csv" u ($3/($4==b)/($2>30)/($1>9)):(-lfac(b)*$5):6 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}"

set title "phonon number deviation; L=10,12,15,18; Nt=32"
set ylabel "<m{/Symbol w}_0x^2>-1/2 - <n_p>_{non-interacting}"
set output "quantum_avg.eps"
plot for [i=1:5] "summary.csv" u ($3/($4==word(betas,i))/($2>30)/($1>9)):($16-phonon_num(word(betas,i))):17 w yerrorbars t "{/Symbol b} = ".word(betas,i)."eV^{-1}"

unset logscale y
set title "phonon number; Nt=32"
set ylabel "<m{/Symbol w}_0x^2>-1/2"
set output "x2_avg.eps"
plot for [b in betas] "summary.csv" u ($3/($4==b)/($2>30)):16:17 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}",\
for [i=1:5] phonon_num(word(betas,i)) lc i t "non-interacting, {/Symbol b} = ".word(betas,i)."eV^{-1}",\
for [i=1:5] 1./(0.006*word(betas,i)) lc i dt 3 t "classical, {/Symbol b} = ".word(betas,i)."eV^{-1}"

set logscale y
set title "electron number; L=10,12,15,18; Nt=32"
set ylabel "<n>"
set output "n_avg.eps"
plot for [b in betas] "summary.csv" u ($3/($4==b)/($2>30)/($1>9)):8:9 w yerrorbars t "{/Symbol b} = ".b."eV^{-1}"

set title "conductivity; L=10,12,15,18; Nt=32"
set ylabel "{/Symbol s}"
set output "cond.eps"
plot for [b in betas] "summary.csv" u ($3/($4==b)/($2>30)/($1>9)):($11):($12) w yerrorbars t "{/Symbol b} = ".b."eV^{-1}"

unset logscale y
mfac(x) = 25.8
set ylabel "{/Symbol m}_{pentacene} [cm^2/(Vs)]"
#mfac(x) = pi/(.08*x)**2 # this goes to 1 at zero hopping
#mfac(x) = 395./x # this appears to be consistent with Alessandro's paper's formulas
#set ylabel "L^2"
set title "mobility; Nt=32"
set output "mobility.eps"
plot for [L in "10 12 15 18"] for [b in betas] "summary.csv" u ($3/($4==b)/($2>30)/($1==L)):(mfac(b)*$11/$8):(mfac(b)*sqrt(($12/$8)**2+($11*$9/$8**2)**2)) w yerrorbars t "{/Symbol b} = ".b."eV^{-1}".", L = ".L

set title "mobility; Nt=32; {/Symbol m}_0 < 0"
set xlabel "{/Symbol b} [eV^{-1}]"
set output "mobility_of_beta.eps"
plot for [L in "10 12 15 18"] "summary.csv" u ($4/($3<0)/($2>30)/($1==L)):(mfac($4)*$11/$8):(mfac($4)*sqrt(($12/$8)**2+($11*$9/$8**2)**2)) w yerrorbars t "L = ".L
