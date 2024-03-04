dyn.load("organic_hmc.so")

simulate_organic_hmc = function(L, Nt, hop=rep(.08,3), dHop=rep(.5,3), mu0=0, mu=mu0-2*sum(hop), omega=0.006, beta=40, steps=1, traj.length=1, therm=0, meas=20, meas.freq=1, obs="/dev/null", corr="/dev/null"){
	stopifnot(length(L) == 2)
	stopifnot(all(L >= 1) & Nt >= 1)
	stopifnot(length(hop) == 3 & length(dHop) == 3)
	stopifnot(missing(mu0) | missing(mu))

	if(!missing(mu)) mu0 = 2*sum(hop) + mu

	dim = L[[1]]*L[[2]]*Nt*3
	cc.dim = 3^2*Nt*meas/meas.freq

	mu = mu*beta
	omega = omega*beta
	hop = hop*beta
	mu0 = mu0*beta

	t0 = sqrt(hop[[1]]*hop[[2]] + hop[[1]]*hop[[3]] + hop[[2]]*hop[[3]])
	if(missing(dHop)) dHop = dHop * sqrt(2*tanh(omega/2))

	free.case = free_el_num_matsubara(L, Nt, hop, mu)
	n0 = free.case$n0
	n = free.case$n
	cc.corr = free.case$cc.corr
	b0 = free_phon_num(omega)

	if(traj.length == 0){
		traj.length = 1 / sqrt(omega/Nt)
	}

	results = list(L = L, Nt = Nt, t = hop, kappa = dHop, mu = mu, omega = omega, mu0 = mu0, t0 = t0, n0 = n0, n = n, cc.corr = cc.corr, b0 = b0, steps = steps, traj.length = traj.length, therm = therm, meas = meas/meas.freq, dim = dim)

	x = .Call("simulate_hmc", L[[1]], L[[2]], Nt, hop/Nt, dHop*hop/Nt, mu/Nt, omega/Nt, steps, traj.length, therm, meas, meas.freq, obs, corr)

	results$greens = array(x[1:cc.dim]*(Nt/beta)^2, dim = c(3, 3, Nt, meas/meas.freq))
	results$res = t(matrix(x[-(1:cc.dim)], ncol=meas/meas.freq))

	return(results)
}

free_el_num_approx = function(mu0, t0, hop, mu){
	if(missing(mu0)) mu0 = 2*sum(hop) + mu
	if(missing(t0)) t0 = sqrt(hop[[1]]*hop[[2]] + hop[[1]]*hop[[3]] + hop[[2]]*hop[[3]])
	if(mu0 < 0) return(exp(mu0) / t0 / (4*pi))
	else return(mu0 / t0 / (4*pi))
}

dispersion = function(hop, k, mu){
	e = mu
	e = e + 2*hop[[1]] * cos(k[[1]])
	e = e + 2*hop[[2]] * cos(k[[1]] + k[[2]])
	e = e + 2*hop[[3]] * cos(k[[2]])
	return(-e)
}

current = function(hop, k){
	mom = c(k[[1]], k[[1]] + k[[2]], k[[2]])
	return(2*hop * sin(mom))
}

cur.cur = function(hop, k, i){
	nn = length(hop)
	cur = current(hop, k)
	return(cur[[(i-1)%/%nn+1]] * cur[[(i-1)%%nn+1]])
}

free_el_num_matsubara = function(L, Nt, hop, mu, beta=1){
	hop = hop*beta
	mu = mu*beta

	p1 = 2*pi * c(0:(L[[1]]-1)) / L[[1]]
	p2 = 2*pi * c(0:(L[[2]]-1)) / L[[2]]
	#frec = pi * (2*c(0:(Nt-1)) + 1) / Nt

	#n = sapply(p1, function(k1){
	#			   sapply(p2, function(k2){
	#						  mean(sapply(frec, function(w){
	#										  1/(1i*(w-pi) + dispersion(hop, c(k1, k2), mu) / Nt)}
	#										 ))}
	#			   )}
	#)

	n.ferm = sapply(p1, function(k1){
				   sapply(p2, function(k2){
							  1/(1 + exp(dispersion(hop, c(k1, k2), mu)))}
				   )}
	)

	current.corr = sapply(1:length(hop)^2, function(i){
							  sapply(p1, function(k1){
										 sapply(p2, function(k2){
													#cur.cur(hop, c(k1, k2), i) * exp(-dispersion(hop, c(k1, k2), mu))}
													cur.cur(hop, c(k1, k2), i)/(1 + exp(dispersion(hop, c(k1, k2), mu)))/(1 + exp(-dispersion(hop, c(k1, k2), mu)))}
													#cur.cur(hop, c(k1, k2), i)/(1 + exp(dispersion(hop, c(k1, k2), mu)))^2}
										 )}
							  )}
	)

	z = mean(
			 sapply(p1, function(k1){
						sapply(p2, function(k2){
								   exp(-dispersion(hop, c(k1, k2), mu))}
						)}
			 )
	)

	return(list(p1 = p1, p2 = p2, rho = n.ferm, n = mean(n.ferm), n0 = free_el_num_approx(hop=hop, mu=mu), cc.corr = matrix(apply(current.corr, 2, mean), length(hop))))
}

free_phon_num = function(omega, beta=1){
	1/(exp(beta*omega) - 1)
}
