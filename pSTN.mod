TITLE  STN ion channels for single compartment model

:
: Na+, K, CaT, CaL, A and AHP current 
: 


NEURON {
	SUFFIX stn
	NONSPECIFIC_CURRENT ilk
	USEION ca READ cai, cao WRITE ica, cai
	USEION k READ ki, ko WRITE ik
	USEION na READ nai, nao WRITE ina
	RANGE ina, ik, ica
	RANGE gnabar, ena, m_inf, h_inf, tau_h, tau_m		 : fast sodium
	RANGE gkdrbar, ek, n_inf, tau_n, ikD                   : delayed K rectifier
	RANGE gl, el, ilk                                      : leak
	RANGE gcatbar, eca, p_inf, tau_p, q_inf, tau_q	       : T-type ca current
	RANGE gcalbar, eca, c_inf, d1_inf, d2_inf, tau_c, tau_d1, tau_d2, icaT, icaL  : L-type ca current
	RANGE gkabar, ek, a_inf, tau_a, b_inf, tau_b, ikA      : A-type K current
	RANGE gkcabar, ek, r_inf, ikAHP                        : ca dependent AHP K current
      RANGE kca, vol, caGain                                 : ca dynamics
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	FARADAY = (faraday) (coulomb)  :units are really coulombs/mole
}

PARAMETER {
	R = 8.31441 (Gas constant)
	T 		(Absolute temp)
	celsius		(degC)

:Fast Na channel
	gnabar   = 49e-3 (S/cm2) 
	theta_m = -40 (mV)
	theta_h = -45.5 (mV) 
	k_m = -8 (mV)    
	k_h = 6.4 (mV)   
	tau_m0 = 0.2 (ms)
	tau_m1 = 3 (ms)
	tau_h0 = 0 (ms)
	tau_h1 = 24.5 (ms) 
	tht_m = -53 (mV)
	tht_h1 = -50 (mV)
	tht_h2 = -50 (mV)
	sig_m = -0.7 (mV)
	sig_h1 = -15 (mV)
	sig_h2 = 16 (mV)

: Delayed rectifier K
	gkdrbar  = 57e-3	(S/cm2)  
	theta_n = -41 (mV)
	k_n = -14 (mV)     
	tau_n0 = 0 (ms)
	tau_n1 = 11 (ms) 
	tht_n1 = -40 (mV)
	tht_n2 = -40 (mV)
	sig_n1 = -40 (mV)
	sig_n2 = 50 (mV)

:Leakage current
	gl	= 0.35e-3	(S/cm2)
	el	= -60	(mV)

:Ca dynamics
	kca   = 2        (1/ms)
      area
      vol = 3.355e-11  (L) :~20um radius sphere
      caGain = .1

:T-type ca current
	gcatbar   = 5e-3 (S/cm2)  
	theta_p = -56 (mV)
	theta_q = -85 (mV) 
	k_p = -6.7 (mV)    
	k_q = 5.8 (mV)  
	tau_p0 = 5 (ms)
	tau_p1 = 0.33 (ms)
	tau_q0 = 0 (ms)
	tau_q1 = 400 (ms) 
	tht_p1 = -27 (mV)
	tht_p2 = -102 (mV)
	tht_q1 = -50 (mV)
	tht_q2 = -50 (mV)
	sig_p1 = -10 (mV)
	sig_p2 = 15 (mV)
	sig_q1 = -15 (mV)
	sig_q2 = 16 (mV)

:Ca L current
	gcalbar   = 15e-3 (S/cm2) 
	theta_c = -30.6 (mV)
	theta_d1 = -60 (mV)
	theta_d2 = 0.1e-3 (mM)
	k_c = -5 (mV)
	k_d1 = 7.5 (mV)
	k_d2 = 0.02e-3 (mM)
	tau_c0 = 45 (ms)
	tau_c1 = 10 (ms)
	tau_d10 = 400 (ms)
	tau_d11 = 500 (ms)
	tht_c1 = -27 (mV)
	tht_c2 = -50 (mV)
	tht_d11 = -40 (mV)
	tht_d12 = -20 (mV)
	sig_c1 = -20 (mV)
	sig_c2 = 15 (mV)
	sig_d11 = -15 (mV)
	sig_d12 = 20 (mV)

	tau_d2 = 130 (ms)

:A current
	gkabar  = 5e-3	(S/cm2)  
	theta_a = -45 (mV)
	theta_b = -90 (mV) 
	k_a = -14.7 (mV)    
	k_b = 7.5 (mV)   
	tau_a0 = 1 (ms)
	tau_a1 = 1 (ms)
	tau_b0 = 0 (ms)
	tau_b1 = 200 (ms) 
	tht_a = -40 (mV)
	tht_b1 = -60 (mV)
	tht_b2 = -40 (mV)
	sig_a = -0.5 (mV)
	sig_b1 = -30 (mV)
	sig_b2 = 10 (mV)

:AHP current (Ca dependent K current)
	gkcabar   = 1e-3 (S/cm2) 
	theta_r = 0.17e-3 (mM)
	k_r = -0.08e-3 (mM)
	tau_r = 2 (ms)
	power_r = 2
	
}

ASSIGNED {
	v	(mV)
	ina	(mA/cm2)
	ik	(mA/cm2) 
	ikD	(mA/cm2)   
	ikA	(mA/cm2) 
	ikAHP	(mA/cm2)  
	ica	(mA/cm2) 
	icaT	(mA/cm2) 
	icaL 	(mA/cm2)
	ilk	(mA/cm2)

:Fast Na
	h_inf
	tau_h	(ms)
	m_inf
	tau_m	(ms)
	ena           (mV)   := 60  

:Delayed rectifier
	n_inf
	tau_n	(ms)
	ek         (mV) := -90

:ca T current
	p_inf
	q_inf
	tau_p	(ms)
	tau_q	(ms)
	eca           (mV)   :calc from Nernst

:ca L current
	c_inf
	tau_c	(ms)
	d1_inf
	tau_d1	(ms)
	d2_inf
	:tau_d2	(ms)  :in PARAMETERS

:A current
	a_inf
	tau_a	(ms)
	b_inf
	tau_b	(ms)

:AHP (Ca dependent K current)
	r_inf
}

STATE {
	m h n
	p q 
	c d1 d2
	cai (mM) <1e-10>
	cao (mM) <1e-10>
	nai (mM) <1e-10>
	nao (mM) <1e-10>
	ki (mM) <1e-10>
	ko (mM) <1e-10>
	a b 
      r
}


BREAKPOINT {
	SOLVE states METHOD cnexp

	T = 273 + celsius - 9.5
	ena = -(R*T)/FARADAY*log(nai/nao)*1000
	ek = (R*T)/FARADAY*log(ko/ki)*1000
	eca = -(R*T)/FARADAY*log(cai/cao)*1000/2
	:printf("%f %f %f\n", ena, ek, eca)

	ina   = gnabar * m*m*m*h * (v - ena)
	ikD   = gkdrbar * n^4 * (v - ek)
	ikA   = gkabar * a*a*b * (v - ek)
	ikAHP   = gkcabar * (v - ek)*r^(power_r)
	ik=ikD+ikA+ikAHP
	icaT   = gcatbar * p*p*q * (v - eca)
	icaL   = gcalbar * c*c*d1*d2 * (v - eca)
	ica=icaT+icaL
	ilk = gl * (v - el)

}

DERIVATIVE states {   
	evaluate_fct(v)
	h' = (h_inf - h)/tau_h
	m' = (m_inf - m)/tau_m
	n' = (n_inf - n)/tau_n
	p' = (p_inf - p)/tau_p
	q' = (q_inf - q)/tau_q

      evaluate_fct2(cai)
	c' = (c_inf - c)/tau_c
	d1' = (d1_inf - d1)/tau_d1
	d2' = (d2_inf - d2)/tau_d2

      :(Ica mA/cm2)*(area um2)*(1e-8 cm2/um2)*(1e-3 A/mA)*(1/(2*F) mol/C)*(1e-3 sec/msec)*(1e3 mMol/mol)(1/volume 1/L)=(mM/msec)
	cai' = caGain*(-ica*area*1e-11/(2*FARADAY*vol) - kca*cai)
:	cai' = -ica*area*somaAreaFrac*1e-11/(2*FARADAY*vol*shellVolFrac) + (5e-6 - cai)/kca

	a' = (a_inf - a)/tau_a
	b' = (b_inf - b)/tau_b

	r' = (r_inf - r)/tau_r
}

UNITSOFF

INITIAL {

	evaluate_fct(v)
	m = m_inf 
	h = h_inf   
	n = n_inf   
	p = p_inf 
	q = q_inf  

	evaluate_fct2(cai)
	c = c_inf 
	d1 = d1_inf  
	d2 = d2_inf   

	a = a_inf 
	b = b_inf   

	r = r_inf 
}

PROCEDURE evaluate_fct(v(mV)) { 
:Fast Na current
	h_inf = 1/(1+exp((v-theta_h)/k_h))
	m_inf = 1/(1+exp((v-theta_m)/k_m))
	tau_h = tau_h0 + tau_h1/(exp(-(v-tht_h1)/sig_h1) + exp(-(v-tht_h2)/sig_h2)) 
	tau_m = tau_m0 + tau_m1/(1+exp(-(v-tht_m)/sig_m)) 

:Delayed rectifier K
	n_inf = 1/(1+exp((v-theta_n)/k_n))
	tau_n = tau_n0 + tau_n1/(exp(-(v-tht_n1)/sig_n1) + exp(-(v-tht_n2)/sig_n2)) 

:Ca T current
	p_inf = 1/(1+exp((v-theta_p)/k_p))
	q_inf = 1/(1+exp((v-theta_q)/k_q))
	tau_p = tau_p0 + tau_p1/(exp(-(v-tht_p1)/sig_p1) + exp(-(v-tht_p2)/sig_p2)) 
	tau_q = tau_q0 + tau_q1/(exp(-(v-tht_q1)/sig_q1) + exp(-(v-tht_q2)/sig_q2))

:Ca L current
	c_inf = 1/(1+exp((v-theta_c)/k_c))
	d1_inf = 1/(1+exp((v-theta_d1)/k_d1))
	tau_c = tau_c0 + tau_c1/(exp(-(v-tht_c1)/sig_c1) + exp(-(v-tht_c2)/sig_c2))  
	tau_d1 = tau_d10 + tau_d11/(exp(-(v-tht_d11)/sig_d11) + exp(-(v-tht_d12)/sig_d12))  

:A current
	a_inf = 1/(1+exp((v-theta_a)/k_a))
	b_inf = 1/(1+exp((v-theta_b)/k_b))
	tau_a = tau_a0 + tau_a1/(1+exp(-(v-tht_a)/sig_a))
	tau_b = tau_b0 + tau_b1/(exp(-(v-tht_b1)/sig_b1) + exp(-(v-tht_b2)/sig_b2))  

}

PROCEDURE evaluate_fct2(cai(mM)) { 
:Ca L current
	d2_inf = 1/(1+exp((cai-theta_d2)/k_d2))

:AHP current
	r_inf = 1/(1+exp((cai-theta_r)/k_r))
}

UNITSON