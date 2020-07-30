TITLE Sodium and Potassium Currents for Thalamic Neuron

COMMENT
  
  Model Reference: 
  
  Rubin, J.E. and Terman, D., 2004. "High frequency stimulation 
  of the subthalamic nucleus eliminates pathological thalamic 
  rhythmicity in a computational model."
  Journal of computational neuroscience, 16(3), pp.211-235.

  
  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
  
  Edits: 
  
ENDCOMMENT


UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX thalamic_i_na_k
	USEION na WRITE ina				: Using na and k ions, treat the reversal potential as a parameter and write to ina and k so the total na and k currents can be tracked
	USEION k WRITE ik
	RANGE g_Na, i_na				: Sodium current, specific conductance and equilibrium potential
	RANGE g_K, i_k					: Potassium current, specific conductance and equilibrium potential
}

PARAMETER {
	ena = 50 (mV)
	ek = -90 (mV)
	i_na = 0.0 (mA/cm2)
	g_Na = 0.3 (S/cm2)
	i_k = 0.0 (mA/cm2)
	g_K = 0.5 (S/cm2)
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	ik (mA/cm2)
	h_inf
	tau_h (ms)
	m_inf
}

STATE {
	h 
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = g_Na*m_inf*m_inf*m_inf*h*(v - ena)
	ik = g_K*(0.75*(1-h))*(0.75*(1-h))*(0.75*(1-h))*(0.75*(1-h))*(v - ek)
	i_na = ina 							: Record i_na (just this sodium current) to check it is working
	i_k = ik							: Record i_k (just this potassium current) to check it is working
}

UNITSOFF

INITIAL {
	settables(v)
	h = h_inf
}

DERIVATIVE states {
	settables(v)
	h' = (h_inf - h)/tau_h
}

PROCEDURE settables(v) {
	TABLE h_inf, m_inf, tau_h FROM -100 TO 100 WITH 400
	
	h_inf = 1/(1+exp((v+41)/4))
	tau_h = 1/((0.128*exp(-(v+46)/18))+(4/(1+exp(-(v+23)/5))))
	m_inf = 1/(1+exp(-(v+37)/7))
}

UNITSON 