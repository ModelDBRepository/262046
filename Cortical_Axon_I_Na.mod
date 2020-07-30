TITLE Sodium Current for Cortical Neuron Axon

COMMENT
  
  Model Reference: 
  
  Foust, A.J., Yu, Y., Popovic, M., Zecevic, D. and McCormick, D.A., 
  2011. "Somatic membrane potential and Kv1 channels control spike 
  repolarization in cortical axon collaterals and presynaptic boutons." 
  Journal of Neuroscience, 31(43), pp.15490-15498.
  
  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
  
  Edits: 
  
ENDCOMMENT


UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cortical_axon_i_na
	USEION na WRITE ina				: Using na ion, treat the reversal potential as a parameter and write to ina so the total na current can be tracked
	RANGE g_Na, i_Na				: Sodium current, specific conductance and equilibrium potential
}

PARAMETER {
	ena = 60 (mV)
	i_Na = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_Na = 0.4 (S/cm2)				: Default is assuming AIS
	Q_s = 3.209						: Temperature rescaling - Q_10 = 2.3 => Q_s = (Q_10)^((37-23)/10) = 3.209
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	alpha_m
	beta_m
	alpha_h
	beta_h
	m_inf
	tau_m (ms)
	h_inf 
	tau_h (ms)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = g_Na*m*m*m*h*(v - ena)
	i_Na = ina 						: Record i_Na (just this sodium current) to check it is working
}

UNITSOFF

INITIAL {
	settables(v)
	m = m_inf
	h = h_inf
}

DERIVATIVE states {
	settables(v)
	m' = (m_inf-m)/tau_m
	h' = (h_inf-h)/tau_h
}

PROCEDURE settables(v) {
	TABLE alpha_m, beta_m, alpha_h, beta_h, m_inf, tau_m, h_inf, tau_h FROM -100 TO 100 WITH 400
	
	alpha_m = Q_s*0.182*vtrap(-(v+30),8)
	beta_m = Q_s*0.124*vtrap((v+30),8)
	alpha_h = Q_s*0.028*vtrap(-(v+45),6)
	beta_h = Q_s*0.0091*vtrap((v+70),6)
	
	m_inf = alpha_m/(alpha_m+beta_m)
	tau_m = 1/(alpha_m+beta_m)
	h_inf = 1/(1+exp((v+60)/6.2))
	tau_h = 1/(alpha_h+beta_h)
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}

UNITSON 