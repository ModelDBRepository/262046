TITLE Slow non-inactivating Potassium Current for Cortical Neuron Soma

COMMENT
  
  Model Reference: 
  
  Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
  Bal, T., Frégnac, Y., Markram, H. and Destexhe, A., 2008. 
  "Minimal Hodgkin–Huxley type models for different classes of 
  cortical and thalamic neurons." 
  Biological cybernetics, 99(4-5), pp.427-441.
  
  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
  
  Edits: 
  
ENDCOMMENT


UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cortical_soma_i_m
	USEION k WRITE ik				: Using k ion, treat the reversal potential as a parameter and write to ik so the total k current can be tracked
	RANGE g_M, i_M					: Potassium current, specific conductance and equilibrium potential
}

PARAMETER {
	ek = -100 (mV)
	i_M = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_M = 7e-5 (S/cm2)
	tau_max = 1000 (ms)
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	p_inf
	tau_p (ms)
}

STATE {
	p
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_M*p*(v - ek)
	i_M = ik 						: Record i_M (just this potassium current) to check it is working
}

UNITSOFF

INITIAL {
	settables(v)					
	p = p_inf
}

DERIVATIVE states {
	settables(v)
	p' = (p_inf - p)/tau_p
}

PROCEDURE settables(v) {
	TABLE p_inf, tau_p DEPEND tau_max FROM -100 TO 100 WITH 400
	
	p_inf = 1/(1+exp(-(v+35)/10))
	tau_p = tau_max/(3.3*exp((v+35)/20)+exp(-(v+35)/20))
}

UNITSON 