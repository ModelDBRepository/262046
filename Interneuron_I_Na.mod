TITLE Sodium Current for Cortical Interneuron

COMMENT
  
  Model Reference: 
  
  Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
  Bal, T., Frégnac, Y., Markram, H. and Destexhe, A., 2008. 
  "Minimal Hodgkin–Huxley type models for different classes of 
  cortical and thalamic neurons." 
  Biological cybernetics, 99(4-5), pp.427-441.
  
  Original Code Link:
  https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=123623
  
  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
  
  Edits: 
  
ENDCOMMENT


UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX interneuron_i_na
	USEION na WRITE ina				: Using na ion, treat the reversal potential as a parameter and write to ik so the total k current can be tracked
	RANGE g_Na, i_Na				: Sodium current, specific conductance and equilibrium potential
}

PARAMETER {
	ena = 50 (mV)
	i_Na = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_Na = 0.05 (S/cm2)
	V_T = -55(mV)
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	alpha_m
	beta_m
	alpha_h
	beta_h
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
	m = 0
	h = 0
}

DERIVATIVE states {
	settables(v)
	m' = alpha_m*(1-m)-beta_m*m
	h' = alpha_h*(1-h)-beta_h*h
}

PROCEDURE settables(v) {
	TABLE alpha_m, beta_m, alpha_h, beta_h DEPEND V_T FROM -100 TO 100 WITH 400
	
	alpha_m = 0.32*vtrap(-(v-V_T-13),4)
	beta_m = 0.28*vtrap((v-V_T-40),5)
	alpha_h = 0.128*exp(-(v-V_T-17)/18)
	beta_h = 4/(1+exp(-(v-V_T-40)/5))
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}

UNITSON 