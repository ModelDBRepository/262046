TITLE Delayed-rectifier Potassium Current for Cortical Neuron Soma

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
	SUFFIX cortical_soma_i_k
	USEION k WRITE ik				: Using k ion, treat the reversal potential as a parameter and write to ik so the total k current can be tracked
	RANGE g_K, i_K					: Potassium current, specific conductance and equilibrium potential
}

PARAMETER {
	ek = -100 (mV)
	i_K = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_K = 5e-3 (S/cm2)
	V_T = -55(mV)
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	alpha_n
	beta_n
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_K*n*n*n*n*(v - ek)
	i_K = ik 						: Record i_K (just this potassium current) to check it is working
}

UNITSOFF

INITIAL {
	settables(v)					: ** Need to double check these intials are correct
	n = 0
}

DERIVATIVE states {
	settables(v)
	n' = alpha_n*(1-n)-beta_n*n
}

PROCEDURE settables(v) {
	TABLE alpha_n, beta_n DEPEND V_T FROM -100 TO 100 WITH 400
	alpha_n = 0.032 * vtrap(-(v-V_T-15), 5)
	beta_n = 0.5*exp(-(v-V_T-10)/40)
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}


UNITSON 