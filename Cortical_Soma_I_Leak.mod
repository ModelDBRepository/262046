TITLE Passive Leak Current for Cortical Neuron Soma

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
	SUFFIX cortical_soma_i_leak
	NONSPECIFIC_CURRENT i_l			: Declare ASSIGNED variables as RANGE variables so that they can be accessed outside of mod file
	RANGE i_l, g_l, e_l				: leak current, specific conductance and equilibrium potential
}

PARAMETER {
	g_l = 0.1e-3 (S/cm2)
	e_l = -70 (mV)
}

ASSIGNED {
	v (mV)	
	i_l (mA/cm2)
}

BREAKPOINT {
	i_l = g_l*(v - e_l)
}