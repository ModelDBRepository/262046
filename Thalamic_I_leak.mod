TITLE Passive Leak Current for Thalamic Neuron

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
	SUFFIX thalamic_i_leak
	NONSPECIFIC_CURRENT i_l			: Declare ASSIGNED variables as RANGE variables so that they can be accessed outside of mod file
	RANGE i_l, g_l, e_l				: leak current, specific conductance and equilibrium potential
}

PARAMETER {
	g_l = 0.005 (S/cm2)
	e_l = -70 (mV)
}

ASSIGNED {
	v (mV)	
	i_l (mA/cm2)
}

BREAKPOINT {
	i_l = g_l*(v - e_l)
}