TITLE Potassium D-Current for Cortical Neuron Axon

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
	SUFFIX cortical_axon_i_kd
	USEION k WRITE ik				: Using k ion, treat the reversal potential as a parameter and write to ik so the total k current can be tracked
	RANGE g_Kd, i_Kd				: Potassium current, specific conductance and equilibrium potential
}

PARAMETER {
	ek = -90 (mV)
	i_Kd = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_Kd = 0.6e-3 (S/cm2)
	tau_m = 1 (ms)
	tau_h = 1500 (ms)
	V_half_m = -43 (mV)
	V_half_h = -67 (mV)
	q_m = 8
	q_h = 7.3
	Q_s = 3.209						: Temperature rescaling - Q_10 = 2.3 => Q_s = (Q_10)^((37-23)/10) = 3.209
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	m_inf
	h_inf
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_Kd*m*h*(v - ek)
	i_Kd = ik 						: Record i_Kv (just this potassium current) to check it is working
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
	TABLE m_inf, h_inf FROM -100 TO 100 WITH 400
	
	m_inf = 1-(1/(1+exp((v-V_half_m)/q_m)))
	h_inf = 1-(1/(1+exp((v-V_half_h)/q_h)))
	
}

UNITSON 