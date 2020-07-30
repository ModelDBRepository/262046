COMMENT
This is a simplified version of xtra.mod for applying extracellular stimulation to neuron sections.
ENDCOMMENT

NEURON {
	SUFFIX xtra
	RANGE rx, er
	RANGE x, y, z
	GLOBAL is
	POINTER im, ex
}

PARAMETER {
	: default transfer resistance between stim electrodes and axon
	rx = 1 (megohm) : mV/nA
	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)
}

ASSIGNED {
	v (millivolts)
	is (milliamp)
	ex (millivolts)
	im (milliamp/cm2)
	er (microvolts)
	area (micron2)
}

INITIAL {
	ex = is*rx*(1e6)
	er = (10)*rx*im*area
}

BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  ex = is*rx*(1e6)
}
AFTER SOLVE { : after each solution step
  er = (10)*rx*im*area
}