

:used to set nai, nao, ki, ko so that NEURON's default values are not used

NEURON {
	SUFFIX myions
	USEION k WRITE ko, ki
	USEION na WRITE nao, nai
}

UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
}

ASSIGNED {
	nai (mM)
	nao (mM)
	ki (mM)
	ko (mM)
}

BREAKPOINT {
	nai = nai
	nao = nao
	ki = ki
	ko = ko
}