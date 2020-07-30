TITLE simple AMPA receptors (discrete connections)

COMMENT
-----------------------------------------------------------------------------

	Simple model for glutamate AMPA receptors
	=========================================

  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS

    Whole-cell recorded postsynaptic currents mediated by AMPA/Kainate
    receptors (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994) were used
    to estimate the parameters of the present model; the fit was performed
    using a simplex algorithm (see Destexhe et al., J. Computational Neurosci.
    1: 195-230, 1994).

  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)

    The simplified model was obtained from a detailed synaptic model that 
    included the release of transmitter in adjacent terminals, its lateral 
    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe
    and Sejnowski, 1995).  Short pulses of transmitter with first-order
    kinetics were found to be the best fast alternative to represent the more
    detailed models.

  - ANALYTIC EXPRESSION

    The first-order model can be solved analytically, leading to a very fast
    mechanism for simulating synapses, since no differential equation must be
    solved (see references below).



References

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
   computing synaptic conductances based on a kinetic model of receptor binding
   Neural Computation 6: 10-14, 1994.  

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a 
   common kinetic formalism, Journal of Computational Neuroscience 1: 
   195-230, 1994.

See also:

   http://www.cnl.salk.edu/~alain
   http://cns.fmed.ulaval.ca

-----------------------------------------------------------------------------
ENDCOMMENT



:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS AMPA_S
	RANGE C, R, R0, R1, g, Cmax
	NONSPECIFIC_CURRENT i
	:GLOBAL Cdur, Alpha, Beta, Erev, blockTime, Rinf, Rtau
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 0.5	(mM)		: max transmitter concentration	(set = 1 to match ~/netcon/ampa.hoc)
	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Alpha	= 0.94	(/ms mM)	: forward (binding) rate
	Beta	= 0.18	(/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		: reversal potential
	blockTime = 4	(ms)		: time window following dbs event during which non-dbs events are blocked
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	R					: fraction of open channels
	R0					: open channels at start of time period
	Rinf				: steady state channels open
	Rtau		(ms)	: time constant of channel binding
	on					: rising phase of PSC
	gmax				: max conductance
	tLast
	nspike
	collisionBlock
}

INITIAL {
	R = 0
	C = 0
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	on = 0
	R0 = 0
	nspike = 0
	collisionBlock = 0
}

BREAKPOINT {
	SOLVE release
	i = R*(v - Erev)
}

PROCEDURE release() {


	if (on) {				: transmitter being released?

	   R = gmax*Rinf + (R0 - gmax*Rinf) * exptable (- (t - tLast) / Rtau)
				
	} else {				: no release occuring

  	   R = R0 * exptable (- Beta * (t - tLast))
	}

}


: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first

NET_RECEIVE(weight, ncType, ncPrb) {LOCAL ok, tmp

	:ncType 0=presyn cell, 1=dbs activated axon
	:MOVED TO dbsStim.mod 4/11/07		ncPrb probability that incoming event causes PSP

	INITIAL {
	}

	: flag is an implicit argument of NET_RECEIVE and  normally 0
        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
		ok = 0

		if (ncType == 1) {
			collisionBlock = collisionBlock + 1
			net_send(blockTime, -1)

:			tmp = scop_random()
:			if (tmp <= ncPrb) {
:				ok = 1
:			}
			ok = 1

		}
		else 
		if (collisionBlock == 0) {
			ok = 1
		} 

		if (ok) {
			if (!on) {
				on = 1
				tLast = t
				R0 = R
				gmax = weight	:weight not additive from separate sources as in original ampa.mod
			}

			nspike = nspike + 1
			: come again in Cdur with flag = current value of nspike
			net_send(Cdur, nspike)
		}
      }
	else
	if (flag == nspike) { : if this associated with last spike then turn off
		if (on) {
			on = 0
			tLast = t
			R0 = R
			gmax = 0
		}
	} 
	else
	if (flag == -1) {
		collisionBlock = collisionBlock - 1
	}


}


FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000

	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}
