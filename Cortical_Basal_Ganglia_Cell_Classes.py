# -*- coding: utf-8 -*-
""" ------------------------------------------------------------------------------------
	Cortical Basal Ganglia Neurons: file containing classes for defining network neurons
	------------------------------------------------------------------------------------
	
								Model References 
	------------------------------------------------------------------------------------
	Cortical Pyramical Cell Soma:
	Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
	Bal, T., Frégnac, Y., Markram, H. and Destexhe, A., 2008. 
	"Minimal Hodgkin–Huxley type models for different classes of 
	cortical and thalamic neurons." 
	Biological cybernetics, 99(4-5), pp.427-441.
	
	Cortical Pyramidal Cell Axon: 
	Foust, A.J., Yu, Y., Popovic, M., Zecevic, D. and McCormick, D.A., 
	2011. "Somatic membrane potential and Kv1 channels control spike 
	repolarization in cortical axon collaterals and presynaptic boutons." 
	Journal of Neuroscience, 31(43), pp.15490-15498.
	
	Cortical Interneurons:
	Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
	Bal, T., Frégnac, Y., Markram, H. and Destexhe, A., 2008. 
	"Minimal Hodgkin–Huxley type models for different classes of 
	cortical and thalamic neurons." 
	Biological cybernetics, 99(4-5), pp.427-441.
  
	STN Neurons:
	Otsuka, T., Abe, T., Tsukagawa, T. and Song, W.J., 2004. 
	"Conductance-based model of the voltage-dependent generation 
	of a plateau potential in subthalamic neurons."
	Journal of neurophysiology, 92(1), pp.255-264.
  
	GP Neurons:
	Terman, D., Rubin, J.E., Yew, A.C. and Wilson, C.J., 2002. 
	"Activity patterns in a model for the subthalamopallidal 
	network of the basal ganglia." 
	Journal of Neuroscience, 22(7), pp.2963-2976.
	
	*Note: The NEURON implementations of the STN and GP neurons are 
		   from the following publication - 
	Hahn, P.J. and McIntyre, C.C., 2010. 
	"Modeling shifts in the rate and pattern of subthalamopallidal
	network activity during deep brain stimulation." 
	Journal of computational neuroscience, 28(3), pp.425-441.
	
	Thalamic Neurons:
	Rubin, J.E. and Terman, D., 2004. High frequency stimulation of 
	the subthalamic nucleus eliminates pathological thalamic rhythmicity 
	in a computational model. Journal of computational neuroscience, 
	16(3), pp.211-235.
	
	Model Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
	
	Edits: 	14-02-20: Added additional getters and setters for model 
					  neuron parameters
			16-01-19: Created classes for cell models so they can be 
					  utilized in PyNN.

	Created on Tues Jan 15 12:51:26 2019

"""

from math import pi
from neuron import h
from nrnutils import Mechanism, Section
from pyNN.neuron import NativeCellType
from pyNN.parameters import Sequence
import numpy as np
from scipy import signal

# Import global variables for GPe DBS
import Global_Variables as GV

try:
	reduce
except NameError:
	from functools import reduce

def _new_property(obj_hierarchy, attr_name):
	"""
	Returns a new property, mapping attr_name to obj_hierarchy.attr_name.

	For example, suppose that an object of class A has an attribute b which
	itself has an attribute c which itself has an attribute d. Then placing
		e = _new_property('b.c', 'd')
	in the class definition of A makes A.e an alias for A.b.c.d
	"""

	def set(self, value):
		obj = reduce(getattr, [self] + obj_hierarchy.split('.'))
		setattr(obj, attr_name, value)

	def get(self):
		obj = reduce(getattr, [self] + obj_hierarchy.split('.'))
		return getattr(obj, attr_name)
	return property(fset=set, fget=get)

class Cortical_Neuron(object):
	
	def __init__(self, **parameters):
		
		# Create cortical pyramidal neuron soma compartment using Pospischil (2008) single compartment model
		self.soma = Section(L=parameters['soma_L'], diam=parameters['soma_diam'], nseg=parameters['soma_nseg'], Ra=parameters['soma_Ra'], cm=parameters['soma_cm'],
							mechanisms=(Mechanism('cortical_soma_i_leak'), Mechanism('cortical_soma_i_na'), Mechanism('cortical_soma_i_k'), Mechanism('cortical_soma_i_m')))
		
		# Create cortical pyramidal neuron axon compartments using Foust (2011) axon model
		self.ais = Section(L=parameters['ais_L'], diam=parameters['ais_diam'], nseg=parameters['ais_nseg'], Ra=parameters['ais_Ra'], cm=parameters['ais_cm'],
							mechanisms=(Mechanism('cortical_axon_i_leak', g_l=3.3e-5), Mechanism('cortical_axon_i_na', g_Na=4000e-4), 
										Mechanism('cortical_axon_i_kv', g_Kv=20e-4), Mechanism('cortical_axon_i_kd', g_Kd=0.015)), parent=self.soma)

		
		# Use loop to create myelin and node sections of axon
		self.myelin = []
		self.node = []
		for i in np.arange(parameters['num_axon_compartments']):
			if i==0: 
				self.myelin.append(Section(L=parameters['myelin_L'], diam=parameters['myelin_diam'], nseg=11, Ra=parameters['myelin_Ra'], cm=parameters['myelin_cm'],
							mechanisms=(Mechanism('cortical_axon_i_leak', g_l=0), Mechanism('cortical_axon_i_na', g_Na=10e-4)),
							parent=self.ais))
			else:
				self.myelin.append(Section(L=parameters['myelin_L'], diam=parameters['myelin_diam'], nseg=11, Ra=parameters['myelin_Ra'], cm=parameters['myelin_cm'],
							mechanisms=(Mechanism('cortical_axon_i_leak', g_l=0), Mechanism('cortical_axon_i_na', g_Na=10e-4)),
							parent=self.node[i-1]))
			
			self.node.append(Section(L=parameters['node_L'], diam=parameters['node_diam'], nseg=parameters['node_nseg'], Ra=parameters['node_Ra'], cm=parameters['node_cm'],
							mechanisms=(Mechanism('cortical_axon_i_leak', g_l=0.02), Mechanism('cortical_axon_i_na', g_Na=2800e-4),
										Mechanism('cortical_axon_i_kv', g_Kv=5e-4), Mechanism('cortical_axon_i_kd', g_Kd=0.0072)),
										parent=self.myelin[i]))
										
		self.collateral = Section(L=parameters['collateral_L'], diam=parameters['collateral_diam'], nseg=parameters['collateral_nseg'], Ra=parameters['collateral_Ra'], cm=parameters['collateral_cm'],
							mechanisms=(Mechanism('cortical_axon_i_leak'), Mechanism('cortical_axon_i_na', g_Na=1333.33333e-4), 
										Mechanism('cortical_axon_i_kv', g_Kv=10e-4), Mechanism('cortical_axon_i_kd', g_Kd=0.006)), parent=self.node[-1])
		
		middle_index = int((parameters['num_axon_compartments']/2.0))
		self.middle_node = self.node[middle_index]
		self.middle_myelin = self.myelin[middle_index]
		
		# Add extracellular and xtra mechanisms to collateral 
		self.collateral.insert('extracellular')
		self.collateral.insert('xtra')
		
		# Assign default rx values to the segments rx_xtra
		#  - these values are updated in the main run file 
		# 	 where rx is calculated as the transfer resistance 
		#    for each collateral segments to the stimulation 
		#    electrode in the homogenous extracellular medium
		for seg in self.collateral :
			seg.xtra.rx = seg.x*3e-1
		
		# Setting pointers to couple extracellular and xtra mechanisms for simulating extracellular DBS 
		for seg in self.collateral:
			h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
			h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
		
		# Add bias current to neuron model - current amplitude is in terms of original model paper, nA 
		self.stim = h.IClamp(0.5, sec=self.soma)
		self.stim.delay = 0
		self.stim.dur = 1e12
		self.stim.amp = parameters['soma_bias_current_amp']
		
		# insert synaptic noise 
		self.noise = h.SynNoise(0.5, sec=self.soma)
		self.noise.f0 = 0
		self.noise.f1 = 0.3
		
		# Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
		self.AMPA = h.AMPA_S(0.5, sec=self.soma)
		self.GABAa = h.GABAa_S(0.5, sec=self.soma)
		
		# needed for PyNN
		self.source = {'soma': self.soma(0.5)._ref_v, 'middle_axon_node': self.middle_node(0.5)._ref_v, 'collateral': self.collateral(0.5)._ref_v}
		self.source_section = {'soma': self.soma, 'middle_axon_node': self.middle_node, 'collateral': self.collateral}
		self.rec = h.NetCon(self.source['collateral'], None, sec=self.source_section['collateral'])		# Needed to clear the simulator
		self.spike_times = h.Vector(0)
		self.traces = {}
		self.recording_time = False
		self.parameter_names = ()
		
	def soma_area(self):
		"""Membrane area in µm²"""
		return pi * self.soma.L * self.soma.diam
	
	def memb_init(self):
		for seg in self.soma:
			seg.v = self.v_init
	
	def _set_collateral_rx(self, sequence_values):
		rx_values = sequence_values.value
		for ii, seg in enumerate(self.collateral):
			seg.xtra.rx = rx_values[ii]
		
	def _get_collateral_rx(self):
		print("Getter Working!")
		rx_values = np.zeros((1,self.collateral.nseg))
		for i, seg in enumerate(self.collateral):
			rx_values[0,i] = seg.xtra.rx
		print(Sequence(rx_values.flatten()))
		
	collateral_rx = property(fget=_get_collateral_rx, fset=_set_collateral_rx)
	
	
class Cortical_Neuron_Type(NativeCellType):
	
	default_parameters = {'soma_L': 35, 'soma_diam': 25, 'soma_nseg': 1, 'soma_Ra': 150, 'soma_cm': 1, 'soma_bias_current_amp': 0.12,
							'ais_L': 20, 'ais_diam': 1.2, 'ais_nseg': 5, 'ais_Ra': 150, 'ais_cm': 0.8,
							'myelin_L': 500, 'myelin_L_0': 80, 'myelin_diam': 1.4, 'myelin_Ra': 150, 'myelin_cm': 0.04,
							'node_L': 2, 'node_diam': 1.2, 'node_nseg': 1,'node_Ra': 150, 'node_cm': 0.8,
							'collateral_L': 500, 'collateral_diam': 0.5, 'collateral_nseg': 11, 'collateral_Ra': 150, 'collateral_cm': 0.8,
							'num_axon_compartments': 10}
	
	# Define initial vector of transfer resistances for the collateral segments
	initial_collateral_rx = np.zeros((1,default_parameters['collateral_nseg'])).flatten()
	initial_collateral_rx_Sequence = Sequence(initial_collateral_rx)
	default_parameters['collateral_rx'] = initial_collateral_rx_Sequence
	
	default_initial_values = {'v': -68.0}
	recordable = ['soma(0.5).v','collateral(0.5).v', 'collateral(0.5).i_membrane_', 'ais(0.5).v', 'middle_node(0.5).v' , 'middle_myelin(0.5).v', 'AMPA.i', 'GABAa.i']
	units = {'soma(0.5).v' : 'mV', 'collateral(0.5).v': 'mV', 'collateral(0.5).i_membrane_': 'nA', 'ais(0.5).v': 'mV', 'middle_node(0.5).v': 'mV' , 'middle_myelin(0.5).v': 'mV', 'AMPA.i': 'nA', 'GABAa.i': 'nA'}    
	receptor_types = ['AMPA', 'GABAa']
	model = Cortical_Neuron

class Interneuron(object):
	
	def __init__(self, **parameters):
		
		# Create single compartment Destexhe Interneuron cell section, i.e. soma section
		self.soma = Section(L=parameters['L'], diam=parameters['diam'], nseg=parameters['nseg'], Ra=parameters['Ra'], cm=parameters['cm'],
							mechanisms=(Mechanism('interneuron_i_leak'), Mechanism('interneuron_i_na'), Mechanism('interneuron_i_k')))
		
		# Add bias current to neuron model - current amplitude is in terms of original model paper, nA 
		self.stim = h.IClamp(0.5, sec=self.soma)
		self.stim.delay = 0
		self.stim.dur = 1e12
		self.stim.amp = parameters['bias_current_amp']	# nA
		
		# insert synaptic noise 
		self.noise = h.SynNoise(0.5, sec=self.soma)
		self.noise.f0 = 0
		self.noise.f1 = 0.3
		
		# Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
		self.AMPA = h.AMPA_S(0.5, sec=self.soma)
		self.GABAa = h.GABAa_S(0.5, sec=self.soma)
		
		# needed for PyNN
		self.source_section = self.soma
		self.source = self.soma(0.5)._ref_v
		self.rec = h.NetCon(self.source, None, sec=self.source_section)		# Needed to clear the simulator
		self.spike_times = h.Vector(0)
		self.parameter_names = ('L', 'diam', 'nseg', 'Ra', 'cm', 'bias_current_amp')
		self.traces = {}
		self.recording_time = False

	L = _new_property('soma', 'L')
	diam = _new_property('soma', 'diam')
	nseg = _new_property('soma', 'nseg')
	Ra = _new_property('soma', 'Ra')
	cm = _new_property('soma', 'cm')
	bias_current_amp = _new_property('stim', 'amp')
	
	def area(self):
		"""Membrane area in µm²"""
		return pi * self.soma.L * self.soma.diam
	
	def memb_init(self):
		for seg in self.soma:
			seg.v = self.v_init

class Interneuron_Type(NativeCellType):
	default_parameters = {'L': 35, 'diam': 25, 'nseg': 1, 'Ra': 150, 'cm': 1, 'bias_current_amp': 0.25}
	default_initial_values = {'v': -68.0}
	recordable = ['soma(0.5).v']
	units = {'soma(0.5).v' : 'mV'}    
	receptor_types = ['AMPA', 'GABAa']
	model = Interneuron

class STN_Neuron(object):
	
	def __init__(self, **parameters):
		
		# Create single compartment Otsuka STN cell section, i.e. soma section
		self.soma = Section(L=parameters['L'], diam=parameters['diam'], nseg=parameters['nseg'], Ra=parameters['Ra'], cm=parameters['cm'],
							mechanisms=[Mechanism('myions'), Mechanism('stn', gnabar=49e-3, gkdrbar=57e-3, gkabar=5e-3, gkcabar=1.0e-3, gcalbar=15e-3,
														gcatbar=5e-3, kca=2, gl=0.35e-3)])
		# Initialize ion concentrations
		h("cai0_ca_ion = 5e-6 ")
		h("cao0_ca_ion = 2")
		h("ki0_k_ion = 105") 
		h("ko0_k_ion = 3")
		h("nao0_na_ion = 108")
		h("nai0_na_ion = 10")
		
		# Add bias current to neuron model
		self.stim = h.IClamp(0.5, sec=self.soma)
		self.stim.delay = 0
		self.stim.dur = 1e12
		self.stim.amp = parameters['bias_current']			# bias current density (nA)
		
		# Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
		self.AMPA = h.AMPA_S(0.5, sec=self.soma)
		self.GABAa = h.GABAa_S(0.5, sec=self.soma)
		
		
		# needed for PyNN
		self.source_section = self.soma
		self.source = self.soma(0.5)._ref_v
		self.rec = h.NetCon(self.source, None, sec=self.soma)		# Needed to clear the simulator
		self.spike_times = h.Vector(0)
		self.parameter_names = ('L', 'diam', 'nseg', 'Ra', 'cm', 'bias_current')
		self.traces = {}
		self.recording_time = False

	L = _new_property('soma', 'L')
	diam = _new_property('soma', 'diam')
	nseg = _new_property('soma', 'nseg')
	Ra = _new_property('soma', 'Ra')
	cm = _new_property('soma', 'cm')
	bias_current_amp = _new_property('stim', 'amp')
	
	def area(self):
		"""Membrane area in µm²"""
		return pi * self.soma.L * self.soma.diam
	
	def memb_init(self):
		for seg in self.soma:
			seg.v = self.v_init
	

class STN_Neuron_Type(NativeCellType):
	default_parameters = {'L': 60, 'diam': 60, 'nseg': 1, 'Ra': 200, 'cm': 1, 'bias_current': 0.0, 'num_AMPA_Synapses': 5, 'num_GABAa_Synapses': 5}
	default_initial_values = {'v': -68.0}
	recordable = ['soma(0.5).v', 'AMPA.i', 'GABAa.i']
	units = {'soma(0.5).v' : 'mV', 'AMPA.i': 'nA', 'GABAa.i': 'nA'}    
	receptor_types = ['AMPA', 'GABAa']
	model = STN_Neuron

class GP_Neuron(object):
	
	def __init__(self, **parameters):
		
		# Create single compartment Rubin and Terman GP cell section, i.e. soma section
		self.soma = Section(L=parameters['L'], diam=parameters['diam'], nseg=parameters['nseg'], Ra=parameters['Ra'], cm=parameters['cm'],
							mechanisms=[Mechanism('myions'), Mechanism('GPeA', gnabar=0.04, gkdrbar=0.0042, gkcabar=0.1e-3,
														gcatbar=6.7e-5, kca=2, gl=4e-5)])
		
		# Initialize ion concentrations
		h("cai0_ca_ion = 5e-6 ")
		h("cao0_ca_ion = 2")
		h("ki0_k_ion = 105") 
		h("ko0_k_ion = 3")
		h("nao0_na_ion = 108")
		h("nai0_na_ion = 10")
		
		# insert current source
		self.stim = h.IClamp(0.5, sec=self.soma)
		self.stim.delay = 0
		self.stim.dur = 1e12
		self.stim.amp = parameters['bias_current']
		
		# Add DBS stimulation current to neuron model
		self.DBS_stim = h.IClamp(0.5, sec=self.soma)
		self.DBS_stim.delay = 0
		self.DBS_stim.dur = 1e9
		self.DBS_stim.amp = 0
		
		# Append the DBS stimulation iclamps to global list
		GV.GPe_stimulation_iclamps.append(self.DBS_stim)
		
		# Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
		self.AMPA = h.AMPA_S(0.5, sec=self.soma)
		self.GABAa = h.GABAa_S(0.5, sec=self.soma)
		
		# needed for PyNN
		self.source_section = self.soma
		self.source = self.soma(0.5)._ref_v
		self.rec = h.NetCon(self.source, None, sec=self.soma)		# Needed to clear the simulator
		self.spike_times = h.Vector(0)
		self.parameter_names = ('L', 'diam', 'nseg', 'Ra', 'cm', 'bias_current_density')
		self.traces = {}
		self.recording_time = False
	
	L = _new_property('soma', 'L')
	diam = _new_property('soma', 'diam')
	nseg = _new_property('soma', 'nseg')
	Ra = _new_property('soma', 'Ra')
	cm = _new_property('soma', 'cm')
	bias_current = _new_property('stim', 'amp')
	
	def area(self):
		"""Membrane area in µm²"""
		return pi * self.soma.L * self.soma.diam
	
	def memb_init(self):
		for seg in self.soma:
			seg.v = self.v_init
	

class GP_Neuron_Type(NativeCellType):
	default_parameters = {'L': 60, 'diam': 60, 'nseg': 1, 'Ra': 200, 'cm': 1.0, 'bias_current': 0.03}
	default_initial_values = {'v': -68.0}
	recordable = ['soma(0.5).v']
	units = {'soma(0.5).v' : 'mV'}    
	receptor_types = ['AMPA', 'GABAa']
	model = GP_Neuron

class Thalamic_Neuron(object):
	
	def __init__(self, **parameters):
		
		# Create single compartment Rubin and Terman Thalamic cell section, i.e. soma section
		self.soma = Section(L=parameters['L'], diam=parameters['diam'], nseg=parameters['nseg'], Ra=parameters['Ra'], cm=parameters['cm'],
							mechanisms=(Mechanism('thalamic_i_leak'), Mechanism('thalamic_i_na_k'), Mechanism('thalamic_i_t')))
		
		"""
		# Note: Thalamic current has no bias current in original paper, i.e. bias_current_density = 0
		#       Can be added in if required.
		# insert current source
		self.stim = h.IClamp(0.5, sec=self.soma)
		self.stim.delay = 0
		self.stim.dur = 1e12
		self.stim.amp = parameters['bias_current_density']*(self.area())*(0.001) 	# (0.001 or 1e-3) is conversion factor so pA -> nA
		"""
		
		# insert synaptic noise 
		self.noise = h.SynNoise(0.5, sec=self.soma)
		self.noise.f0 = 0
		self.noise.f1 = 0.3
		
		# Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
		self.AMPA = h.AMPA_S(0.5, sec=self.soma)
		self.GABAa = h.GABAa_S(0.5, sec=self.soma)
		
		# needed for PyNN
		self.source_section = self.soma
		self.source = self.soma(0.5)._ref_v
		self.rec = h.NetCon(self.source, None, sec=self.soma)		# Needed to clear the simulator
		self.spike_times = h.Vector(0)
		self.parameter_names = ('L', 'diam', 'nseg', 'Ra', 'cm', 'bias_current_density')
		self.traces = {}
		self.recording_time = False

	L = _new_property('soma', 'L')
	diam = _new_property('soma', 'diam')
	nseg = _new_property('soma', 'nseg')
	Ra = _new_property('soma', 'Ra')
	cm = _new_property('soma', 'cm')
	
	def area(self):
		"""Membrane area in µm²"""
		return pi * self.soma.L * self.soma.diam
	
	def memb_init(self):
		for seg in self.soma:
			seg.v = self.v_init
	

class Thalamic_Neuron_Type(NativeCellType):
	default_parameters = {'L': 100, 'diam': 100, 'nseg': 1, 'Ra': 150, 'cm': 100, 'bias_current_amplitude': 0}
	default_initial_values = {'v': -68.0}
	recordable = ['soma(0.5).v']
	units = {'soma(0.5).v' : 'mV'}    
	receptor_types = ['AMPA', 'GABAa']
	model = Thalamic_Neuron