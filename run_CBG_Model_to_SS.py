# -*- coding: utf-8 -*-
"""
Created on Wed April 03 14:27:26 2019

Description: Cortico-Basal Ganglia Network Model implemented in PyNN using the simulator Neuron.
			This version of the model runs an initial run of the model integrating the model to
			the steady state. The steady state can then be loaded in subsequent simulations to 
			test DBS controllers. Full documentation of the model and controllers used is given in:
			 
			https://www.frontiersin.org/articles/10.3389/fnins.2020.00166/
			 
@author: John Fleming, john.fleming@ucdconnect.ie
"""

import neuron
h = neuron.h

from pyNN.neuron import setup, run, reset, run_until, run_to_steady_state, run_from_steady_state, end, simulator, Population, SpikeSourcePoisson, SpikeSourceArray, Projection, OneToOneConnector, AllToAllConnector, FromFileConnector, FixedNumberPreConnector, StaticSynapse, NativeCellType, SpikeSourcePoisson, SpikeSourceArray, NoisyCurrentSource, StepCurrentSource
from pyNN.random import RandomDistribution, NumpyRNG
from pyNN import space
from Cortical_Basal_Ganglia_Cell_Classes import Cortical_Neuron_Type, Interneuron_Type, STN_Neuron_Type, GP_Neuron_Type, Thalamic_Neuron_Type
from Electrode_Distances import distances_to_electrode, collateral_distances_to_electrode
from pyNN.parameters import Sequence
from Controllers import Constant_Controller, ON_OFF_Controller, Dual_Threshold_Controller, standard_PID_Controller
import random
import neo.io
import quantities as pq
import numpy as np
import math
from scipy import signal
import os
import sys

# Import global variables for GPe DBS
import Global_Variables as GV

def generate_poisson_spike_times(pop_size, start_time, duration, fr, timestep, random_seed):
    """ generate_population_spike_times generates (N = pop_size) poisson distributed spiketrains
        with firing rate fr.
    
	Example inputs:
		pop_size = 10			
		start_time = 0.0		# ms
		end_time = 6000.0		# ms
		timestep = 1  			# ms
		fr = 1					# Hz
	"""
	
    # Convert to sec for calculating the spikes matrix
    dt = float(timestep)/1000.0                         # sec
    tSim = float(((start_time+duration) - start_time)/1000.0)         # sec
    nBins = int(np.floor(tSim/dt))

    spikeMat = np.where(np.random.uniform(0,1,(pop_size, nBins)) < fr*dt)

    # Create time vector - ms
    tVec = np.arange(start_time, start_time+duration, timestep)

    # Make array of spike times
    for neuron_index in np.arange(pop_size):
        neuron_spike_times = tVec[spikeMat[1][np.where(spikeMat[0][:]==neuron_index)]]
        if neuron_index == 0:
            spike_times = Sequence(neuron_spike_times)
        else:
            spike_times = np.vstack((spike_times, Sequence(neuron_spike_times)))

    return spike_times
	
def generate_DBS_Signal(start_time, stop_time, dt, amplitude, frequency, pulse_width, offset):
	"""Generate monophasic square-wave DBS signal
	
	Example inputs:
		start_time = 0				# ms
		stop_time = 12000			# ms
		dt = 0.01					# ms 
		amplitude = -1.0			# mA - (amplitude<0 = cathodic stimulation, amplitude>0 = anodic stimulation)
		frequency = 130.0			# Hz
		pulse_width	= 0.06			# ms
		offset = 0					# mA
	"""
	
	times = np.round(np.arange(start_time, stop_time, dt), 2)
	tmp = np.arange(0, stop_time - start_time, dt)/1000.0
		
	# Calculate the duty cycle of the DBS signal
	T = (1.0/frequency)*1000.0	# time is in ms, so *1000 is conversion to ms
	duty_cycle = ((pulse_width)/T)
	#DBS_Signal = offset + amplitude * (1.0+signal.square(2.0 * np.pi * frequency * tmp, duty=duty_cycle))/2.0
	DBS_Signal = offset + 1.0 * (1.0+signal.square(2.0 * np.pi * frequency * tmp, duty=duty_cycle))/2.0 		# Need to initially set value > 0 to find last spike time, then can scale by amplitude
	DBS_Signal[-1] = 0.0
	
	# Calculate the time for the first pulse of the next segment
	last_pulse_index = np.where(np.diff(DBS_Signal)<0)[0][-1]
	next_pulse_time = times[last_pulse_index] + T - pulse_width
	
	# Rescale amplitude
	DBS_Signal *= amplitude
	
	return DBS_Signal, times, next_pulse_time
	
def make_beta_cheby1_filter(Fs, N, rp, low, high):
	"""Calculate bandpass filter coefficients (1st Order Chebyshev Filter)"""
	nyq = 0.5*Fs
	lowcut = low/nyq
	highcut = high/nyq

	b, a = signal.cheby1(N, rp, [lowcut, highcut], 'band')
	
	return b, a
	
def calculate_avg_beta_power(lfp_signal, tail_length, beta_b, beta_a):
	"""Calculate the average power in the beta-band for the current LFP signal window, i.e. beta Average Rectified Value (ARV)
	
	Exaqmple inputs:
		lfp_signal 			- window of LFP signal												# samples
		tail_length 		- tail length which will be discarded due to filtering artifact		# samples
		beta_b, beta_a 		- filter coefficients for filtering the beta-band from the signal	
	"""
	
	lfp_beta_signal = signal.filtfilt(beta_b, beta_a, lfp_signal)
	lfp_beta_signal_rectified = np.absolute(lfp_beta_signal)
	avg_beta_power = np.mean(lfp_beta_signal_rectified[-2*tail_length:-tail_length])
	
	return avg_beta_power
	
if __name__ == '__main__':
	
	# Setup simulation
	setup(timestep=0.01, rngseed=3695)
	steady_state_duration = 6000.0					# Duration of simulation steady state
	simulation_duration = steady_state_duration		# Total simulation time
	rec_sampling_interval = 0.5						# Fs = 2000Hz
	Pop_size = 100
	
	# Make beta band filter centred on 25Hz (cutoff frequencies are 21-29 Hz) for biomarker estimation
	beta_b, beta_a = make_beta_cheby1_filter(Fs=(1.0/rec_sampling_interval)*1000, N=4, rp=0.5, low=21, high=29)
	
	# Use CVode to calculate i_membrane_ for fast LFP calculation
	cvode = h.CVode()
	cvode.active(0)
	
	# Get the second spatial derivative (the segment current) for the collateral
	cvode.use_fast_imem(1)
	
	# Set initial values for cell membrane voltages
	v_init = -68

	# Create random distribution for cell membrane noise current
	r_init = RandomDistribution('uniform',(0, Pop_size))
	
	# Create Spaces for STN Population
	STN_Electrode_space=space.Space(axes='xy')
	STN_space = space.RandomStructure(boundary=space.Sphere(2000))				# Sphere with radius 2000um
	
	# Generate poisson distributed spike time striatal input
	striatal_spike_times =  np.load('Striatal_Spike_Times.npy')					# Load spike times from file
	
	# Generate the cortico-basal ganglia neuron populations
	Cortical_Pop = Population(Pop_size, Cortical_Neuron_Type(soma_bias_current_amp=0.245), structure=STN_space, label='Cortical Neurons') # Better than above (ibias=0.2575)
	Interneuron_Pop = Population(Pop_size, Interneuron_Type(bias_current_amp=0.070), initial_values={'v': v_init}, label='Interneurons')
	STN_Pop = Population(Pop_size, STN_Neuron_Type(bias_current=-0.125), structure=STN_space, initial_values={'v': v_init}, label='STN Neurons')
	GPe_Pop = Population(Pop_size, GP_Neuron_Type(bias_current=-0.009), initial_values={'v': v_init}, label='GPe Neurons')					# GPe/i have the same parameters, but different bias currents 
	GPi_Pop = Population(Pop_size, GP_Neuron_Type(bias_current=0.006), initial_values={'v': v_init}, label='GPi Neurons')					# GPe/i have the same parameters, but different bias currents 
	Striatal_Pop = Population(Pop_size, SpikeSourceArray(spike_times=striatal_spike_times[0][0]), label='Striatal Neuron Spike Source')	
	Thalamic_Pop = Population(Pop_size, Thalamic_Neuron_Type(), initial_values={'v': v_init}, label='Thalamic Neurons')
	
	# Load Cortical Bias currents for beta burst modulation - turn bursts off when running model to steady state
	burst_times_script = "burst_times_1.txt"
	burst_level_script = "burst_level_1.txt"
	modulation_times = np.loadtxt(burst_times_script, delimiter=',')
	modulation_signal = np.loadtxt(burst_level_script, delimiter=',')
	modulation_signal = 0.00*modulation_signal 								# Scale the modulation signal - turn off
	cortical_modulation_current = StepCurrentSource(times=modulation_times, amplitudes=modulation_signal)
	Cortical_Pop.inject(cortical_modulation_current)
	
	# Generate Noisy current sources for cortical pyramidal and interneuron populations
	Cortical_Pop_Membrane_Noise = [NoisyCurrentSource(mean=0,stdev=0.005,start=0.0, stop=simulation_duration, dt=1.0) for count in range(Pop_size)]
	Interneuron_Pop_Membrane_Noise = [NoisyCurrentSource(mean=0,stdev=0.005,start=0.0, stop=simulation_duration, dt=1.0) for count in range(Pop_size)]
	
	# Inject each membrane noise current into each cortical and interneuron in network
	for Cortical_Neuron, Cortical_Neuron_Membrane_Noise in zip(Cortical_Pop, Cortical_Pop_Membrane_Noise):
		Cortical_Neuron.inject(Cortical_Neuron_Membrane_Noise)
	
	for Interneuron, Interneuron_Membrane_Noise in zip(Interneuron_Pop, Interneuron_Pop_Membrane_Noise):
		Interneuron.inject(Interneuron_Membrane_Noise)
	
	# Update the spike times for the striatal populations
	for i in range(0,Pop_size):
		Striatal_Pop[i].spike_times=striatal_spike_times[i][0]
		
	# Load cortical positions - Comment/Remove to generate new positions
	Cortical_Neuron_xy_Positions = np.loadtxt('cortical_xy_pos.txt', delimiter=',')
	Cortical_Neuron_x_Positions = Cortical_Neuron_xy_Positions[0,:]
	Cortical_Neuron_y_Positions = Cortical_Neuron_xy_Positions[1,:]
	
	# Set cortical xy positions to those loaded in
	for cell_id, Cortical_cell in enumerate(Cortical_Pop):
		Cortical_cell.position[0] = Cortical_Neuron_x_Positions[cell_id]
		Cortical_cell.position[1] = Cortical_Neuron_y_Positions[cell_id]
	
	# Load STN positions - Comment/Remove to generate new positions
	STN_Neuron_xy_Positions = np.loadtxt('STN_xy_pos.txt', delimiter=',')
	STN_Neuron_x_Positions = STN_Neuron_xy_Positions[0,:]
	STN_Neuron_y_Positions = STN_Neuron_xy_Positions[1,:]
	
	# Set STN xy positions to those loaded in
	for cell_id, STN_cell in enumerate(STN_Pop):
		STN_cell.position[0] = STN_Neuron_x_Positions[cell_id]
		STN_cell.position[1] = STN_Neuron_y_Positions[cell_id]
		STN_cell.position[2] = 500
	
	"""
	# Position Check - 	
	# 1) Make sure cells are bounded in 4mm space in x, y coordinates
	# 2) Make sure no cells are placed inside the stimulating/recording electrode -0.5mm<x<0.5mm, -1.5mm<y<2mm
	for Cortical_cell in Cortical_Pop:
		#while(((np.abs(Cortical_cell.position[0])>2000) or ((np.abs(Cortical_cell.position[1])>2000))) or ((np.abs(Cortical_cell.position[0])<500) and (Cortical_cell.position[1]>-1500 and Cortical_cell.position[1]<2000))):
		while(((np.abs(Cortical_cell.position[0])>3000) or ((np.abs(Cortical_cell.position[1])>3000))) or ((np.abs(Cortical_cell.position[0])<700) and Cortical_cell.position[1]>-1500)):
			Cortical_cell.position = STN_space.generate_positions(1).flatten()
	
	#np.savetxt('cortical_xy_pos.txt', Cortical_Pop.positions, delimiter=',')	# Save the generated cortical xy positions to a textfile
	
	for STN_cell in STN_Pop:
		#while(((np.abs(STN_cell.position[0])>2000) or ((np.abs(STN_cell.position[1])>2000))) or ((np.abs(STN_cell.position[0])<500) and (STN_cell.position[1]>-1500 and STN_cell.position[1]<2000))):
		while(((np.abs(STN_cell.position[0])>2000) or ((np.abs(STN_cell.position[1])>2000))) or ((np.abs(STN_cell.position[0])<635) and STN_cell.position[1]>-1500)):
			STN_cell.position = STN_space.generate_positions(1).flatten()
	
	#np.savetxt('STN_xy_pos.txt', STN_Pop.positions, delimiter=',')	# Save the generated STN xy positions to a textfile
	"""
	
	# Assign Positions for recording and stimulating electrode point sources
	recording_electrode_1_position = np.array([0,-1500,250])
	recording_electrode_2_position = np.array([0,1500,250])
	stimulating_electrode_position = np.array([0,0,250])
	
	# Calculate STN cell distances to each recording electrode - only using xy coordinates for distance calculations
	STN_recording_electrode_1_distances = distances_to_electrode(recording_electrode_1_position, STN_Pop)
	STN_recording_electrode_2_distances = distances_to_electrode(recording_electrode_2_position, STN_Pop)
	
	# Calculate Cortical Collateral distances from the stimulating electrode - uses xyz coordinates for distance calculation - these distances need to be in um for xtra
	Cortical_Collateral_stimulating_electrode_distances = collateral_distances_to_electrode(stimulating_electrode_position, Cortical_Pop, L=500, nseg=11)
	#np.savetxt('cortical_collateral_electrode_distances.txt', Cortical_Collateral_stimulating_electrode_distances, delimiter=',')	# Save the generated cortical collateral stimulation electrode distances to a textfile
	
	# Synaptic Connections
	# Add variability to Cortical connections - cortical interneuron connection weights are random from uniform distribution
	gCtxInt_max_weight = 2.5e-3									# Ctx -> Int max coupling value
	gIntCtx_max_weight = 6.0e-3									# Int -> Ctx max coupling value
	gCtxInt = RandomDistribution('uniform',(0, gCtxInt_max_weight), rng=NumpyRNG(seed=3695))
	gIntCtx = RandomDistribution('uniform',(0, gIntCtx_max_weight), rng=NumpyRNG(seed=3695))
	
	# Define other synaptic connection weights and delays
	syn_CorticalAxon_Interneuron = StaticSynapse(weight=gCtxInt, delay=2)	
	syn_Interneuron_CorticalSoma = StaticSynapse(weight=gIntCtx, delay=2)
	syn_CorticalSpikeSourceCorticalAxon = StaticSynapse(weight=0.25, delay=0)
	syn_CorticalCollateralSTN = StaticSynapse(weight=0.12, delay=1)
	syn_STNGPe = StaticSynapse(weight=0.111111, delay=4)
	syn_GPeGPe = StaticSynapse(weight=0.015, delay=4)
	syn_GPeSTN = StaticSynapse(weight=0.111111, delay=3)
	syn_StriatalGPe = StaticSynapse(weight=0.01, delay=1)
	syn_STNGPi = StaticSynapse(weight=0.111111, delay=2)
	syn_GPeGPi = StaticSynapse(weight=0.111111, delay=2)
	syn_GPiThalamic = StaticSynapse(weight=3.0, delay=2)
	syn_ThalamicCortical = StaticSynapse(weight=5, delay=2)
	syn_CorticalThalamic = StaticSynapse(weight=0.0, delay=2)
	
	"""
	# Create new Synaptic Projections
	prj_CorticalAxon_Interneuron = Projection(Cortical_Pop, Interneuron_Pop,  FixedNumberPreConnector(n=10, allow_self_connections=False), syn_CorticalAxon_Interneuron, source='middle_axon_node', receptor_type='AMPA')
	prj_Interneuron_CorticalSoma = Projection(Interneuron_Pop, Cortical_Pop,  FixedNumberPreConnector(n=10, allow_self_connections=False), syn_Interneuron_CorticalSoma, receptor_type='GABAa')
	prj_CorticalSTN = Projection(Cortical_Pop, STN_Pop, FixedNumberPreConnector(n=5, allow_self_connections=False), syn_CorticalCollateralSTN, source='collateral(0.5)', receptor_type='AMPA')
	prj_STNGPe = Projection(STN_Pop, GPe_Pop, FixedNumberPreConnector(n=1, allow_self_connections=False), syn_STNGPe, source='soma(0.5)', receptor_type='AMPA')
	prj_GPeGPe = Projection(GPe_Pop, GPe_Pop, FixedNumberPreConnector(n=1, allow_self_connections=False), syn_GPeGPe, source='soma(0.5)', receptor_type='GABAa')
	prj_GPeSTN = Projection(GPe_Pop, STN_Pop, FixedNumberPreConnector(n=2, allow_self_connections=False), syn_GPeSTN, source='soma(0.5)', receptor_type='GABAa')	
	prj_StriatalGPe = Projection(Striatal_Pop, GPe_Pop, FixedNumberPreConnector(n=1, allow_self_connections=False), syn_StriatalGPe, source='soma(0.5)', receptor_type='GABAa')	
	prj_STNGPi = Projection(STN_Pop, GPi_Pop, FixedNumberPreConnector(n=1,allow_self_connections=False), syn_STNGPi, source='soma(0.5)', receptor_type='AMPA')
	prj_GPeGPi = Projection(GPe_Pop, GPi_Pop, FixedNumberPreConnector(n=1,allow_self_connections=False), syn_GPeGPi, source='soma(0.5)', receptor_type='GABAa')
	prj_GPiThalamic = Projection(GPi_Pop, Thalamic_Pop, FixedNumberPreConnector(n=1,allow_self_connections=False), syn_GPiThalamic, source='soma(0.5)', receptor_type='GABAa')
	prj_ThalamicCortical = Projection(Thalamic_Pop, Cortical_Pop, FixedNumberPreConnector(n=1,allow_self_connections=False), syn_ThalamicCortical, source='soma(0.5)', receptor_type='AMPA')
	prj_CorticalThalamic = Projection(Cortical_Pop, Thalamic_Pop, FixedNumberPreConnector(n=1,allow_self_connections=False), syn_CorticalThalamic, source='soma(0.5)', receptor_type='AMPA')
	"""
	
	# Load network topology from file
	prj_CorticalAxon_Interneuron = Projection(Cortical_Pop, Interneuron_Pop,  FromFileConnector("CorticalAxonInterneuron_Connections.txt"), syn_CorticalAxon_Interneuron, source='middle_axon_node', receptor_type='AMPA')
	prj_Interneuron_CorticalSoma = Projection(Interneuron_Pop, Cortical_Pop,  FromFileConnector("InterneuronCortical_Connections.txt"), syn_Interneuron_CorticalSoma, receptor_type='GABAa')
	prj_CorticalSTN = Projection(Cortical_Pop, STN_Pop, FromFileConnector("CorticalSTN_Connections.txt"), syn_CorticalCollateralSTN, source='collateral(0.5)', receptor_type='AMPA')
	prj_STNGPe = Projection(STN_Pop, GPe_Pop, FromFileConnector("STNGPe_Connections.txt"), syn_STNGPe, source='soma(0.5)', receptor_type='AMPA')
	prj_GPeGPe = Projection(GPe_Pop, GPe_Pop, FromFileConnector("GPeGPe_Connections.txt"), syn_GPeGPe, source='soma(0.5)', receptor_type='GABAa')
	prj_GPeSTN = Projection(GPe_Pop, STN_Pop, FromFileConnector("GPeSTN_Connections.txt"), syn_GPeSTN, source='soma(0.5)', receptor_type='GABAa')	
	prj_StriatalGPe = Projection(Striatal_Pop, GPe_Pop, FromFileConnector("StriatalGPe_Connections.txt"), syn_StriatalGPe, source='soma(0.5)', receptor_type='GABAa')	
	prj_STNGPi = Projection(STN_Pop, GPi_Pop, FromFileConnector("STNGPi_Connections.txt"), syn_STNGPi, source='soma(0.5)', receptor_type='AMPA')
	prj_GPeGPi = Projection(GPe_Pop, GPi_Pop, FromFileConnector("GPeGPi_Connections.txt"), syn_GPeGPi, source='soma(0.5)', receptor_type='GABAa')
	prj_GPiThalamic = Projection(GPi_Pop, Thalamic_Pop, FromFileConnector("GPiThalamic_Connections.txt"), syn_GPiThalamic, source='soma(0.5)', receptor_type='GABAa')
	prj_ThalamicCortical = Projection(Thalamic_Pop, Cortical_Pop, FromFileConnector("ThalamicCorticalSoma_Connections.txt"), syn_ThalamicCortical, source='soma(0.5)', receptor_type='AMPA')
	prj_CorticalThalamic = Projection(Cortical_Pop, Thalamic_Pop, FromFileConnector("CorticalSomaThalamic_Connections.txt"), syn_CorticalThalamic, source='soma(0.5)', receptor_type='AMPA')

	"""
	# Save the network topology so it can be reloaded 
	#prj_CorticalSpikeSourceCorticalSoma.saveConnections(file="CorticalSpikeSourceCorticalSoma_Connections.txt")
	prj_CorticalAxon_Interneuron.saveConnections(file="CorticalAxonInterneuron_Connections.txt")
	prj_Interneuron_CorticalSoma.saveConnections(file="InterneuronCortical_Connections.txt")
	prj_CorticalSTN.saveConnections(file="CorticalSTN_Connections.txt")
	prj_STNGPe.saveConnections(file="STNGPe_Connections.txt")
	prj_GPeGPe.saveConnections(file="GPeGPe_Connections.txt")
	prj_GPeSTN.saveConnections(file="GPeSTN_Connections.txt")
	prj_StriatalGPe.saveConnections(file="StriatalGPe_Connections.txt")
	prj_STNGPi.saveConnections(file="STNGPi_Connections.txt")
	prj_GPeGPi.saveConnections(file="GPeGPi_Connections.txt")
	prj_GPiThalamic.saveConnections(file="GPiThalamic_Connections.txt")
	prj_ThalamicCortical.saveConnections(file="ThalamicCorticalSoma_Connections.txt")
	prj_CorticalThalamic.saveConnections(file="CorticalSomaThalamic_Connections.txt")
	"""
	
	# Define state variables to record from each population
	Cortical_Pop.record('soma(0.5).v', sampling_interval=rec_sampling_interval)
	Cortical_Pop.record('collateral(0.5).v', sampling_interval=rec_sampling_interval)
	Interneuron_Pop.record('soma(0.5).v', sampling_interval=rec_sampling_interval)
	STN_Pop.record('soma(0.5).v', sampling_interval=rec_sampling_interval)
	STN_Pop.record('AMPA.i', sampling_interval=rec_sampling_interval)
	STN_Pop.record('GABAa.i', sampling_interval=rec_sampling_interval)
	Striatal_Pop.record('spikes')
	GPe_Pop.record('soma(0.5).v', sampling_interval=rec_sampling_interval)
	GPi_Pop.record('soma(0.5).v', sampling_interval=rec_sampling_interval)
	Thalamic_Pop.record('soma(0.5).v', sampling_interval=rec_sampling_interval)
		
	# Conductivity and resistivity values for homogenous, isotropic medium
	sigma = 0.27							# Latikka et al. 2001 - Conductivity of Brain tissue S/m
	rho = (1/(sigma*1e-2))					# rho needs units of ohm cm for xtra mechanism (S/m -> S/cm)
	
	# Calculate transfer resistances for each collateral segment for xtra - units are Mohms
	collateral_rx = (rho/(4*math.pi))*(1/Cortical_Collateral_stimulating_electrode_distances)*(0.01)
	
	# Convert ndarray to array of Sequence objects - needed to set cortical collateral transfer resistances
	collateral_rx_seq = np.ndarray(shape=(1, Pop_size), dtype=Sequence).flatten()
	for ii in range(0, Pop_size):
		collateral_rx_seq[ii] = Sequence(collateral_rx[ii,:].flatten())
	
	# Assign transfer resistances values to collaterals
	for ii, cortical_neuron in enumerate(Cortical_Pop):
		cortical_neuron.collateral_rx = collateral_rx_seq[ii]
	
	# Initialise STN LFP list
	STN_LFP = []
	STN_LFP_AMPA = []
	STN_LFP_GABAa = []
	
	# Variables for writing simulation data
	last_write_time = 0
	
	DBS_Signal_neuron = h.Vector([0, 0])
	DBS_times_neuron = h.Vector([0, steady_state_duration+10])
	
	# Play DBS signal to global variable is_xtra
	DBS_Signal_neuron.play(h._ref_is_xtra, DBS_times_neuron, 1)
	
	# Make new GPe DBS vector for each GPe neuron - each GPe neuron needs a pointer to it's own DBS signal
	GPe_DBS_Signal_neuron = []
	GPe_DBS_times_neuron = []
	updated_GPe_DBS_signal = []
	for i in range(0, Pop_size):
	
		# Neuron vector of GPe DBS signals
		GPe_DBS_Signal_neuron.append(h.Vector([0, 0]))
		GPe_DBS_times_neuron.append(h.Vector([0, steady_state_duration+10]))
		
		# Play the stimulation into eacb GPe neuron
		GPe_DBS_Signal_neuron[i].play(GV.GPe_stimulation_iclamps[i]._ref_amp, GPe_DBS_times_neuron[i], 1)
	
	# Run the model to the steady state
	run_to_steady_state(steady_state_duration)
	
	# Calculate the LFP and biomarkers, etc. 
	STN_AMPA_i = np.array(STN_Pop.get_data('AMPA.i').segments[0].analogsignals[0])
	STN_GABAa_i = np.array(STN_Pop.get_data('GABAa.i').segments[0].analogsignals[0])
	STN_Syn_i = STN_AMPA_i + STN_GABAa_i
			
	# STN LFP Calculation - Syn_i is in units of nA -> LFP units are mV
	STN_LFP_1 = (1/(4*math.pi*sigma))*np.sum((1/(STN_recording_electrode_1_distances*1e-6))*STN_Syn_i.transpose(),axis=0)*1e-6 
	STN_LFP_2 = (1/(4*math.pi*sigma))*np.sum((1/(STN_recording_electrode_2_distances*1e-6))*STN_Syn_i.transpose(),axis=0)*1e-6 
	STN_LFP = STN_LFP_1 - STN_LFP_2
	
	# STN LFP AMPA and GABAa Contributions
	STN_LFP_AMPA_1 = (1/(4*math.pi*sigma))*np.sum((1/(STN_recording_electrode_1_distances*1e-6))*STN_AMPA_i.transpose(),axis=0)*1e-6 
	STN_LFP_AMPA_2 = (1/(4*math.pi*sigma))*np.sum((1/(STN_recording_electrode_2_distances*1e-6))*STN_AMPA_i.transpose(),axis=0)*1e-6 
	STN_LFP_AMPA = STN_LFP_AMPA_1 - STN_LFP_AMPA_2
	STN_LFP_GABAa_1 = (1/(4*math.pi*sigma))*np.sum((1/(STN_recording_electrode_1_distances*1e-6))*STN_GABAa_i.transpose(),axis=0)*1e-6 
	STN_LFP_GABAa_2 = (1/(4*math.pi*sigma))*np.sum((1/(STN_recording_electrode_2_distances*1e-6))*STN_GABAa_i.transpose(),axis=0)*1e-6 
	STN_LFP_GABAa = STN_LFP_GABAa_1 - STN_LFP_GABAa_2
	
	"""
	# Simulation Label for writing model output data - uncomment to write the specified variables to file
	simulation_label = "Steady_State_Simulation"
	
	# Write population membrane voltage data to file
	Cortical_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/Cortical_Pop/Cortical_Collateral_v.mat", 'collateral(0.5).v', clear=False)
	Cortical_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/Cortical_Pop/Cortical_Soma_v.mat", 'soma(0.5).v', clear=True)
	Interneuron_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/Interneuron_Pop/Interneuron_Soma_v.mat", 'soma(0.5).v', clear=True)
	STN_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/STN_Pop/STN_Soma_v.mat", 'soma(0.5).v', clear=True)
	GPe_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/GPe_Pop/GPe_Soma_v.mat", 'soma(0.5).v', clear=True)
	GPi_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/GPi_Pop/GPi_Soma_v.mat", 'soma(0.5).v', clear=True)
	Thalamic_Pop.write_data("/Simulation_Output_Results/"+simulation_label+"/Thalamic_Pop/Thalamic_Soma_v.mat", 'soma(0.5).v', clear=True)
	
	# Write the STN LFP to .mat file
	STN_LFP_Block = neo.Block(name='STN_LFP')
	STN_LFP_seg = neo.Segment(name='segment_0')
	STN_LFP_Block.segments.append(STN_LFP_seg)
	STN_LFP_signal = neo.AnalogSignal(STN_LFP, units='mV', t_start=0*pq.ms, sampling_rate=pq.Quantity(simulator.state.dt, '1/ms'))
	STN_LFP_seg.analogsignals.append(STN_LFP_signal)
	
	w = neo.io.NeoMatlabIO(filename="/Simulation_Output_Results/"+simulation_label+"/STN_LFP.mat") 
	w.write_block(STN_LFP_Block)
	
	# Write LFP AMPA and GABAa conmponents to file
	STN_LFP_AMPA_Block = neo.Block(name='STN_LFP_AMPA')
	STN_LFP_AMPA_seg = neo.Segment(name='segment_0')
	STN_LFP_AMPA_Block.segments.append(STN_LFP_AMPA_seg)
	STN_LFP_AMPA_signal = neo.AnalogSignal(STN_LFP_AMPA, units='mV', t_start=0*pq.ms, sampling_rate=pq.Quantity(simulator.state.dt, '1/ms'))
	STN_LFP_AMPA_seg.analogsignals.append(STN_LFP_AMPA_signal)
	w = neo.io.NeoMatlabIO(filename="/Simulation_Output_Results/"+simulation_label+"/STN_LFP_AMPA.mat") 
	w.write_block(STN_LFP_AMPA_Block)
	
	STN_LFP_GABAa_Block = neo.Block(name='STN_LFP_GABAa')
	STN_LFP_GABAa_seg = neo.Segment(name='segment_0')
	STN_LFP_GABAa_Block.segments.append(STN_LFP_GABAa_seg)
	STN_LFP_GABAa_signal = neo.AnalogSignal(STN_LFP_GABAa, units='mV', t_start=0*pq.ms, sampling_rate=pq.Quantity(simulator.state.dt, '1/ms'))
	STN_LFP_GABAa_seg.analogsignals.append(STN_LFP_GABAa_signal)
	w = neo.io.NeoMatlabIO(filename="/Simulation_Output_Results/"+simulation_label+"/STN_LFP_GABAa.mat") 
	w.write_block(STN_LFP_GABAa_Block)
	"""
	
	print("Steady State Simulation Done!")
	
	end()
