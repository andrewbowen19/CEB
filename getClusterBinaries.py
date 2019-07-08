"""Class of getClusterBinaries that runs our HS-period finding code, uses calculated sigma values
Takes in these params: mass, Rhm, age, metallicity, velocity dispersion, and number of binaries requested
should calculate hard-soft binary for each binary drawn with cosmic and return # of binaries requested on input
output of this should be a numpy array of arrays
"""

# Importing needed mdoules
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import independent
from cosmic.sample.sampler import multidim
from cosmic.evolve import Evolve
import numpy as np
import pandas as pd

# BSE dictionary copied from cosmic's documentation: https://cosmic-popsynth.github.io
BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, \
'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, \
'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, \
'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'ppsn': 1,\
 'natal_kick_array' : [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,\
  'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 0, \
  'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsnp' : 2.5, 'ecsn_mlow' : 1.6, 'aic' : 1, 'sigmadiv' :-20.0}

# Class that samples hard binaries for every cluster and takes in outputs from Aaron
class getClusterBinaries:

	def __init__(self, mass, Rhm, age, Z, sigma, Nbin):
		# Input data from clusters
		self.mass = mass
		self.Rhm = Rhm
		self.age = age
		self.Z = Z
		self.sigma = sigma
		self.Nbin = Nbin

	@staticmethod
	def get_Phs(row):
	# """Function to calculate the hard-soft period cutoff given by Eq 1 in Geller, Leigh 2015
    	# Masses are in solar masses (converted later to kg), m1 and m2 are binary component masses,
    	# and m3 is the mass of the incoming (disrupting) object, velocity dispersions are given in km/s"""
		G = 1.334 * (10 ** 11) # Gravitational Constant in units of km^3 M_sun ^ -1 s ^ -2 (consistent with cosmic output units)

		const = (np.pi*G/np.sqrt(2))
		m1 = row['mass1_binary']
		m2 = row['mass2_binary']
		m3 = row['m3']
		sigma = row['sigma']

		Phs = const*(np.sqrt(m1*m2/0.5)**3)*((np.sqrt(3)*0.8)**(-3))/(np.sqrt(m1+m2))
		return Phs

	# New sampler function - only pulls initial binaries
	@classmethod#Setting this up as a class method so it doesn't need the self instance
	def Initial_Binary_Sample(cls, age, Nbin, Z, random_seed):
#	"""Creates and evolves a set of binaries with given 
#	age (to evolve to), number of binaries, metallicity, and velocity dispersion. """

		
		# Initial (input) binares -- using sampler method from cosmic #1234 - random seed
		InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('multidim',\
		 [11], [11],random_seed, 1, 'delta_burst', age, Z, Nbin)# porb_hi = [3])



		# Inclination and omega values
		inc = np.arccos(2.*np.random.uniform(0,1,Nbin) - 1.)
		omega = np.random.uniform(0,2*np.pi,Nbin)
		OMEGA = np.random.uniform(0,2*np.pi,Nbin)

		# Input binary params
		p_i = InitialBinaries['porb']
		m1_i = InitialBinaries['mass1_binary']#In Msun units
		m2_i = InitialBinaries['mass2_binary']
		ecc_i = InitialBinaries['ecc']
		tphysf_i = InitialBinaries['tphysf']#Time to evolve to (Myr)
		Z_i = InitialBinaries['metallicity']
		
# Applying HS formula to every binary/row in Initial Binary Table
		InitialBinaries['m3'] = 0.5 * np.ones(len(InitialBinaries))#Setting up m3 values (0.5 Msun)
		InitialBinaries['sigma'] = 0.8 * np.ones(len(InitialBinaries))#Sample velocity dispersion, will use cluster sigma values

		InitialBinaries['HS_Cutoff'] = InitialBinaries.apply(cls.get_Phs, axis = 1) #args = (row['mass1_binary'], row['mass2_binary'], 0.5, 1))
		InitialBinaries['HS_Cutoff'] = InitialBinaries['HS_Cutoff'] / (24 *3600) #converting from days to seconds
		HS_Cutoff = InitialBinaries['HS_Cutoff']

		 # Throwing out soft binaries
		InitialBinaries = InitialBinaries.loc[np.where(p_i < HS_Cutoff)]#Throwing out soft binaries

		return(InitialBinaries)

	# Evolving hard binaries
	def EvolveBinaries(cls):
		bpp, bcm, initC  = Evolve.evolve(initialbinarytable = cls.InitialBinaries, BSEDict = BSEDict)






"""def get_Phs(m1, m2, m3, sigma):
    '''Function to calculate the hard-soft period cutoff given by Eq 1 in Geller, Leigh 2015
    Masses are in solar masses (converted later to kg), m1 and m2 are binary component masses,
    and m3 is the mass of the incoming (disrupting) object, 
    velocity dispersions are given in km/s'''
    G = 6.67408*10**-20#Gravitational constant in km^3 kg^-1 s^-2
    Msun = 1.989*10**30#Solar mass in kg
    const = (np.pi*G/np.sqrt(2)) * Msun#Last term used to convert from solar masses to kg
    Phs = const*(np.sqrt(m1*m2/m3)**3)*((np.sqrt(3)*sigma)**(-3))/(np.sqrt(m1+m2))
    return Phs

"""

# Test variable
MyBinaries = getClusterBinaries.Initial_Binary_Sample(13000, 10, 0.01, 12)
print(MyBinaries)





	# # Evolving input binaries
	# 	bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

	# 	# Making returned variables global
	# 	global p_f, m1_f, m2_f, ecc_f, tphysf_f, r1_f, r2_f, lum1_f, lum2_f, Teff1, Teff2
	# 	# Evolved Binary Params (made global for plotting later, can )
	# 	p_f = bcm['porb']
	# 	m1_f = bcm['mass_1']
	# 	m2_f = bcm['mass_2']
	# 	ecc_f = bcm['ecc']
	# 	tphysf_f = bcm['tphys']
	# 	r1_f = bcm['rad_1']
	# 	r2_f = bcm['rad_2']
	# 	lum1_f = bcm['lumin_1']
	# 	lum2_f = bcm['lumin_2']
	# 	Teff1 = bcm['teff_1']#Effective temperature - just in case
	# 	Teff2 = bcm['teff_2']





