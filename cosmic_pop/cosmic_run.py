"""Creating a function to run on cosmic (can input our datatables)
	Want to take data from our files (GC and OC files) and be abel to imput those params (real cluster data)
	into cosmic, and evolve a binary population for each cluster (will likely use the same # of binaries 
	for each cluster)

	will also need to define a hard-soft boundary for each population of binaries
"""

from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import independent
from cosmic.sample.sampler import multidim
from cosmic.evolve import Evolve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, \
'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, \
'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, \
'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'ppsn': 1,\
 'natal_kick_array' : [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,\
  'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 0, \
  'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsnp' : 2.5, 'ecsn_mlow' : 1.6, 'aic' : 1, 'sigmadiv' :-20.0}

# Creating function to run cosmic like in bin_cosmic_test.py, inpust are cluster gae and number of binaries
def cosmic_run(cluster_age, Nbin):
	n_grid = Nbin
	m1 = np.random.uniform(0.5, 2, 1000)
	m2 = np.random.uniform(0.5, 2, 1000)
	metal = np.ones(n_grid)
	porb = np.random.uniform(0,1000, 1000)
	ecc = np.random.random(1000)
	k1 = np.ones(n_grid)
	k2 = np.ones(n_grid)
	t_p = np.ones(n_grid)*cluster_age

	# Initial grid of binaries (input parameters)
	binary_grid = InitialBinaryTable.MultipleBinary(m1=m1, m2=m2, \
	porb=porb, ecc=ecc, tphysf=t_p,\
	 kstar1=k1, kstar2=k2, metallicity=metal)
	# Setting up variables drawn from initial binary table (binary_grid)
	p_i = binary_grid['porb']
	m1_i = binary_grid['mass1_binary']
	m2_i = binary_grid['mass2_binary']
	ecc_i = binary_grid['ecc']
	tphysf_i = binary_grid['tphysf']
	Z_i = binary_grid['metallicity']

	# Using this class to evolve the binaries taken from the initial binary grid
	bpp, bcm, initC  = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

	# Setting up variables from evolved binaries (bcm - certain times)
	p_f = bcm['porb']
	m1_f = bcm['mass_1']
	m2_f = bcm['mass_2']
	ecc_f = bcm['ecc']
	tphysf_f = bcm['tphys']
	r1_f = bcm['rad_1']
	r2_f = bcm['rad_2']

	return bcm


# Now we read in in the data files for OCs and GCs

# Names in GC data file
names_gc = ['ID_x', 'Name', 'RA', 'DEC', 'L','B','R_Sun','R_gc','X','Y', 'Z', 'key_0','[Fe/H]_x', 'wt', 'E(B-V)_x',\
 'V_HB','(m-M)V_x', 'V_t', 'M_V,t', 'U-B', 'B-V', 'V-R', 'V-I', 'spt', 'ellip', 'ID_y', 'v_r', '+/-', 'v_LSR' ,'sig_v' ,'+/-.1', 'c', 'r_c', 'r_h', 'mu_V',\
  'rho_', 'lg(tc)', 'lg(th)', 'Mcl[Msun]', 'rh[pc]', '[Fe/H]_y', 'age[Gyr]', '(m-M)V_y', 'E(B-V)_y', 'log10(rho[Msun]/pc^3)',\
 'rc', 'sigma0[km/s]', 'esigma0[km/s]', 'fb', 'efb', '[M/H]', 'Rgc[kpc]','Rsun[kpc]']

# names in OC datafile
names_oc = ['Cluster_name', 'RA', 'DEC', 'l', 'b', 'Dist Mod', 'EB-V', 'Age', 'ST' ,'Z', 'Diam', 'Fe/H', 'MRV',\
 'pm RA', 'pm Dec', 'logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]', 'rtK[pc]', 'Rhm[pc]',\
  '[Fe/H]K]', 'deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt' ,'Rgc' ,'z' ,'Diam[pc]', 'd[pc]']
path = '/Users/andrewbowen/ceb_project/cosmic_pop/'

# Globular cluster read in
GCs = pd.read_csv(path + 'gc_data.txt', sep = ' ', header = 0, names = names_gc)

# Galactic coords (kiloparsecs: 10^3pc)
xGX = GCs['X']
yGY = GCs['Y']
zGZ = GCs['Z']

#Distance from sun (kiloparsecs)
gc_dist_kpc = GCs['R_Sun']

# Open Cluster file read in
OCs = pd.read_csv(path + 'oc_data.txt', sep = ' ', header = 0)

oc_dist = OCs['d[pc]']
oc_dist_kpc = oc_dist / 1000
oc_age = OCs['Age']#Webda age for OCs (Check units!!!!!)


# Adding Binary # column to dataframes - for now using 40k binaries per cluster
GCs['Nbin'] = pd.Series()
OCs['Nbin'] = pd.Series()
gc_nbin = GCs['Nbin']
oc_nbin = OCs['Nbin']

nbins = np.ones(len(GCs))*10
GCs['Nbin'] = nbins
nbins = np.ones(len(OCs))*10
OCs['Nbin'] = nbins


# InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('multidim', [11], [11], 2, 1, 'delta_burst', 10000.0, 0.02, 10)

# ############################################# Sampler Method ########################################################
# Using InitialBinaryTable.sampler method (from cosmic's documentation)
# Eventually will use cosmic_run/cosmic_sampler function on each cluster (inputting age/Nbin/Metallicity/velo. dispersion)
# Cosmic output period is in years (*365)
# tphys = np.ones(200)*13000

# sigma - velocity dispersion of a cluster 

# Function with same purpose as above - with sampler class from InitialBinaries - different sampling
def cosmic_sampler(age, Nbin, Z, sigma):
	"""Creates and evolves a set of binaries with given 
	age (to evolve to), number of binaries, metallicity, and velocity dispersion. 
	Will later loop through globular and open clusters and apply this for loop"""
	n_grid = Nbin

	# Initial (input) binares -- using sampler method from cosmic
	InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('multidim',\
	 [11], [11], 2, 1, 'delta_burst', age, Z, Nbin)

	# Inclination and omega values
	inc = np.arccos(2.*np.random.uniform(0,1,Nbin) - 1.)
	omega = np.random.uniform(0,2*np.pi,Nbin)
	OMEGA = np.random.uniform(0,2*np.pi,Nbin)
	print(InitialBinaries)
# Making Input variables global (for plotting later)
	global p_i, m1_i, m2_i, ecc_i, tphysf_i, Z_i
	# Input binary params (for plotting later)
	p_i = InitialBinaries['porb']
	m1_i = InitialBinaries['mass1_binary']
	m2_i = InitialBinaries['mass2_binary']
	ecc_i = InitialBinaries['ecc']
	tphysf_i = InitialBinaries['tphysf']
	Z_i = InitialBinaries['metallicity']
	
	# Evolving input binaries
	bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

	# Making returned variables global
	global p_f, m1_f, m2_f, ecc_f, tphysf_f, r1_f, r2_f, lum1_f, lum2_f, Teff1, Teff2
	# Evolved Binary Params (made global for plotting later, can )
	p_f = bcm['porb']
	m1_f = bcm['mass_1']
	m2_f = bcm['mass_2']
	ecc_f = bcm['ecc']
	tphysf_f = bcm['tphys']
	r1_f = bcm['rad_1']
	r2_f = bcm['rad_2']
	lum1_f = bcm['lumin_1']
	lum2_f = bcm['lumin_2']
	Teff1 = bcm['teff_1']#Effective temperature - just in case
	Teff2 = bcm['teff_2']


	return bcm


# Define from cluster params: xGX, yGX, zGX , dist_kpc


# Get these from Evolve output: m1, m2, logp, ecc, rad1, rad2, Lum1, Lum2, xGX, yGX, zGX, dist_kpc, inc, OMEGA, omega (check email from Aaron if you need definitions)

# bpp: contains binary parameters at important stages in the binary’s evolution
# bcm: contains several binary parameters at user specified time steps during the binary’s evolution.


# ################################################## Cluster Loops ############################################################
# Want to loop through cluster dataframes and use cosmic run for each cluster

print('Going into OC for loop...')

# Looping through open clusters 
for index, row in OCs.iterrows():
	"""For loop that will run our cosmic_sampler function for every open cluster in our data frame.
	Each OC data source (Webda, Pisk, Solaris) has age and metallicity data, so we're using that
	for the input parameters of our function"""

	# Goes through webda age - metallicity combos to perfomr cosmic_sampler on each cluster
	if np.isfinite(row['Age']):

		oc_age = row['Age'] * (10**9)
		oc_age = int(oc_age)

		if row['Fe/H'] != -99.99:
			print(index)
			OC_run = cosmic_sampler(oc_age, int(row['Nbin']), row['Fe/H'], 1)
			print('Webda')
			print(OC_run)

	# Goes through Solaris age-metallicity combos to perform cosmic_sampler on each Solaris cluster
	elif np.isfinite(row['t']):
		print(index)
		OC_run = cosmic_sampler(row['t'], int(row['Nbin']), row['[FeH]'], 1)
		print('Solaris')
		print(OC_run)

	# Goes through Piskunov age-metallicity combos ... "          "
	elif np.isfinite(row['log(t[yr])K']):
		print(index)
		OC_run = cosmic_sampler(row['log(t[yr])K'], int(row['Nbin']), row['[Fe/H]K'], 1)
		print('Piskunov')
		print(OC_run)
	# else:
	# 	pass


print('Done with Open Clusters')

# Globular Cluster for loop - applying cosmic_sampler to each cluster (if there is no velocity dispersion given we just assume 1km/s, can change this value later)
for index, row in GCs.iterrows():
	sigma_v = row['sig_v']
	Z = row['[Fe/H]_x']
	Nbin = row['Nbin']
	gc_age = (row['age[Gyr]'] * (10**9))
	gc_age = int(gc_age)

	# Verifying there is a velodicty dispersion parameter - if not we'll give the function one to use (1 in this case)
	if np.isfinite(row['sig_v']):
		print(row['ID_x'])
		GC_run = cosmic_sampler(gc_age, int(Nbin), Z, int(sigma_v))
	# else:
	# 	GC_run = cosmic_sampler(gc_age, int(Nbin), Z, 1)
	print(GC_run)

print('Done with GC loop!')





# ######################################################## Plotting ####################################################################

# Creating function instance to use for plotting (test function with arbitrary params)
EvolvedBinaries = cosmic_sampler(10000, 100, 0.02, 8)

# Histograms (aligned)
f, axarr = plt.subplots(2,2, figsize = (10,6), sharex = False)
f.subplots_adjust(wspace=1, hspace = 1)

# Period Histogram - initial
axarr[0,0].hist(p_i, bins = 50, range = (0,100))
axarr[0,0].set_xlabel('Period (log-days)')
axarr[0,0].set_title('Initial Binary Population periods', fontsize = 12)

# Mass 1 Histogram - initial
axarr[0,1].hist(m1_i, bins = 50, range = (0,15))
axarr[0,1].set_xlabel('M1 (Solar Masses)', fontsize = 12)
axarr[0,1].set_title('Initial Binary Population mass 1')

# Mass 2 Histogram - initial
axarr[1,0].hist(m2_i, bins = 50, range = (0,15))
axarr[1,0].set_xlabel('M2 (Solar Masses)', fontsize = 12)
axarr[1,0].set_title('Initial Binary Population mass 2')

# Eccentricity Histogram - initial
axarr[1,1].hist(ecc_i, bins = 50, range = (0,1))
axarr[1,1].set_xlabel('Eccentricity', fontsize = 12)
axarr[1,1].set_title('Initial Binary Population Eccentricity')


f, axarr = plt.subplots(2,2, figsize = (10,6), sharex = False)
f.subplots_adjust(wspace=2, hspace =0.5)

# Period Histogram - evolved
axarr[0,0].hist(EvolvedBinaries['porb'], bins = 50, range = (0,3), color = 'orange')
axarr[0,0].set_xlabel('Period (log-seconds)', fontsize = 12)
axarr[0,0].set_title('Evolved Binary Population periods')

# Mass 1 Histogram - evolved
axarr[0,1].hist(EvolvedBinaries['mass_1'], bins = 50, range = (0.5,2), color = 'orange')
axarr[0,1].set_xlabel('M 1 (Solar Masses)', fontsize = 12)
axarr[0,1].set_title('Evolved Binary Population mass 1')

# Mass 2 Histogram - evolved
axarr[1,0].hist(EvolvedBinaries['mass_2'], bins = 50, range = (0.5,2), color = 'orange')
axarr[1,0].set_xlabel('M 2 (Solar Masses))', fontsize = 12)
axarr[1,0].set_title('Evolved Binary Population mass 2')

# Eccentricity Histogram - evolved
axarr[1,1].hist(EvolvedBinaries['ecc'], bins = 50, range = (0,1), color = 'orange')
axarr[1,1].set_xlabel('Eccentricity', fontsize = 12)
axarr[1,1].set_title('Evolved Binary Population Eccentricity')

# Radii Histograms - evolved
f, ax = plt.subplots(figsize = (6,5))
ax.hist(EvolvedBinaries['rad_1'], bins = 50, range = (0,5), color = 'orange')
ax.set_xlabel('Radii 1')
ax.set_title('Evolved Binary Population Radius 1')

f, ax = plt.subplots(figsize = (6,5))
ax.hist(EvolvedBinaries['rad_2'], bins = 50, range = (0,5), color = 'orange')
ax.set_xlabel('Radii 2')
ax.set_title('Evolved Binary Population Radius 2')

# Luminosity Histograms - evolved
f, ax = plt.subplots(figsize = (6,5))
ax.hist(EvolvedBinaries['lumin_1'], bins = 50, range = (0,5), color = 'orange')
ax.set_xlabel('Luminosity 1')
ax.set_title('Evolved Binary Population Luminosity 1')

f, ax = plt.subplots(figsize = (6,5))
ax.hist(EvolvedBinaries['lumin_2'], bins = 50, range = (0,5), color = 'orange')
ax.set_xlabel('Luminosity 2')
ax.set_title('Evolved Binary Population Luminosity 2')


# plt.show()
