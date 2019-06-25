# test script for cosmic runs - let's generate some binaries!

# https://cosmic-popsynth.github.io   --- documentation

# Importing needed libraries
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import independent
from cosmic.evolve import Evolve
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('MacOSX')
# ##############################################################################################################################

# BSEDict pulled from cosmic's documentation: likely won't change these parameters
BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, \
'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, \
'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, \
'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'ppsn': 1,\
 'natal_kick_array' : [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,\
  'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 0, \
  'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsnp' : 2.5, 'ecsn_mlow' : 1.6, 'aic' : 1, 'sigmadiv' :-20.0}



# Setting up input parameters
n_grid = 1000
m1 = np.random.uniform(0.5, 2, 1000)
m2 = np.random.uniform(0.5, 2, 1000)
metal = np.ones(n_grid)
porb = np.random.uniform(0,1000, 1000)
ecc = np.random.random(1000)
k1 = np.ones(n_grid)
k2 = np.ones(n_grid)
t_p = np.ones(n_grid)*13500.0

# Generating BSE test - gives random binary properties
binary_grid = InitialBinaryTable.MultipleBinary(m1=m1, m2=m2, \
	porb=porb, ecc=ecc, tphysf=t_p,\
	 kstar1=k1, kstar2=k2, metallicity=metal)

 # Can print initial binaries if needed
# print('Initial Binaries:')
# print(binary_grid)
# print(' ')

# Setting up variables drawn from initial binary table (binary_grid)
p_i = binary_grid['porb']
m1_i = binary_grid['mass1_binary']
m2_i = binary_grid['mass2_binary']
ecc_i = binary_grid['ecc']
tphysf_i = binary_grid['tphysf']
Z_i = binary_grid['metallicity']

# Using this class to evolve the binaries taken from the initial binary grid
bpp, bcm, initC  = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

# Printing Evolved binaries
# print('Evolved Binaries:')
# print(bpp, bcm, initC)
print(np.size(bcm))

# Setting up variables from evolved binaries (bcm - certain times)
p_f = bcm['porb']
m1_f = bcm['mass_1']
m2_f = bcm['mass_2']
ecc_f = bcm['ecc']
tphysf_f = bcm['tphys']
r1_f = bcm['rad_1']
r2_f = bcm['rad_2']
print(bcm['radc_1'])
# Z_f = bpp['metallicity']


# ################################################## Plotting #################################################################################

# Histograms of initial and final distributions for: period, mass1, mass1, ecc, etc.

# Histograms (aligned)
f, axarr = plt.subplots(2,2, figsize = (10,6), sharex = False)
f.subplots_adjust(wspace=1, hspace = 1)

# Period Histogram - initial
axarr[0,0].hist(p_i, bins = 50, range = (0,100))
axarr[0,0].set_xlabel('Period (days)')
axarr[0,0].set_title('Initial Binary Population periods', fontsize = 12)

# Mass 1 Histogram - initial
axarr[0,1].hist(m1_i, bins = 50, range = (0.5,2))
axarr[0,1].set_xlabel('M1 (Solar Masses)', fontsize = 12)
axarr[0,1].set_title('Initial Binary Population mass 1')

# Mass 2 Histogram - initial
axarr[1,0].hist(m2_i, bins = 50, range = (0.5,2))
axarr[1,0].set_xlabel('M2 (Solar Masses)', fontsize = 12)
axarr[1,0].set_title('Initial Binary Population mass 2')

# Eccentricity Histogram - initial
axarr[1,1].hist(ecc_i, bins = 50, range = (0,1))
axarr[1,1].set_xlabel('Eccentricity', fontsize = 12)
axarr[1,1].set_title('Initial Binary Population Eccentricity')



#### Evolved Histrograms ####

# Evolved Binaries Histograms (aligned)
f, axarr = plt.subplots(2,2, figsize = (10,6), sharex = False)
f.subplots_adjust(wspace=1, hspace =1)

# Period Histogram - evolved
axarr[0,0].hist(p_f, bins = 50, range = (0,3), color = 'orange')
axarr[0,0].set_xlabel('Period (log-seconds)', fontsize = 12)
axarr[0,0].set_title('Evolved Binary Population periods')

# Mass 1 Histogram - evolved
axarr[0,1].hist(m1_f, bins = 50, range = (0.5,2), color = 'orange')
axarr[0,1].set_xlabel('M 1 (Solar Masses)', fontsize = 12)
axarr[0,1].set_title('Evolved Binary Population mass 1')

# Mass 2 Histogram - evolved
axarr[1,0].hist(m2_f, bins = 50, range = (0.5,2), color = 'orange')
axarr[1,0].set_xlabel('M 2 (Solar Masses))', fontsize = 12)
axarr[1,0].set_title('Evolved Binary Population mass 2')

# Eccentricity Histogram - evolved
axarr[1,1].hist(ecc_f, bins = 50, range = (0,1), color = 'orange')
axarr[1,1].set_xlabel('Eccentricity', fontsize = 12)
axarr[1,1].set_title('Evolved Binary Population Eccentricity')


# Radius 1 histogram (evolved stars)
f,ax = plt.subplots(figsize = (6,5))
ax.hist(r1_f, bins = 50, range = (0,2), color = 'red')
ax.set_title('R_1 (Solar radii)')

# Radius 2 histogram (evolved stars)
f,ax = plt.subplots(figsize = (6,5))
ax.hist(r2_f, bins = 50, range = (0,2), color = 'red')
ax.set_title('R_2 (Solar radii)')


plt.show()


