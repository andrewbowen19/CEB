# Trying to get central velocity dispersion for both Globular and open clusters

# Check units

import pandas as pd
import numpy as np
import math

# Names in GC data file
names_gc = ['ID_x', 'Name', 'RA', 'DEC', 'L','B','R_Sun','R_gc','X','Y', 'Z', 'key_0','[Fe/H]_x', 'wt', 'E(B-V)_x',\
 'V_HB','(m-M)V_x', 'V_t', 'M_V,t', 'U-B', 'B-V', 'V-R', 'V-I', 'spt', 'ellip', 'ID_y', 'v_r', '+/-', 'v_LSR' ,'sig_v' ,'+/-.1', 'c', 'r_c', 'r_h', 'mu_V',\
  'rho_', 'lg(tc)', 'lg(th)', 'Mcl[Msun]', 'rh[pc]', '[Fe/H]_y', 'age[Gyr]', '(m-M)V_y', 'E(B-V)_y', 'log10(rho[Msun]/pc^3)',\
 'rc', 'sigma0[km/s]', 'esigma0[km/s]', 'fb', 'efb', '[M/H]', 'Rgc[kpc]','Rsun[kpc]']

# names in OC datafile
names_oc = ['Cluster_name', 'RA', 'DEC', 'l', 'b', 'Dist Mod', 'EB-V', 'Age', 'ST' ,'Z', 'Diam', 'Fe/H', 'MRV',\
 'pm RA', 'pm Dec', 'logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]', 'rtK[pc]', 'Rhm[pc]',\
  '[Fe/H]K]', 'deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt' ,'Rgc' ,'z' ,'Diam[pc]', 'd[pc]']
path = '/Users/andrewbowen/ceb_project/data/data_files/'

# Globular cluster read in
GCs = pd.read_csv(path + 'gc_data.txt', sep = ' ', header = 0, names = names_gc)

# Open Cluster Read-in
OCs = pd.read_csv(path + 'oc_data.txt', sep = ' ', header = 0)

# Globular Cluster core radius
gc_rc = GCs['r_c']

gc_distance = GCs['Rsun[kpc]']
a = math.sqrt(2) * gc_rc

# Squared z velocity dispersion (from core radius) - will need cluster mass maybe 

# Velocity Dispersions for Globular Clusters
for index, cluster in GCs.iterrows():
	rc = cluster['rh[pc]'] * (3.086 * (10**13))#core radius (parsecs)

	a = math.sqrt(2) * rc
	mass = cluster['Mcl[Msun]']
	dist = cluster['Rsun[kpc]']

	if mass != -9.99:
		sigma_z = math.sqrt(((3/64)*(math.pi)*(6.6741 * (10**(-20) * (mass * 2 * (10 **30)))/rc)))
		# print(sigma_z, "km/s")

print(OCs.columns)
"""Make histograms of age, metallicity, size, mass comparing all clusters with ones that we 
can get a velocity dispersion for"""



# Diam column in arcmin (i think),