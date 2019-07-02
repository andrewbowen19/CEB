"""Hard-Soft Period cutoff script to be run on quest. Going to sample about 60k binaries for the clusters to find the 
	P_hs avg value for each cluster. Copied from hs_cutofff.ipynb in our cosmic repo"""

# Importing needed libraries -- make sure to run in cosmic environment
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import multidim
from cosmic.evolve import Evolve
import numpy as np
import matplotlib.pyplot as plt

# Dictionary neeeded for evolving the binaries
BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0,\
           'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, \
           'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, 'mxns': 3.0, \
           'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, \
           'merger': 0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, \
           'ppsn': 1, 'natal_kick_array' :[-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], \
           'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, \
           'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],\
           'cekickflag' : 0, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsnp' : 2.5, \
           'ecsn_mlow' : 1.6, 'aic' : 1, 'sigmadiv' :-20.0}

# Final kstar types to evolve up to (10- max evolution time, 11 - binary disruption, 12 - begin simbiotic phase): https://cosmic-popsynth.github.io/examples/index.html
final_kstar1 = [11,12]
final_kstar2 = [10]

# Sampling initial binaries - pulled from cosmic's documentation
InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('multidim', [0,14], [0,14], \
                            2,1, SFH_model='delta_burst', \
                            component_age=13000.0, met=0.02, size=60000)

# Initial Binary Period Histogram
f,ax = plt.subplots(figsize = (8,5))

ax.hist(np.log10(p_i), bins = 50, color = '#CF0A2C')
ax.set_xlabel('Input Periods (log-days)')
ax.grid(None)
f.savefig('/projects/p30137/abowen/CEB/cosmic/plots/InputBinaryHist-60k.png')

# Evolving the binaries
bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

# Defining initial and final periods
p_i = InitialBinaries['porb']
p_f = bcm['porb']


# Evolved Binary Histogram - Period

f,ax = plt.subplots(figsize = (8,5))
ax.hist(np.log10(p_f), bins = 50, color = 'black')
ax.set_xlabel('Evolved periods (log-days)')
ax.grid(None)
f.savefig('/projects/p30137/abowen/CEB/cosmic/plots/EvolvedBinaryHist-60k.png')

# Making copy of bcm for later
EvolvedBinaries = bcm

# Using hard-soft cutoff from Aaron's paper: arXiv:1506.08830
m1_f = bcm['mass_1']
m2_f = bcm['mass_2']


import pandas as pd
GC_sigma = pd.read_csv('/projects/p30137/abowen/CEB/cosmic/gc-sigma.txt', names = ['index','sigma_v'])
OC_sigma = pd.read_csv('/projects/p30137/abowen/CEB/cosmic/oc-sigma.txt', names = ['index','sigma_v'])

# Pulling only sigma values (no indices)
gc_sigma = GC_sigma['sigma_v']
oc_sigma = OC_sigma['sigma_v']

# Defining function that takes in mass of stars in binary (m1,m2), incoming/disrupting mass (m3) and velocity dispersion of a cluster
def get_Phs(m1, m2, m3, sigma):
    '''Function to calculate the hard-soft period cutoff given by Eq 1 in Geller, Leigh 2015
    Masses are in solar masses (converted later to kg), m1 and m2 are binary component masses,
    and m3 is the mass of the incoming (disrupting) object, 
    velocity dispersions are given in km/s'''
    G = 6.67408*10**-20#Gravitational constant in km^3 kg^-1 s^-2
    Msun = 1.989*10**30#Solar mass in kg
    const = (np.pi*G/np.sqrt(2)) * Msun#Last term used to convert from solar masses to kg
    Phs = const*(np.sqrt(m1*m2/m3)**3)*((np.sqrt(3)*sigma)**(-3))/(np.sqrt(m1+m2))
    return Phs

# Reading in OC and GC data files
names_gc = ['ID_x', 'Name', 'RA', 'DEC', 'L','B','R_Sun','R_gc','X','Y', 'Z', 'key_0','[Fe/H]_x', 'wt', 'E(B-V)_x',\
 'V_HB','(m-M)V_x', 'V_t', 'M_V,t', 'U-B', 'B-V', 'V-R', 'V-I', 'spt', 'ellip', 'ID_y', 'v_r', '+/-', 'v_LSR' ,'sig_v' ,'+/-.1', 'c', 'r_c', 'r_h', 'mu_V',\
  'rho_', 'lg(tc)', 'lg(th)', 'Mcl[Msun]', 'rh[pc]', '[Fe/H]_y', 'age[Gyr]', '(m-M)V_y', 'E(B-V)_y', 'log10(rho[Msun]/pc^3)',\
 'rc', 'sigma0[km/s]', 'esigma0[km/s]', 'fb', 'efb', '[M/H]', 'Rgc[kpc]','Rsun[kpc]']

# names in OC datafile
names_oc = ['Cluster_name', 'RA', 'DEC', 'l', 'b', 'Dist Mod', 'EB-V', 'Age', 'ST' ,'Z', 'Diam', 'Fe/H', 'MRV',\
 'pm RA', 'pm Dec', 'logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]', 'rtK[pc]', 'Rhm[pc]',\
  '[Fe/H]K]', 'deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt' ,'Rgc' ,'z' ,'Diam[pc]', 'd[pc]']
path = '/projects/p30137/abowen/CEB/cosmic/'

# Globular cluster read in
GCs = pd.read_csv(path + 'gc_data.txt', sep = ' ', header = 0, names = names_gc)
# Open Cluster file read in
OCs = pd.read_csv(path + 'oc_data.txt', sep = ' ', header = 0)

# Adding velocity dispersion columns to dataframe
GCs['sigma_v'] = gc_sigma
OCs['sigma_v'] = oc_sigma

""" Going to loop through clusters and binary pop (1k binaries for now - speed) 
    and get a period cutoff (get_phs) for each binary (using cluster velo dispersion)"""


# Globular Cluster for loop to find average period cutoff for each population of binaries in each cluster
gc_phs = []#going to be a list of the average hard-soft period cutoff for every GC
for index, cluster in GCs.iterrows():
    sig = cluster['sigma_v']
    bin_p = []#list of binary periods for each cluster
    for index, binary in EvolvedBinaries.iterrows():
        
        m1 = binary['mass_1']# * (1.989*10**30)#converting to kilograms from solar masses
        m2 = binary['mass_2']# * (1.989*10**30)
        m3 = 0.5#Average stellar mass given by Kroupa '93 IMF
#         Checking for calculated cluster sigma value, if not we'll feed it a given sigma value
        if (cluster['sigma_v'] != 0):
            period_cutoff = get_Phs(m1, m2, m3, sig)
            print(period_cutoff)
#         else:
#             period_cutoff = get_Phs(m1,m2,m3,8)
        if np.isfinite(period_cutoff):#only will calculate average cutoff if the value is finite
            bin_p.append(period_cutoff)#Appending cutoffs to list
    cluster_cutoff = np.mean(bin_p)#finding average of cutoffs (one average for each cluster)
    gc_phs.append(cluster_cutoff)

gc_p = np.asarray(gc_phs)#turning this list of avg period cutoffs into an array
gc_p = gc_p[np.where(gc_p != 0)] /(24 * 3600)#converting units into days (from seconds) 

# Histogram of average period cutoff for Globular Clusters
f,ax = plt.subplots(figsize = (8,5))
ax.hist(np.log10(gc_p), bins = 15, color = '#0E3386')
ax.set_xlabel('Average Period cutoff (log-days)')
ax.set_title('Globular Clusters')
ax.grid(None)

f.savefig('/projects/p30137/abowen/CEB/cosmic/plots/gc-avg-hsperiod-60k.png')

# OPEN CLUSTERS
oc_phs = []#going to be a list of the average period cutoff for every GC
for index, cluster in OCs.iterrows():
    sigma = cluster['sigma_v']
    bin_p = []#list of binary periods for each cluster
    for index, binary in EvolvedBinaries.iterrows():
        
        m1 = binary['mass_1']# * (1.989 * 10 ** 30)#converting to kilograms from solar masses
        m2 = binary['mass_2']# * (1.989 * 10 **30)
        m3 = 0.5#Average stellar mass given by Kroupa '93 IMF
#         Checking for calculated cluster sigma value, if not we'll feed it a given sigma value
        if (cluster['sigma_v'] != 0):
            period_cutoff = get_Phs(m1, m2, m3, sigma)
            print(period_cutoff)
        if np.isfinite(period_cutoff) and period_cutoff != 0:#only will calculate average cutoff if the value is finite
            bin_p.append(period_cutoff)#Appending cutoffs to list - each cluster will have a list of cutoffs for binaries 'in' that cluster
    cluster_cutoff = np.mean(bin_p)#finding average of cutoffs (one average for each cluster)
    oc_phs.append(cluster_cutoff)

oc_p = np.asarray(oc_phs)#turning this list of avg period cutoffs for OCs into an array
oc_p = oc_p[np.where(oc_p != 0) and np.isfinite(oc_p)] /(24 * 3600 * np.sqrt(3)**3)

# Histogram of average period cutoff for Open Clusters
f,ax = plt.subplots(figsize = (8,5))
ax.hist(np.log10(oc_p),bins = 35, color = '#0A162B')
ax.set_xlabel('Average Period cutoff (log-days)')
ax.set_title('Open Clusters')
ax.grid(None)
f.savefig('/projects/p30137/abowen/CEB/cosmic/plots/oc-avg-hsperiod-60k.png')




