# Andrew's plotting script for period-eccentricity tidal circularization plots
# Want scatter plot of eccentricity vs period colored by Nrec
# 2D histogram (CIERA REU Python sequence part 5) - Nrec and %rec in each ecc-p bin

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.cm as cm

datAll = pd.read_csv('/Users/andrewbowen/ceb_project/testing/GlobularClusters/pecc/all-ecc-p.csv', header = 0, names = ['e', 'p'])
datObs = pd.read_csv('/Users/andrewbowen/ceb_project/testing/GlobularClusters/pecc/obs-ecc-p.csv', header = 0, names = ['e', 'p'])
datRec = pd.read_csv('/Users/andrewbowen/ceb_project/testing/GlobularClusters/pecc/rec-ecc-p.csv', header = 0, names = ['e', 'p'])

# Test 2D hist - 'all' dataframe
# plt.hist2d(datAll['p'], datAll['e'], bins = [10,10])
# plt.xlabel('period (days)')
# plt.ylabel('eccentricity')
# plt.show()


# Function copied from Gaia_M67_SOLUTIONS (/CIERA_REU/PythonTutorialSequence/Part5/)
# Want to convert from plotting proper motion to ecc and p
# Can call this for each scenario (GC/OC, baseline/colossus, crowding/noCrowding)
def plotPecc(p, ecc, xmin, xmax, Nx, ymin, ymax, Ny, plot_title, norm = None):
	f = plt.figure(figsize=(8, 8)) 
	gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[2])
	ax3 = plt.subplot(gs[3])

	#histograms
	pBins = np.linspace(xmin,xmax,Nx)
	eccBins = np.linspace(ymin,ymax,Ny)
	hx1D, x1D, im = ax1.hist(p, bins=pBins, histtype='step', fill=False)
	hy1D, y1D, im = ax3.hist(ecc, bins=eccBins, histtype='step', fill=False, orientation="horizontal")

	#heatmap
	#ax2.hexbin(r['p']*np.cos(r['dec']*np.pi/180.), r['ecc'], gridsize=30, cmap=cm.Blues, bins='log', extent=(xmin, xmax, ymin, ymax))
	h2D, x2D, y2D, im = ax2.hist2d(p, ecc, bins=[Nx, Ny], \
									range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)
	ax1.set_xlim(xmin, xmax)
	ax2.set_xlim(xmin, xmax)
	ax2.set_ylim(ymin, ymax)
	ax3.set_ylim(ymin, ymax)
	ax2.set_xlabel(r'period (days)', fontsize=16)
	ax2.set_ylabel(r'binary eccentricity', fontsize=16)
	plt.setp(ax1.get_yticklabels()[0], visible=False)
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.setp(ax3.get_yticklabels(), visible=False)
	plt.setp(ax3.get_xticklabels()[0], visible=False)
	plt.title(plot_title, fontsize = 22) # Can add titles to distinguish between plots
	f.subplots_adjust(hspace=0., wspace=0.)
	# plt.show()
	# f.savefig('./plots/contour_plots/local/'+ plot_title + '-p-ecc-hist2D.pdf')

	return x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, (ax1, ax2, ax3)

# #################################################################################################################


# New percent rec function - takes in rec and obs 2D arrays (period and ecc)
def pRecHist(recHist, obsHist):

	# Dividing 
	pRecHist2D = np.divide(obsHist, recHist)
	pRec2D = np.nan_to_num(pRecHist2D)
	print('Rec',recHist)
	print('Obs')
	# f,ax 
	plt.hexbin(np.log10(pRec2D[0][0:]), pRec2D[1][0:])

	plt.title(r'% recovered')
	plt.xlabel('period (log-days)')
	plt.ylabel('Eccentricity')
	plt.show()



# Function to make p-ecc percent recovered histogram (Nrec/Nobs in each bin)
# hist needs to be 2D histogram already - basically just plotting this
# def percentRecHist(recHist, obsHist, xmin, xmax, Nx, ymin, ymax, Ny, plot_title, norm = None):
	# f = plt.figure(figsize=(8, 8)) 
	# gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
	# ax1 = plt.subplot(gs[0])
	# ax2 = plt.subplot(gs[2])
	# ax3 = plt.subplot(gs[3])

	# #histograms
	# pBins = np.linspace(xmin,xmax,Nx)
	# eccBins = np.linspace(ymin,ymax,Ny)
	# hx1D, x1D, im = ax1.hist(recHist[0], bins=pBins, histtype='step', fill=False)
	# hy1D, y1D, im = ax3.hist(recHist[1], bins=eccBins, histtype='step', fill=False, orientation="horizontal")

	# # dividing to get percent recovered
	# pRecHist2D = np.divide(obsHist, recHist)
	# pRec2D = np.nan_to_num(pRecHist2D)

	# plt.hexbin(pRec2D[0][0:], pRec2D[1][0:])
	# for a in pRec2D:
	# 	print('Percent recovered array (10000 bins): ', pRec2D)

	# # Test plotting - can just index 2D percent recovered hist from above
	# f = plt.figure(figsize=(8, 8)) 
	# gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
	# ax1 = plt.subplot(gs[0])
	# ax2 = plt.subplot(gs[2])
	# ax3 = plt.subplot(gs[3])

	# # Making heatmap
	# h2D, x2D, y2D, im = ax2.hist2d(pRec2D[0][0:], pRec2D[1][0:], bins=[Nx, Ny], \
	# 						range=[[xmin, xmax], [ymin, ymax]], norm = None, cmap = cm.Blues)

	# ax1.set_xlim(xmin, xmax)
	# ax2.set_xlim(xmin, xmax)
	# ax2.set_ylim(ymin, ymax)
	# ax3.set_ylim(ymin, ymax)
	# ax2.set_xlabel(r'period (days) - foo', fontsize=16)
	# ax2.set_ylabel(r'binary eccentricity', fontsize=16)
	# plt.setp(ax1.get_yticklabels()[0], visible=False)
	# plt.setp(ax1.get_xticklabels(), visible=False)
	# plt.setp(ax3.get_yticklabels(), visible=False)
	# plt.setp(ax3.get_xticklabels()[0], visible=False)
	# plt.title(plot_title, fontsize = 22) # Can add titles to distinguish between plots
	# f.subplots_adjust(hspace=0., wspace=0.)
	# plt.title('Foo foo I like poo')
	# plt.show()



	# return h2D

# #################################################################################################################

# Setting up bins: will likely change
xmin, xmax, Nx = 0, 50, 100
ymin, ymax, Ny = 0, 1, 100

# Testing with local files - baseline GCs no crowding
# All
x2DAll, y2DAll, h2DAll, x1DAll, hx1DAll, y1DAll, hy1DAll, axAll = plotPecc(datAll['p'], datAll['e'], xmin, xmax, Nx, ymin, ymax, Ny, 'All',norm=mpl.colors.LogNorm())
# Obs
x2DObs, y2DObs, h2DObs, x1DObs, hx1DObs, y1DObs, hy1DObs, axObs = plotPecc(datObs['p'], datObs['e'], xmin, xmax, Nx, ymin, ymax, Ny, 'Obs',norm=mpl.colors.LogNorm())
# Rec
x2DRec, y2DRec, h2DRec, x1DRec, hx1DRec, y1DRec, hy1DRec, axRec = plotPecc(datRec['p'], datRec['e'], xmin, xmax, Nx, ymin, ymax, Ny, 'Rec',norm=mpl.colors.LogNorm())
print(np.max(h2DRec.flatten()))

# Going through every directory to find ecc-p files
import os
n = 0
for root, dirs, files in os.walk(".", topdown=True):
	for name in files:
		# Only reading in p-ecc output files (distinguishes from other csv files in directory)
		if '-ecc-p.csv' in name:
			df = pd.read_csv(os.path.join(root,name), header = 0, names = ['e', 'p'])

			n += 1
			# Making the histograms and saving them in 'plots/contour_plots' directory
			x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, ax = plotPecc(df['p'],df['e'], 
				xmin, xmax, Nx, ymin, ymax, Ny, name[0:3], norm=mpl.colors.LogNorm())

			# Only want matching scenario obs-rec dataframes

			if ('obs' in name) and ('.csv' in name):
				obsHist2D = h2D
			elif 'rec' in name and '.csv' in name:
				recHist2D = h2D
			# Making 2D histgram of percent recovered
				# pRecHist = percentRecHist(recHist2D, obsHist2D,xmin, xmax, Nx, ymin, ymax, Ny, r'% rec', norm=mpl.colors.LogNorm())
				pRecHist = pRecHist(recHist2D, obsHist2D)
		else:
			pass


# Creating percent recovered (rec/obs) maps
# Want arrays of percent recovered in each bin: for each bin we should divide Nrec by Nobs (if it's zero it will be zero)
pRecHist2D = np.divide(h2DRec,h2DObs)
print('percent recovered hist:' , pRecHist2D)

# Period hists (just p)
print('2D Rec: ', h2DRec)
print('2D Obs', h2DObs)


pRec2D = np.divide(h2DObs, h2DRec)
pRec2D = np.nan_to_num(pRec2D)

print(pRec2D.shape)
print(pRec2D[0][0:].shape)




