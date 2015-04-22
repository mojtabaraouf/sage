# Make plots to see how the disc module has performed with SAGE

from pylab import *
import sys
import os
sys.path.insert(0, '../../../Swinburne Shared/6-MonthProject')
from galprops import galplot as gp
from galprops import galread as gr
from galprops import galcalc as gc

fsize = 28
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman'})

outdir = 'results/millennium/'

G = gr.sagesnap('model_z0.000', 0, 7, 'results/millennium/', disc=True)

DiscBinEdge = np.append(0, np.array([0.5*200*1.2**i for i in range(30)])) / 0.73

# Just get the galaxies of interest
filt = (G.Vvir>=200) * (G.Vvir<=235) * np.isfinite(G.Vvir) * (G.StellarMass/0.73 > 1.0) * (G.ColdGas/0.73 > 10**-0.5)
#randgals = np.array(np.random.random(10) * len(G.Vvir[filt]), dtype=int)
#randgals = np.unique(randgals)
#N = len(randgals)
N = len(G.Vvir[filt])
print N

Sigma_star_arr = np.zeros((N,30))
Sigma_gas_arr = np.zeros((N,30))

for i in xrange(N):
	radius_bins = DiscBinEdge / G.Vvir[filt][i]
	Sigma_star = (G.DiscStars[filt][i]*1e10/0.73) / (np.pi*(radius_bins[1:]**2 - radius_bins[:-1]**2)*1e6) # Solar masses per pc^2
	Sigma_star_arr[i,:] = Sigma_star
	Sigma_gas_arr[i,:] = (G.DiscGas[filt][i]*1e10/0.73) / (np.pi*(radius_bins[1:]**2 - radius_bins[:-1]**2)*1e6)

#Sigma_star_arr[np.where(Sigma_star_arr==0)] = np.min(Sigma_star_arr[Sigma_star_arr>0])
#Sigma_star_arr = np.sort(np.log10(Sigma_star_arr), axis=0)
Sigma_star_arr = np.sort(Sigma_star_arr, axis=0)
Sigma_star_low = Sigma_star_arr[int(N*0.16), :]
Sigma_star_av = np.mean(Sigma_star_arr, axis=0)
Sigma_star_med = np.median(Sigma_star_arr, axis=0)
Sigma_star_high = Sigma_star_arr[int(N*0.84), :]
Sigma_star_std = np.std(Sigma_star_arr, axis=0)

#Sigma_gas_arr[np.where(Sigma_gas_arr==0)] = np.min(Sigma_gas_arr[Sigma_gas_arr>0])
#Sigma_gas_arr = np.sort(np.log10(Sigma_gas_arr), axis=0)
Sigma_gas_arr = np.sort(Sigma_gas_arr, axis=0)
Sigma_gas_low = Sigma_gas_arr[int(N*0.16), :]
Sigma_gas_av = np.mean(Sigma_gas_arr, axis=0)
Sigma_gas_med = np.median(Sigma_gas_arr, axis=0)
Sigma_gas_high = Sigma_gas_arr[int(N*0.84), :]
Sigma_gas_std = np.std(Sigma_gas_arr, axis=0)

rad_plot = (DiscBinEdge[1:]+DiscBinEdge[:-1])/(200+235)

plt.figure()
plt.plot(rad_plot, Sigma_star_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_star_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_star_high, Sigma_star_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
#plt.fill_between(rad_plot, Sigma_star_av+Sigma_star_std, Sigma_star_av-Sigma_star_std, color='b', alpha=0.3)
#plt.plot((radius_bins[1:]+radius_bins[:-1])/2, Sigma_star, lw=1)
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e4])
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm star}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
gp.savepng(outdir+'StarSurface')

plt.figure()
plt.plot(rad_plot, Sigma_gas_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_gas_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_gas_high, Sigma_gas_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
#plt.fill_between(rad_plot, Sigma_gas_av+Sigma_gas_std, Sigma_gas_av-Sigma_gas_std, color='b', alpha=0.3)
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e3])
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm gas}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
gp.savepng(outdir+'GasSurface')