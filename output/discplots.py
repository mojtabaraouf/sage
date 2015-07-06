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

indir = 'results/millennium/'
outdir = indir+'plots/'
if not os.path.exists(outdir): os.makedirs(outdir)

G = gr.sagesnap('model_z0.000', 0, 7, indir, disc=True) #Commit 4 - good results/

h = 0.678;

DiscBinEdge = np.append(0, np.array([0.5*200*1.2**i for i in range(30)])) / h
#DiscBinEdge = np.append(0, np.array([0.2*200*1.2**i for i in range(30)])) / h
#DiscBinEdge = np.append(0, np.array([0.5*200*1.15**i for i in range(30)])) / h

# Just get the galaxies of interest
filt = (G.Vvir>=200) * (G.Vvir<=235) * np.isfinite(G.Vvir) * (G.StellarMass/h > 1.0) * (G.ColdGas/h > 10**-0.5)# * (G.StellarMass/h < 10**0.8) * (G.ColdGas/h < 10**0.2)
filt2 = (G.Vvir>=175) * (G.Vvir<200) * np.isfinite(G.Vvir) * (G.StellarMass/h > 1.0) * (G.ColdGas/h > 10**-0.8) # Maybe 1.4 for stellar mass
#filt = (G.Vvir>=200) * (G.Vvir<=235) * np.isfinite(G.Vvir) * ((G.ClassicalBulgeMass+G.SecularBulgeMass) <= 0.15*G.StellarMass)
#filt2 = (G.Vvir>=175) * (G.Vvir<200) * np.isfinite(G.Vvir) * ((G.ClassicalBulgeMass+G.SecularBulgeMass) <= 0.15*G.StellarMass)


#randgals = np.array(np.random.random(10) * len(G.Vvir[filt]), dtype=int)
#randgals = np.unique(randgals)
#N = len(randgals)

N = len(G.Vvir[filt])
N2 = len(G.Vvir[filt2])

print N, N2

Sigma_star_arr = np.zeros((N,30))
Sigma_gas_arr = np.zeros((N,30))
Sigma_HI_arr = np.zeros((N,30))
Sigma_H2_arr = np.zeros((N,30))

Sigma_star_arr2 = np.zeros((N2,30))
Sigma_gas_arr2 = np.zeros((N2,30))

for i in xrange(N):
	radius_bins = DiscBinEdge / G.Vvir[filt][i]
	Sigma_star = (G.DiscStars[filt][i]*1e10/h) / (np.pi*(radius_bins[1:]**2 - radius_bins[:-1]**2)*1e6) # Solar masses per pc^2
	Sigma_star_arr[i,:] = Sigma_star
	Sigma_gas = (G.DiscGas[filt][i]*1e10/h) / (np.pi*(radius_bins[1:]**2 - radius_bins[:-1]**2)*1e6)
	Sigma_gas_arr[i,:] = Sigma_gas
	#
	H2_HI = 1.306e-3 * (Sigma_gas**2 + 0.1*Sigma_gas * np.sqrt(Sigma_star*Sigma_star[0]))**0.92
	H2_gas = 0.75 * 1.0/(gc.divide(np.ones(len(H2_HI)),H2_HI) + 1) * (1 - gc.divide(G.DiscGasMetals[filt][i], G.DiscGas[filt][i]))
	#H2_gas = 1.0/(gc.divide(np.ones(len(H2_HI)),H2_HI) + 1)
	Sigma_H2_arr[i,:] = Sigma_gas_arr[i,:] * H2_gas
	Sigma_HI_arr[i,:] = gc.divide(Sigma_H2_arr[i,:], H2_HI)

for i in xrange(N2):
	radius_bins = DiscBinEdge / G.Vvir[filt2][i]
	Sigma_star = (G.DiscStars[filt2][i]*1e10/h) / (np.pi*(radius_bins[1:]**2 - radius_bins[:-1]**2)*1e6) # Solar masses per pc^2
	Sigma_star_arr2[i,:] = Sigma_star
	Sigma_gas_arr2[i,:] = (G.DiscGas[filt2][i]*1e10/h) / (np.pi*(radius_bins[1:]**2 - radius_bins[:-1]**2)*1e6)

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

Sigma_HI_arr = np.sort(Sigma_HI_arr, axis=0)
Sigma_HI_low = Sigma_HI_arr[int(N*0.16), :]
Sigma_HI_av = np.mean(Sigma_HI_arr, axis=0)
Sigma_HI_med = np.median(Sigma_HI_arr, axis=0)
Sigma_HI_high = Sigma_HI_arr[int(N*0.84), :]

Sigma_H2_arr = np.sort(Sigma_H2_arr, axis=0)
Sigma_H2_low = Sigma_H2_arr[int(N*0.16), :]
Sigma_H2_av = np.mean(Sigma_H2_arr, axis=0)
Sigma_H2_med = np.median(Sigma_H2_arr, axis=0)
Sigma_H2_high = Sigma_H2_arr[int(N*0.84), :]


rad_plot = (DiscBinEdge[1:]+DiscBinEdge[:-1])/(200+235)

gp.figure()
plt.plot(rad_plot, Sigma_star_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_star_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_star_high, Sigma_star_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
#plt.fill_between(rad_plot, Sigma_star_av+Sigma_star_std, Sigma_star_av-Sigma_star_std, color='b', alpha=0.3)
#plt.plot((radius_bins[1:]+radius_bins[:-1])/2, Sigma_star, lw=1)
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e4])
plt.errorbar([0,0], [0,1], [1,1], ecolor='k', color='k', label=r'Leroy et al. (2008)')
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm star}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$V_{\rm vir} \in$ [200,235] km s$^{-1}$')
gp.Leroygals()
gp.savepng(outdir+'StarSurface')

gp.figure()
plt.plot(rad_plot, Sigma_gas_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_gas_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_gas_high, Sigma_gas_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
#plt.fill_between(rad_plot, Sigma_gas_av+Sigma_gas_std, Sigma_gas_av-Sigma_gas_std, color='b', alpha=0.3)
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e3])
plt.errorbar([0,0], [0,1], [1,1], ecolor='k', color='k', label=r'Leroy et al. (2008)')
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm gas}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$V_{\rm vir} \in$ [200,235] km s$^{-1}$')
gp.Leroygals(HI=True, H2=True)
gp.savepng(outdir+'GasSurface')

gp.figure()
plt.plot(rad_plot, Sigma_HI_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_HI_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_HI_high, Sigma_HI_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
#plt.fill_between(rad_plot, Sigma_gas_av+Sigma_gas_std, Sigma_gas_av-Sigma_gas_std, color='b', alpha=0.3)
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e3])
plt.errorbar([0,0], [0,1], [1,1], ecolor='k', color='k', label=r'Leroy et al. (2008)')
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm HI}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$V_{\rm vir} \in$ [200,235] km s$^{-1}$')
gp.Leroygals(HI=True)
gp.savepng(outdir+'HISurface')

gp.figure()
plt.plot(rad_plot, Sigma_H2_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_H2_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_H2_high, Sigma_H2_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
#plt.fill_between(rad_plot, Sigma_gas_av+Sigma_gas_std, Sigma_gas_av-Sigma_gas_std, color='b', alpha=0.3)
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e3])
plt.errorbar([0,0], [0,1], [1,1], ecolor='k', color='k', label=r'Leroy et al. (2008)')
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm H_2}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$V_{\rm vir} \in$ [200,235] km s$^{-1}$')
gp.Leroygals(H2=True)
gp.savepng(outdir+'H2Surface')






rad_plot = (DiscBinEdge[1:]+DiscBinEdge[:-1])/(200+175)


# Lower vel range
Sigma_star_arr = np.sort(Sigma_star_arr2, axis=0)
Sigma_star_low = Sigma_star_arr[int(N2*0.16), :]
Sigma_star_av = np.mean(Sigma_star_arr, axis=0)
Sigma_star_med = np.median(Sigma_star_arr, axis=0)
Sigma_star_high = Sigma_star_arr[int(N2*0.84), :]

Sigma_gas_arr = np.sort(Sigma_gas_arr2, axis=0)
Sigma_gas_low = Sigma_gas_arr[int(N2*0.16), :]
Sigma_gas_av = np.mean(Sigma_gas_arr, axis=0)
Sigma_gas_med = np.median(Sigma_gas_arr, axis=0)
Sigma_gas_high = Sigma_gas_arr[int(N2*0.84), :]

gp.figure()
plt.plot(rad_plot, Sigma_star_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_star_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_star_high, Sigma_star_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e4])
plt.errorbar([0,0], [0,1], [1,1], ecolor='k', color='k', label=r'Leroy et al. (2008)')
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm star}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$V_{\rm vir} \in$ [175,200] km s$^{-1}$')
gp.Leroygals(HighVvir=False, LowVvir=True)
gp.savepng(outdir+'StarSurface2')

gp.figure()
plt.plot(rad_plot, Sigma_gas_av, 'r-', lw=2, label=r'Mean')
plt.plot(rad_plot, Sigma_gas_med, 'r--', lw=2, label=r'Median')
plt.fill_between(rad_plot, Sigma_gas_high, Sigma_gas_low, color='b', alpha=0.3)
plt.plot([-1,-1],[0,1], 'b-', lw=5, alpha=0.5, label=r'$16^{\rm th}$--$84^{\rm th}$ percentile')
plt.yscale('log', nonposy='clip')
plt.axis([0,25,1e0,1e3])
plt.errorbar([0,0], [0,1], [1,1], ecolor='k', color='k', label=r'Leroy et al. (2008)')
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\Sigma_{\rm gas}$ [M$_{\bigodot}$ pc$^{-2}$]')
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$V_{\rm vir} \in$ [175,200] km s$^{-1}$')
gp.Leroygals(HI=True, H2=True, HighVvir=False, LowVvir=True)
gp.savepng(outdir+'GasSurface2')



### OTHER PLOTS FOR INTEGRATED PROPERTIES

# Stellar mass function
gp.figure()
gp.massfunction(G.StellarMass*1e10/h, 62.5/h, extra=2, h=h, step=False, fsize=fsize, binwidth=0.1, label=r'SAGE Disc')
#gp.massfunction(G.StellarMass*1e10*h, 500*(8./512)**(1./3), extra=2, h=1, step=False, fsize=fsize, binwidth=0.1, label=r'SAGE Disc')
plt.axis([8,11.8,5e-6,0.2])
gp.savepng(outdir+'SMF', xpixplot=768, ypixplot=512)

# Black hole -- bulge relation
filt = (G.BlackHoleMass > 1e-5) * ((G.ClassicalBulgeMass+G.SecularBulgeMass) > 0.01)
gp.figure()
gp.bhbulge(G.BlackHoleMass[filt]*1e10/h, (G.ClassicalBulgeMass[filt]+G.SecularBulgeMass[filt])*1e10/h, h=h, fsize=fsize, Nbins=100, extra=2, label=r'SAGE Disc')
gp.savepng(outdir+'BHbulge', xpixplot=768, ypixplot=512)

# Stellar mass -- gas metallicity relation
filt = (G.Type==0) * (G.ColdGas > 0) * (G.MetalsColdGas > 0) * (G.StellarMass > 0.01)
gp.figure()
gp.massmet(G.StellarMass[filt]*1e10/h, np.log10(G.MetalsColdGas[filt]/G.ColdGas[filt]/0.02)+9, h=h, fsize=fsize, Nbins=100, extra=True, label=r'SAGE Disc')
gp.savepng(outdir+'MassMet', xpixplot=768, ypixplot=512)

# Baryonic Tully-Fisher
filt = (G.Type==0) * (G.StellarMass+G.ColdGas > 0) * (G.ClassicalBulgeMass+G.SecularBulgeMass > 0.1*G.StellarMass) * (G.ClassicalBulgeMass+G.SecularBulgeMass < 0.5*G.StellarMass)
gp.figure()
gp.btf((G.StellarMass[filt]+G.ColdGas[filt])*1e10/h, G.Vmax[filt], h=h, fsize=fsize, extra=True, label=r'SAGE Disc')
gp.savepng(outdir+'BTF', xpixplot=768, ypixplot=512)


### BH mass history plot
zstr = ['0.000', '0.089', '0.208', '0.362', '0.564', '0.828', '1.173', '1.630', '2.239', '3.060', '4.179', '5.724', '7.883', '10.944', '15.343']
if os.path.isfile(indir+'model_z'+zstr[1]+'_0'):
    z = np.array(zstr, dtype='f4')
    SFRD = np.zeros(len(z))
    gp.figure()
    for s in xrange(len(z)):
        M = gr.sagesnap('model_z'+zstr[s], 0, 7, indir, disc=True)
        SFRD[s] = np.log10(np.sum(M.SfrDisk + M.SfrBulge)*h**3 / (62.5**3))
        gp.point_with_spread(z[s], M.BlackHoleMass[M.BlackHoleMass>0]*1e10/h, mean=False, alpha=0.8)
    plt.xlabel(r'Redshift')
    plt.yscale('log')
    plt.ylabel(r'$\log_{10} (M_{\rm BH}\ [\mathrm{M}_{\bigodot}])$')
    plt.axis([0, 12, 1e3, 1e8])
    gp.savepng(outdir+'BHhistory', xpixplot=768, ypixplot=512)



### Madau plot
if os.path.isfile(indir+'model_z'+zstr[1]+'_0'):
	gp.figure()
	gp.SFRD_obs(h)
	plt.plot(z, SFRD, 'k-', lw=2, label=r'\textsc{sage} Disc')
	plt.axis([0,7.5,-2.8,-0.5])
	plt.xlabel(r'Redshift')
	plt.ylabel(r'$\log_{10} \left(\bar{\rho}_{\rm SFR}\ [\mathrm{M}_{\bigodot}\ \mathrm{yr}^{-1}\ \mathrm{cMpc}^{-3}] \right)$')
	plt.legend(fontsize=fsize-4, loc='best', frameon=False, numpoints=1)
	gp.savepng(outdir+'SFRD', xpixplot=768, ypixplot=512)



	