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

h = 0.73;
boxlen = 62.5/h

#angle_limit = 180 # Number of degrees between stellar and gas discs to use stellar info on HI/H2 ratio
#cos_angle_limit = np.cos(angle_limit*np.pi/180)

DiscBinEdge = np.append(0, np.array([0.5*200*1.2**i for i in range(30)])) / h
#DiscBinEdge = np.append(0, np.array([0.2*200*1.2**i for i in range(30)])) / h
#DiscBinEdge = np.append(0, np.array([0.5*200*1.15**i for i in range(30)])) / h

# Just get the galaxies of interest
filt = (G.Vvir>=200) * (G.Vvir<=235) * np.isfinite(G.Vvir) * (G.StellarMass/h > 1.0) * (G.ColdGas/h > 10**-0.5) * (G.Type==0)# * (G.StellarMass/h < 10**0.8) * (G.ColdGas/h < 10**0.2)
filt2 = (G.Vvir>=175) * (G.Vvir<200) * np.isfinite(G.Vvir) * (G.StellarMass/h > 1.0) * (G.ColdGas/h > 10**-0.8) * (G.Type==0)# Maybe 1.4 for stellar mass
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
    cos_angle = G.SpinGas[filt][i,0]*G.SpinStars[filt][i,0] + G.SpinGas[filt][i,1]*G.SpinStars[filt][i,1] + G.SpinGas[filt][i,2]*G.SpinStars[filt][i,2]
    #if cos_angle >= cos_angle_limit:
    H2_HI = 1.306e-3 * (Sigma_gas**2 + 0.1*Sigma_gas * np.sqrt(Sigma_star*Sigma_star[0]))**0.92
    #else:
    #H2_HI = 1.306e-3 * Sigma_gas**1.84
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
"""
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
    """
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
"""
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
    """


### OTHER PLOTS FOR INTEGRATED PROPERTIES

# Stellar mass function
gp.figure()
gp.massfunction(G.StellarMass*1e10/h, boxlen, extra=2, h=h, step=False, fsize=fsize, binwidth=0.1, label=r'\textsc{dark sage}')
red = (G.SfrDisk+G.SfrBulge)/(G.StellarMass*1e10/h) < 1e-11
gp.massfunction(G.StellarMass[red]*1e10/h, boxlen, extra=0, h=h, step=False, fsize=fsize, binwidth=0.1, colour='r', lw=1)
gp.massfunction(G.StellarMass[True-red]*1e10/h, boxlen, extra=0, h=h, step=False, fsize=fsize, binwidth=0.1, colour='b', lw=1)
#gp.massfunction(G.StellarMass*1e10*h, 500*(8./512)**(1./3), extra=2, h=1, step=False, fsize=fsize, binwidth=0.1, label=r'\textsc{dark sage}')
plt.axis([8,11.8,5e-6,0.2])
gp.savepng(outdir+'SMF', xpixplot=768, ypixplot=512)

pcs = [0.68, 0.95]

# Black hole -- bulge relation
filt = (G.BlackHoleMass > 1e-5) * ((G.ClassicalBulgeMass+G.SecularBulgeMass) > 0.01)
gp.figure()
gp.bhbulge(G.BlackHoleMass[filt]*1e10/h, (G.ClassicalBulgeMass[filt]+G.SecularBulgeMass[filt])*1e10/h, h=h, fsize=fsize, Nbins=100, extra=2, label=r'\textsc{dark sage}', pcs=pcs)
gp.savepng(outdir+'BHbulge', xpixplot=768, ypixplot=512)

# Stellar mass -- gas metallicity relation
filt = (G.Type==0) * (G.ColdGas > 0) * (G.MetalsColdGas > 0) * (G.StellarMass > 0.01)
gp.figure()
gp.massmet(G.StellarMass[filt]*1e10/h, np.log10(G.MetalsColdGas[filt]/G.ColdGas[filt]/0.02)+9, h=h, fsize=fsize, Nbins=100, extra=True, label=r'\textsc{dark sage}', pcs=pcs)
gp.savepng(outdir+'MassMet', xpixplot=768, ypixplot=512)

# Baryonic Tully-Fisher
filt = (G.Type==0) * (G.StellarMass+G.ColdGas > 0) * (G.ClassicalBulgeMass+G.SecularBulgeMass > 0.1*G.StellarMass) * (G.ClassicalBulgeMass+G.SecularBulgeMass < 0.5*G.StellarMass)
gp.figure()
gp.btf((G.StellarMass[filt]+G.ColdGas[filt])*1e10/h, G.Vmax[filt], h=h, fsize=fsize, extra=True, label=r'\textsc{dark sage}', pcs=pcs)
gp.savepng(outdir+'BTF', xpixplot=768, ypixplot=512)

# Quiescent satellites
filt = (G.Type!=0) * (G.StellarMass > 0.0)
gp.figure()
gp.quiescenthalo(G.StellarMass[filt]*1e10/h, (G.SfrBulge[filt]+G.SfrDisk[filt]), G.CentralMvir[filt]*1e10/h, h=h, Nbins=15, label=r'\textsc{dark sage}', extra=True)
gp.savepng(outdir+'QuiescentSat', xpixplot=768, ypixplot=512)

# HI and H2 mass functions
filt = (G.ColdGas > 0.0) * (G.Vvir > 0.0)
Sigma_gas, Sigma_star = np.zeros((len(G[filt]), 30)), np.zeros((len(G[filt]), 30))
DiscHI, DiscH2 = np.zeros((len(G[filt]), 30)), np.zeros((len(G[filt]), 30))
cos_angle = G.SpinGas[filt][:,0]*G.SpinStars[filt][:,0] + G.SpinGas[filt][:,1]*G.SpinStars[filt][:,1] + G.SpinGas[filt][:,2]*G.SpinStars[filt][:,2]
#angle_filt = (cos_angle < cos_angle_limit)
for i in xrange(30):
    area = np.pi * (DiscBinEdge[i+1]**2 - DiscBinEdge[i]**2) * 1e6 / G.Vvir[filt]**2
    Sigma_gas[:,i] = (G.DiscGas[filt][:,i]*1e10/h) / area
    Sigma_star[:,i] = (G.DiscStars[filt][:,i]*1e10/h) / area
    H2_HI = 1.306e-3 * (Sigma_gas[:,i]**2 + 0.1*Sigma_gas[:,i] * np.sqrt(Sigma_star[:,i]*Sigma_star[:,0]))**0.92
    #H2_HI[angle_filt] = 1.306e-3 * Sigma_gas[angle_filt][:,i]**1.84
    DiscH2[:,i] = 0.75 * 1.0/(gc.divide(np.ones(len(H2_HI)),H2_HI) + 1) * (G.DiscGas[filt][:,i] - G.DiscGasMetals[filt][:,i]) * 1e10/h
DiscHI = 0.75*G.DiscGas[filt]*1e10/h - G.DiscGasMetals[filt]*1e10/h - DiscH2

#HImass = np.sum(DiscHI,axis=1)
H2mass = np.sum(DiscH2,axis=1)
HImass = 0.75*G.ColdGas[filt]*1e10/h - G.MetalsColdGas[filt]*1e10/h - H2mass

gp.figure()
gp.massfunction(HImass, boxlen, extra=0, h=h, step=False, fsize=fsize, binwidth=0.1, label=r'H\,\textsc{i} \textsc{dark sage}')
gp.massfunction(H2mass, boxlen, extra=0, h=h, step=False, fsize=fsize, binwidth=0.1, label=r'H$_2$ \textsc{dark sage}', colour='c', ls='--')
gp.massfunction_HI_H2_obs()
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
gp.savepng(outdir+'GMFs', xpixplot=768, ypixplot=512)

# Quiescent fraction of galaxies
gp.figure()
gp.quiescent(G.StellarMass*1e10/h, G.SfrDisk+G.SfrBulge, label=r'\textsc{dark sage}', extra=True, h=h)
gp.savepng(outdir+'QuiescentGals', xpixplot=768, ypixplot=512)

# Gas fraction
gp.figure()
SFR = G.SfrDisk + G.SfrBulge
filt = (G.StellarMass > 0.0) * (G.ColdGas > 0.0) * (SFR/(G.StellarMass*1e10/h) > 1e-11)
gp.contour(np.log10(G.StellarMass[filt]*1e10/h), np.log10(G.ColdGas[filt]/G.StellarMass[filt]), range=[[8.5,12],[-1.5,1.5]], Nlevels=5, Nbins=50, pcs=pcs)
gp.gas_fraction_obs(h)
plt.axis([8.5, 12, -1.45, 1.5])
gp.savepng(outdir+'Gasfraction_SFgals', xpixplot=768, ypixplot=512)

plt.close("all")

# Spread of angles between stellar and gas discs
filt = (G.ColdGas>0.0) * (G.StellarMass - G.SecularBulgeMass - G.ClassicalBulgeMass > 0.0)
cos_theta = G.SpinGas[filt][:,0]*G.SpinStars[filt][:,0] + G.SpinGas[filt][:,1]*G.SpinStars[filt][:,1] + G.SpinGas[filt][:,2]*G.SpinStars[filt][:,2]
cos_theta[cos_theta>1.0] = 1.0
cos_theta[cos_theta<-1.0] = -1.0
theta = np.arccos(cos_theta)*180/np.pi
gp.figure()
plt.hist(theta, bins=30)
plt.xlabel(r'Angle between gas and stellar discs [degrees]')
plt.ylabel(r'Number of galaxies')
ymin, ymax = plt.gca().get_ylim()
plt.text(100, 0.2*ymax, str(100.0*len(theta[theta<=10])/len(theta))[:4]+r' $\% \leq 10^{\circ}$' )
plt.axis([0,180,9e-1,ymax])
plt.yscale('log', nonposy='clip')
gp.savepng(outdir+'SpinOffset', xpixplot=768, ypixplot=512)

# Mass distributions for offset galaxies
gp.figure()
N_gas, edges = np.histogram(np.log10(G.ColdGas[filt][theta>20]*1e10/h), range=[8,12], bins=30)
N_gas_tot, edges = np.histogram(np.log10(G.ColdGas[filt]*1e10/h), range=[8,12], bins=30)
N_stars, edges = np.histogram(np.log10(G.StellarMass[filt][theta>20]*1e10/h), range=[8,12], bins=30)
N_stars_tot, edges = np.histogram(np.log10(G.StellarMass[filt]*1e10/h), range=[8,12], bins=30)
N_stars_disc, edges = np.histogram(np.log10((G.StellarMass[filt][theta>20]-G.SecularBulgeMass[filt][theta>20]-G.ClassicalBulgeMass[filt][theta>20])*1e10/h), range=[8,12], bins=30)
N_stars_disc_tot, edges = np.histogram(np.log10((G.StellarMass[filt]-G.SecularBulgeMass[filt]-G.ClassicalBulgeMass[filt])*1e10/h), range=[8,12], bins=30)
print 'N_stars_disc', N_stars_disc
print 'N_stars_disc_tot', N_stars_disc_tot
y_gas = (1.0*N_gas)/N_gas_tot
y_stars = (1.0*N_stars)/N_stars_tot
y_stars_disc = (1.0*N_stars_disc)/N_stars_tot
y_stars_disc2 = (1.0*N_stars_disc)/N_stars_disc_tot
plt.step(edges, np.append(y_gas[0],y_gas), color='b', label=r'Cold Gas')
plt.step(edges, np.append(y_stars[0],y_stars), color='r', label=r'Stars')
#plt.step(edges, np.append(y_stars_disc[0],y_stars_disc), color='m', label=r'Disc Stars')
plt.step(edges, np.append(y_stars_disc2[0],y_stars_disc2), color='m', ls='-', label=r'Disc Stars')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
plt.xlabel(r'$\log_{10}(M_{\rm stars}\ \mathrm{or}\ M_{\rm gas}\ [\mathrm{M}_{\bigodot}])$')
plt.ylabel('Fraction of galaxies with offset $> 20^{\circ}$')
gp.savepng(outdir+'MassDistOffset', xpixplot=768, ypixplot=512)

# As above but absolute number
gp.figure()
plt.step(edges, np.append(N_gas[0],N_gas), color='b', label=r'Cold Gas')
plt.step(edges, np.append(N_stars[0],N_stars), color='r', label=r'Stars')
plt.step(edges, np.append(N_stars_disc[0],N_stars_disc), color='m', label=r'Disc Stars')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
plt.xlabel(r'$\log_{10}(M_{\rm stars}\ \mathrm{or}\ M_{\rm gas}\ [\mathrm{M}_{\bigodot}])$')
plt.ylabel('Number of galaxies with offset $> 20^{\circ}$')
gp.savepng(outdir+'MassDistOffsetAbs', xpixplot=768, ypixplot=512)


# Offset distribution based on relative masses of the components
discfrac = G.ColdGas[filt] / (G.StellarMass[filt] - G.SecularBulgeMass[filt] - G.ClassicalBulgeMass[filt])
gos_tot, edges = np.histogram(discfrac, bins=25, range=[0,1]) # gas on stars
gos_off, edges = np.histogram(discfrac[theta>20], bins=25, range=[0,1])
gos_frac = (1.0*gos_off)/gos_tot
sog_tot, edges = np.histogram(1/discfrac, bins=25, range=[0,1]) # stars on gas
sog_off, edges = np.histogram(1/discfrac[theta>20], bins=25, range=[0,1])
sog_frac = (1.0*sog_off)/sog_tot
gp.figure()
plt.step(edges, np.append(gos_frac[0],gos_frac), color='r', label=r'$M_{\rm cold} / M_{\rm star}$')
plt.step(edges, np.append(sog_frac[0],sog_frac), color='b', label=r'$M_{\rm star} / M_{\rm cold}$')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
plt.xlabel(r'Disc mass ratio')
plt.ylabel('Fraction of galaxies with offset $> 20^{\circ}$')
gp.savepng(outdir+'MassRatioOffset', xpixplot=768, ypixplot=512)

# As above but absolute number
gp.figure()
plt.step(edges, np.append(gos_off[0],gos_off), color='r', label=r'$M_{\rm cold} / M_{\rm star}$')
plt.step(edges, np.append(sog_off[0],sog_off), color='b', label=r'$M_{\rm star} / M_{\rm cold}$')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
plt.xlabel(r'Disc mass ratio')
plt.ylabel('Number of galaxies with offset $> 20^{\circ}$')
gp.savepng(outdir+'MassRatioOffsetAbs', xpixplot=768, ypixplot=512)

plt.close("all")

"""
    # Cumulative mass profiles of z=0 galaxies
    #filt = (G.ColdGas*1e10/h >= 1e8) * (G.ColdGas*1e10/h <= 1e12) * ((G.StellarMass-G.SecularBulgeMass-G.ClassicalBulgeMass)*1e10/h >= 1e8) * ((G.StellarMass-G.SecularBulgeMass-G.ClassicalBulgeMass)*1e10/h <= 1e12) * (G.Vvir>100) * (G.Vvir<300)
    filt = (G.Vvir>=200) * (G.Vvir<=235) * np.isfinite(G.Vvir) * (G.StellarMass/h > 1.0) * (G.ColdGas/h > 10**-0.5) * (G.Type==0)
    gals = np.random.random_integers(0, len(G.Type[filt])-1, 20)
    gals = np.unique(gals)
    gp.figure()
    for gal in gals:
    xplot = DiscBinEdge[1:]/G.Vvir[filt][gal]/(G.Rvir[filt][gal]*1e3/h)
    try:
    xcut = np.where(xplot >= G.Rvir[filt][gal]/h)[0][0] + 1
    except:
    continue
    if xcut<29: continue
    xplot = xplot[:xcut]
    xplot /= np.max(xplot)
    #yplot = np.cumsum((G.DiscGas[filt][gal,:]+G.DiscStars[filt][gal,:])/(G.ColdGas[filt][gal]+G.StellarMass[filt][gal]-G.SecularBulgeMass[filt][gal]-G.ClassicalBulgeMass[filt][gal])) # Haven't included bulge masses yet which may cause an issue!
    yplot = np.cumsum((G.DiscGas[filt][gal,:xcut]+G.DiscStars[filt][gal,:xcut]))
    yplot /= np.max(yplot)
    plt.plot(xplot, yplot, lw=2)
    plt.ylabel(r'$M_{\rm bary}(<r) / M_{\rm bary}(<r_{\rm vir})$')
    plt.xlabel(r'$r / r_{\rm vir}$')
    plt.axis([0.1,0.7,0.101,1])
    gp.savepng(outdir+'CumBary0', xpixplot=768, ypixplot=512)
    """



# How were the stars formed?
filt = (G.StellarMass>0.0)
InSituFrac = G.StarsInSitu[filt] / G.StellarMass[filt]
InstabFrac = G.StarsInstability[filt] / G.StellarMass[filt]
MergeFrac = G.StarsMergeBurst[filt] / G.StellarMass[filt]
print 'Summed fractions', InSituFrac+InstabFrac+MergeFrac
logM = np.log10(G.StellarMass[filt]*1e10/h)
#
Nbins = 25
range = [8.5, 12.0]
InSitu_bin, xedge = np.histogram(logM, bins=Nbins, range=range, weights=InSituFrac)
Instab_bin, xedge = np.histogram(logM, bins=Nbins, range=range, weights=InstabFrac)
Merge_bin, xedge = np.histogram(logM, bins=Nbins, range=range, weights=MergeFrac)
N_bin, xedge = np.histogram(logM, bins=Nbins, range=range)
xplot = (xedge[1:]+xedge[:-1])/2.
#
gp.figure()
plt.plot(xplot, InSitu_bin/N_bin, 'r-', label=r'In situ', lw=2)
plt.plot(xplot, Instab_bin/N_bin, 'g-.', label=r'Instability', lw=3)
plt.plot(xplot, Merge_bin/N_bin, 'b--', label=r'Merger Starburst', lw=2)
#
plt.xlabel(r'$\log_{10}(M_{\rm stars}\ [\mathrm{M}_{\bigodot}])$')
plt.ylabel(r'Fraction of star formation channel')
plt.legend(loc='best', frameon=False, fontsize=fsize-6)
plt.axis([8.51,12,0,1])
gp.savepng(outdir+'StarFractions', xpixplot=768, ypixplot=512)



# Baryon fraction breakdown
filt = (G.Mvir/h>1.0)
Mvir = G.Mvir[filt]
Ngal = len(G.Mvir[filt])
#
tot = (G.HotGas[filt] + G.ColdGas[filt] + G.EjectedMass[filt] + G.StellarMass[filt] + G.IntraClusterStars[filt] + G.BlackHoleMass[filt]) / Mvir
hot = G.HotGas[filt] / Mvir
cold = G.ColdGas[filt] / Mvir
ejected = G.EjectedMass[filt] / Mvir
stars = G.StellarMass[filt] / Mvir
ics = G.IntraClusterStars[filt] / Mvir
#
range = [11.0, np.log10(np.max(Mvir)*1e10/h)]
Nbins = 15
#
tot_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range, weights=tot)
hot_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range, weights=hot)
cold_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range, weights=cold)
ejected_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range, weights=ejected)
stars_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range, weights=stars)
ics_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range, weights=ics)
N_bin, xedge = np.histogram(np.log10(Mvir*1e10/h), bins=Nbins, range=range)
xplot = (xedge[1:]+xedge[:-1])/2.
#
tot_low, tot_high = np.zeros(Nbins), np.zeros(Nbins)
for i in xrange(Nbins):
    filt = (np.log10(Mvir*1e10/h) >= xedge[i]) * (np.log10(Mvir*1e10/h) < xedge[i+1])
    binvals = np.sort(tot[filt])
    if len(binvals)==0: continue
    i_low = int(len(binvals)*0.16)
    i_high = int(len(binvals)*0.84)
    tot_low[i] = binvals[i_low]
    tot_high[i] = binvals[i_high]
#
gp.figure()
plt.fill_between(xplot, tot_high, tot_low, color='k', alpha=0.2)
plt.plot(xplot, tot_bin/N_bin, 'k-', lw=1, label=r'Total')
plt.plot(xplot, stars_bin/N_bin, 'm--', lw=2, label=r'Stars')
plt.plot(xplot, ics_bin/N_bin, 'y-', lw=2, label=r'Intra-cluster stars', dashes=[6,3,2,3])
plt.plot(xplot, hot_bin/N_bin, 'r-', lw=2, label=r'Hot gas')
plt.plot(xplot, cold_bin/N_bin, 'b:', lw=3, label=r'Cold gas')
plt.plot(xplot, ejected_bin/N_bin, 'g-', lw=2, label=r'Ejected gas', dashes=[5,3,2,3,2,3])
plt.axis([np.min(xplot), np.max(xplot), 0, np.max(tot_high)])
plt.xlabel(r'$\log_{10}(M_{\rm vir}\ [\mathrm{M}_{\bigodot}])$')
plt.ylabel(r'Baryon fraction')
plt.legend(loc='best', frameon=False, fontsize=fsize-6)
gp.savepng(outdir+'baryon_frac', xpixplot=768, ypixplot=512)


# z=2 Stellar mass function
G = gr.sagesnap('model_z2.070', 0, 7, indir, disc=True)
gp.figure()
gp.massfunction(G.StellarMass*1e10/h, boxlen, extra=0, h=h, step=False, fsize=fsize, binwidth=0.1, label=r'\textsc{dark sage}')
gp.smf_nifty_obs(h, z=2, haxes=False)
plt.legend(fontsize=fsize-4, loc='best', frameon=False, title=r'$z=2$')
plt.axis([8,11.8,5e-6,0.2])
gp.savepng(outdir+'SMF_z2', xpixplot=768, ypixplot=512)
"""
    # Cumulative mass profiles of z=2 galaxies
    filt = (G.ColdGas*1e10/h >= 1e8) * (G.ColdGas*1e10/h <= 1e12) * ((G.StellarMass-G.SecularBulgeMass-G.ClassicalBulgeMass)*1e10/h >= 1e8) * ((G.StellarMass-G.SecularBulgeMass-G.ClassicalBulgeMass)*1e10/h <= 1e12) * (G.Vvir>100) * (G.Vvir<300)
    gals = np.random.random_integers(0, len(G.Type[filt])-1, 50)
    gals = np.unique(gals)
    gp.figure()
    for gal in gals:
    xplot = DiscBinEdge[1:]/G.Vvir[filt][gal]/(G.Rvir[filt][gal]*1e3/h)
    try:
    xcut = np.where(xplot >= G.Rvir[filt][gal]/h)[0][0] + 1
    except:
    continue
    if xcut<10: continue
    xplot = xplot[:xcut]
    xplot /= np.max(xplot)
    #yplot = np.cumsum((G.DiscGas[filt][gal,:]+G.DiscStars[filt][gal,:])/(G.ColdGas[filt][gal]+G.StellarMass[filt][gal]-G.SecularBulgeMass[filt][gal]-G.ClassicalBulgeMass[filt][gal])) # Haven't included bulge masses yet which may cause an issue!
    yplot = np.cumsum((G.DiscGas[filt][gal,:xcut]+G.DiscStars[filt][gal,:xcut]))
    yplot /= np.max(yplot)
    plt.plot(xplot, yplot)
    plt.ylabel(r'$M_{\rm bary}(<r) / M_{\rm bary}(<r_{\rm vir})$')
    plt.xlabel(r'$r / r_{\rm vir}$')
    plt.axis([0.201,1,0.3,1])
    gp.savepng(outdir+'CumBary2', xpixplot=768, ypixplot=512)
    """


# SFR function
G = gr.sagesnap('model_z0.144', 0, 7, indir, disc=True)
gp.figure()
gp.massfunction(G.SfrDisk+G.SfrBulge, boxlen, label=r'\textsc{dark sage}', step=False, range=[-0.5,2.5], Nbins=6)
gp.SFR_function_obs(h)
gp.savepng(outdir+'SFRF', xpixplot=768, ypixplot=512)

### BH mass history plot
#zstr = ['0.000', '0.089', '0.208', '0.362', '0.564', '0.828', '1.173', '1.630', '2.239', '3.060', '4.179', '5.724', '7.883', '10.944', '15.343']
zstr = ['0.000', '0.144', '0.280', '0.509', '0.755', '1.078', '1.504', '2.070', '2.619', '3.060', '3.576', '4.179', '4.888', '5.289', '5.724', '6.197', '6.712', '7.272']
if os.path.isfile(indir+'model_z'+zstr[1]+'_0'):
    z = np.array(zstr, dtype='f4')
    SFRD = np.zeros(len(z))
    SMsum = np.zeros(len(z))
    t = np.zeros(len(z))
    gp.figure()
    for s in xrange(len(z)):
        M = gr.sagesnap('model_z'+zstr[s], 0, 7, indir, disc=True)
        SFRD[s] = np.log10(np.sum(M.SfrDisk + M.SfrBulge) / (boxlen**3))
        SMsum[s] = np.sum(M.StellarMass*1e10/h) + np.sum(M.IntraClusterStars*1e10/h)
        t[s] = gc.z2t(z[s], 100*h, 0, 0.25, 0.75)*1e9 # Time in years since Big Bang
        gp.point_with_spread(z[s], M.BlackHoleMass[M.BlackHoleMass>0]*1e10/h, mean=False, alpha=0.8)
    plt.xlabel(r'Redshift')
    plt.yscale('log')
    plt.ylabel(r'$\log_{10} (M_{\rm BH}\ [\mathrm{M}_{\bigodot}])$')
    plt.axis([0, 8, 1e3, 1e8])
    gp.savepng(outdir+'BHhistory', xpixplot=768, ypixplot=512)



### Madau plot
if os.path.isfile(indir+'model_z'+zstr[1]+'_0'):
    gp.figure()
    gp.SFRD_obs(h)
    plt.plot(z, SFRD, 'k-', lw=2, label=r'\textsc{dark sage}')
    plt.axis([0,7.5,-2.8,-0.5])
    plt.xlabel(r'Redshift')
    plt.ylabel(r'$\log_{10} \left(\bar{\rho}_{\rm SFR}\ [\mathrm{M}_{\bigodot}\ \mathrm{yr}^{-1}\ \mathrm{cMpc}^{-3}] \right)$')
    plt.legend(fontsize=fsize-4, loc='best', frameon=False, numpoints=1)
    gp.savepng(outdir+'SFRD', xpixplot=768, ypixplot=512)

### Star formation history build-up versus actual total stellar masses
SFRH = np.append(0, 10**SFRD[::-1] * (boxlen)**3)
time = np.append(0, t[::-1])
SM_integrated = 0.57 * np.cumsum(np.diff(time)*(SFRH[1:]+SFRH[:-1])/2)[::-1]
#
gp.figure()
plt.plot(z, SMsum, 'k-', label=r'Actual stellar mass', lw=2)
plt.plot(z, SM_integrated, 'b--', label=r'Integrated SFR', lw=2)
plt.xlabel(r'Redshift')
plt.ylabel(r'$M_{\rm stars}\ [\mathrm{M}_{\bigodot}]$')
plt.yscale('log')
plt.legend(fontsize=fsize-4, loc='best', frameon=False)
gp.savepng(outdir+'StellarHistory', xpixplot=768, ypixplot=512)


plt.close("all")