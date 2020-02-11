# Plots that require multiple snapshots of Dark Sage. Potentially usable for calibration. First version just does the Madau plot.

from pylab import *
import os
import routines as r
import random

# Warnings are annoying
import warnings
warnings.filterwarnings("ignore")


###### USER NEEDS TO SET THESE THINGS ######
indir = 'results/millennium/' # directory where the Dark Sage data are
sim = 0 # which simulation Dark Sage has been run on -- if it's new, you will need to set its defaults below.
#   0 = Mini Millennium, 1 = Full Millennium, 2 = SMDPL

Nannuli = 30 # number of annuli used for discs in Dark Sage
###### ============================== ######



##### SET PLOTTING DEFAULTS #####
fsize = 26
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman', 'legend.numpoints': 1, 'legend.columnspacing': 1, 'legend.fontsize': fsize-4, 'xtick.top': True, 'ytick.right': True})

NpartMed = 100 # minimum number of particles for finding relevant medians for minima on plots

outdir = indir+'plots/' # where the plots will be saved
if not os.path.exists(outdir): os.makedirs(outdir)
######  =================== #####



##### SEARCH DIRECTORY FOR SNAPSHOTS AND FILENUMBERS PRESENT #####
filenames = os.listdir(indir)
redshiftstr = np.array([], dtype=str)
filenumbers = np.array([], dtype=np.int32)
fpre = None
for f in filenames:
    w1 = f.find('_z')
    if w1==-1: continue
    w2 = f.find('_', w1+3)
    zstr = f[w1+2:w2]
    fno = int(f[w2+1:])
    if zstr not in redshiftstr: redshiftstr = np.append(redshiftstr, zstr)
    if fno not in filenumbers: filenumbers = np.append(filenumbers, fno)
    if fpre is None: fpre = f[:w1+2]
redshifts = redshiftstr.astype(np.float32)
args = np.argsort(redshifts) # want redshifts to be listed in order
redshiftstr = redshiftstr[args]
redshifts = redshifts[args]
##### ====================================================== #####



##### SIMULATION DEFAULTS #####
if sim==0:
    h = 0.73
    vol = (62.5/h)**3 * len(filenumbers)/8. # comoving volume of the (part of the) simulation
elif sim==1:
    h = 0.73
    vol = (500.0/h)**3 * len(filenumbers)/512.
elif sim==2:
    h = 0.6777
    vol = (400.0/h)**3 * len(filenumbers)/1000.
# add here 'elif sim==3:' etc for a new simulation
else:
    print('Please specify a valid simulation.  You may need to add its defaults to this code.')
    quit()
######  ================= #####



##### READ DARK SAGE DATA AND BUILD RELEVANT ARRAYS #####
Nsnap = len(redshifts)
SFRD = np.zeros(Nsnap)
SFRD_resolved = np.zeros(Nsnap)

for i in range(Nsnap):
    G = r.darksage_snap(indir+fpre+redshiftstr[i], filenumbers, Nannuli=Nannuli)
    SFRD[i] = np.sum(G['SfrFromH2']+G['SfrInstab']+G['SfrMergeBurst']) / vol
    SFRD_resolved[i] = np.sum((G['SfrFromH2']+G['SfrInstab']+G['SfrMergeBurst'])[G['LenMax']>=100]) / vol
##### ============================================= #####




##### PLOT 1: MADAU-LILLY DIAGRAM (UNIVERSAL STAR FORMATION RATE DENSITY HISTORY) #####
try:
    fig, ax  = plt.subplots(1, 1)
    plt.plot(1+redshifts, np.log10(SFRD), 'k--', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p}\!\geq\!20$')
    plt.plot(1+redshifts, np.log10(SFRD_resolved), 'k.-', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p,max}\!\geq\!100$')
    r.SFRD_obs(h, plus=1)
    plt.xlabel(r'Redshift')
    plt.ylabel(r'$\log_{10}\left( \bar{\rho}_{\rm SFR}~[{\rm M}_{\odot}\, {\rm yr}^{-1}\, {\rm cMpc}^{-3}] \right)$')
    plt.xscale('log')
    plt.axis([1,9,-2.8,0.5])
    plt.minorticks_off()
    plt.xticks(range(1,10), (str(i) for i in range(9)))
    plt.legend(loc='best', frameon=False, ncol=2)
    fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
    r.savepng(outdir+'H1-SFRDH', xsize=768, ysize=400)
except Exception as excptn:
    print('Unexpected issue with plot H1: {0}'.format(excptn))

##### =========================================================================== #####

