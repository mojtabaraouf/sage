# Routines used for reading and plotting data (e.g. used by plot_constraints_z0.py).  Many of these are copied from the arhstevens/Dirty-AstroPy GitHub repository.

from pylab import *


def galdtype_darksage(Nannuli=30):
    floattype = np.float32
    Galdesc_full = [
                    ('Type'                         , np.int32),
                    ('GalaxyIndex'                  , np.int64),
                    ('HaloIndex'                    , np.int32),
                    ('SimulationHaloIndex'          , np.int32),
                    ('TreeIndex'                    , np.int32),
                    ('SnapNum'                      , np.int32),
                    ('CentralGalaxyIndex'           , np.int64),
                    ('CentralMvir'                  , floattype),
                    ('mergeType'                    , np.int32),
                    ('mergeIntoID'                  , np.int32),
                    ('mergeIntoSnapNum'             , np.int32),
                    ('dT'                           , floattype),
                    ('Pos'                          , (floattype, 3)),
                    ('Vel'                          , (floattype, 3)),
                    ('Spin'                         , (floattype, 3)),
                    ('Len'                          , np.int32),
                    ('LenMax'                       , np.int32),
                    ('Mvir'                         , floattype),
                    ('Rvir'                         , floattype),
                    ('Vvir'                         , floattype),
                    ('Vmax'                         , floattype),
                    ('VelDisp'                      , floattype),
                    ('DiscRadii'                    , (floattype, Nannuli+1)), 
                    ('ColdGas'                      , floattype),
                    ('StellarMass'                  , floattype),
                    ('MergerBulgeMass'              , floattype),
                    ('InstabilityBulgeMass'          , floattype),
                    ('HotGas'                       , floattype),
                    ('EjectedMass'                  , floattype),
                    ('BlackHoleMass'                , floattype),
                    ('IntraClusterStars'            , floattype),
                    ('DiscGas'                      , (floattype, Nannuli)),
                    ('DiscStars'                    , (floattype, Nannuli)),
                    ('SpinStars'                    , (floattype, 3)),
                    ('SpinGas'                      , (floattype, 3)),
                    ('SpinClassicalBulge'           , (floattype, 3)),
                    ('StarsInSitu'                  , floattype),
                    ('StarsInstability'             , floattype),
                    ('StarsMergeBurst'              , floattype),
                    ('DiscHI'                       , (floattype, Nannuli)),
                    ('DiscH2'                       , (floattype, Nannuli)),
                    ('DiscSFR'                      , (floattype, Nannuli)), 
                    ('MetalsColdGas'                , floattype),
                    ('MetalsStellarMass'            , floattype),
                    ('ClassicalMetalsBulgeMass'     , floattype),
                    ('SecularMetalsBulgeMass'       , floattype),
                    ('MetalsHotGas'                 , floattype),
                    ('MetalsEjectedMass'            , floattype),
                    ('MetalsIntraClusterStars'      , floattype),
                    ('DiscGasMetals'                , (floattype, Nannuli)),
                    ('DiscStarsMetals'              , (floattype, Nannuli)),
                    ('SfrFromH2'                    , floattype),
                    ('SfrInstab'                    , floattype),
                    ('SfrMergeBurst'                , floattype),
                    ('SfrDiskZ'                     , floattype),
                    ('SfrBulgeZ'                    , floattype),
                    ('DiskScaleRadius'              , floattype),
                    ('CoolScaleRadius'              , floattype), 
                    ('StellarDiscScaleRadius'       , floattype),
                    ('Cooling'                      , floattype),
                    ('Heating'                      , floattype),
                    ('LastMajorMerger'              , floattype),
                    ('LastMinorMerger'              , floattype),
                    ('OutflowRate'                  , floattype),
                    ('infallMvir'                   , floattype),
                    ('infallVvir'                   , floattype),
                    ('infallVmax'                   , floattype)
                    ]
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    return Galdesc



def darksage_out_single(fname, fields=[], Nannuli=30):
    # Read a single Dark Sage output file, returning all the galaxy data
    # fname is the full name for the file to read, including its path
    # fields is the list of fields you want to read in.  If empty, will read all fields.
    
    Galdesc = galdtype_darksage(Nannuli)
    if len(fields)==0: fields=list(Galdesc.names)
    
    fin = open(fname, 'rb')  # Open the file
    Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
    NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
    GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
    G = np.fromfile(fin, Galdesc, NtotGals) # Read all the galaxy data
    G = G[fields]
    return G 



def darksage_snap(fpre, filelist, fields=[], Nannuli=30):
    # Read full Dark Sage snapshot, going through each file and compiling into 1 array
    # fpre is the name of the file up until the _ before the file number
    # filelist contains all the file numbers you want to read in
    
    Galdesc = galdtype_darksage()
    Glist = []
    Ngal = np.array([],dtype=np.int32)
    G = darksage_out_single(fpre+'_'+str(filelist[0]), fields, Nannuli)
    
    for i in filelist[1:]:
        G1 = darksage_out_single(fpre+'_'+str(i), fields)
        G = np.append(G, G1)
    return G



def massfunction(mass, Lbox, range=[8,12.5], c='k', lw=2, ls='-', label='', ax=None):
    masslog = np.log10(mass[(mass>0)*np.isfinite(mass)])
    N, edges = np.histogram(masslog, bins=np.arange(range[0],range[1]+0.1,0.1))
    binwidth = edges[1]-edges[0]
    x = edges[:-1] + binwidth/2
    y = N/(binwidth*Lbox**3)
    
    if ax is None: ax = plt.gca()
    
    if len(label)>0:
        ax.plot(x, y, c+ls, linewidth=lw, label=label)
    else:
        ax.plot(x, y, c+ls, linewidth=lw)


def schechter(phistar, Mstar, alpha, Mlog=False, range=[7,12], Npoints=2000, logM=None):
    if Mlog: Mstar = 10**Mstar
    if logM is None: logM = np.linspace(range[0],range[1],Npoints)
    M = 10**logM
    Phi = np.log(10.) * (phistar) * (M/Mstar)**(alpha+1) * np.exp(-M/Mstar)
    return Phi, logM


def stellar_massfunction_obsdata(h=0.678, ax=None):
    B = np.array([
                  [7.05, 1.3531e-01, 6.0741e-02],
                  [7.15, 1.3474e-01, 6.0109e-02],
                  [7.25, 2.0971e-01, 7.7965e-02],
                  [7.35, 1.7161e-01, 3.1841e-02],
                  [7.45, 2.1648e-01, 5.7832e-02],
                  [7.55, 2.1645e-01, 3.9988e-02],
                  [7.65, 2.0837e-01, 4.8713e-02],
                  [7.75, 2.0402e-01, 7.0061e-02],
                  [7.85, 1.5536e-01, 3.9182e-02],
                  [7.95, 1.5232e-01, 2.6824e-02],
                  [8.05, 1.5067e-01, 4.8824e-02],
                  [8.15, 1.3032e-01, 2.1892e-02],
                  [8.25, 1.2545e-01, 3.5526e-02],
                  [8.35, 9.8472e-02, 2.7181e-02],
                  [8.45, 8.7194e-02, 2.8345e-02],
                  [8.55, 7.0758e-02, 2.0808e-02],
                  [8.65, 5.8190e-02, 1.3359e-02],
                  [8.75, 5.6057e-02, 1.3512e-02],
                  [8.85, 5.1380e-02, 1.2815e-02],
                  [8.95, 4.4206e-02, 9.6866e-03],
                  [9.05, 4.1149e-02, 1.0169e-02],
                  [9.15, 3.4959e-02, 6.7898e-03],
                  [9.25, 3.3111e-02, 8.3704e-03],
                  [9.35, 3.0138e-02, 4.7741e-03],
                  [9.45, 2.6692e-02, 5.5029e-03],
                  [9.55, 2.4656e-02, 4.4359e-03],
                  [9.65, 2.2885e-02, 3.7915e-03],
                  [9.75, 2.1849e-02, 3.9812e-03],
                  [9.85, 2.0383e-02, 3.2930e-03],
                  [9.95, 1.9929e-02, 2.9370e-03],
                  [10.05, 1.8865e-02, 2.4624e-03],
                  [10.15, 1.8136e-02, 2.5208e-03],
                  [10.25, 1.7657e-02, 2.4217e-03],
                  [10.35, 1.6616e-02, 2.2784e-03],
                  [10.45, 1.6114e-02, 2.1783e-03],
                  [10.55, 1.4366e-02, 1.8819e-03],
                  [10.65, 1.2588e-02, 1.8249e-03],
                  [10.75, 1.1372e-02, 1.4436e-03],
                  [10.85, 9.1213e-03, 1.5816e-03],
                  [10.95, 6.1125e-03, 9.6735e-04],
                  [11.05, 4.3923e-03, 9.6254e-04],
                  [11.15, 2.5463e-03, 5.0038e-04],
                  [11.25, 1.4298e-03, 4.2816e-04],
                  [11.35, 6.4867e-04, 1.6439e-04],
                  [11.45, 2.8294e-04, 9.9799e-05],
                  [11.55, 1.0617e-04, 4.9085e-05],
                  [11.65, 3.2702e-05, 2.4546e-05],
                  [11.75, 1.2571e-05, 1.2571e-05],
                  [11.85, 8.4589e-06, 8.4589e-06],
                  [11.95, 7.4764e-06, 7.4764e-06]
                  ], dtype=np.float32)
    if ax is None: ax = plt.gca()
    ax.fill_between(B[:,0]+np.log10(0.7**2)-np.log10(h**2), (B[:,1]+B[:,2])*h**3, (B[:,1]-B[:,2])*h**3, facecolor='purple', alpha=0.2)
    ax.plot([1,1], [1,2], color='purple', linewidth=8, alpha=0.3, label=r'Baldry et al.~(2008)') # Just for the legend



def HIH2_massfunction_obsdata(h=0.678, HI=True, H2=True, K=True, OR=False, ax=None, Z=True, M=False, B=False):
    Zwaan = np.array([[6.933,   -0.333],
                      [7.057,   -0.490],
                      [7.209,   -0.698],
                      [7.365,   -0.667],
                      [7.528,   -0.823],
                      [7.647,   -0.958],
                      [7.809,   -0.917],
                      [7.971,   -0.948],
                      [8.112,   -0.927],
                      [8.263,   -0.917],
                      [8.404,   -1.062],
                      [8.566,   -1.177],
                      [8.707,   -1.177],
                      [8.853,   -1.312],
                      [9.010,   -1.344],
                      [9.161,   -1.448],
                      [9.302,   -1.604],
                      [9.448,   -1.792],
                      [9.599,   -2.021],
                      [9.740,   -2.406],
                      [9.897,   -2.615],
                      [10.053,  -3.031],
                      [10.178,  -3.677],
                      [10.335,  -4.448],
                      [10.492,  -5.083]])
        
    Martin_data = np.array([[6.302,    -0.504],
                              [6.500,    -0.666],
                              [6.703,    -0.726],
                              [6.904,    -0.871],
                              [7.106,    -1.135],
                              [7.306,    -1.047],
                              [7.504,    -1.237],
                              [7.703,    -1.245],
                              [7.902,    -1.254],
                              [8.106,    -1.414],
                              [8.306,    -1.399],
                              [8.504,    -1.476],
                              [8.705,    -1.591],
                              [8.906,    -1.630],
                              [9.104,    -1.695],
                              [9.309,    -1.790],
                              [9.506,    -1.981],
                              [9.707,    -2.141],
                              [9.905,    -2.317],
                              [10.108,    -2.578],
                              [10.306,    -3.042],
                              [10.509,    -3.780],
                              [10.703,    -4.534],
                              [10.907,    -5.437]])     
                    
    Martin_mid = Martin_data[:,1] + 3*np.log10(h/0.7)
    Martin_x = Martin_data[:,0] + 2*np.log10(0.7/h)
    Martin_high = np.array([-0.206, -0.418, -0.571, -0.725, -1.003, -0.944, -1.144, -1.189, -1.189, -1.358, -1.344, -1.417, -1.528, -1.586, -1.651, -1.753, -1.925, -2.095, -2.281, -2.537, -3.003, -3.729, -4.451, -5.222]) + 3*np.log10(h/0.7)
    Martin_low = np.array([-0.806, -0.910, -0.885, -1.019, -1.268, -1.173, -1.313, -1.314, -1.320, -1.459, -1.443, -1.530, -1.647, -1.669, -1.736, -1.838, -2.021, -2.191, -2.359, -2.621, -3.098, -3.824, -4.618, -5.663]) + 3*np.log10(h/0.7)
      
    Keres_high = np.array([-1.051, -1.821, -1.028, -1.341, -1.343, -1.614, -1.854, -2.791,  -3.54 , -5.021]) + 3*np.log10(h)
    Keres_mid = np.array([-1.271, -1.999, -1.244, -1.477, -1.464, -1.713, -1.929, -2.878,   -3.721, -5.22 ]) + 3*np.log10(h)
    Keres_low = np.array([-1.706, -2.302, -1.71 , -1.676, -1.638, -1.82 , -2.033, -2.977,   -4.097, -5.584]) + 3*np.log10(h)
    Keres_M = np.array([  6.953,   7.353,   7.759,   8.154,   8.553,   8.96 ,   9.365,  9.753,  10.155,  10.558]) - 2*np.log10(h)
      
    ObrRaw_high = np.array([-0.905, -1.122, -1.033, -1.1  , -1.242, -1.418, -1.707, -2.175, -2.984, -4.868]) + 3*np.log10(h)
    ObrRaw_mid = np.array([-1.116, -1.308, -1.252, -1.253, -1.373, -1.509, -1.806, -2.261,  -3.198, -5.067]) + 3*np.log10(h)
    ObrRaw_low = np.array([-1.563, -1.602, -1.73 , -1.448, -1.537, -1.621, -1.918, -2.369,  -3.556, -5.413]) + 3*np.log10(h)
    ObrRaw_M = np.array([ 7.301,  7.586,  7.862,  8.133,  8.41 ,  8.686,  8.966,  9.242,    9.514,  9.788]) - 2*np.log10(h)
      
    HI_x = Zwaan[:,0] - 2*np.log10(h)
    HI_y = 10**Zwaan[:,1] * h**3

    Boselli_const_XCO = np.array([[7.39189, -3.06989, -3.32527, -2.86828],
                                [7.78378, -2.45161, -2.54570, -2.37097],
                                [8.18919, -1.91398, -1.96774, -1.84677],
                                [8.62162, -2.12903, -2.20968, -2.03495],
                                [9.01351, -2.41129, -2.51882, -2.31720],
                                [9.41892, -2.62634, -2.80108, -2.53226],
                                [9.81081, -2.73387, -2.85484, -2.54570],
                                [10.2297, -3.64785, -5.97312, -3.36559]])

    Boselli_var_XCO = np.array([[7.59030, -3.19086, -3.58065, -2.98925],
                              [7.98113, -2.55914, -2.72043, -2.45161],
                              [8.37197, -2.22312, -2.30376, -2.14247],
                              [8.78976, -1.94086, -1.99462, -1.90054],
                              [9.18059, -1.98118, -2.06183, -1.90054],
                              [9.59838, -2.72043, -2.92204, -2.62634], 
                              [9.98922, -3.67473, -5.98656, -3.31183]])
    
    if ax is None: ax = plt.gca()         
    if HI and Z: ax.plot(HI_x, HI_y, '-', color='g', lw=8, alpha=0.4, label=r'Zwaan et al.~(2005)')
    if HI and M: ax.fill_between(Martin_x, 10**Martin_high, 10**Martin_low, color='c', alpha=0.4)
    if HI and M: ax.plot([0,1], [1,1], 'c-', lw=8, alpha=0.4, label=r'Martin et al.~(2010)')

    if H2 and K: ax.fill_between(Keres_M, 10**Keres_high, 10**Keres_low, color='teal', alpha=0.4)
    if H2 and K: ax.plot([0,1], [1,1], '-', color='teal', lw=8, alpha=0.4, label=r'Keres et al.~(2003)')

    if H2 and OR: ax.fill_between(ObrRaw_M, 10**ObrRaw_high, 10**ObrRaw_low, color='darkcyan', alpha=0.4)
    if H2 and OR: ax.plot([0,1], [1,1], '-', color='darkcyan', lw=8, alpha=0.4, label=r'Obreschkow \& Rawlings (2009)')
                        
    if H2 and B: ax.fill_between(Boselli_const_XCO[:,0]+2*np.log10(0.7/h), 10**Boselli_const_XCO[:,2]*(h/0.7)**3/0.4, 10**Boselli_const_XCO[:,3]*(h/0.7)**3/0.4, color='orange', alpha=0.4)
    if H2 and B: ax.plot([0,1], [1,1], '-', color='orange', lw=8, alpha=0.4, label=r'Boselli et al.~(2014), const.~$X_{\rm CO}$')
    if H2 and B: ax.fill_between(Boselli_var_XCO[:,0]+2*np.log10(0.7/h), 10**Boselli_var_XCO[:,2]*(h/0.7)**3/0.4, 10**Boselli_var_XCO[:,3]*(h/0.7)**3/0.4, color='violet', alpha=0.4)
    if H2 and B: ax.plot([0,1], [1,1], '-', color='violet', lw=8, alpha=0.4, label=r'Boselli et al.~(2014), var.~$X_{\rm CO}$')

    ax.set_xlabel(r'$\log_{10}(M_{\mathrm{H}\,\huge\textsc{i}}\ \mathrm{or}\ M_{\mathrm{H}_2}\ [\mathrm{M}_{\bigodot}])$')
    ax.set_ylabel(r'$\Phi\ [\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1}]$')
    ax.axis([8,11.5,1e-6,1e-1])



def savepng(filename, xsize=1024, ysize=None, fig=None, transparent=False):
    # Save a figure as a PNG with a normalised size / aspect ratio
    xpix = 2560
    ypix = 1440
    ss = 27
    
    if ysize==None: ysize = int(xsize*9./16)
    
    mydpi = np.sqrt(xpix**2 + ypix**2)/ss 
    xinplot = xsize*(9./7.)/mydpi
    yinplot = ysize*(9./7.)/mydpi
    if fig is None: fig = plt.gcf()
    fig.set_size_inches(xinplot,yinplot)
    fig.set_dpi(mydpi)
    
    filename = str(filename)
    if filename[-4:] != '.png':
        filename = filename+'.png'
    fig.savefig(filename, dpi=mydpi, bbox_inches='tight', transparent=transparent)

