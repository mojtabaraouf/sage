"""
Run this test once you have installed Dark Sage to ensure everything is working correctly.
This will fetch pre-made Dark Sage output, run the code on your machine, then compare the outputs.
This script should run straight out of the box with "python test.py"
"""

from __future__ import print_function
import os
import sys
try: # Python 2
    from urllib import urlretrieve
except ImportError: # Python 3
    from urllib.request import urlretrieve
import filecmp
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    xrange
except NameError:
    xrange = range

def galdtype():
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
        ('DiscRadii'                    , (floattype, 31)),
        ('ColdGas'                      , floattype),
        ('StellarMass'                  , floattype),
        ('ClassicalBulgeMass'           , floattype),
        ('SecularBulgeMass'             , floattype),
        ('HotGas'                       , floattype),
        ('EjectedMass'                  , floattype),
        ('BlackHoleMass'                , floattype),
        ('IntraClusterStars'            , floattype),
        ('DiscGas'                      , (floattype, 30)),
        ('DiscStars'                    , (floattype, 30)),
        ('SpinStars'                    , (floattype, 3)),
        ('SpinGas'                      , (floattype, 3)),
        ('SpinClassicalBulge'           , (floattype, 3)),
        ('StarsInSitu'                  , floattype),
        ('StarsInstability'             , floattype),
        ('StarsMergeBurst'              , floattype),
        ('DiscHI'                       , (floattype, 30)),
        ('DiscH2'                       , (floattype, 30)),
        ('DiscSFR'                      , (floattype, 30)),
        ('MetalsColdGas'                , floattype),
        ('MetalsStellarMass'            , floattype),
        ('ClassicalMetalsBulgeMass'     , floattype),
        ('SecularMetalsBulgeMass'       , floattype),
        ('MetalsHotGas'                 , floattype),
        ('MetalsEjectedMass'            , floattype),
        ('MetalsIntraClusterStars'      , floattype),
        ('DiscGasMetals'                , (floattype, 30)),
        ('DiscStarsMetals'              , (floattype, 30)),
        ('SfrDisk'                      , floattype),
        ('SfrBulge'                     , floattype),
        ('SfrDiskZ'                     , floattype),
        ('SfrBulgeZ'                    , floattype),
        ('DiskScaleRadius'              , floattype),
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


def read_darksage(fname):
    Galdesc = galdtype()
    with open(fname, 'rb') as fin:
        Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
        NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
        GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
        G = np.fromfile(fin, Galdesc, NtotGals) # Read all the galaxy data
    return G




###============= Main code =============###

# Make directory for input and output of this test
dir = 'test/'
if not os.path.exists(dir):
    os.makedirs(dir)

# Fetch the input for the test
if not os.path.isfile(dir+'model_to_test_against_z2.239_0'):
    zip = urlretrieve('https://github.com/arhstevens/DarkSageTest/archive/master.zip')
    subprocess.call(['unzip', '-j', zip[0], '-d', dir])
    subprocess.call(['rm', zip[0]])

# Delete any old data produced by this script
if os.path.isfile(dir+'model_z2.239_0'):
    subprocess.call(['rm', dir+'model_z2.239_0'])

# Run Dark Sage
subprocess.call(['./darksage', dir+'test.par'])


# Read produced and fetched Dark Sage output
G_test = read_darksage(dir+'model_to_test_against_z2.239_0')
G_out = read_darksage(dir+'model_z2.239_0')

# Switch off unhelpful warnings
import warnings
warnings.filterwarnings("ignore")

# Build a histogram of stellar masses to gauge if differences are physical or numerical
plt.figure()
h = 0.73
mmin, mmax = 8.0, 12.0
plt.hist(np.log10(G_test['StellarMass']*1e10/h), bins=np.arange(mmin,mmax,0.2), histtype='step', lw=2, color='k', label='Expected', log=True)
plt.hist(np.log10(G_out['StellarMass']*1e10/h), bins=np.arange(mmin,mmax,0.2), histtype='step', lw=2, color='b', ls='dashed', label='Result', log=True)
plt.xlabel('log Stellar Mass [solar]')
plt.ylabel('Number of galaxies')
plt.legend(loc='best', frameon=False)
figname = 'SMF_test.png'
plt.savefig(dir+figname, bbox_inches='tight')

# Compare output from installed Dark Sage to fetched data
success = True
for field in G_out.dtype.names:
    if not bool(np.allclose(G_out[field], G_test[field])):
        success = False
        diff_abs = G_out[field] - G_test[field]
        diff_rel = diff_abs / G_test[field]
        diff_rel = diff_rel[np.isfinite(diff_rel)] # Get rid of divide-by-0 entries
        print('Field {0} differed by max(abs,rel)=({1},{2}), mean(abs,rel)=({3},{4}), min(abs,rel)=({5},{6})'.format(field,np.max(diff_abs),np.max(diff_rel),np.mean(diff_abs),np.mean(diff_rel),np.min(diff_abs),np.min(diff_rel)))

# Declare success or not
if success:
    print('Success! Dark Sage output matches what is expected!')
else:
    print('Uh oh! The Dark Sage output did not match what was expected!')
    print('This can happen if {0}test.py or any of the Dark Sage codebase was modified from the main repository.'.format(dir))
    print('If you recently updated your local repository for Dark Sage, try deleting the `{0}\' directory and running this again.'.format(dir))
    print('See {0} to check if the difference is significant.'.format(dir+figname))
    sys.exit(1)