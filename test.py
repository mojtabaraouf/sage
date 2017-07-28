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

import subprocess
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    xrange
except NameError: # Python 3
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
        ('MergerBulgeMass'              , floattype),
        ('InstabilityBulgeMass'         , floattype),
        ('HotGas'                       , floattype),
        ('EjectedMass'                  , floattype),
        ('BlackHoleMass'                , floattype),
        ('IntraClusterStars'            , floattype),
        ('DiscGas'                      , (floattype, 30)),
        ('DiscStars'                    , (floattype, 30)),
        ('SpinStars'                    , (floattype, 3)),
        ('SpinGas'                      , (floattype, 3)),
        ('SpinMergerBulge'              , (floattype, 3)),
        ('StarsInSitu'                  , floattype),
        ('StarsInstability'             , floattype),
        ('StarsMergeBurst'              , floattype),
        ('DiscHI'                       , (floattype, 30)),
        ('DiscH2'                       , (floattype, 30)),
        ('DiscSFR'                      , (floattype, 30)),
        ('MetalsColdGas'                , floattype),
        ('MetalsStellarMass'            , floattype),
        ('MetalsMergerBulgeMass'        , floattype),
        ('MetalsInstabilityBulgeMass'   , floattype),
        ('MetalsHotGas'                 , floattype),
        ('MetalsEjectedMass'            , floattype),
        ('MetalsIntraClusterStars'      , floattype),
        ('DiscGasMetals'                , (floattype, 30)),
        ('DiscStarsMetals'              , (floattype, 30)),
        ('SfrDisk'                      , floattype),
        ('SfrBulge'                     , floattype),
        ('SfrDiskZ'                     , floattype),
        ('SfrBulgeZ'                    , floattype),
        ('CoolingScaleRadius'           , floattype),
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
SMfigname = dir+'SMF_test.png'
figdir = dir+'differences/'

# Fetch the input for the test
if not os.path.isfile(dir+'model_to_test_against_z2.239_0'):
    zip = urlretrieve('https://github.com/arhstevens/DarkSageTest/archive/master.zip')
    subprocess.call(['unzip', '-j', zip[0], '-d', dir])
    subprocess.call(['rm', zip[0]])

# Delete any old data produced by this script
if os.path.isfile(dir+'model_z2.239_0'):
    subprocess.call(['rm', dir+'model_z2.239_0'])
if os.path.exists(figdir):
    subprocess.call(['rm', '-rf', figdir])
if os.path.isfile(SMfigname):
    subprocess.call(['rm', SMfigname])

# Run Dark Sage
subprocess.call(['./darksage', dir+'test.par'])

# Check the produced and fetched data are the same size
size_test = os.path.getsize(dir+'model_to_test_against_z2.239_0')
size_out = os.path.getsize(dir+'model_z2.239_0')
if size_test!=size_out:
    print('\nUh oh! The Dark Sage output did not match what was expected!')
    print('The size of your output file is different.')
    print('Please report this issue if you cannot find a fast solution.')
    sys.exit(1)

# Read produced and fetched Dark Sage output
G_test = read_darksage(dir+'model_to_test_against_z2.239_0')
G_out = read_darksage(dir+'model_z2.239_0')

# Check that the same number of galaxies is present
if len(G_test)!=len(G_out):
    print('\nUh oh! The Dark Sage output did not match what was expected!')
    print('The number of galaxies is different to what was expected!')
    sys.exit(1)

# Check fields that are not modified by Dark Sage, and hence should be identical across machines and compilers
halo_fields = ['Type', 'GalaxyIndex', 'HaloIndex', 'SimulationHaloIndex',
               'TreeIndex', 'SnapNum', 'CentralGalaxyIndex', 'CentralMvir',
               'dT', 'Pos', 'Vel', 'Spin',
               'Len', 'LenMax', 'Mvir', 'Rvir', 'Vvir', 'Vmax', 'VelDisp',
               'CoolingScaleRadius', 'infallMvir', 'infallVvir', 'infallVmax']
for field in halo_fields:
    if not bool(np.allclose(G_out[field], G_test[field])):
        print('\nUh oh! The Dark Sage output did not match what was expected!')
        print('The properties that don\'t match should not be affected by your compiler.')
        print('This error was sprung by the property {0}, but is likely not limited to it.'.format(field))
        print('Please report this issue if you cannot find a fast solution.')
        sys.exit(1)

# Travis test should pass if this point is reached

# Reduce galaxies to those that are reasonably well resolved
f = (G_test['LenMax']>=50)
G_test, G_out = G_test[f], G_out[f]

# Switch off unhelpful warnings
import warnings
warnings.filterwarnings("ignore")

# Build a histogram of stellar masses as a sanity check
fig = plt.figure()
plt.clf()
fig.subplots_adjust(left=0, bottom=0)
plt.subplot(111)
h = 0.73
mmin, mmax = 8.0, 12.0
plt.hist(np.log10(G_test['StellarMass']*1e10/h), bins=np.arange(mmin,mmax,0.2), histtype='step', lw=2, color='k', label='Expected', log=True)
plt.hist(np.log10(G_out['StellarMass']*1e10/h), bins=np.arange(mmin,mmax,0.2), histtype='step', lw=2, color='b', ls='dashed', label='Result', log=True)
plt.xlabel('log Stellar Mass [solar]')
plt.ylabel('Number of galaxies')
plt.axis([mmin, mmax, 1, 1e3])
plt.legend(loc='best', frameon=False)
plt.savefig(SMfigname, bbox_inches='tight')

# Compare galaxy properties from installed Dark Sage to fetched data
great_success = True
fields_to_check = np.array(G_out.dtype.names) # returns all fields
for field in halo_fields+['mergeType', 'mergeIntoID', 'mergeIntoSnapNum']: # reduce fields to relevant galaxy evolution ones
    fields_to_check = np.delete(fields_to_check, np.where(fields_to_check==field)[0])
for field in fields_to_check:
    cut = (G_test[field]>np.percentile(G_test[field],2)) * np.isfinite(G_test[field]) * (G_test[field]>0) # cut out the lowest values, as these are the most likely to cause scientifically inconsequential problems
    field_test, field_out = G_test[field][cut], G_out[field][cut]
    diff_abs = abs(field_out - field_test)
    diff_rel = diff_abs / field_test
    frac_bad = (1.0*len(diff_rel[diff_rel>=0.02])) / (1.0*len(diff_rel)) # fraction with 2% of greater difference

    # If there are too many (>2%) galaxies with differences, plot where those differences are
    if frac_bad>=0.02:
        great_success = False
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        figname = figdir+field+'.png'
        plt.clf()
        plt.scatter(field_test[diff_rel<0.01], diff_rel[diff_rel<0.01], c='k')
        plt.scatter(field_test[diff_rel>=0.01], diff_rel[diff_rel>=0.01], c='r')
        plt.xlabel(field+' -- test output [internal units]')
        plt.ylabel('Fractional difference to your output ('+str(round(100*frac_bad,2))+'% red)')
        if(np.log10(np.max(field_test))>np.log10(np.min(field_test))+1): plt.xscale('log')
        plt.xlim(np.min(field_test), np.max(field_test))
        plt.ylim(-0.1,2)
        plt.savefig(figname, bbox_inches='tight')
        #print('Some differences in {0} of galaxies were found -- see {1}'.format(field,figname))

# Declare concerns
if great_success:
    print('\nSuccess! Dark Sage output matches what is expected!\n')
else:
    print('\nDark Sage galaxy properties differed more than they ideally should.')
    print('This might be due to your C compiler being different to that used for the test case.')
    print('This will also happen if {0}test.py or any of the Dark Sage codebase was modified from the main repository.'.format(dir))
    print('If you recently pulled updates for Dark Sage, try `rm -rf {0}\' then run this again.'.format(dir))
    print('See {0} and plots in {1} to check if the differences are significant.\n'.format(SMfigname, figdir))