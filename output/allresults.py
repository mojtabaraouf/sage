#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import h5py as h5
import numpy as np
import pylab as plt
from random import sample, seed
from os.path import getsize as getFileSize
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
# ================================================================================
# Basic variables
# ================================================================================

# Set up some basic attributes of the run

whichsimulation = 0
whichimf = 1        # 0=Slapeter; 1=Chabrier
dilute = 7500       # Number of galaxies to plot in scatter plots
sSFRcut = -11.0     # Divide quiescent from star forming galaxies (when plotmags=0)


matplotlib.rcdefaults()
plt.rc('axes', color_cycle=[
    'k',
    'b',
    'r',
    'g',
    'm',
    '0.5',
    ], labelsize='x-large')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

OutputDir = '' # set in main below

#Mac book Air
#OutputDir = '/Volumes/09158118409/sage-master/output/new_SAM/allresult/model/'

#IMAC
#OutputDir = '/Users/mraoufhajarzarrin/Desktop/sage-master/output/new_SAM/allresult/model/'



OutputFormat = '.png'
TRANSPARENT = False

OutputList = []

class Results:

    """ The following methods of this class generate the figures and plot them.
    """

    def __init__(self):
        """Here we set up some of the variables which will be global to this
        class."""

        self.Hubble_h = 0.73

        if whichsimulation == 0:    # Mini-Millennium
          self.BoxSize = 62.5       # Mpc/h
          self.MaxTreeFiles = 8     # FilesPerSnapshot

        elif whichsimulation == 1:  # Millennium
          self.BoxSize = 500        # Mpc/h
          self.MaxTreeFiles = 512   # FilesPerSnapshot

        elif whichsimulation == 2:  # Bolshoi
          self.BoxSize = 250.0      # Mpc/h
          self.MaxTreeFiles = 12987 # FilesPerSnapshot

        elif whichsimulation == 3:  # GiggleZ MR
          self.BoxSize = 125.0      # Mpc/h
          self.MaxTreeFiles = 8 # FilesPerSnapshot

        else:
          print "Please pick a valid simulation!"
          exit(1)


    def read_gals(self, model_name, first_file, last_file):

        # The input galaxy structure:
        Galdesc_full = [
            ('SnapNum'                      , np.int32),                    
            ('Type'                         , np.int32),                    
            ('GalaxyIndex'                  , np.int64),                    
            ('CentralGalaxyIndex'           , np.int64),                    
            ('SAGEHaloIndex'                , np.int32),                    
            ('SAGETreeIndex'                , np.int32),                    
            ('SimulationFOFHaloIndex'       , np.int32),                    
            ('mergeType'                    , np.int32),                    
            ('mergeIntoID'                  , np.int32),                    
            ('mergeIntoSnapNum'             , np.int32),                    
            ('dT'                           , np.float32),                    
            ('Pos'                          , (np.float32, 3)),             
            ('Vel'                          , (np.float32, 3)),             
            ('Spin'                         , (np.float32, 3)),             
            ('Len'                          , np.int32),                    
            ('Mvir'                         , np.float32),                  
            ('CentralMvir'                  , np.float32),                  
            ('Rvir'                         , np.float32),                  
            ('Vvir'                         , np.float32),                  
            ('Vmax'                         , np.float32),                  
            ('VelDisp'                      , np.float32),                  
            ('ColdGas'                      , np.float32),                  
            ('StellarMass'                  , np.float32),                  
            ('BulgeMass'                    , np.float32),                  
            ('HotGas'                       , np.float32),                  
            ('EjectedMass'                  , np.float32),                  
            ('BlackHoleMass'                , np.float32),                  
            ('IntraClusterStars'            , np.float32),                  
            ('MetalsColdGas'                , np.float32),                  
            ('MetalsStellarMass'            , np.float32),                  
            ('MetalsBulgeMass'              , np.float32),                  
            ('MetalsHotGas'                 , np.float32),                  
            ('MetalsEjectedMass'            , np.float32),                  
            ('MetalsIntraClusterStars'      , np.float32),                  
            ('SfrDisk'                      , np.float32),                  
            ('SfrBulge'                     , np.float32),                  
            ('SfrDiskZ'                     , np.float32),                  
            ('SfrBulgeZ'                    , np.float32),                  
            ('DiskRadius'                   , np.float32),                  
            ('Cooling'                      , np.float32),                  
            ('Heating'                      , np.float32),
            ('r_heat'                       , np.float32),
            ('QuasarModeBHaccretionMass'    , np.float32),
            ('TimeSinceMajorMerger'         , np.float32),
            ('TimeSinceMinorMerger'         , np.float32),
            ('OutflowRate'                  , np.float32),
            ('infallMvir'                   , np.float32),
            ('infallVvir'                   , np.float32),
            ('infallVmax'                   , np.float32),
            ('Qjet'                         , np.float32),
            ('Rcocoon'                      , np.float32),
            ('Rshocked'                     , np.float32),
            ('t_AGN_returne'                , np.float32),
            ('t_AGN_on'                     , np.float32),
            ('Tshocked'                     , np.float32),
            ('Mshocked'                     , np.float32),
            ('RadioLuminosity'              , (np.float32, 7)),
            ('RadioAGNaccretionRate'        , np.float32),
            ('rho_zero_Makino'              , np.float32),
            ('rho_zero_Capelo'              , np.float32),
            ('rho_zero_iso'                 , np.float32),
            ('b_gas'                        , np.float32),
            ('Rs'                           , np.float32),
            ('concentration'                , np.float32),
            ('Temp_Gas'                     , np.float32),
            ('Lx_bol'                       , np.float32),
            ('R_index'                      , np.float32),
            ('Q_index'                      , np.float32),
            ('R_cool'                       , np.float32),
            ('fcool'                        , np.float32),
            ('t_static'                     , np.float32),
            ('t_AGN_off'                    , np.float32),
            ('time_to_next_on'              , np.float32),
            ('delta'                        , np.float32)
            ]
        names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
        formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
        Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)


        # Initialize variables.
        TotNTrees = 0
        TotNGals = 0
        FileIndexRanges = []

        print "Determining array storage requirements."
        
        # Read each file and determine the total number of galaxies to be read in
        goodfiles = 0
        for fnr in xrange(first_file,last_file+1):
            fname = model_name+'_'+str(fnr)  # Complete filename
        
            if not os.path.isfile(fname):
              # print "File\t%s  \tdoes not exist!  Skipping..." % (fname)
              continue

            if getFileSize(fname) == 0:
                print "File\t%s  \tis empty!  Skipping..." % (fname)
                continue
        
            fin = open(fname, 'rb')  # Open the file
            Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
            NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
            TotNTrees = TotNTrees + Ntrees  # Update total sim trees number
            TotNGals = TotNGals + NtotGals  # Update total sim gals number
            goodfiles = goodfiles + 1  # Update number of files read for volume calculation
            fin.close()

        print
        print "Input files contain:\t%d trees ;\t%d galaxies ." % (TotNTrees, TotNGals)
        print

        # Initialize the storage array
        G = np.empty(TotNGals, dtype=Galdesc)

        offset = 0  # Offset index for storage array

        # Open each file in turn and read in the preamble variables and structure.
        print "Reading in files."
        for fnr in xrange(first_file,last_file+1):
            fname = model_name+'_'+str(fnr)  # Complete filename
        
            if not os.path.isfile(fname):
              continue

            if getFileSize(fname) == 0:
                continue
        
            fin = open(fname, 'rb')  # Open the file
            Ntrees = np.fromfile(fin, np.dtype(np.int32), 1)  # Read number of trees in file
            NtotGals = np.fromfile(fin, np.dtype(np.int32), 1)[0]  # Read number of gals in file.
            GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
            print ":   Reading N=", NtotGals, "   \tgalaxies from file: ", fname
            GG = np.fromfile(fin, Galdesc, NtotGals)  # Read in the galaxy structures
        
            FileIndexRanges.append((offset,offset+NtotGals))
        
            # Slice the file array into the global array
            # N.B. the copy() part is required otherwise we simply point to
            # the GG data which changes from file to file
            # NOTE THE WAY PYTHON WORKS WITH THESE INDICES!
            G[offset:offset+NtotGals]=GG[0:NtotGals].copy()
            
            del(GG)
            offset = offset + NtotGals  # Update the offset position for the global array
        
            fin.close()  # Close the file


        print
        print "Total galaxies considered:", TotNGals

        # Convert the Galaxy array into a recarray
        G = G.view(np.recarray)

        w = np.where(G.StellarMass > 1.0)[0]
        print "Galaxies more massive than 10^10Msun/h:", len(w)

        print

        # Calculate the volume given the first_file and last_file
        self.volume = self.BoxSize**3.0 * goodfiles / self.MaxTreeFiles
        print self.volume,self.BoxSize,goodfiles,self.MaxTreeFiles
        return G

# --------------------------------------------------------

    def StellarMassFunction(self, G):

        print 'Plotting the stellar mass function'

        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        binwidth = 0.1  # mass function histogram bin width

        # calculate all
        w = np.where(G.StellarMass > 0.0)[0]
        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        sSFR = (G.SfrDisk[w] + G.SfrBulge[w]) / (G.StellarMass[w] * 1.0e10 / self.Hubble_h)

        mi = np.floor(min(mass)) - 2
        ma = np.floor(max(mass)) + 2
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        
        # additionally calculate red
        w = np.where(sSFR < 10.0**sSFRcut)[0]
        massRED = mass[w]
        (countsRED, binedges) = np.histogram(massRED, range=(mi, ma), bins=NB)

        # additionally calculate blue
        w = np.where(sSFR > 10.0**sSFRcut)[0]
        massBLU = mass[w]
        (countsBLU, binedges) = np.histogram(massBLU, range=(mi, ma), bins=NB)

        # Baldry+ 2008 modified data used for the MCMC fitting
        Baldry = np.array([
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
            [11.95, 7.4764e-06, 7.4764e-06],
            ], dtype=np.float32)
            
        No_jet_mMS = np.array( [
             [  8.05000000e+00,   1.09600000e+03],
             [  8.15000000e+00,   1.25100000e+03],   
             [  8.25000000e+00,   1.39700000e+03],   
             [  8.35000000e+00,   1.66400000e+03],   
             [  8.45000000e+00,   1.75400000e+03],   
             [  8.55000000e+00,   1.70800000e+03],   
             [  8.65000000e+00,   1.64200000e+03],   
             [  8.75000000e+00,   1.49600000e+03],   
             [  8.85000000e+00,   1.18800000e+03],   
             [  8.95000000e+00,   1.06400000e+03],   
             [  9.05000000e+00,   9.80000000e+02],   
             [  9.15000000e+00,   8.50000000e+02],   
             [  9.25000000e+00,   7.98000000e+02],   
             [  9.35000000e+00,   7.07000000e+02],   
             [  9.45000000e+00,   6.67000000e+02],   
             [  9.55000000e+00,   6.26000000e+02],   
             [  9.65000000e+00,   6.31000000e+02],   
             [  9.75000000e+00,   5.92000000e+02],   
             [  9.85000000e+00,   5.21000000e+02],   
             [  9.95000000e+00,   5.72000000e+02],   
             [  1.00500000e+01,   5.01000000e+02],   
             [  1.01500000e+01,   4.63000000e+02],   
             [  1.02500000e+01,   4.30000000e+02],   
             [  1.03500000e+01,   4.27000000e+02],   
             [  1.04500000e+01,   4.19000000e+02],   
             [  1.05500000e+01,   3.11000000e+02],   
             [  1.06500000e+01,   3.19000000e+02],   
             [  1.07500000e+01,   2.69000000e+02],   
             [  1.08500000e+01,   2.39000000e+02],   
             [  1.09500000e+01,   2.00000000e+02],   
             [  1.10500000e+01,   1.94000000e+02],   
             [  1.11500000e+01,   1.69000000e+02],   
             [  1.12500000e+01,   1.39000000e+02],   
             [  1.13500000e+01,   1.36000000e+02],   
             [  1.14500000e+01,   1.28000000e+02],   
             [  1.15500000e+01,   6.20000000e+01],   
             [  1.16500000e+01,   4.80000000e+01],   
             [  1.17500000e+01,   2.60000000e+01],   
             [  1.18500000e+01,   1.30000000e+01],   
             [  1.19500000e+01,   5.00000000e+00],   
             [  1.20500000e+01,   4.00000000e+00],   
             [  1.21500000e+01,   0.00000000e+00],   
             [  1.22500000e+01,   0.00000000e+00],   
             [  1.23500000e+01,   0.00000000e+00],   
             ], dtype=np.float32)
             
        No_jet_Iso_mMS = np.array( [
             [  8.02848837e+00,   1.65100000e+03],
             [  8.17965116e+00,   1.90000000e+03],
             [  8.33081395e+00,   2.39800000e+03],
             [  8.48197674e+00,   2.63400000e+03],
             [  8.63313953e+00,   2.54800000e+03],
             [  8.78430233e+00,   2.19300000e+03],
             [  8.93546512e+00,   1.73000000e+03],
             [  9.08662791e+00,   1.44300000e+03],
             [  9.23779070e+00,   1.24700000e+03],
             [  9.38895349e+00,   1.05400000e+03],
             [  9.54011628e+00,   9.72000000e+02],
             [  9.69127907e+00,   8.97000000e+02],
             [  9.84244186e+00,   8.15000000e+02],
             [  9.99360465e+00,   8.35000000e+02],
             [  1.01447674e+01,   6.89000000e+02],
             [  1.02959302e+01,   6.58000000e+02],
             [  1.04470930e+01,   6.13000000e+02],
             [  1.05982558e+01,   4.68000000e+02],
             [  1.07494186e+01,   4.06000000e+02],
             [  1.09005814e+01,   3.37000000e+02],
             [  1.10517442e+01,   2.94000000e+02],
             [  1.12029070e+01,   2.33000000e+02],
             [  1.13540698e+01,   2.22000000e+02],
             [  1.15052326e+01,   1.16000000e+02],
             [  1.16563953e+01,   6.20000000e+01],
             [  1.18075581e+01,   2.70000000e+01],
             [  1.19587209e+01,   1.50000000e+01],
             [  1.21098837e+01,   5.00000000e+00],
             [  1.22610465e+01,   2.00000000e+00],
             [  1.24122093e+01,   0.00000000e+00],
             [  1.25633721e+01,   0.00000000e+00],
             [  1.27145349e+01,   0.00000000e+00],
             [  1.28656977e+01,   0.00000000e+00],
             ], dtype=np.float32)

        heat_mMS = np.array( [
             [  8.05000000e+00,   1.09600000e+03],
             [  8.15000000e+00,   1.25100000e+03],
             [  8.25000000e+00,   1.40000000e+03],
             [  8.35000000e+00,   1.66300000e+03],
             [  8.45000000e+00,   1.75800000e+03],
             [  8.55000000e+00,   1.70500000e+03],
             [  8.65000000e+00,   1.64200000e+03],
             [  8.75000000e+00,   1.49600000e+03],
             [  8.85000000e+00,   1.19000000e+03],
             [  8.95000000e+00,   1.06000000e+03],
             [  9.05000000e+00,   9.82000000e+02],
             [  9.15000000e+00,   8.50000000e+02],
             [  9.25000000e+00,   7.99000000e+02],
             [  9.35000000e+00,   7.05000000e+02],
             [  9.45000000e+00,   6.68000000e+02],
             [  9.55000000e+00,   6.26000000e+02],
             [  9.65000000e+00,   6.33000000e+02],
             [  9.75000000e+00,   5.96000000e+02],
             [  9.85000000e+00,   5.24000000e+02],
             [  9.95000000e+00,   5.74000000e+02],
             [  1.00500000e+01,   5.08000000e+02],
             [  1.01500000e+01,   4.61000000e+02],
             [  1.02500000e+01,   4.57000000e+02],
             [  1.03500000e+01,   4.54000000e+02],
             [  1.04500000e+01,   4.48000000e+02],
             [  1.05500000e+01,   3.64000000e+02],
             [  1.06500000e+01,   3.32000000e+02],
             [  1.07500000e+01,   3.44000000e+02],
             [  1.08500000e+01,   2.81000000e+02],
             [  1.09500000e+01,   2.49000000e+02],
             [  1.10500000e+01,   1.85000000e+02],
             [  1.11500000e+01,   1.50000000e+02],
             [  1.12500000e+01,   1.14000000e+02],
             [  1.13500000e+01,   7.10000000e+01],
             [  1.14500000e+01,   4.40000000e+01],
             [  1.15500000e+01,   3.00000000e+01],
             [  1.16500000e+01,   9.00000000e+00],
             [  1.17500000e+01,   5.00000000e+00],
             [  1.18500000e+01,   1.00000000e+00],
             [  1.19500000e+01,   0.00000000e+00],
               ], dtype=np.float32)
        Heat_MS = np.array( [
             [  8.05000000e+00,   1.73850000e+04],
             [  8.15000000e+00,   2.01670000e+04],
             [  8.25000000e+00,   2.35420000e+04],
             [  8.35000000e+00,   2.61030000e+04],
             [  8.45000000e+00,   2.79610000e+04],
             [  8.55000000e+00,   2.78380000e+04],
             [  8.65000000e+00,   2.63030000e+04],
             [  8.75000000e+00,   2.33800000e+04],
             [  8.85000000e+00,   2.03780000e+04],
             [  8.95000000e+00,   1.72640000e+04],
             [  9.05000000e+00,   1.53250000e+04],
             [  9.15000000e+00,   1.39060000e+04],
             [  9.25000000e+00,   1.25580000e+04],
             [  9.35000000e+00,   1.18230000e+04],
             [  9.45000000e+00,   1.12180000e+04],
             [  9.55000000e+00,   1.06720000e+04],
             [  9.65000000e+00,   1.01820000e+04],
             [  9.75000000e+00,   9.52300000e+03],
             [  9.85000000e+00,   9.33000000e+03],
             [  9.95000000e+00,   8.72800000e+03],
             [  1.00500000e+01,   8.36300000e+03],
             [  1.01500000e+01,   7.90000000e+03],
             [  1.02500000e+01,   7.49200000e+03],
             [  1.03500000e+01,   7.13200000e+03],
             [  1.04500000e+01,   6.69200000e+03],
             [  1.05500000e+01,   6.37200000e+03],
             [  1.06500000e+01,   5.98000000e+03],
             [  1.07500000e+01,   5.67800000e+03],
             [  1.08500000e+01,   4.96300000e+03],
             [  1.09500000e+01,   3.66100000e+03],
             [  1.10500000e+01,   2.94100000e+03],
             [  1.11500000e+01,   2.41700000e+03],
             [  1.12500000e+01,   1.68000000e+03],
             [  1.13500000e+01,   9.64000000e+02],
             [  1.14500000e+01,   6.61000000e+02],
             [  1.15500000e+01,   3.87000000e+02],
             [  1.16500000e+01,   1.58000000e+02],
             [  1.17500000e+01,   6.00000000e+01],
             [  1.18500000e+01,   2.20000000e+01],
             [  1.19500000e+01,   7.00000000e+00],
             [  1.20500000e+01,   4.00000000e+00],
             [  1.21500000e+01,   3.00000000e+00],
             [  1.22500000e+01,   0.00000000e+00],
                ], dtype=np.float32)
        No_jet_MS = np.array( [
             [  8.05000000e+00,   1.73850000e+04],
             [  8.15000000e+00,   2.01590000e+04],
             [  8.25000000e+00,   2.35370000e+04],
             [  8.35000000e+00,   2.60930000e+04],
             [  8.45000000e+00,   2.79360000e+04],
             [  8.55000000e+00,   2.78350000e+04],
             [  8.65000000e+00,   2.62780000e+04],
             [  8.75000000e+00,   2.33530000e+04],
             [  8.85000000e+00,   2.03960000e+04],
             [  8.95000000e+00,   1.72780000e+04],
             [  9.05000000e+00,   1.53440000e+04],
             [  9.15000000e+00,   1.39090000e+04],
             [  9.25000000e+00,   1.25570000e+04],
             [  9.35000000e+00,   1.18130000e+04],
             [  9.45000000e+00,   1.11980000e+04],
             [  9.55000000e+00,   1.06730000e+04],
             [  9.65000000e+00,   1.01610000e+04],
             [  9.75000000e+00,   9.50900000e+03],
             [  9.85000000e+00,   9.27800000e+03],
             [  9.95000000e+00,   8.64600000e+03],
             [  1.00500000e+01,   8.29700000e+03],
             [  1.01500000e+01,   7.72100000e+03],
             [  1.02500000e+01,   7.29700000e+03],
             [  1.03500000e+01,   6.77400000e+03],
             [  1.04500000e+01,   6.14900000e+03],
             [  1.05500000e+01,   5.57900000e+03],
             [  1.06500000e+01,   5.24800000e+03],
             [  1.07500000e+01,   4.56000000e+03],
             [  1.08500000e+01,   4.09700000e+03],
             [  1.09500000e+01,   3.61900000e+03],
             [  1.10500000e+01,   3.08300000e+03],
             [  1.11500000e+01,   2.51300000e+03],
             [  1.12500000e+01,   2.25000000e+03],
             [  1.13500000e+01,   2.17900000e+03],
             [  1.14500000e+01,   1.84400000e+03],
             [  1.15500000e+01,   9.75000000e+02],
             [  1.16500000e+01,   5.78000000e+02],
             [  1.17500000e+01,   3.50000000e+02],
             [  1.18500000e+01,   2.14000000e+02],
             [  1.19500000e+01,   1.26000000e+02],
             [  1.20500000e+01,   6.10000000e+01],
             [  1.21500000e+01,   1.20000000e+01],
             [  1.22500000e+01,   6.00000000e+00],
             [  1.23500000e+01,   2.00000000e+00],
             [  1.24500000e+01,   2.00000000e+00],
             [  1.25500000e+01,   0.00000000e+00],
                ], dtype=np.float32)
        No_jet_MS_Iso = np.array( [
             [  8.05000000e+00,   1.73320000e+04],
             [  8.15000000e+00,   1.98940000e+04],
             [  8.25000000e+00,   2.33320000e+04],
             [  8.35000000e+00,   2.57660000e+04],
             [  8.45000000e+00,   2.78220000e+04],
             [  8.55000000e+00,   2.79800000e+04],
             [  8.65000000e+00,   2.69330000e+04],
             [  8.75000000e+00,   2.42900000e+04],
             [  8.85000000e+00,   2.15860000e+04],
             [  8.95000000e+00,   1.81000000e+04],
             [  9.05000000e+00,   1.59400000e+04],
             [  9.15000000e+00,   1.42940000e+04],
             [  9.25000000e+00,   1.28100000e+04],
             [  9.35000000e+00,   1.19140000e+04],
             [  9.45000000e+00,   1.12530000e+04],
             [  9.55000000e+00,   1.07270000e+04],
             [  9.65000000e+00,   1.01720000e+04],
             [  9.75000000e+00,   9.49200000e+03],
             [  9.85000000e+00,   9.29200000e+03],
             [  9.95000000e+00,   8.63500000e+03],
             [  1.00500000e+01,   8.29900000e+03],
             [  1.01500000e+01,   7.70100000e+03],
             [  1.02500000e+01,   7.27000000e+03],
             [  1.03500000e+01,   6.77400000e+03],
             [  1.04500000e+01,   6.09300000e+03],
             [  1.05500000e+01,   5.57200000e+03],
             [  1.06500000e+01,   5.21100000e+03],
             [  1.07500000e+01,   4.53700000e+03],
             [  1.08500000e+01,   4.09400000e+03],
             [  1.09500000e+01,   3.64500000e+03],
             [  1.10500000e+01,   3.13300000e+03],
             [  1.11500000e+01,   2.58900000e+03],
             [  1.12500000e+01,   2.37100000e+03],
             [  1.13500000e+01,   2.25100000e+03],
             [  1.14500000e+01,   1.78500000e+03],
             [  1.15500000e+01,   8.31000000e+02],
             [  1.16500000e+01,   5.15000000e+02],
             [  1.17500000e+01,   3.28000000e+02],
             [  1.18500000e+01,   2.14000000e+02],
             [  1.19500000e+01,   1.19000000e+02],
             [  1.20500000e+01,   8.70000000e+01],
             [  1.21500000e+01,   7.00000000e+01],
             [  1.22500000e+01,   3.90000000e+01],
             [  1.23500000e+01,   1.60000000e+01],
             [  1.24500000e+01,   1.40000000e+01],
             [  1.25500000e+01,   5.00000000e+00],
             [  1.26500000e+01,   4.00000000e+00],
             [  1.27500000e+01,   2.00000000e+00],
             [  1.28500000e+01,   2.00000000e+00],
             [  1.29500000e+01,   0.00000000e+00],
             [  1.30500000e+01,   0.00000000e+00],
             ], dtype=np.float32)

              
        # Finally plot the data
        # plt.errorbar(
        #     Baldry[:, 0],
        #     Baldry[:, 1],
        #     yerr=Baldry[:, 2],
        #     color='g',
        #     linestyle=':',
        #     lw = 1.5,
        #     label='Baldry et al. 2008',
        #     )

        Baldry_xval = np.log10(10 ** Baldry[:, 0]  /self.Hubble_h/self.Hubble_h)
        if(whichimf == 1):  Baldry_xval = Baldry_xval - 0.26  # convert back to Chabrier IMF
        Baldry_yvalU = (Baldry[:, 1]+Baldry[:, 2]) * self.Hubble_h*self.Hubble_h*self.Hubble_h
        Baldry_yvalL = (Baldry[:, 1]-Baldry[:, 2]) * self.Hubble_h*self.Hubble_h*self.Hubble_h

        plt.fill_between(Baldry_xval, Baldry_yvalU, Baldry_yvalL, 
            facecolor='red', alpha=0.35, label='Baldry et al. 2008 (z=0.1)')
            
            
            
        No_jet_Iso_mMS_xval=No_jet_Iso_mMS[:,0]
        No_jet_Iso_mMS_yval=No_jet_Iso_mMS[:,1]
        No_jet_xval = No_jet_mMS[:,0]
        No_jet_yval = No_jet_mMS[:,1]

        heat_mMS_xval = heat_mMS[:,0]
        heat_mMS_yval = heat_mMS[:,1]
        
        No_jet_MS_xval = No_jet_MS[:,0]
        No_jet_MS_yval = No_jet_MS[:,1]
        No_jet_MS_Iso_xval = No_jet_MS_Iso[:,0]
        No_jet_MS_Iso_yval = No_jet_MS_Iso[:,1]
        heat_MS_xval = Heat_MS[:,0]
        heat_MS_yval = Heat_MS[:,1]


        # This next line is just to get the shaded region to appear correctly in the legend

#        plt.plot(xaxeshisto, counts / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, label='Baldry et al. 2008', color='red')
        x = [1,1,1]
        y = [2,2,2]
        plt.plot(x,y, label='Baldry et al. 2008', color='red',lw= 4, alpha= 0.35)
#        print np.column_stack((xaxeshisto, counts.astype(int)))
        # # Cole et al. 2001 SMF (h=1.0 converted to h=0.73)
        # M = np.arange(7.0, 13.0, 0.01)
        # Mstar = np.log10(7.07*1.0e10 /self.Hubble_h/self.Hubble_h)
        # alpha = -1.18
        # phistar = 0.009 *self.Hubble_h*self.Hubble_h*self.Hubble_h
        # xval = 10.0 ** (M-Mstar)
        # yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)      
        # plt.plot(M, yval, 'g--', lw=1.5, label='Cole et al. 2001')  # Plot the SMF
        
        # Overplot the model histograms
        
############################################################
#################mili-Millennium plots######################
############################################################
#        plt.plot( No_jet_Iso_mMS_xval, No_jet_Iso_mMS_yval   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / 0.15, color ='grey', lw = 5, alpha=0.45, label='No-AGN (Isothermal)')
#        plt.plot(No_jet_xval, No_jet_yval   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / 0.1, 'b--', lw = 3, alpha=0.85, label='No-AGN (Makino+98)')
#        plt.plot( heat_mMS_xval, heat_mMS_yval   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / 0.1, 'm--', lw = 4, alpha=0.85, label='Jet-model (+Heat)')
#
#        plt.plot(xaxeshisto, counts    / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'b-', lw = 3, alpha=0.45, label='Jet-model (+Heat, +Uplift)')

############################################################
#################full-millennium plots######################
############################################################

        plt.plot(No_jet_MS_Iso_xval, No_jet_MS_Iso_yval   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / 0.1, 'k:', label='W/O-AGN (Iso)')

        plt.plot(No_jet_MS_xval, No_jet_MS_yval   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / 0.1, 'b--', label='W/O-AGN (Makino+98)')
        plt.plot( heat_MS_xval, heat_MS_yval   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / 0.1, 'm--', label='Jet-model (+Heat)')#

        plt.plot(xaxeshisto, counts    / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'b-', label='Jet-model (+Heat, +Uplift)')

#
#        print  np.column_stack((xaxeshisto, counts))

#        plt.plot(xaxeshisto, countsRED / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'r:', lw=2, label='Model - Red')
#        plt.plot(xaxeshisto, countsBLU / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'b:', lw=2, label='Model - Blue')

        plt.yscale('log', nonposy='clip')
        plt.axis([8.4, 12.7, 2.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} (M_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

#        plt.text(12.2, 0.03, whichsimulation, size = 'large')

        leg = plt.legend(loc='lower left', numpoints=1,
                         labelspacing=0.3)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')

        outputFile = OutputDir + '1_StellarMassFunction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
   
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def Metallicity(self, G):
    
        print 'Plotting the metallicities'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        w = np.where((G.Type == 0) & (G.ColdGas / (G.StellarMass + G.ColdGas) > 0.1) & (G.StellarMass > 0.01))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        mass = np.log10(G.StellarMass[w] * 1e10 / self.Hubble_h)
        Z = np.log10((G.MetalsColdGas[w] / G.ColdGas[w]) / 0.02) + 9.0
        
        plt.scatter(mass, Z, marker='o', s=10, color = 'grey', alpha=0.05)
        total_bins = 30
        X = mass
        Y = Z
        bins = np.linspace(X.min(),X.max(), total_bins)
#        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2, running_median, running_std, color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)
#        print  np.column_stack((bins-delta/2, running_median, running_std))
        # overplot Tremonti et al. 2003 (h=0.7)
#        mass , p2.5, p16, p50, p84, p97.5
        Tremonti2004 = np.array([
            [8.57  ,     8.18  ,     8.25 ,     8.44  ,     8.64  ,     8.77  ],
            [8.67  ,     8.11  ,     8.28 ,     8.48  ,     8.65  ,     8.84  ],
            [8.76  ,     8.13  ,     8.32 ,     8.57  ,     8.70  ,     8.88  ],
            [8.86  ,     8.14  ,     8.37 ,     8.61  ,     8.73  ,     8.89  ],
            [8.96  ,     8.21  ,     8.46 ,     8.63  ,     8.75  ,     8.95  ],
            [9.06  ,     8.26  ,     8.56 ,     8.66  ,     8.82  ,     8.97  ],
            [9.16  ,     8.37  ,     8.59 ,     8.68  ,     8.82  ,     8.95  ],
            [9.26  ,     8.39  ,     8.60 ,     8.71  ,     8.86  ,     9.04  ],
            [9.36  ,     8.46  ,     8.63 ,     8.74  ,     8.88  ,     9.03  ],
            [9.46  ,     8.53  ,     8.66 ,     8.78  ,     8.92  ,     9.07  ],
            [9.57  ,     8.59  ,     8.69 ,     8.82  ,     8.94  ,     9.08  ],
            [9.66  ,     8.60  ,     8.72 ,     8.84  ,     8.96  ,     9.09  ],
            [9.76  ,     8.63  ,     8.76 ,     8.87  ,     8.99  ,     9.10  ],
            [9.86  ,     8.67  ,     8.80 ,     8.90  ,     9.01  ,     9.12  ],
            [9.96  ,     8.71  ,     8.83 ,     8.94  ,     9.05  ,     9.14  ],
            [10.06  ,     8.74  ,     8.85 ,     8.97  ,     9.06  ,     9.15  ],
            [10.16  ,     8.77  ,     8.88 ,     8.99  ,     9.09  ,     9.16  ],
            [10.26  ,     8.80  ,     8.92 ,     9.01  ,     9.10  ,     9.17  ],
            [10.36  ,     8.82  ,     8.94 ,     9.03  ,     9.11  ,     9.18  ],
            [10.46  ,     8.85  ,     8.96 ,     9.05  ,     9.12  ,     9.21  ],
            [10.56  ,     8.87  ,     8.98 ,     9.07  ,     9.14  ,     9.21  ],
            [10.66  ,     8.89  ,     9.00 ,     9.08  ,     9.15  ,     9.23  ],
            [10.76  ,     8.91  ,     9.01 ,     9.09  ,     9.15  ,     9.24  ],
            [10.86  ,     8.93  ,     9.02 ,     9.10  ,     9.16  ,     9.25  ],
            [10.95  ,     8.93  ,     9.03 ,     9.11  ,     9.17  ,     9.26  ],
            [11.05  ,     8.92  ,     9.03 ,     9.11  ,     9.17  ,     9.27  ],
            [11.15  ,     8.94  ,     9.04 ,     9.12  ,     9.18  ,     9.29  ],
            [11.25  ,     8.93  ,     9.03 ,     9.12  ,     9.18  ,     9.29  ],
            ], dtype=np.float32)
        
        Tremonti2004_xval = np.log10(10 ** Tremonti2004[:, 0]  /self.Hubble_h)
        if(whichimf == 1):  Tremonti2004_xval = Tremonti2004_xval - 0.26  # convert back to Chabrier IMF
        Tremonti2004_yval = Tremonti2004[:, 3]
        Tremonti2004_yvalU = Tremonti2004[:, 5]
        Tremonti2004_yvalL = Tremonti2004[:, 1]
        plt.fill_between(Tremonti2004_xval, Tremonti2004_yvalU, Tremonti2004_yvalL,
                         facecolor='red', alpha=0.15, label='Tremonti et al. 2004')
        plt.scatter(Tremonti2004_xval, Tremonti2004_yval, color='red', s= 50, lw=2.0, alpha=0.6, marker='*', label='Tremonti et al. (2004)')

        C06_radiomode = np.array([
            [  8.1696498  ,  8.40445518         , 0.0],
            [  8.23550807  , 8.41744995         , 0.0],
            [  8.30136635   ,8.42498684         , 0.0],
            [  8.36722463  , 8.43815804         , 0.0],
            [  8.4330829   , 8.45020866         , 0.0],
            [  8.49894118,   8.48432159         , 0.0],
            [  8.56479945,   8.49309158         , 0.0],
            [  8.63065773,   8.50572205         , 0.0],
            [  8.69651601,   8.51731014         , 0.0],
            [  8.76237428,   8.54400826         , 0.0],
            [  8.82823256,   8.53828716         , 0.0],
            [  8.89409084,   8.57867527  , 0.32811964],
            [  8.95994911,   8.58334637  , 0.39820516],
            [  9.02580739,   8.59026527  , 0.33167547],
            [  9.09166567,   8.63316917  , 0.33505383],
            [  9.15752394,   8.59598255  , 0.31745344],
            [  9.22338222,   8.64529228  , 0.19593245],
            [  9.2892405,    8.66448784  , 0.24415597],
            [  9.35509877 ,  8.68904495  , 0.19811964],
            [  9.42095705,   8.66328144  , 0.21577145],
            [  9.48681533,   8.72275257  , 0.14804168],
            [  9.5526736 ,   8.74805069  , 0.20016955],
            [  9.61853188,   8.77426147  , 0.17289273],
            [  9.68439016,   8.76579189  , 0.17022833],
            [  9.75024843,   8.79179001  , 0.19856688],
            [  9.81610671,   8.8088789   , 0.16956642],
            [  9.88196499,   8.85399151  , 0.16543479],
            [  9.94782326,   8.86523056  , 0.15174147],
            [  10.01368154,   8.86518764  , 0.14398499],
            [  10.07953981 ,  8.89381886  , 0.15836425],
            [  10.14539809 ,  8.91288757  , 0.1228577 ],
            [  10.21125637 ,  8.93867111  , 0.14206758],
            [  10.27711464 ,  8.92370796  , 0.18610322],
            [  10.34297292 ,  8.95251179  , 0.12795138],
            [  10.4088312  ,  8.98543072  , 0.12930936],
            [  10.47468947 ,  9.0313797   , 0.19201504],
            [  10.54054775 ,  9.03530788  , 0.13325645],
            [  10.60640603 ,  9.01967525  , 0.09557171],
            [  10.6722643  ,  9.0589819   , 0.13003948],
            [  10.73812258 ,  9.06770897  , 0.10807342],
            [  10.80398086 ,  9.10868645  , 0.11382841],
            [  10.86983913 ,  9.10411835  , 0.12160621],
            [  10.93569741 ,  9.13024998  , 0.05859008],
            [  11.00155569 ,  9.15722942  , 0.06320351],
            [  11.06741396 ,  9.18078995  , 0.04941788],
            [  11.13327224 ,  9.08304596  , 0.10199235],
            [  11.19913052 ,  9.1738081   , 0.06388161],
            [  11.33084707 ,  8.98505878   , 0.04602838],
            ], dtype=np.float32)
        WO_AGN = np.array([
             [  8.19350724,   8.42535496 ,        0.0],
             [  8.30715862,   8.46208572  ,       0.0],
             [  8.42080999,   8.48798752  ,       0.0],
             [  8.53446137,   8.51355362  ,       0.0],
             [  8.64811274,   8.53928375  ,       0.0],
             [  8.76176412,   8.55736923  ,       0.0],
             [  8.87541549,   8.58065796   ,      0.0],
             [  8.98906686,   8.60317039   ,      0.0],
             [  9.10271824,   8.62510967   ,      0.0],
             [  9.21636961,   8.64912796   ,      0.0],
             [  9.33002099,   8.67246246   ,      0.0],
             [  9.44367236,   8.70321083   ,      0.0],
             [  9.55732374,   8.73611259   ,      0.0],
             [  9.67097511,   8.77640629   ,      0.0],
             [  9.78462648,   8.81467915    ,     0.0],
             [  9.89827786,   8.85155678   ,0.1271788 ],
             [ 10.01192923,   8.88812637     ,    0.0],
             [ 10.12558061,   8.92440796   ,0.13707912],
             [ 10.23923198,   8.95521259   ,0.12232325],
             [ 10.35288336,   8.98885345   ,0.13806941],
             [ 10.46653473,   9.01888561   ,0.13000216],
             [ 10.5801861 ,   9.04424286   ,      0.0],
             [ 10.69383748,   9.06837082   ,0.14887953],
             [ 10.80748885,   9.09390259   ,0.14143373],
             [ 10.92114023,   9.11397934   ,0.14185718],
             [ 11.0347916 ,   9.11975288   ,0.15042779],
             [ 11.14844298,   9.09342766   ,0.1593359 ],
             [ 11.26209435,   9.0403862    ,0.14703996],
             [ 11.37574572,   9.02881432   ,0.08471279],
            ], dtype=np.float32)
        plt.plot( WO_AGN[:,0], WO_AGN[:,1],'k-', label='W/O AGN Jet')
#        plt.fill_between(WO_AGN[:,0], WO_AGN[:,1]+WO_AGN[:,2],WO_AGN[:,1]-WO_AGN[:,2],
#                             facecolor='blue', alpha=0.15)

#        plt.errorbar(C06_radiomode[:,0], C06_radiomode[:,1], C06_radiomode[:,2], color='k', lw=2.0, alpha=0.6, marker='o', markersize=3, ls='none', label='C16-Radio mode(Median$~\pm~\sigma$)', mew=1)
        plt.plot(C06_radiomode[:,0],C06_radiomode[:,1],'m-', label='C16(Radio mode)')
#        plt.fill_between(C06_radiomode[:,0], C06_radiomode[:,1]+C06_radiomode[:,2],C06_radiomode[:,1]-C06_radiomode[:,2],
#                         facecolor='red', alpha=0.15)



#        w = np.arange(7.0, 13.0, 0.1)
#        Zobs = -1.492 + 1.847*w - 0.08026*w*w
#        if(whichimf == 0):
#            # Conversion from Kroupa IMF to Slapeter IMF
#            plt.plot(np.log10((10**w *1.5)), Zobs, 'b-', lw=2.0, label='Tremonti et al. 2003')
#        elif(whichimf == 1):
#            # Conversion from Kroupa IMF to Slapeter IMF to Chabrier IMF
#            plt.plot(np.log10((10**w *1.5 /1.8)), Zobs, 'b-', lw=2.0, label='Tremonti et al. 2003')

        plt.ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
        plt.xlabel(r'$\log_{10} (M_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([8.5, 12.0, 8.0, 9.5])
            
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '2_Metallicity' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)
    

# ---------------------------------------------------------

    def BlackHoleBulgeRelationship(self, G):
    
        print 'Plotting the black hole-bulge relationship'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
    
        w = np.where((G.BulgeMass > 0.01) & (G.BlackHoleMass > 0.00001))[0]
        if(len(w) > dilute): w = sample(w, dilute)
    
        bh = np.log10(G.BlackHoleMass[w] * 1.0e10 / self.Hubble_h)
        bulge = np.log10(G.BulgeMass[w] * 1.0e10 / self.Hubble_h)
                    
        plt.scatter(bulge, bh, marker='o', s=10, color = 'grey', alpha=0.05)
        
        total_bins = 30
        X = bulge
        Y = bh
        bins = np.linspace(X.min(),X.max(), total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)
#        print  np.column_stack((bins-delta/2 , running_median, running_std))

#        BHmass, errorU, errorL, bulgemass
        Haring_Rix2004 = np.array([
                                 [3.0e9, 1.0e9, 1.0e9, 6e11  ],
                                 [1.4e7, 1.3e7, 0.7e7, 2.3e10],
                                 [1.0e8, 0.6e8, 0.5e8, 6.8e10],
                                 [4.3e8, 3.2e8, 1.7e8, 3.6e11],
                                 [5.2e8, 1.0e8, 1.1e8, 3.6e11],
                                 [5.3e8, 2.0e8, 4.0e8, 5.6e11],
                                 [3.3e8, 2.3e8, 1.3e8, 2.9e11],
                                 [1.4e7, 0.4e7, 0.5e7, 6.2e9 ],
                                 [3.7e7, 1.7e7, 1.5e7, 1.3e11],
                                 [2.5e9, 0.5e9, 0.4e9, 2.9e11],
#                                 group2
                                 [4.5e7, 4.0e7, 2.5e7, 3.7e10],
                                 [2.5e6, 0.5e6, 0.5e6, 8.0e8 ],
                                 [4.4e7, 0.5e7, 0.5e7, 6.9e10],
                                 [1.4e7, 1.6e7, 0.8e7, 7.6e10],
                                 [1.0e9, 1.0e9, 0.6e9, 1.2e11],
                                 [2.1e8, 1.0e8, 0.6e8, 6.8e10],
                                 [1.0e8, 0.9e8, 0.1e8, 1.6e10],
                                 [1.6e7, 0.1e7, 0.2e7, 2.0e10],
                                 [1.9e8, 1.0e8, 0.6e8, 9.7e10],
                                 [3.1e8, 1.3e8, 1.1e8, 1.3e11],
                                 [3.0e8, 1.7e8, 1.0e8, 1.2e10],
                                 [1.1e8, 0.4e8, 1.0e8, 9.2e10],
                                 [5.6e7, 0.3e7, 0.7e7, 4.4e10],
                                 [1.0e9, 1.0e9, 0.7e9, 2.7e11],
                                 [2.0e9, 0.5e9, 1.0e9, 4.9e11],
                                 [1.7e8, 0.2e8, 0.1e8, 1.1e11],
                                 [2.4e8, 0.4e8, 1.4e8, 3.7e10],
                                 [1.3e7, 0.6e7, 0.5e7, 1.5e10],
                                 [3.5e6, 1.1e6, 1.4e6, 7.0e9 ],
                                 [3.7e6, 1.5e6, 1.5e6, 1.1e10],
                                 ], dtype=np.float32)
                                 
        Haring_Rix2004_xval = np.log10(Haring_Rix2004[:, 3])
        Haring_Rix2004_yval = np.log10(Haring_Rix2004[:, 0])
        Haring_Rix2004_yvalU = np.log10(Haring_Rix2004[:, 0] + Haring_Rix2004[:, 1]) - Haring_Rix2004_yval
        Haring_Rix2004_yvalL = Haring_Rix2004_yval - np.log10(Haring_Rix2004[:, 0] - Haring_Rix2004[:, 2])
        plt.errorbar(Haring_Rix2004_xval, Haring_Rix2004_yval, yerr=[Haring_Rix2004_yvalL ,Haring_Rix2004_yvalU], color='red', lw=2.0, alpha=0.3, marker='*', markersize=12, ls='none', label='Haring \& Rix (2004)', mew=1)
        Scott2013 = np.array([
                                 [ 39,     4,     5 , 69, 59, 32 ],
                                 [ 11,     2,     2 , 37, 32, 17 ],
                                 [ 0.45,  0.17,  0.10 , 1.4, 2.0, 0.8 ],
                                 [ 25,     7,     7 , 55, 80, 33 ],
                                 [ 24,    10,    10 , 27, 23, 12 ],
                                 [ 0.044, 0.044, 0.022 , 2.4, 3.5, 1.4 ],
                                 [ 1.4,   0.9,   0.3 , 0.46, 0.68, 0.28 ],
                                 [ 0.73, 0.0, 0.0 , 1.0, 1.5, 0.6 ],
                                 [  9.0,   0.9,   0.8 , 19, 16, 9 ],
                                 [ 58,   3.5,   3.5 , 23,   19,   10 ],
                                 [ 0.10,  0.10,  0.05 , 0.61, 0.89, 0.36 ],
                                 [ 8.3,   2.7,   1.3 , 4.6, 6.6, 2.7 ],
                                 [ 0.39,  0.26,  0.09 , 11, 9, 5 ],
                                 [ 0.42,  0.04,  0.04 , 1.9, 2.7, 1.1 ],
                                 [ 0.084, 0.003, 0.003 , 4.5, 6.6, 2.7 ],
                                 [ 0.66,  0.03,  0.03 , 1.4, 2.1, 0.8 ],
                                 [ 0.73,  0.69,  0.35 , 0.66, 0.97, 0.40 ],
                                 [ 15,     2,     2 , 4.7, 6.9, 2.8 ],
                                 [ 4.7,   0.6,   0.6 , 26, 22, 12 ],
                                 [ 0.083, 0.004, 0.004 , 2.0, 2.9, 1.2 ],
                                 [ 0.14,  0.02,  0.13 , 0.39, 0.57, 0.23 ],
                                 [ 0.15,  0.09,   0.1 , 0.35, 0.52, 0.21 ],
                                 [ 0.40,  0.04,  0.05 , 0.30, 0.45, 0.18 ],
                                 [ 0.12, 0.005, 0.005 , 3.5, 5.1, 2.1 ],
                                 [  1.7,   0.2,   0.2 , 6.7, 5.7, 3.1 ],
                                 [ 0.024, 0.024, 0.012 , 0.88, 1.28, 0.52 ],
                                 [ 8.8,  10.0,   2.7 , 1.9, 2.7, 1.1 ],
                                 [ 0.14,  0.10,  0.06 , 0.93, 1.37, 0.56 ],
                                 [  2.0,   0.5,   0.5 , 1.24, 1.8, 0.7 ],
                                 [ 0.073, 0.015, 0.015 , 0.86, 1.26, 0.51 ],
                                 [ 0.77,  0.04,  0.06 , 2.0, 1.7, 0.9 ],
                                 [  4.0,   1.0,   1.0 , 5.4, 4.7, 2.5 ],
                                 [ 0.17,  0.01,  0.02 , 1.2, 1.7, 0.7 ],
                                 [ 0.34,  0.02,  0.02 , 4.9, 7.1, 2.9 ],
                                 [ 2.4,   0.3,   0.3 , 2.0, 2.9, 1.2 ],
                                 [ 0.058, 0.008, 0.008 , 0.66, 0.97, 0.40 ],
                                 [ 3.1,   1.4,   0.6 , 5.1, 7.4, 3.0 ],
                                 [  1.3,   0.5,   0.5 , 2.6, 3.8, 1.5 ],
                                 [  2.0,   1.1,   0.6 , 3.2, 2.7, 1.5 ],
                                 [ 97,    30,    26 , 100, 86, 46 ],
                                 [ 8.1,   2.0,   1.9 , 1.4, 2.1, 0.9 ],
                                 [  1.8,   0.6,   0.3 , 0.88, 1.30, 0.53 ],
                                 [ 0.65,  0.07,  0.07 , 1.3, 1.9, 0.8 ],
                                 [ 0.39,  0.01,  0.01 , 0.56, 0.82, 0.34 ],
                                 [  5.0,   1.0,   1.0 , 29, 25, 13 ],
                                 [ 3.3,   0.9,   2.5 , 6.1, 5.2, 2.8 ],
                                 [  4.5,   2.3,   1.5 , 0.65, 0.96, 0.39 ],
                                 [ 0.075, 0.002, 0.002 , 3.3, 4.9, 2.0 ],
                                 [ 0.68,  0.13,  0.13 , 2.0, 3.0, 1.2 ],
                                 [  1.2,   0.4,   0.9 , 6.9, 5.9, 3.2 ],
                                 [ 0.13,  0.08,  0.08 , 1.4, 1.2, 0.6 ],
                                 [ 4.7,   0.5,   0.5 , 7.7, 6.6, 3.6 ],
                                 [ 0.59,  0.03,  0.09 , 0.90, 1.3, 0.5 ],
                                 [ 6.4,   0.4,   0.4 , 3.9, 5.7, 2.3 ],
                                 [ 0.79,  0.38,  0.33 , 1.8, 2.7, 1.1 ],
                                 [ 3.9,   0.4,   0.4 , 8.4, 7.2, 3.9 ],
                                 [ 47,    10,    10 , 27, 23, 12 ],
                                 [  1.8,   0.2,   0.1 , 6.0, 5.2, 2.8 ],
                                 [ 0.060, 0.014, 0.014 , 0.43, 0.64, 0.26 ],
                                 [ 0.016, 0.004, 0.004 , 1.0, 1.5, 0.6 ],
                                 [ 210,   160,   160 , 122, 105, 57 ],
                                 [ 0.014, 0.014, 0.007 , 0.30, 0.45, 0.18 ],
                                 [ 7.4,   4.7,   3.0 , 29, 25, 13 ],
                                 [ 1.6,   0.3,   0.4 , 11, 10, 5 ],
                                 [ 6.8,   0.7,   0.7 , 20, 17, 9 ],
                                 [ 2.6,   0.4,   1.5 , 2.8, 2.4, 1.3 ],
                                 [ 11,     1,     1 , 24, 20, 11 ],
                                 [ 37,    18,    11 , 78, 67, 36 ],
                                 [ 5.9,   2.0,   2.0 , 96, 83, 44 ],
                                 [ 0.31, 0.004, 0.004 , 3.6, 5.2, 2.1 ],
                                 [ 0.10, 0.001, 0.001 , 2.6, 3.8, 1.5 ],
                                 [ 3.7,   2.6,   1.5 , 55, 48, 26 ],
                                 [ 0.55,  0.26,  0.19 , 1.4, 2.0, 0.8 ],
                                 [ 13,     5,     4 , 64, 55, 30 ],
                                 [ 0.11, 0.005, 0.005, 1.2, 1.8, 0.7 ],
                                 ], dtype=np.float32)
        Scott2013_xval = np.log10(Scott2013[:, 3] * 1e10)
        Scott2013_xvalU = np.log10(Scott2013[:, 3]* 1e10 + Scott2013[:, 4] * 1e10) - Scott2013_xval
        Scott2013_xvalL = Scott2013_xval - np.log10(Scott2013[:, 3]*1e10 - Scott2013[:, 5]*1e10)
                                 
        Scott2013_yval = np.log10(Scott2013[:, 0] * 1e8)
        Scott2013_yvalU = np.log10(Scott2013[:, 0]* 1e8 + Scott2013[:, 1] * 1e8) - Scott2013_yval
        Scott2013_yvalL = Scott2013_yval - np.log10(Scott2013[:, 0]*1e8 - Scott2013[:, 2]*1e8)
                                 
        plt.errorbar(Scott2013_xval, Scott2013_yval, xerr=[Scott2013_xvalL ,Scott2013_xvalU], yerr=[Scott2013_yvalL ,Scott2013_yvalU], color='g', lw=2.0, alpha=0.3, marker='o', markersize=12, ls='none', label='Scott et al. (2013)', mew=1)
            
            
        C16_radiomod = np.array([
             [  8.20096016  , 5.78055906  , 0.34268743],
             [  8.32373238  , 5.81182098  , 0.35518801],
             [  8.44650459  , 5.81877995  , 0.36261836],
             [  8.56927681  , 5.93450165  , 0.38017115],
             [  8.69204903  , 6.01289749  , 0.41153997],
             [  8.81482124  , 6.09399128  , 0.40203154],
             [  8.93759346  , 6.1391654   , 0.43416524],
             [  9.06036568  , 6.26994991  , 0.49449337],
             [  9.18313789  , 6.24147987  , 0.49529985],
             [  9.30591011  , 6.3880167   , 0.49670631],
             [  9.42868233  , 6.41732264  , 0.49717754],
             [  9.55145454  , 6.57046556  , 0.45355967],
             [  9.67422676  , 6.65723515  , 0.4715611 ],
             [  9.79699898  , 6.83294773  , 0.42552239],
             [  9.91977119   ,6.95776367  , 0.37278837],
             [ 10.04254341   ,7.12588692  , 0.32109106],
             [ 10.16531563   ,7.19048214  , 0.36343625],
             [ 10.28808784  , 7.4281044   , 0.31336313],
             [ 10.41086006  , 7.57510424  , 0.24604779],
             [ 10.53363228  , 7.64705086  , 0.23507479],
             [ 10.6564045   , 7.74794483  , 0.22545446],
             [ 10.77917671  , 7.91944838  , 0.21864554],
             [ 10.90194893  , 7.94540882  , 0.1898195 ],
             [ 11.02472115  , 7.98879862  , 0.22896363],
             [ 11.14749336  , 8.19942093  , 0.17171529],
             [ 11.27026558   ,8.05909729  , 0.17347971],
             [ 11.3930378    ,8.32024097  , 0.01881277],
             [ 11.51581001   ,8.60910416  , 0.15775472],
              ], dtype=np.float32)
            
            
        WO_AGN = np.array([
             [  8.20437264 ,  5.66519308,   0.29150358],
             [  8.33903763 ,  5.73355293,   0.30672476],
             [  8.47370263 ,  5.82126808,   0.30891222],
             [  8.60836762  , 5.88087845,   0.31300241],
             [  8.74303262 ,  5.94557285,   0.34254357],
             [  8.87769762  , 6.02951717,   0.35799038],
             [  9.01236261  , 6.03450775,   0.41259292],
             [  9.14702761  , 6.16672707,   0.40676951],
             [  9.2816926   , 6.24311447,   0.41857803],
             [  9.4163576   , 6.37715626,   0.41059718],
             [  9.5510226   , 6.48192883,   0.42418519],
             [  9.68568759  , 6.60839558,   0.40825772],
             [  9.82035259  , 6.74312592,   0.36236078],
             [  9.95501758  , 6.89757776,   0.33002424],
             [ 10.08968258  , 7.11632252,   0.33210322],
             [ 10.22434757 ,  7.30670309,   0.3038148 ],
             [ 10.35901257 ,  7.46141768,   0.27826899],
             [ 10.49367757 ,  7.61683369,   0.19997555],
             [ 10.62834256 ,  7.79985046,   0.2084181 ],
             [ 10.76300756  , 7.91753674,   0.19088672],
             [ 10.89767255  , 8.05400467,   0.17636988],
             [ 11.03233755 ,  8.18591022,   0.15330727],
             [ 11.16700255 ,  8.29750443,   0.14088827],
             [ 11.30166754  , 8.35466194,   0.16252218],
             [ 11.43633254  , 8.43300915,   0.1569801 ],
             [ 11.57099753  , 8.46488857,   0.11263688],
             [ 11.70566253  , 8.64614296,   0.1248166 ],
             [ 11.84032753  , 8.71203327,   0.12789141],
             [ 11.97499252  , 8.73795414,   0.17633615],
             ], dtype=np.float32)
        plt.plot( WO_AGN[:,0], WO_AGN[:,1],'b-', alpha=0.15,lw=4, label='W/O AGN Jet')
        plt.fill_between(WO_AGN[:,0], WO_AGN[:,1]+WO_AGN[:,2],WO_AGN[:,1]-WO_AGN[:,2],
                         facecolor='blue', alpha=0.15)
             

        plt.plot(C16_radiomod[:,0],C16_radiomod[:,1],'m-',lw=4,alpha=.15, label='C06(Radio mode)')
#        plt.errorbar(C16_radiomod[:,0],C16_radiomod[:,1], C16_radiomod[:,2],color='m', lw=2.0, alpha=0.6, marker='o', markersize=8, ls='none', label='C06(Radio mode)', mew=1)
        plt.fill_between(C16_radiomod[:,0], C16_radiomod[:,1]+C16_radiomod[:,2],C16_radiomod[:,1]-C16_radiomod[:,2],
                         facecolor='red', alpha=0.15)
        # overplot Haring & Rix 2004
#        w = 10. ** np.arange(20)
#        BHdata = 10. ** (8.2 + 1.12 * np.log10(w / 1.0e11))
#        plt.plot(np.log10(w), np.log10(BHdata), 'b-', label="Haring \& Rix 2004")
             
        plt.ylabel(r'$\log\ M_{\mathrm{BH}}\ (M_{\odot})$')  # Set the y...
        plt.xlabel(r'$\log\ M_{\mathrm{bulge}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([8.0, 12.0, 6.0, 10.0])
            
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '3_BlackHoleBulgeRelationship' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def Lradio_Qjet(self, G):
    
        print 'Plotting the Radio Luminosity -- Qjet relation'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        w = np.where((np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)
        w1 = np.where((G.Type == 0) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0)&(G.RadioLuminosity[:,5]>0)&(G.RadioLuminosity[:,5]<1e40))[0]
#        if(len(w1) > dilute): w1 = sample(w1, dilute)

        w2 = np.where((G.Type == 1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
        if(len(w2) > dilute): w2 = sample(w2, dilute)


        Lradio1400_0 = np.log10(G.RadioLuminosity[w,0])
        Lradio1400_1 = np.log10(G.RadioLuminosity[w,1])
        Lradio1400_2 = np.log10(G.RadioLuminosity[w,2])
        Lradio1400_3 = np.log10(G.RadioLuminosity[w,3])
        Lradio1400_4 = np.log10(G.RadioLuminosity[w,4])
        Lradio1400_5 = np.log10(G.RadioLuminosity[w,5])
        Lradio1400_6 = np.log10(G.RadioLuminosity[w,6])
        Q_jet_1      = np.log10(G.Qjet[w])


        Lradio1400_5_S = np.log10(G.RadioLuminosity[w2,5])
        Q_jet_S      = np.log10(G.Qjet[w2])
        Lradio1400_5_E = np.log10(G.RadioLuminosity[w1,5])
        Q_jet_E      = np.log10(G.Qjet[w1])
        

        Heckman2014_1 = np.array([
                                 [20.01386138613861, 33.269085411942555],
                                 [21.275247524752473, 34.354497354497354],
                                 [22.728382838283828, 35.600151171579746],
                                 [24.768316831683165, 37.3567649281935],
                                 [26.547194719471946, 38.88662131519274],
                                  ], dtype=np.float32)

        Heckman2014_2 = np.array([
                                    [20.01735733221337,  34.033570935094446],
                                    [22.775669620493755,  35.910935356386275],
                                    [24.769305876740944,  37.25960257094492],
                                    [26.916561848624763,  38.71132532014027],
                                     ], dtype=np.float32)


        # For fw = 10,20 convert to f= 5
        Heckman2014_xval_1 = Heckman2014_1[:, 0]
        Heckman2014_yval_1 = Heckman2014_1[:, 1]
        Heckman2014_xval_2 = Heckman2014_2[:, 0]
        Heckman2014_yval_2 = Heckman2014_2[:, 1]
                                
 
        plt.scatter(Lradio1400_5_E, Q_jet_E, marker='o',s=10, color='grey', alpha=0.05)
#        plt.scatter(Lradio1400_5_S, Q_jet_S, marker='o', s=10, color='b', alpha=0.15,label='Satellite Galaxies (Jet-model)')

        total_bins = 30
        X = Lradio1400_5_E
        Y = Q_jet_E
        bins = np.linspace(X.min(),X.max(), total_bins)
#        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)
        
        plt.plot(Heckman2014_xval_2, Heckman2014_yval_2, 'k-', label=(r'Heckman \& Best(2014), (Eq.2)'))
        plt.plot(Heckman2014_xval_1, Heckman2014_yval_1, 'k--', label=(r'Heckman \& Best(2014), (Eq.1)'))


        plt.ylabel(r'$\log_{10} (Q_{jet} ~ [W])$')  # Set the y...
        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
                                
        plt.axis([21., 26.0, 34, 38.0])
        #        plt.axis([14.0, 25.0, 25, 37.0])
        leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
                                
        outputFile = OutputDir + '4_Lradio_Qjet' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
    
        # Add this plot to our output list
        OutputList.append(outputFile)
    
    
# ---------------------------------------------------------
    def Rshocked_Rvir(self, G):
        
        print 'Plotting the R shocked -- R virial relation'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        w = np.where((G.Type == 0) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0)&(G.StellarMass>0) )[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        R_shocked = G.Rshocked[w]* 1000.0
        R_vir      = G.Rvir[w]* 1000.0
        
        w1 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>12.5))[0]
        if(len(w1) > dilute): w1 = sample(w1, dilute)

        R_vir_Hi      = G.Rvir[w1]* 1000.0
        R_shocked_Hi      = G.Rshocked[w1]* 1000.0
        w2 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) <12.5))[0]
        if(len(w2) > dilute): w2 = sample(w2, dilute)

        R_vir_Lo      = G.Rvir[w2]* 1000.0
        R_shocked_Lo      = G.Rshocked[w2]* 1000.0
        
        plt.scatter(R_vir, R_shocked, marker='o', s=10, color = 'grey', alpha=0.05)
#        plt.scatter(R_vir_Hi, R_shocked_Hi, marker='o', s=30, c='r', alpha=0.25, label='Satellite-Galaxies-Hi')
#        plt.scatter(R_vir_Lo, R_shocked_Lo, marker='o', s=30, c='b', alpha=0.45, label='Satellite-Galaxies-Lo')


        total_bins = 50
        X = R_vir
        Y = R_shocked
        bins = np.linspace(X.min(),X.max(), total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)


        m_0 = 1.0
        b_0 = 1.0
        x = [10,100,300,500, 1000]
        y = [10,100,300,500, 1000]
        plt.plot(x, y, 'k--', lw = 2,label='1:1')

        plt.ylabel(r'$R_{\mathrm{Shock}}~ [kpc]$')  # Set the y...
        plt.xlabel(r'$R_{\mathrm{vir}}~ [kpc]$')  # and the x-axis labels

        plt.axis([100.0, 700, 0.0, 700.0])
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '5_Rshocked_Rvir' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)

# ---------------------------------------------------------

    def Lradio_Rshock(self, G):
    
        print 'Plotting the Radio Luminosity -- Rshock relation'
        
        seed(2222)


        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        w = np.where((G.Type == 0) & (G.Rshocked >0.0001)& (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)
        mass = np.log10(G.StellarMass[w] * 1e10 / self.Hubble_h)
        Lradio1400 = np.log10(G.RadioLuminosity[w,5])
        R_shocked      = np.log10(G.Rshocked[w]* 1000.0/self.Hubble_h)
        
        w1 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>13)& (G.Rshocked >0.0001))[0]
#        if(len(w1) > dilute): w1 = sample(w1, dilute)

        Lradio1400_Hi = np.log10(G.RadioLuminosity[w1,5])
        R_shocked_Hi      = np.log10(G.Rshocked[w1]* 1000.0/self.Hubble_h)
        w2 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)<13)& (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>11)& (G.Rshocked >0.001))[0]
#        if(len(w2) > dilute): w2 = sample(w2, dilute)

        Lradio1400_Lo = np.log10(G.RadioLuminosity[w2,5])
        R_shocked_Lo      = np.log10(G.Rshocked[w2]* 1000.0/self.Hubble_h)
        
        total_bins = 50
        X = R_shocked
        Y = Lradio1400
        bins = np.linspace(X.min(),X.max(), total_bins)
        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=10, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)

        # Shabala et al. (2008)
        Shabala2008_LS = np.array([
                                   [51.05240355303475, 1.3186490373585515E22],
                                   [99.30812142284744, 3.114864039219891E23],
                                   [139.18268120338283, 1.2706630032922422E24],
                                   [181.09948998145737, 3.2958443369251625E24],
                                   [276.07295608175383, 1.4548869997215123E25],
                                   [394.42704247506873, 4.971831440193331E25],
                                   [502.5800391819623, 1.0973531662830234E26],
                                   ], dtype=np.float32)
        Shabala2008_LS_2 = np.array([
                                    [136.89181681379173, 1.1744199801610529E22],
                                    [235.3088586574712, 9.544703947366687E22],
                                    [390.11800635348067, 5.773782468179141E23],
                                    [936.5523719538558, 1.073259900813307E25],
                                    ], dtype=np.float32)

                                    
        Shabala2008_LS_xval = np.log10(Shabala2008_LS[:,0]/2.0)
        Shabala2008_LS_yval = np.log10(Shabala2008_LS[:,1])
        
        Shabala2008_LS_2_xval = np.log10(Shabala2008_LS_2[:,0]/2.0)
        Shabala2008_LS_2_yval = np.log10(Shabala2008_LS_2[:,1])


        plt.scatter(R_shocked, Lradio1400, marker='o', s=10, color = 'grey', alpha=0.05)
#        plt.scatter(R_shocked_Hi, Lradio1400_Hi, marker='o', s=30, c='r', alpha=0.25, label='Satellite-Galaxies-Hi')
#        plt.scatter(R_shocked_Lo, Lradio1400_Lo, marker='o', s=30, c='b', alpha=0.45, label='Satellite-Galaxies-Lo')
        plt.plot(Shabala2008_LS_xval, Shabala2008_LS_yval, 'k--', label='Shabala et al. (2008)~(Low-$t_{on}$)')
        plt.plot(Shabala2008_LS_2_xval, Shabala2008_LS_2_yval, 'k-', label='Shabala et al. (2008)~(High-$t_{on}$)')

        plt.ylabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # Set the y...
        plt.xlabel(r'$\log_{10} (R_{\mathrm{Shock}}~ [kpc])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        #        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

        plt.axis([1, 3, 22, 26.0])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
    
        outputFile = OutputDir + '6_Lradio_Rshock' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def Lradio_Mass(self, G):
    
        print 'Plotting the Radio Luminosity -- Mass'
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        w = np.where((G.Type == 0)&(np.log10(G.CentralMvir* 1e10/self.Hubble_h) > 11) &(G.StellarMass>0.0) & (G.RadioLuminosity[:,5]>0)& (G.RadioLuminosity[:,5]< 1e40))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        Lradio1400 = np.log10(G.RadioLuminosity[w,5])
        mass      = np.log10(G.StellarMass[w] * 1e10/self.Hubble_h)
        plt.scatter( Lradio1400, mass, marker='o', s=10, color = 'grey', alpha=0.15)

        total_bins = 30
        X = Lradio1400
        Y = mass
        bins = np.linspace(X.min(),X.max(), total_bins)
#        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)


#        High extinsion
        Best2012_HE = np.array([
        [23.216549052858767,  10.954998965825974],
        [23.719985739887846,  11.19047171353235],
        [24.21784325438884,  11.166462965721834],
        [24.715522710192126,  11.224002817581828],
        [25.21641531613127,  11.217296596771932],
        [25.714046210471547,  11.297076975814793],
        [26.214830902048444,  11.339794148744602],
        ], dtype=np.float32)
#        Low extinsion
        Best2012_LE = np.array([
        [23.277412753165166,  11.22691441522883],
        [23.77804096891681,  11.341295509081206],
        [24.27562869751219,  11.440845245619949],
        [24.776337849035514,  11.51815879416755],
        [25.277003834813943,  11.615241700211032],
        [25.774969263677182,  11.541809558660812],
        [26.272902318231743,  11.483204435232503],
        ], dtype=np.float32)
            
        Best2012_HE_xval = Best2012_HE[:,0]
        Best2012_HE_yval = Best2012_HE[:,1]
        if(whichimf == 1):  Best2012_HE_yval = Best2012_HE_yval - 0.26
        
        Best2012_LE_xval = Best2012_LE[:,0]

        Best2012_LE_yval = Best2012_LE[:,1]
        if(whichimf == 1):  Best2012_LE_yval = Best2012_LE_yval - 0.26
        
        plt.errorbar(Best2012_HE_xval, Best2012_HE_yval, yerr=0.0, color='red', alpha=0.6, marker='o', markersize=17, ls='none', label=r'Best \& Heckman (2012)-HERGs', mew=1, lw=2.0)
        
        plt.errorbar(Best2012_LE_xval, Best2012_LE_yval, yerr=0.0, color='orange', alpha=0.6, marker='d', markersize=17, ls='none', label=r'Best \& Heckman (2012)-LERGs', mew=1, lw=2.0)



        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # Set the y...
        plt.ylabel(r'$\log_{10} (M_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        plt.axis([ 23.0, 26.5,10.5, 12.0])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
    
        outputFile = OutputDir + '7_Lradio_Mass' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)

# ---------------------------------------------------------

    def Temp_hist(self, G):
    
        print 'Plotting the Temp_hist relation'
        
        seed(2222)
#        plt.figure(figsize=(16,6))  # New figure
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        w = np.where((G.Type == 0) & (np.log10(G.CentralMvir*1e10/self.Hubble_h) > 11)&(G.StellarMass>0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        Temp = np.log10(G.Temp_Gas[w])
        Tvir = np.log10(35.9 * G.Vvir[w] * G.Vvir[w])
        
        plt.scatter(mass,Temp/Tvir, marker='o', s=10, color = 'grey', alpha=0.05)

        total_bins = 70
        X = mass
        Y = Temp/Tvir
        bins = np.linspace(X.min(),X.max(), total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)

        x=[9,11,12,13]
        y=[1,1,1,1]
        plt.plot(x, y, 'k--', lw=2)


        plt.ylabel(r'$T_{\mathrm{new-hot}}/T_{\mathrm{vir}} $')  # Set the y...
        plt.xlabel(r'$\log_{10} (M_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        plt.axis([10, 12.0, 0.96, 1.15])
        
        
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')




        outputFile = OutputDir + '8_Temp_hist' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)
# ---------------------------------------------------------

    def RadioLF(self, G):
    
        print 'Plotting the Radio Luminosity function'
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        binwidth = 0.5  # Radio Luminosity function histogram bin width
        w1 = np.where((G.RadioLuminosity[:,5] > 0)& (G.RadioLuminosity[:,5] < 1e30)&(np.log10(G.CentralMvir * 1e10/self.Hubble_h) > 11) & (G.fcool>0))[0]
#        if(len(w1) > dilute): w1 = sample(w1, dilute)

        Lradio1400_5 = np.log10(G.RadioLuminosity[w1,5])
        delta_duty  =  0.03 * (G.StellarMass[w1]* 1e10 / 1e11)**1.5
        mi_5 = np.floor(min(Lradio1400_5)) - 3.0
        ma_5 = np.floor(max(Lradio1400_5)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        
        
        (counts_5, binedges_5) = np.histogram(Lradio1400_5, range=(mi_5, ma_5), bins=NB_5 )
        
        # Set the x-axis values to be the centre of the bins
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 12)
        
        Lradio1400_5_12 = Lradio1400_5[w]
        (counts_5_12, binedges_5) = np.histogram(Lradio1400_5_12, range=(mi_5, ma_5), bins=NB_5 )
        
        
#             Lradio, logp_all, Uy_all, Ly_all, logp_AGN, Uy_AGN, Ly_AGN,
        Best2012 = np.array([
                            [  22.0, -3.09, 0.03, 0.03, -4.44, 0.15, 0.23],
                            [  22.3, -3.49, 0.02, 0.03, -4.45, 0.05, 0.05],
                            [  22.6, -3.87, 0.02, 0.02, -4.45, 0.04, 0.04],
                            [  22.9, -4.22, 0.02, 0.02, -4.48, 0.02, 0.03],
                            [  23.2, -4.56, 0.02, 0.02, -4.69, 0.02, 0.02],
                            [  23.5, -4.75, 0.01, 0.01, -4.80, 0.01, 0.01],
                            [  23.8, -4.90, 0.02, 0.02, -4.91, 0.02, 0.02],
                            [  24.1, -5.08, 0.01, 0.01, -5.09, 0.01, 0.02],
                            [  24.4, -5.25, 0.02, 0.02, -5.26, 0.02, 0.02],
                            [  24.7, -5.54, 0.02, 0.02, -5.54, 0.02, 0.02],
                            [  25.0, -5.82, 0.03, 0.03, -5.82, 0.03, 0.03],
                            [  25.3, -6.32, 0.05, 0.06, -6.32, 0.05, 0.06],
                            [  25.6, -6.58, 0.07, 0.08, -6.58, 0.07, 0.08],
                            [  25.9, -7.18, 0.12, 0.17, -7.18, 0.12, 0.17],
                            [  26.2, -7.78, 0.21, 0.43, -7.78, 0.21, 0.43],
                             ], dtype=np.float32)
                             
                             
#        xplot= np.log10(10**(Best2012[:,0]+0.15))
#        yplot = (10**Best2012[:,1] )
#        yerr2 = 10**(Best2012[:,1]+Best2012[:,2]) -  yplot
#        yerr1 = yplot - 10**(Best2012[:,1]-Best2012[:,3])
##        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='c', lw=2.0, alpha=0.6, marker='s', markersize=3, ls='none', label='Best 2012 (All-Radio)', mew=1)
#        plt.fill_between(xplot, yplot+yerr2, yplot-yerr1,facecolor='red', alpha=0.35)
#        0.15 is base on the bin center in Best 2012
        xplot= np.log10(10**(Best2012[:,0]+0.15))
        yplot = (10**Best2012[:,4])
        yerr2 = 10**(Best2012[:,4]+Best2012[:,5]) - yplot
        yerr1 = yplot - 10**(Best2012[:,4]-Best2012[:,6])
#        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='r', lw=2.0, alpha=0.6, marker='o', markersize=3, ls='none', label='Best 2012 (All-Radio)', mew=1)
        plt.fill_between(xplot, yplot+yerr2, yplot-yerr1,facecolor='red', alpha=0.35,label= 'Best \& Heckman (2012)')
        plt.plot(xaxeshisto_5, counts_5_12   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, color='red',lw = 4, alpha=0.35, label='Best \& Heckman (2012)')

        plt.plot(xaxeshisto_5, counts_5_12   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5,'c', lw=5.0, alpha=0.8, marker='s', markersize=13, label='Jet model (0 $< f_{cool} <$ 1)')
        plt.plot(xaxeshisto_5, counts_5_12  / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'b--',lw=5,alpha=.8)


        plt.yscale('log', nonposy='clip')
        plt.axis([22, 26.5, 1.0e-8, 1.0e-4])
                              
        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
                              
        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # Set the y...
   
        leg = plt.legend(loc='upper right', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
                              
        outputFile = OutputDir + '9_RadioLF' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
    
        # Add this plot to our output list
        OutputList.append(outputFile)
    
    
    # ---------------------------------------------------------
    
    def Density_profile(self, G):
        
        print 'Plotting Density profile of hot gas'
        
        seed(2222)
        
#        plt.figure(figsize=(8,7))  # New figure
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        w = np.where((G.Type == 0) & (G.CentralMvir > 3500) & (G.rho_zero_Makino > 0))[0]
        if(len(w) > dilute): w = sample(w, dilute)
        # Makino density profile 
        r_plot_i  = np.arange(-3,1,0.1)
        r_plot_Makino  = np.zeros(len(r_plot_i))
        rho_gas_Makino = np.zeros(len(r_plot_i))
        j=2
        for i in range(len(r_plot_i)):
            r_plot_Makino[i]=  G.Rvir[w][j] * (10 ** r_plot_i[i])
            rho_gas_Makino[i] = G.rho_zero_Makino[w][j]  * np.exp(-13.5 * G.b_gas[w][j]) * ((1+r_plot_Makino[i] / (G.Rs[w][j])) ** (13.5 * G.b_gas[w][j] /(r_plot_Makino[i]/(G.Rs[w][j]))))
        # Capelo density profile
        w1 = np.where((G.Type == 0) & (G.CentralMvir > 3500) & (G.rho_zero_Makino > 0))[0]
        if(len(w1) > dilute): w1 = sample(w1, dilute)
        
        r_plot_Capelo  = np.zeros(len(r_plot_i))
        rho_gas_Capelo = np.zeros(len(r_plot_i))
        for i in range(len(r_plot_i)):
            r_plot_Capelo[i]=  G.Rvir[w1][j] * (10 ** r_plot_i[i])
            rho_gas_Capelo[i] = G.rho_zero_Capelo[w1][j]* np.exp(-13.5 * G.b_gas[w1][j]) * ((1+r_plot_Capelo[i] / (G.Rs[w][j])) ** (13.5 * G.b_gas[w1][j] /(r_plot_Capelo[i]/(G.Rs[w1][j]))))

        # Isothermal density profile
        w2 = np.where((G.Type == 0) & (G.CentralMvir > 3500) & (G.rho_zero_Makino > 0))[0]
        if(len(w2) > dilute): w2 = sample(w2, dilute)
        
        r_plot_iso  = np.zeros(len(r_plot_i))
        rho_gas_iso = np.zeros(len(r_plot_i))
        for i in range(len(r_plot_i)):
            
            r_plot_iso[i]=  G.Rvir[w1][j] * (10 ** r_plot_i[i])
            rho_gas_iso[i] = G.rho_zero_iso[w1][j]*(1e10 / self.Hubble_h) * (1.989e30)/((3.085678e22)**3)  /((r_plot_iso[i])** 2.0)
            #            rho_gas_iso[i]*=(1e10 / self.Hubble_h) * (1.989e30)/(3.085678e22**3)
            print r_plot_iso[i], np.log10(rho_gas_iso[i]), np.log10(rho_gas_Makino[i]), np.log10(rho_gas_Capelo[i])

        r_plot_beta  = np.zeros(len(r_plot_i))
        rho_gas_beta = np.zeros(len(r_plot_i))
        for i in range(len(r_plot_i)):
            
            r_plot_beta[i]=  G.Rvir[w1][j] * (10 ** r_plot_i[i])
            rho_gas_beta[i] = G.rho_zero_Makino[w1][j]  /(1+(r_plot_beta[i]/(0.22*G.Rs[w1][j]))** 2.0)**(3*0.9*G.b_gas[w1][j]/2)
            print r_plot_iso[i], np.log10(rho_gas_iso[i]), np.log10(rho_gas_Makino[i]), np.log10(rho_gas_Capelo[i]),rho_gas_beta[i]

#        plt.plot(r_plot_Capelo/G.Rvir[w1][j], rho_gas_Capelo,  'g-', lw = 4, alpha=0.5, label='Capelo')
        plt.plot(r_plot_Makino/G.Rvir[w1][j], rho_gas_Makino,  'r-', lw = 4, alpha=0.5, label='Makino et al. (1998)')
        plt.plot(r_plot_beta/G.Rvir[w1][j], rho_gas_beta,  'b--', lw = 4, alpha=0.65, label=(r'$\beta$ - model ( $\beta_{eff}$ = 0.9 b , $r_0$ = 0.22 $r_s$)'))
        plt.plot(r_plot_iso/G.Rvir[w1][j], rho_gas_iso,  'k--', lw = 4, alpha=0.35, label='Isothermal')


        plt.xscale('log', nonposy='clip')
        plt.yscale('log', nonposy='clip')
        
        #        plt.ylabel(r'$\rho  [kg / m^{3}]$ ')  # Set the y...
        plt.ylabel(r'$\rho ~[kg /m^{3}]$ ')  # Set the y...
        plt.xlabel(r'$r/r_{vir}$')  # and the x-axis labels
        
        # Set the x and y axis minor ticks
#        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
        #        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        
        plt.axis([1e-2, 9,1e-27,2e-21])
        
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')
        
        outputFile = OutputDir + '10_Density_profile' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)
# ---------------------------------------------------------

    def cooling_Temp(self, G):
    
        print 'Plotting the b_gas -- Rshock relation'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        w = np.where((G.Type == 0)&(np.log10(G.CentralMvir * 1e10 /self.Hubble_h) > 11) & (G.Cooling>39) & (G.Temp_Gas>1e4))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        E_cooling = G.Cooling[w]-40.0 #- 2*np.log10(self.Hubble_h)
        temp_x      = G.Temp_Gas[w] * 8.617328149741e-8  # [K_b T] in [kev]
#        temp_x =35.9*(G.Vvir[w]*G.Vvir[w]) / 11604.5 / 1.0e3

        plt.scatter(np.log10(temp_x), E_cooling, marker='o', s=10, color = 'grey', alpha=0.1)

#        m_0, b_0  = np.polyfit(np.log10(temp_x), E_cooling, 1)
#        plt.plot(np.log10(temp_x), m_0*(np.log10(temp_x))+b_0, 'b-', lw = 4,label='Linear fit', alpha=0.55)

        total_bins = 30
        X = np.log10(temp_x)
        Y = E_cooling
        bins = np.linspace(X.min(),X.max(), total_bins)
        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=10, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)


        w = np.where((G.Type == 1)&(np.log10(G.CentralMvir* 1e10 /self.Hubble_h) > 11))[0]
        if(len(w) > dilute): w = sample(w, dilute)

        E_cooling = G.Cooling[w]-40.0 #- 3*np.log10(self.Hubble_h)
        temp_x      = G.Temp_Gas[w] * 8.617328149741e-8  # [K_b T] in [kev]


#        plt.scatter(np.log10(temp_x), E_cooling, marker='o', s=20, c='r', alpha=0.25, label='Satellite-Galaxies Jet-model (Heat, Uplift) ')
        #Peres et al.  1998
        Obs = np.array([
                [6.2, 2.7, 3.7, 0.2, 6.2, 4.8, 1.1, 1.2, 46.0, 58.0            ],
                [0.0, 0.0, 0.0, 0.0, 5.1, 0.0, 0.0, 0.0, -1.0, 0.0             ],
                [0.0, 0.0, 0.0, 0.0, 2.4, 0.3, 0.02, 0.02, -1.0, 131.0         ],
                [3.6, 0.3, 0.7, 0.1, 3.6, 0.5, 0.03, 0.04, -1.0, -1.0          ],
                [5.8, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0             ],
                [0.0, 0.0, 0.0, 0.0, 7.8, 0.9, 1.6, 0.9, -1.0, 0.0             ],
                [4.1, 7.8, 0.3, 0.7, 4.1, 7.0, 1.2, 1.2, 690.0, -1.0           ],
                [5.5, 6.8, 999.9, 0.0, 5.5, 13.7, 0.4, 0.4, 42000.0, 21200.0   ],
                [3.0, 3.8, 1.3, 0.3, 3.0, 5.0, 0.3, 0.3, -1.0, 24.1            ],
                [0.0, 0.0, 0.0, 0.0, 5.5, 0.4, 1.0, 0.4, -1.0, -1.0            ],
                [6.8, 14.5, 0.5, 0.3, 6.8, 17.3, 1.7, 2.7, -1.0, 35.5          ],
                [6.2, 0.1, 0.03, 0.03, 6.2, 0.0, 0.7, 0.0, -1.0, -1.0          ],
                [4.7, 1.9, 0.2, 0.1, 4.7, 2.0, 0.1, 0.2, 44.0, -1.0            ],
                [0.0, 0.0, 0.0, 0.0, 5.2, 0.0, 0.1, 0.0, 1922.0, -1.0          ],
                [4.3, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0             ],
                [8.5, 34.9, 1.8, 1.8, 8.5, 43.6, 4.6, 2.8, 480.0, 2370.0       ],
                [6.6, 4.5, 0.6, 1.0, 6.6, 4.1, 2.3, 0.6, -1.0, -1.0            ],
                [8.7, 0.0, 0.0, 0.0, 8.7, 0.0, 0.7, 0.0, -1.0, 9.0             ],
                [3.8, 4.7, 999.9, 0.0, 3.8, 4.7, 1.0, 1.2, 14000.0, 40800.0    ],
                [3.3, 0.1, 0.01, 0.01, 3.3, 0.15, 0.04, 0.1, -1.0, 0.0         ],
                [0.0, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, -1.0, 0.0             ],
                [2.4, 0.1, 999.9, 0.0, 2.4, 0.3, 0.01, 0.01, 61000.0, 22400.0  ],
                [3.6, 0.4, 999.9, 0.0, 3.6, 0.5, 0.1, 0.1, 1530.0, -1.0        ],
                [8.0, 0.01, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 84.0, 207.0          ],
                [4.7, 0.2, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 112.0, 99.0           ],
                [0.0, 0.0, 0.0, 0.0, 4.4, 0.0, 0.3, 0.0, 440.0, 1160.0         ],
                [5.5, 5.4, 4.1, 0.9, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0            ],
                [0.0, 0.0, 0.0, 0.0, 7.0, 3.3, 1.0, 1.0, -1.0, -1.0            ],
                [10.1, 18.0, 1.0, 1.0, 10.1, 26.1, 9.5, 1.5, -1.0, 0.0         ],
                [4.6, 0.02, 0.03, 999.9, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0         ],
                [6.5, 0.8, 0.4, 0.6, 6.5, 1.0, 0.9, 0.2, -1.0, 4.5             ],
                [0.0, 0.0, 0.0, 0.0, 3.8, 0.5, 0.3, 0.3, -1.0, 0.0             ],
                [7.6, 1.8, 0.6, 0.7, 7.6, 1.6, 1.2, 0.5, -1.0, 8.4             ],
                [5.1, 10.4, 0.4, 0.4, 5.1, 8.5, 0.6, 0.2, 260.0, 930.0         ],
                [7.8, 16.8, 3.5, 1.1, 7.8, 17.8, 1.2, 2.9, 89.0, 550.0         ],
                [3.4, 1.7, 0.2, 0.1, 3.4, 1.9, 0.4, 0.02, 1030.0, 5400.0       ],
                [3.0, 1.4, 0.4, 0.5, 3.0, 2.3, 0.1, 0.6, -1.0, 126.0           ],
                [8.4, 0.5, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, -1.0, 14.0            ],
                [0.0, 0.0, 0.0, 0.0, 4.1, 0.6, 0.1, 0.3, -1.0, 15.0            ],
                [11.0, 11.2, 2.1, 2.8, 11.0, 13.0, 1.8, 4.6, -1.0, 0.0         ],
                [4.7, 3.0, 0.2, 0.5, 4.7, 2.7, 0.3, 0.1, 480.0, 3700.0         ],
                [9.0, 40.3, 8.8, 4.2, 9.0, 40.1, 4.6, 4.6, -1.0, 70.1          ],
                [0.0, 0.0, 0.0, 0.0, 7.9, 0.7, 1.7, 999.9, -1.0, -1.0          ],
                [0.0, 0.0, 0.0, 0.0, 7.1, 6.6, 0.8, 4.0, -1.0, 2.41            ],
                [7.5, 0.0, 0.26, 0.0, 7.5, 0.0, 0.2, 0.0, -1.0, 0.0            ],
                [9.0, 1.3, 0.1, 0.01, 9.0, 3.8, 0.9, 2.6, -1.0, 29.1           ],
                [0.0, 0.0, 0.0, 0.0, 7.3, 0.0, 0.1, 0.0, -1.0, 0.0             ],
                [9.9, 0.1, 2.6, 999.9, 9.9, 0.5, 1.5, 999.9, -1.0, 0.0         ],
                [7.3, 8.1, 0.3, 0.4, 7.3, 9.3, 0.8, 1.1, 210000.0, -1.0        ],
                [0.0, 0.0, 0.0, 0.0, 6.5, 0.0, 0.2, 0.0, -1.0, -1.0            ],
                [6.0, 8.2, 0.7, 0.9, 6.0, 8.2, 1.3, 1.9, 410.0, 1880.0         ],
                [0.0, 0.0, 0.0, 0.0, 3.3, 1.1, 0.3, 0.2, -1.0, 28.4            ],
                [3.5, 1.5, 0.4, 0.1, 3.5, 1.6, 0.2, 0.2, 117.0, 1290.0         ],
                ], dtype=np.float32)
    
        cut = 10000.0
        
        OtempH = Obs[:, 0]
        OlumH = Obs[:, 1] * 1.0e4
        OlumHerrU = Obs[:, 2] * 1.0e4
        OlumHerrD = Obs[:, 3] * 1.0e4
                
        OtempP = Obs[:, 4]
        OlumP = Obs[:, 5] * 1.0e4
        OlumPerrU = Obs[:, 6] * 1.0e4
        OlumPerrD = Obs[:, 7] * 1.0e4
                
        radioL = Obs[:, 8]
        radioS = Obs[:, 9]
                
        w = np.where((OlumHerrU < cut*OtempH) & ((OlumHerrD < cut*OtempH)))[0]
        xplot= np.log10(OtempH[w])
        yplot = np.log10(OlumH[w])
        yerr2 = np.log10(OlumH[w]+OlumHerrU[w]) -  yplot
        yerr1 = yplot - np.log10(OlumH[w]-OlumHerrD[w])
        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='b', lw=2.0, alpha=0.3, marker='o', markersize=8, ls='none', label='P98 HRI', mew=1)

        w = np.where((OlumPerrU < cut*OtempP) & ((OlumPerrD < cut*OtempP)))[0]
        xplot = np.log10(OtempP[w])
        yplot = np.log10(OlumP[w])
        yerr2 = np.log10(OlumP[w]+OlumPerrU[w])-yplot
        yerr1 = yplot-np.log10(OlumP[w]-OlumPerrD[w])
        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='b', lw=2.0, alpha=0.3, marker='*', markersize=12, ls='none', label='P98 PSPC', mew=1)

        # Data from Ponman et al 1996
        P_temp = np.array([0.89, 0.44, 0.3, 0.61, 0.91, 0.67, 0.82, 1.09, 0.82, 0.64, 0.96, 0.82, 0.54, 0.59, 0.68, 0.75, 0.87])
        P_temp_err = np.array([0.12, 0.08, 0.05, 0.3, 0.18, 0.11, 0.03, 0.21, 0.27, 0.19, 0.04, 0.19, 0.15, 0.0, 0.12, 0.08, 0.05])
        P_loglum = np.array([42.31, 41.8, 41.68, 41.77, 42.35, 42.12, 42.16, 41.58, 41.98, 41.89, 43.04, 41.69, 41.27, 42.43, 41.48, 42.16, 42.78])
        P_loglum_err = np.array([0.08, 0.12, 0.06, 0.11, 0.11, 0.06, 0.02, 0.14, 0.21, 0.11, 0.03, 0.1, 0.26, 0.24, 0.09, 0.04, 0.02])
                
        P_xplot = np.log10(P_temp)
        P_xerr = np.log10(P_temp+P_temp_err) - P_xplot
        plt.errorbar(P_xplot, P_loglum-40, xerr=P_xerr, yerr=P_loglum_err, color='brown', alpha=0.6, marker='>', markersize=8, ls='none', label=r'P96', mew=1, lw=2)

        # Data from Bharadwaj et al (2015)
        B_temp = np.array([1.9, 1.79, 2.1, 1.43, 0.98, 1.98, 3.58, 2.01, 3.25, 2.06, 1.45, 14.3, 0.81, 1.25, 1.40, 1.06, 0.97, 1.05, 1.96, 2.05, 2.13, 0.64, 2.05, 1.43, 1.78, 0.91])
        B_temp_errH = np.array([0.09, 0.07, 0.05, 0.04, 0.02, 0.04, 0.14, 0.03, 0.08, 0.08, 0.02, 0.04, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.12, 0.1, 0.1, 0.01, 0.1, 0.06, 0.29, 0.01])
        B_temp_errL = np.array([0.09, 0.07, 0.07, 0.05, 0.02, 0.05, 0.14, 0.03, 0.08, 0.09, 0.02, 0.04, 0.02, 0.01, 0.01, 0.01, 0.02, 0.01, 0.15, 0.12, 0.05, 0.01, 0.12, 0.1, 0.22, 0.01])
        B_lum = np.array([2.69, 1.22, 2.02, 0.399, 0.358, 3.42, 3.1, 2.86, 6.92, 3.22, 1.67, 0.413, 0.232, 0.929, 2., 1.25, 0.293, 0.445, 0.444, 2.66, 3.63, 0.111, 2.62, 0.602, 1.98, 0.666])
        B_lum_err = np.array([0.38, 0.17, 0.24, 0.053, 0.139, 0.17, 0.25, 0.05, 0.58, 0.42, 0.02, 0.066, 0.032, 0.15, 0.11, 0.1, 0.04, 0.076, 0.064, 0.22, 0.56, 0.012, 0.74, 1.09, 0.33, 0.045])
                
        B_xplot = np.log10(B_temp)
        B_xerr1 = np.log10(B_temp+B_temp_errH) - B_xplot
        B_xerr2 = B_xplot - np.log10(B_temp-B_temp_errL)
        B_yplot = np.log10(B_lum)+3.0
        B_yerr = np.log10(B_lum+B_lum_err)+3.0 - B_yplot
        plt.errorbar(B_xplot, B_yplot, xerr=[B_xerr1, B_xerr2], yerr=B_yerr, color='orange', alpha=0.6, marker='^', markersize=8, ls='none', label=r'B15', mew=1, lw=2.0)
        
        # Anderson et al. 2015
        A_Ltot = np.array([43.82, 43.46, 43.39, 42.98, 42.64, 42.34, 41.80, 41.52, 41.29, 40.97, 40.58, 40.40, 39.96, 40.10, 39.60, 38.96, 39.94, 40.00, 39.60])
        A_Ltot_uncert_m = np.array([0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.04, 0.07, 0.09, 0.27, 0.19, 0.97, 0.86, 0.21, 0.19, 0.97])
        A_Ltot_uncert_b = np.array([0.21, 0.11, 0.09, 0.06, 0.05, 0.06, 0.06, 0.05, 0.07, 0.11, 0.10, 0.19, 0.46, 0.63, 0.78, 0.83, 0.28, 0.47, 0.86])
        A_Ltot_uncert = np.max(np.array([A_Ltot_uncert_m, A_Ltot_uncert_b]), axis=0)
        A_temp = np.array([5.0, 4.0, 3.4, 2.5, 1.9, 1.5, 1.1, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1])
        plt.errorbar(np.log10(A_temp), A_Ltot-40, yerr=A_Ltot_uncert, marker='_', alpha=0.8, color='purple', markersize=15, ls='none', label=r'A15', mew=2, lw=2)


        plt.xlabel(r' $log_{10} \ T_{\mathrm{new-hot}}~ [Kev]$')  # Set the y...
        plt.ylabel(r' $log_{10}$ (Net cooling $[10^{40}erg~ s^{-1}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #        plt.xscale('log', nonposy='clip')
        #        plt.yscale('log', nonposy='clip')
        plt.axis([-0.55, 1.2, -2, 7])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('small')
        
        outputFile = OutputDir + '11_Cooling_Temp' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)

# ---------------------------------------------------------

    def RadioFraction(self, G):
    
        print 'Plotting the quiescent fraction vs stellar mass'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        groupscale = 12.5
        
        w = np.where((np.log10(G.StellarMass * 1e10/self.Hubble_h)>10))[0]
        Lradio1400 = np.log10(G.RadioLuminosity[w,5])
        Q_jet = G.Qjet[w]
        StellarMass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        CentralMvir = np.log10(G.CentralMvir[w] * 1.0e10 / self.Hubble_h)
        Type = G.Type[w]
        f_cool = G.fcool[w]
#        sSFR = (G.SfrDisk[w] + G.SfrBulge[w]) / (G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        on_radio =  G.t_AGN_on[w]/(G.t_AGN_on[w]+G.t_AGN_off[w])
        Radiocut = 22
        MinRange = 9.0
        MaxRange = 12.0
        Interval = 0.1
        Nbins = int((MaxRange-MinRange)/Interval)
        Range = np.arange(MinRange, MaxRange, Interval)
        
        Mass = []
        Fraction = []
        CentralFraction = []
        SatelliteFraction = []
        SatelliteFractionLo = []
        SatelliteFractionHi = []
        fract_on_Fraction = []
        for i in xrange(Nbins-1):
            
            w = np.where((StellarMass >= Range[i]) & (StellarMass < Range[i+1]))[0]
            if len(w) > 0:
                wQ = np.where((StellarMass >= Range[i]) & (StellarMass < Range[i+1])&(Lradio1400 > Radiocut))[0]
                Fraction.append( 1.0 * len(wQ) / len(w))
#                fract_on = np.median(1 - f_cool[wQ])
                fract_on = np.mean(on_radio[wQ])
                fract_on_Fraction.append(fract_on * 1.0 * len(wQ) / len(w))
            else:
                Fraction.append(0.0)
                fract_on_Fraction.append(0.0)

            w = np.where( (Type == 0) &(StellarMass >= Range[i]) & (StellarMass < Range[i+1]))[0]
            if len(w) > 0:
                wQ = np.where( (Type == 0) &(StellarMass >= Range[i]) & (StellarMass < Range[i+1])&(Lradio1400 > Radiocut) )[0]
                CentralFraction.append(1.0*len(wQ) / len(w))
            else:
                CentralFraction.append(0.0)
            
            w = np.where((Type == 1) &(StellarMass >= Range[i]) & (StellarMass < Range[i+1]))[0]
            if len(w) > 0:
                wQ = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1])&(Lradio1400 > Radiocut) )[0]
                SatelliteFraction.append(1.0*len(wQ) / len(w))
                wQ = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (CentralMvir < groupscale)&(Lradio1400 > Radiocut))[0]
                SatelliteFractionLo.append(1.0*len(wQ) / len(w))
                wQ = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (CentralMvir > groupscale)&(Lradio1400 > Radiocut))[0]
                SatelliteFractionHi.append(1.0*len(wQ) / len(w))
            else:
                SatelliteFraction.append(0.0)
                SatelliteFractionLo.append(0.0)
                SatelliteFractionHi.append(0.0)

            Mass.append((Range[i] + Range[i+1]) / 2.0)
        # print '  ', Mass[i], Fraction[i], CentralFraction[i], SatelliteFraction[i]
        
        Mass = np.array(Mass)
        Fraction = np.array(Fraction)
        CentralFraction = np.array(CentralFraction)
        SatelliteFraction = np.array(SatelliteFraction)
        SatelliteFractionLo = np.array(SatelliteFractionLo)
        SatelliteFractionHi = np.array(SatelliteFractionHi)
        fract_on_Fraction = np.array(fract_on_Fraction)
#        w = np.where((Fraction > 0))[0]
#        plt.plot(Mass[w],  Fraction[w], 'b-', label='Jet-model -All')

        w = np.where((fract_on_Fraction > 0))[0]
        plt.plot(Mass[w],  fract_on_Fraction[w], 'c', lw=5.0, alpha=0.8, marker='s', markersize=10, label='Jet model (Active-AGN)')
        plt.plot(Mass[w],  fract_on_Fraction[w], 'b--', lw=5.0)
        Best2007 = np.array([
                    [10.200254141026278, 2.6910517657478307E-4,  0.00009240177381,0.000101242701],
                    [10.399533979171528, 0.0010630129078379156, 0.0001766265079, 0.0002088631324],
                    [10.600430643785689, 0.0023178891610480174,0.0002990764554 , 0.0003625841001],
                    [10.79817970772826, 0.008052101871536859, 0.0006896188308, 0.0009945327869],
                    [10.998982670985853, 0.03856776775500328, 0.00364806146, 0.003916260167],
                    [11.199833441058022, 0.1236436535322261, 0.01107479225, 0.01225214058],
                    [11.400725324990725, 0.28064851443748245, 0.02481196515, 0.02919245295],
                    [11.598593905969732, 0.3573218329016989, 0.1015240761, 0.07206889011],
                    ], dtype=np.float32)


        Best2007_xval = np.log10(10**Best2007[:,0]/self.Hubble_h)
        Best2027_yval = Best2007[:,1]
        yerr2 = Best2007[:,2]
        yerr1 = Best2007[:,3]
        #        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='r', lw=2.0, alpha=0.6, marker='o', markersize=3, ls='none', label='Best 2012 (All-Radio)', mew=1)
        plt.fill_between(Best2007_xval, Best2027_yval+yerr2, Best2027_yval-yerr1,facecolor='red', alpha=0.35)
        plt.plot(Best2007_xval, Best2027_yval, color='r',lw = 4, alpha=0.35, label='Best et al. (2007)')
        
        
        
        plt.xlabel(r'$\log_{10} M_{\mathrm{stellar}}\ (M_{\odot})$')  # Set the x-axis label
        plt.ylabel(r'Fraction~(Radio-loud)')  # Set the y-axis label
        
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        plt.yscale('log', nonposy='clip')
        plt.axis([10.3, 12.0, 0.0001, 1.4])
        
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '12_RadioFraction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)

# ---------------------------------------------------------

    def Rshock_hist(self, G):
    
        print 'Plotting the Rshock_Rheat_hist relation'
        
        seed(2222)
#        plt.figure(figsize=(16,6))  # New figure
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        w = np.where((G.Type == 0) & (np.log10(G.CentralMvir * 1.0e10 / self.Hubble_h) > 11.0) & (G.StellarMass>0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
#        Rshock = G.Rshocked[w] * 1000.0 # Kpc
        Rshock = G.Rshocked[w]/G.Rvir[w]
        plt.scatter(mass, Rshock, marker='o', s=10, color = 'grey', alpha=0.05)

        total_bins = 70
        X = mass
        Y = Rshock
        bins = np.linspace(X.min(),X.max(), total_bins)
#        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)
        
        
        x=[9,11,12,13]
        y=[1,1,1,1]
        plt.plot(x, y, 'k--', lw=2)
        
        plt.ylabel(r'$R_{Shock}/R_{vir}$')  # Set the y...
        plt.xlabel(r'$\log_{10} (M_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
#        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
#        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        
        plt.axis([10, 12.0, 0, 1.3])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '13_Rshock_hist' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)
# ---------------------------------------------------------

    def Lradio_BHmass(self, G):
    
        print 'Plotting the delta-- mass relation'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        w = np.where((G.Type == 0)&(np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>11) &(G.BlackHoleMass>0)&(G.RadioLuminosity[:,5]>0)&(G.RadioLuminosity[:,5]< 1e40))[0]
        #        if(len(w) > dilute): w = sample(w, dilute)
        
        Lradio1400 = np.log10(G.RadioLuminosity[w,5])
        BHmass = np.log10(G.BlackHoleMass[w]* 1e10/self.Hubble_h)
        plt.scatter(Lradio1400, BHmass, marker='o', s=10, color = 'grey', alpha=0.15)
        
        total_bins = 70
        X = Lradio1400
        Y = BHmass
        bins = np.linspace(X.min(),X.max(), total_bins)
#        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=4,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='Jet-model (Median$~\pm~\sigma$)', mew=1)
        
        
#        High excitation
        Best2012_HE = np.array([
        [23.222246074630203,  7.7987147970657595],
        [23.722141243939475,  8.175823264965134],
        [24.219926695597884,  8.183359544900306],
        [24.72082546460808,  8.175376105084753],
        [25.220922972566733,  8.474845075460957],
        [25.7155789197814,  8.504112234948709],
                ], dtype=np.float32)
        Best2012_LE = np.array([
        [23.279100537257293,  8.385824813665208],
        [23.779837435347986,  8.439952971868191],
        [24.277452922540913,  8.512706429722828],
        [24.77814935290173,  8.582362487430446],
        [25.28191053933849,  8.6551323047929],
        [25.78003591992785,  8.532234228889147],
                ], dtype=np.float32)
        
        
        Best2012_HE_xval = Best2012_HE[:,0]
        Best2012_HE_yval = Best2012_HE[:,1]
        if(whichimf == 1):  Best2012_HE_yval = Best2012_HE_yval - 0.26
        
        Best2012_LE_xval = Best2012_LE[:,0]
        Best2012_LE_yval = Best2012_LE[:,1]
        if(whichimf == 1):  Best2012_LE_yval = Best2012_LE_yval - 0.26
        
        plt.errorbar(Best2012_HE_xval, Best2012_HE_yval, yerr=0.0, color='red', alpha=0.6, marker='o', markersize=17, ls='none', label=r'Best \& Heckman (2012)-HERGs', mew=1, lw=2.0)
        
        plt.errorbar(Best2012_LE_xval, Best2012_LE_yval, yerr=0.0, color='orange', alpha=0.6, marker='d', markersize=17, ls='none', label=r'Best \& Heckman (2012)-LERGs', mew=1, lw=2.0)
        
        plt.ylabel(r'$\log_{10} (M_{\mathrm{BH}}\ [M_{\odot}])$')  # Set the y...
        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # and the x-axis labels
 
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        plt.axis([23, 26, 7, 9.0])
        
        
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '14_Lradio_BHmass' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)

# =================================================================


#  'Main' section of code.  This if statement executes if the code is run from the 
#   shell command line, i.e. with 'python allresults.py'
if __name__ == '__main__':

    from optparse import OptionParser
    import os

    parser = OptionParser()
    parser.add_option(
        '-d',
        '--dir_name',
        dest='DirName',
        default='./results/millennium/',
        help='input directory name (default: ./results/millennium/)',
              
        metavar='DIR',
        )
    parser.add_option(
        '-f',
        '--file_base',
        dest='FileName',
        default='model_z0.000',
        help='filename base (default: model_z0.000)',
        metavar='FILE',
        )
    parser.add_option(
        '-n',
        '--file_range',
        type='int',
        nargs=2,
        dest='FileRange',
        default=(0, 7),
        help='first and last filenumbers (default: 0 7)',
        metavar='FIRST LAST',
        )

    (opt, args) = parser.parse_args()

    if opt.DirName[-1] != '/':
        opt.DirName += '/'

    OutputDir = opt.DirName + 'plots/'

    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    res = Results()

    print 'Running allresults...'

    FirstFile = opt.FileRange[0]
    LastFile = opt.FileRange[1]

    fin_base = opt.DirName + opt.FileName
    G = res.read_gals(fin_base, FirstFile, LastFile)

    res.StellarMassFunction(G)
    res.Metallicity(G)
    res.BlackHoleBulgeRelationship(G)
    res.Lradio_Qjet(G)
    res.Rshocked_Rvir(G)
    res.Lradio_Rshock(G)
    res.Lradio_Mass(G)
    res.Temp_hist(G)
    res.RadioLF(G)
    res.Density_profile(G)
    res.cooling_Temp(G)
    res.RadioFraction(G)
    res.Rshock_hist(G)
    res.Lradio_BHmass(G)
