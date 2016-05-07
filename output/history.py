#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import h5py as h5
import numpy as np
import pylab as plt
from random import sample, seed
from os.path import getsize as getFileSize

# ================================================================================
# Basic variables
# ================================================================================

# Set up some basic attributes of the run

whichsimulation = 0
whichimf = 1        # 0=Slapeter; 1=Chabrier


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

#OutputDir = '/Users/mojtabaraouf/sage-master/output/new_SAM' # set in main below

#Mac book Air
#OutputDir = '/Volumes/09158118409/sage-master/output/new_SAM/allresult/model/'

#IMAC
OutputDir = '/Users/mraoufhajarzarrin/Desktop/sage-master/output/new_SAM/allresult/model/'


OutputFormat = '.png'
TRANSPARENT = False

OutputList = []


class Results:

    """ The following methods of this class generate the figures and plot them.
    """

    def __init__(self):
        """Here we set up some of the variables which will be global to this
        class."""

        if whichsimulation == 0:    # Mini-Millennium
          self.Hubble_h = 0.73
          self.BoxSize = 62.5       # Mpc/h
          self.MaxTreeFiles = 8     # FilesPerSnapshot

        elif whichsimulation == 1:  # Full Millennium
          self.Hubble_h = 0.73
          self.BoxSize = 500        # Mpc/h
          self.MaxTreeFiles = 512   # FilesPerSnapshot

        else:
          print "Please pick a valid simulation!"
          exit(1)


        if whichsimulation == 0 or whichsimulation == 1 :
          
#          self.SMFsnaps = [63, 37, 32, 27, 23, 20, 18, 16]
#          self.SMFsnaps = [63, 48, 41]
          self.SMFsnaps = [63, 62, 61, 60, 59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40, 37, 32, 27, 23, 20, 18, 16]
          self.SMFz =[0.000, 0.020, 0.041, 0.064, 0.089, 0.116, 0.144, 0.175, 0.208, 0.242, 0.280, 0.320, 0.362, 0.408, 0.457, 0.509, 0.564, 0.624, 0.687,0.755,0.828,0.905,0.989]

          self.redshift_file = ['_z127.000', '_z79.998', '_z50.000', '_z30.000', '_z19.916', '_z18.244', '_z16.725', '_z15.343', '_z14.086', '_z12.941', '_z11.897', '_z10.944', '_z10.073', '_z9.278', '_z8.550', '_z7.883', '_z7.272', '_z6.712', '_z6.197', '_z5.724', '_z5.289', '_z4.888', '_z4.520', '_z4.179', '_z3.866', '_z3.576', '_z3.308', '_z3.060', '_z2.831', '_z2.619', '_z2.422', '_z2.239', '_z2.070', '_z1.913', '_z1.766', '_z1.630', '_z1.504', '_z1.386', '_z1.276', '_z1.173', '_z1.078', '_z0.989', '_z0.905', '_z0.828', '_z0.755', '_z0.687', '_z0.624', '_z0.564', '_z0.509', '_z0.457', '_z0.408', '_z0.362', '_z0.320', '_z0.280', '_z0.242', '_z0.208', '_z0.175', '_z0.144', '_z0.116', '_z0.089', '_z0.064', '_z0.041', '_z0.020', '_z0.000']

          self.redshift = [127.000, 79.998, 50.000, 30.000, 19.916, 18.244, 16.725, 15.343, 14.086, 12.941, 11.897, 10.944, 10.073, 9.278, 8.550, 7.883, 7.272, 6.712, 6.197, 5.724, 5.289, 4.888, 4.520, 4.179, 3.866, 3.576, 3.308, 3.060, 2.831, 2.619, 2.422, 2.239, 2.070, 1.913, 1.766, 1.630, 1.504, 1.386, 1.276, 1.173, 1.078, 0.989, 0.905, 0.828, 0.755, 0.687, 0.624, 0.564, 0.509, 0.457, 0.408, 0.362, 0.320, 0.280, 0.242, 0.208, 0.175, 0.144, 0.116, 0.089, 0.064, 0.041, 0.020, 0.000]


    def read_gals(self, model_name, first_file, last_file, thissnap):

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
            ('RadioLuminosity'              , (np.float32,7)),
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
        goodfiles = 0
            
        if thissnap in self.SMFsnaps:

            print
            print "Determining array storage requirements."
        
            # Read each file and determine the total number of galaxies to be read in
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

            print "Input files contain:\t%d trees ;\t%d galaxies ." % (TotNTrees, TotNGals)

        # Initialize the storage array
        G = np.empty(TotNGals, dtype=Galdesc)

        if thissnap in self.SMFsnaps:

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


            print "Total galaxies considered:", TotNGals
            print

        # Convert the Galaxy array into a recarray
        G = G.view(np.recarray)

        # Calculate the volume given the first_file and last_file
        self.volume = self.BoxSize**3.0 * goodfiles / self.MaxTreeFiles

        return G


# --------------------------------------------------------

    def StellarMassFunction(self, G_history):

        print 'Plotting the stellar mass function'

        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        binwidth = 0.15  # mass function histogram bin width

        # Marchesini et al. 2009ApJ...701.1765M SMF, h=0.7
        M = np.arange(7.0, 11.8, 0.01)
        Mstar = np.log10(10.0**10.96)
        alpha = -1.18
        phistar = 30.87*1e-4
        xval = 10.0 ** (M-Mstar)
        yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)
        if(whichimf == 0):
            plt.plot(np.log10(10.0**M *1.6), yval, 'r', lw=10, alpha=0.3, label='Marchesini et al. 2009 z=[0.1]')
        elif(whichimf == 1):
            plt.plot(np.log10(10.0**M *1.6 /1.8), yval, 'r', lw=10, alpha=0.3, label='Marchesini et al. 2009 z=[0.1]')
            

        M = np.arange(9.3, 11.8, 0.01)
        Mstar = np.log10(10.0**10.91)
        alpha = -0.99
        phistar = 10.17*1e-4
        xval = 10.0 ** (M-Mstar)
        yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)
        if(whichimf == 0):
            plt.plot(np.log10(10.0**M *1.6), yval, 'b', lw=10, alpha=0.3, label='Marchesini et al. 2009 z=[1.3,2.0]')
        elif(whichimf == 1):
            plt.plot(np.log10(10.0**M *1.6/1.8), yval, 'b', lw=10, alpha=0.3, label='Marchesini et al. 2009 z=[1.3,2.0]')
#
#        M = np.arange(9.7, 11.8, 0.01)
#        Mstar = np.log10(10.0**10.96)
#        alpha = -1.01
#        phistar = 3.95*1e-4
#        xval = 10.0 ** (M-Mstar)
#        yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)
#        if(whichimf == 0):
#            plt.plot(np.log10(10.0**M *1.6), yval, 'g:', lw=10, alpha=0.5, label='... z=[2.0,3.0]')
#        elif(whichimf == 1):
#            plt.plot(np.log10(10.0**M *1.6/1.8), yval, 'g:', lw=10, alpha=0.5, label='... z=[2.0,3.0]')
#
#        M = np.arange(10.0, 11.8, 0.01)
#        Mstar = np.log10(10.0**11.38)
#        alpha = -1.39
#        phistar = 0.53*1e-4
#        xval = 10.0 ** (M-Mstar)
#        yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)
#        if(whichimf == 0):
#            plt.plot(np.log10(10.0**M *1.6), yval, 'r:', lw=10, alpha=0.5, label='... z=[3.0,4.0]')
#        elif(whichimf == 1):
#            plt.plot(np.log10(10.0**M *1.6/1.8), yval, 'r:', lw=10, alpha=0.5, label='... z=[3.0,4.0]')


        ###### z=0
        
        w = np.where(G_history[self.SMFsnaps[0]].StellarMass > 0.0)[0]
        mass = np.log10(G_history[self.SMFsnaps[0]].StellarMass[w] * 1.0e10 /self.Hubble_h)

        mi = np.floor(min(mass)) - 2
        ma = np.floor(max(mass)) + 2
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        # Overplot the model histograms
        plt.plot(xaxeshisto, counts / self.volume *self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'b-', lw = 4, alpha=0.65, label='Jet model (z = 0)')

        ###### z=1.3  SMFsnap[1] to [24]
        
        w = np.where(G_history[self.SMFsnaps[24]].StellarMass > 0.0)[0]
        mass = np.log10(G_history[self.SMFsnaps[24]].StellarMass[w] * 1.0e10 /self.Hubble_h)

        mi = np.floor(min(mass)) - 2
        ma = np.floor(max(mass)) + 2
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        # Overplot the model histograms
        plt.plot(xaxeshisto, counts / self.volume *self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'm--', lw = 4, alpha=0.65, label='Jet model (z = 0.5)')

        ###### z=2     SMFsnap[2] to [25]
        
        w = np.where(G_history[self.SMFsnaps[25]].StellarMass > 0.0)[0]
        mass = np.log10(G_history[self.SMFsnaps[25]].StellarMass[w] * 1.0e10 /self.Hubble_h)

        mi = np.floor(min(mass)) - 2
        ma = np.floor(max(mass)) + 2
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        # Overplot the model histograms
        plt.plot(xaxeshisto, counts / self.volume *self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'r:', lw = 4, alpha=0.65, label='Jet model (z = 1.0)')

#        ###### z=3    SMFsnap[3] to [26]
#        
#        w = np.where(G_history[self.SMFsnaps[26]].StellarMass > 0.0)[0]
#        mass = np.log10(G_history[self.SMFsnaps[26]].StellarMass[w] * 1.0e10 /self.Hubble_h)
#
#        mi = np.floor(min(mass)) - 2
#        ma = np.floor(max(mass)) + 2
#        NB = (ma - mi) / binwidth
#
#        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)
#
#        # Set the x-axis values to be the centre of the bins
#        xaxeshisto = binedges[:-1] + 0.5 * binwidth
#
#        # Overplot the model histograms
#        plt.plot(xaxeshisto, counts / self.volume *self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'r-')
#
        ######

        plt.yscale('log', nonposy='clip')

        plt.axis([8.0, 12.5, 2.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1}$)')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$')  # and the x-axis labels

        leg = plt.legend(loc='lower left', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = OutputDir + 'A_StellarMassFunction_z' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)



# ---------------------------------------------------------

    def PlotHistory_SFRdensity(self, G_history):
    
        print 'Plotting SFR density evolution for all galaxies'

        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        #Somerville et al. (2001)
        ObsSFRdensity = np.array([
            [0, 0.0158489, 0, 0, 0.0251189, 0.01000000],
            [0.150000, 0.0173780, 0, 0.300000, 0.0181970, 0.0165959],
            [0.0425000, 0.0239883, 0.0425000, 0.0425000, 0.0269153, 0.0213796],
            [0.200000, 0.0295121, 0.100000, 0.300000, 0.0323594, 0.0269154],
            [0.350000, 0.0147911, 0.200000, 0.500000, 0.0173780, 0.0125893],
            [0.625000, 0.0275423, 0.500000, 0.750000, 0.0331131, 0.0229087],
            [0.825000, 0.0549541, 0.750000, 1.00000, 0.0776247, 0.0389045],
            [0.625000, 0.0794328, 0.500000, 0.750000, 0.0954993, 0.0660693],
            [0.700000, 0.0323594, 0.575000, 0.825000, 0.0371535, 0.0281838],
            [1.25000, 0.0467735, 1.50000, 1.00000, 0.0660693, 0.0331131],
            [0.750000, 0.0549541, 0.500000, 1.00000, 0.0389045, 0.0776247],
            [1.25000, 0.0741310, 1.00000, 1.50000, 0.0524807, 0.104713],
            [1.75000, 0.0562341, 1.50000, 2.00000, 0.0398107, 0.0794328],
            [2.75000, 0.0794328, 2.00000, 3.50000, 0.0562341, 0.112202],
            [4.00000, 0.0309030, 3.50000, 4.50000, 0.0489779, 0.0194984],
            [0.250000, 0.0398107, 0.00000, 0.500000, 0.0239883, 0.0812831],
            [0.750000, 0.0446684, 0.500000, 1.00000, 0.0323594, 0.0776247],
            [1.25000, 0.0630957, 1.00000, 1.50000, 0.0478630, 0.109648],
            [1.75000, 0.0645654, 1.50000, 2.00000, 0.0489779, 0.112202],
            [2.50000, 0.0831764, 2.00000, 3.00000, 0.0512861, 0.158489],
            [3.50000, 0.0776247, 3.00000, 4.00000, 0.0416869, 0.169824],
            [4.50000, 0.0977237, 4.00000, 5.00000, 0.0416869, 0.269153],
            [5.50000, 0.0426580, 5.00000, 6.00000, 0.0177828, 0.165959],
            [3.00000, 0.120226, 2.00000, 4.00000, 0.173780, 0.0831764],
            [3.04000, 0.128825, 2.69000, 3.39000, 0.151356, 0.109648],
            [4.13000, 0.114815, 3.78000, 4.48000, 0.144544, 0.0912011],
            [0.350000, 0.0346737, 0.200000, 0.500000, 0.0537032, 0.0165959],
            [0.750000, 0.0512861, 0.500000, 1.00000, 0.0575440, 0.0436516],
            [1.50000, 0.0691831, 1.00000, 2.00000, 0.0758578, 0.0630957],
            [2.50000, 0.147911, 2.00000, 3.00000, 0.169824, 0.128825],
            [3.50000, 0.0645654, 3.00000, 4.00000, 0.0776247, 0.0512861],
            ], dtype=np.float32)
            
        ObsRedshift = ObsSFRdensity[:, 0]
        xErrLo = ObsSFRdensity[:, 0]-ObsSFRdensity[:, 2]
        xErrHi = ObsSFRdensity[:, 3]-ObsSFRdensity[:, 0]
            
        ObsSFR = np.log10(ObsSFRdensity[:, 1])
        yErrLo = np.log10(ObsSFRdensity[:, 1])-np.log10(ObsSFRdensity[:, 4])
        yErrHi = np.log10(ObsSFRdensity[:, 5])-np.log10(ObsSFRdensity[:, 1])
            
        # plot observational data (compilation used in Croton et al. 2006)
        plt.errorbar(ObsRedshift, ObsSFR, yerr=[yErrLo, yErrHi], xerr=[xErrLo, xErrHi], color='r', lw=1.0, alpha=0.3, marker='o', ls='none', label='Somerville et al. (2001)')
            
            
        C06SFRdensity = np.array([
             [ 7.272   ,   -1.5079008 ],
             [ 6.197   ,   -1.23429604],
             [ 5.289   ,   -1.0677671 ],
             [ 4.179   ,   -0.93931313],
             [ 3.06    ,   -0.91103558],
             [ 2.07    ,   -0.99985986],
             [ 1.386   ,   -1.15112987],
             [ 1.078   ,   -1.25634226],
             [ 0.989   ,   -1.29094775],
             [ 0.905   ,   -1.32623472],
             [ 0.828   ,   -1.36155264],
             [ 0.755   ,   -1.39545521],
             [ 0.687   ,   -1.42912022],
             [ 0.624   ,   -1.46213665],
             [ 0.564   ,   -1.49477694],
             [ 0.509   ,   -1.5265802 ],
             [ 0.457   ,   -1.55854078],
             [ 0.408   ,   -1.58770231],
             [ 0.362   ,   -1.61646077],
             [ 0.32    ,   -1.64317747],
             [ 0.28    ,   -1.67014686],
             [ 0.242   ,   -1.69629931],
             [ 0.208   ,   -1.72010743],
             [ 0.175   ,   -1.74408391],
             [ 0.144   ,   -1.76571565],
             [ 0.116   ,   -1.78554315],
             [ 0.089   ,   -1.80602209],
             [ 0.064   ,   -1.8242581 ],
             [ 0.041   ,   -1.84199117],
             [ 0.02    ,   -1.85798567],
             [ 0.      ,   -1.87108685],
             ], dtype=np.float32)
        C06SFRdensity_xval = C06SFRdensity[:, 0]
        C06SFRdensity_yval = C06SFRdensity[:, 1]
             
        plt.plot(C06SFRdensity_xval, C06SFRdensity_yval, 'm-', lw=3.0, label='C16 (Radio mode)')
        
        WOAGNSFRdensity = np.array([
             [ 7.272   ,   -1.48943879],
             [ 6.197   ,   -1.19834185],
             [ 5.289   ,   -1.00922241],
             [ 4.179   ,   -0.84202861],
             [ 3.06    ,   -0.76386072],
             [ 2.07    ,   -0.7946669 ],
             [ 1.386   ,   -0.88909161],
             [ 1.078   ,   -0.95986964],
             [ 0.989   ,   -0.9837791 ],
             [ 0.905   ,   -1.00762045],
             [ 0.828   ,   -1.03105964],
             [ 0.755   ,   -1.05505066],
             [ 0.687   ,   -1.07824001],
             [ 0.624   ,   -1.10115936],
             [ 0.564   ,   -1.12401459],
             [ 0.509   ,   -1.14609809],
             [ 0.457   ,   -1.16839575],
             [ 0.408   ,   -1.18952079],
             [ 0.362   ,   -1.20905302],
             [ 0.32    ,   -1.22823874],
             [ 0.28    ,   -1.24810857],
             [ 0.242   ,   -1.26749936],
             [ 0.208   ,   -1.28488588],
             [ 0.175   ,   -1.30199313],
             [ 0.144   ,   -1.31880621],
             [ 0.116   ,   -1.33450582],
             [ 0.089   ,   -1.34964766],
             [ 0.064   ,   -1.36357204],
             [ 0.041   ,   -1.37715987],
             [ 0.02    ,   -1.38956864],
             [ 0.      ,   -1.40131551],
            ], dtype=np.float32)
        WOAGNSFRdensity_xval = WOAGNSFRdensity[:, 0]
        WOAGNSFRdensity_yval = WOAGNSFRdensity[:, 1]
 
        plt.plot(WOAGNSFRdensity_xval, WOAGNSFRdensity_yval, 'k-', lw=3.0, label='W/O AGN Jet')
        
        SFR_density = np.zeros((LastSnap+1-FirstSnap))       
        for snap in xrange(FirstSnap,LastSnap+1):
          SFR_density[snap-FirstSnap] = sum(G_history[snap].SfrDisk+G_history[snap].SfrBulge) / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h

        z = np.array(self.redshift)
        nonzero = np.where(SFR_density > 0.0)[0]
        plt.plot(z[nonzero], np.log10(SFR_density[nonzero]), 'b--', lw=3.0)
        plt.plot(z[nonzero], np.log10(SFR_density[nonzero]),color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, label='Jet-model ')
#        print np.column_stack((z[nonzero], np.log10(SFR_density[nonzero])))
        plt.ylabel(r'$\log_{10} \mathrm{SFR\ density}\ (M_{\odot}\ \mathrm{yr}^{-1}\ \mathrm{Mpc}^{-3})$')  # Set the y...
        plt.xlabel(r'$\mathrm{redshift}$')  # and the x-axis labels
    
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    
        plt.axis([0.0, 8.0, -2.0, -0.4])
        
        leg = plt.legend(loc='upper right', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
    
        outputFile = OutputDir + 'B_History-SFR-density' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
    
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def StellarMassDensityEvolution(self, G_history):

        print 'Plotting stellar mass density evolution'

        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        # SMD observations taken from Marchesini+ 2009, h=0.7
        # Values are (minz, maxz, rho,-err,+err)
        dickenson2003 = np.array(((0.6,1.4,8.26,0.08,0.08),
                         (1.4,2.0,7.86,0.22,0.33),
                         (2.0,2.5,7.58,0.29,0.54),
                         (2.5,3.0,7.52,0.51,0.48)),float)
        drory2005 = np.array(((0.25,0.75,8.3,0.15,0.15),
                    (0.75,1.25,8.16,0.15,0.15),
                    (1.25,1.75,8.0,0.16,0.16),
                    (1.75,2.25,7.85,0.2,0.2),
                    (2.25,3.0,7.75,0.2,0.2),
                    (3.0,4.0,7.58,0.2,0.2)),float)
        # Perez-Gonzalez (2008)
        pg2008 = np.array(((0.2,0.4,8.41,0.06,0.06),
                 (0.4,0.6,8.37,0.04,0.04),
                 (0.6,0.8,8.32,0.05,0.05),
                 (0.8,1.0,8.24,0.05,0.05),
                 (1.0,1.3,8.15,0.05,0.05),
                 (1.3,1.6,7.95,0.07,0.07),
                 (1.6,2.0,7.82,0.07,0.07),
                 (2.0,2.5,7.67,0.08,0.08),
                 (2.5,3.0,7.56,0.18,0.18),
                 (3.0,3.5,7.43,0.14,0.14),
                 (3.5,4.0,7.29,0.13,0.13)),float)
        glazebrook2004 = np.array(((0.8,1.1,7.98,0.14,0.1),
                         (1.1,1.3,7.62,0.14,0.11),
                         (1.3,1.6,7.9,0.14,0.14),
                         (1.6,2.0,7.49,0.14,0.12)),float)
        fontana2006 = np.array(((0.4,0.6,8.26,0.03,0.03),
                      (0.6,0.8,8.17,0.02,0.02),
                      (0.8,1.0,8.09,0.03,0.03),
                      (1.0,1.3,7.98,0.02,0.02),
                      (1.3,1.6,7.87,0.05,0.05),
                      (1.6,2.0,7.74,0.04,0.04),
                      (2.0,3.0,7.48,0.04,0.04),
                      (3.0,4.0,7.07,0.15,0.11)),float)
        rudnick2006 = np.array(((0.0,1.0,8.17,0.27,0.05),
                      (1.0,1.6,7.99,0.32,0.05),
                      (1.6,2.4,7.88,0.34,0.09),
                      (2.4,3.2,7.71,0.43,0.08)),float)
        elsner2008 = np.array(((0.25,0.75,8.37,0.03,0.03),
                     (0.75,1.25,8.17,0.02,0.02),
                     (1.25,1.75,8.02,0.03,0.03),
                     (1.75,2.25,7.9,0.04,0.04),
                     (2.25,3.0,7.73,0.04,0.04),
                     (3.0,4.0,7.39,0.05,0.05)),float)

        obs = (dickenson2003,drory2005,pg2008,glazebrook2004,
               fontana2006,rudnick2006,elsner2008)

        for o in obs:
            xval = ((o[:,1]-o[:,0])/2.)+o[:,0]
            if(whichimf == 0):
                ax.errorbar(xval, np.log10(10**o[:,2] *1.6), xerr=(xval-o[:,0], o[:,1]-xval), yerr=(o[:,3], o[:,4]), alpha=0.3, lw=1.0, marker='o', ls='none')
            elif(whichimf == 1):
                ax.errorbar(xval, np.log10(10**o[:,2] *1.6/1.8), xerr=(xval-o[:,0], o[:,1]-xval), yerr=(o[:,3], o[:,4]), alpha=0.3, lw=1.0, marker='o', ls='none')
                

        smd = np.zeros((LastSnap+1-FirstSnap))       

        for snap in xrange(FirstSnap,LastSnap+1):
          w = np.where((G_history[snap].StellarMass/self.Hubble_h > 0.01) & (G_history[snap].StellarMass/self.Hubble_h < 1000.0))[0]
          if(len(w) > 0):
            smd[snap-FirstSnap] = sum(G_history[snap].StellarMass[w]) *1.0e10/self.Hubble_h / (self.volume /self.Hubble_h/self.Hubble_h/self.Hubble_h)

        z = np.array(self.redshift)
        nonzero = np.where(smd > 0.0)[0]
        plt.plot(z[nonzero], np.log10(smd[nonzero]), 'k-', lw=3.0, label='Jet-model ')

        plt.ylabel(r'$\log_{10}\ \phi\ (M_{\odot}\ \mathrm{Mpc}^{-3})$')  # Set the y...
        plt.xlabel(r'$\mathrm{redshift}$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

        plt.axis([0.0, 4.2, 6.5, 9.0])
        
        leg = plt.legend(loc='lower left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = OutputDir + 'C_History-stellar-mass-density' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)



# ---------------------------------------------------------

    def RadioLF(self, G_history):
    
        print 'Plotting the Radio Luminosity function'
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        binwidth = 0.6  # Radio Luminosity function histogram bin width
        
        snap = 63
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0)& (G_history[snap].fcool>0))[0]
        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        mi_5 = np.floor(min(Lradio1400_5)) - 3.0
        ma_5 = np.floor(max(Lradio1400_5)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        (counts_5, binedges_5) = np.histogram(Lradio1400_5, range=(mi_5, ma_5), bins=NB_5 )
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
#        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 2, alpha=0.55, label='Jet model, 63')

        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5,'c', lw=3.0, alpha=0.8, marker='s', markersize=10, label='Jet model (z = 0)')
        plt.plot(xaxeshisto_5, counts_5  / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'b--',lw=3,alpha=.8)


        snap = 48
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0)& (G_history[snap].fcool>0))[0]
        #        if(len(w1) > dilute): w1 = sample(w1, dilute)
        
        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        mi_5 = np.floor(min(Lradio1400_5)) - 3.0
        ma_5 = np.floor(max(Lradio1400_5)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        (counts_5, binedges_5) = np.histogram(Lradio1400_5, range=(mi_5, ma_5), bins=NB_5 )
        # Set the x-axis values to be the centre of the bins
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'm--',lw = 4, alpha=0.95, label='Jet model (z = 0.5)')
        
        snap = 41
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0)& (G_history[snap].fcool>0))[0]
        #        if(len(w1) > dilute): w1 = sample(w1, dilute)
        
        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        mi_5 = np.floor(min(Lradio1400_5)) - 3.0
        ma_5 = np.floor(max(Lradio1400_5)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        (counts_5, binedges_5) = np.histogram(Lradio1400_5, range=(mi_5, ma_5), bins=NB_5 )
        # Set the x-axis values to be the centre of the bins
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'r:',lw = 4, alpha=0.95, label='Jet model (z = 1.0)')
        
        
        
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
            
            
            
        xplot= np.log10(10**(Best2012[:,0]+0.15))
        yplot = (10**Best2012[:,4])
        yerr2 = 10**(Best2012[:,4]+Best2012[:,5]) -  yplot
        yerr1 = yplot - 10**(Best2012[:,4]-Best2012[:,6])
#        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='r', lw=2.0, alpha=0.6, marker='o', markersize=3, ls='none', label='Best 2012 (All-Radio)', mew=1)
        plt.fill_between(xplot, yplot+yerr2, yplot-yerr1,facecolor='red', alpha=0.35,label= 'Best \& Heckman (2012)(z $\le$ 0.1)')
        plt.plot(xplot, yplot, color='red',lw = 3, alpha=0.35, label='Best \& Heckman (2012)(z $\le$ 0.1)')
        #             Lradio, logp_all, Uy_all, Ly_all, logp_AGN, Uy_AGN, Ly_AGN,
        Best2014_05 = np.array([
                             [  23.3, -4.69, 0.14, 0.22],
                             [  23.9, -5.23, 0.20, 0.39],
                             [  24.5, -5.16, 0.13, 0.18],
                             [  25.1, -5.81, 0.15, 0.23],
                             [  25.7, -6.66, 0.24, 0.58],
                             [  26.3, -7.67, 0.26, 0.70],
                             [  26.9, -8.41, 0.26, 0.80],
                             ], dtype=np.float32)
            
        xplot= Best2014_05[:,0]+0.3
        yplot = (10**Best2014_05[:,1])
        yerr2 = 10**(Best2014_05[:,1]+Best2014_05[:,2]) -  yplot
        yerr1 = yplot - 10**(Best2014_05[:,1]-Best2014_05[:,3])
#        plt.errorbar(xplot, yplot, yerr=0, color='green', lw=4.0, alpha=0.6, marker='o', markersize=1, ls='none', label='Best 2014 (z= 0.5-0.7)', mew=1)
        plt.fill_between(xplot, yplot+yerr2, yplot-yerr1,facecolor='green', alpha=0.15)
        plt.plot(xplot, yplot, color='green',lw = 4, alpha=0.15, label='Best et al. (2014) (z = 0.5-0.7)')
        Best2014_1 = np.array([
                        [  23.9, -6.23, 0.26, 0.70],
                        [  24.5, -5.95, 0.23, 0.54],
                        [  25.1, -5.96, 0.17, 0.27],
                        [  25.7, -6.47, 0.20, 0.37],
                        [  26.3, -7.36, 0.23, 0.55],
                        [  26.9, -8.27, 0.26, 0.78],
                        ], dtype=np.float32)
    
        xplot= Best2014_1[:,0]+0.3
        yplot = (10**Best2014_1[:,1])
        yerr2 = 10**(Best2014_1[:,1]+Best2014_1[:,2]) -  yplot
        yerr1 = yplot - 10**(Best2014_1[:,1]-Best2014_1[:,3])
#        plt.errorbar(xplot, yplot, yerr=0, color='green', lw=4.0, alpha=0.6, marker='o', markersize=1, ls='none', label='Best 2014 (z= 0.5-0.7)', mew=1)
        plt.fill_between(xplot, yplot+yerr2, yplot-yerr1,facecolor='blue', alpha=0.15,label= 'Best \& Heckman (2012)-z=1')
        plt.plot(xplot, yplot, color='blue',lw = 4, alpha=0.15, label='Best et al. (2014) (z = 0.7-1.0)')


        plt.yscale('log', nonposy='clip')
        plt.axis([22, 27.1, 1.0e-8, 1.0e-4])
                             
        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
            
        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # Set the y...
        
        leg = plt.legend(loc='lower left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
                              
        outputFile = OutputDir + 'D_RadioLF_allz' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
    
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def Lradio_Qjet(self, G_history):
    
        print 'Plotting the Radio Luminosity -- Qjet relation'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure


        snap = 63
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)&(G_history[snap].Qjet > 1e0)& (G_history[snap].Qjet < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0)& (G_history[snap].fcool>0))[0]
        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        Q_jet = np.log10 (G_history[snap].Qjet[w1])
#        plt.scatter(Lradio1400_5, Q_jet, marker='o',s=10, color='grey')
        total_bins = 30
        X = Lradio1400_5
        Y = Q_jet
        bins = np.linspace(X.min(),X.max(), total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b--',lw=2,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='c', lw=2.0, marker='s', markersize=3, ls='none', alpha=0.6,label='Jet model (z = 0)', mew=1)

        
        snap = 48
#        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0))[0]
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)&(G_history[snap].Qjet > 1e0)& (G_history[snap].Qjet < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0))[0]
        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        Q_jet = np.log10 (G_history[snap].Qjet[w1])
#        plt.scatter(Lradio1400_5, Q_jet, marker='o',s=10, color='m', alpha=0.05)
        total_bins = 30
        X = Lradio1400_5
        Y = Q_jet
        bins = np.linspace(X.min(),X.max(), total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'m--',lw=2,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='m', lw=2.0, alpha=0.6, marker='s', markersize=3, ls='none', label='Jet model (z = 0.5)', mew=1)
        
        snap = 41
#        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0))[0]
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)&(G_history[snap].Qjet > 1e0)& (G_history[snap].Qjet < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0))[0]
        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        Q_jet = np.log10 (G_history[snap].Qjet[w1])
#        plt.scatter(Lradio1400_5, Q_jet, marker='o',s=10, color='r', alpha=0.05)
        total_bins = 30
        X = Lradio1400_5
        Y = Q_jet
        bins = np.linspace(X.min(),X.max(), total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'r:',lw=2,alpha=.8)
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='r', lw=2.0, alpha=0.6, marker='s', markersize=3, ls='none', label='Jet model (z = 1.0)', mew=1)
                    
        
        plt.ylabel(r'$Log (Q_jet) [W]$')  # Set the y...
        plt.xlabel(r'$Log (L_{1.4 GHz} Radio) [W/Hz]$')  # and the x-axis labels
                                
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([21., 26.0, 34, 38.0])
        #        plt.axis([14.0, 25.0, 25, 37.0])
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
                                
        outputFile = OutputDir + 'E_Lradio_Qjet_allz' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
    
        # Add this plot to our output list
        OutputList.append(outputFile)

# ---------------------------------------------------------
    
    def Lradio_AGE(self, G_history):
        
        print 'Plotting the Radio Luminosity -- AGE relation'
        
        seed(2222)
        
        plt.figure(figsize=(16,6))  # New figure
        plt.figure()
#        ax = plt.subplot(131)  # 1 plot on the figure
        ax = plt.subplot(111)
        w = np.where ((G_history[63].StellarMass > 0) & (np.log10(G_history[63].CentralMvir*1e10) > 11 ))[0]

#        f = open('galaxies_mMS_all.txt', 'w+')
#        f.write('# GalaxyIndex  CentralGalaxyIndex  StellarMass RadioLuminosity BlackHoleMass Temp_Gas Lx_bol Cooling Type'+"\n")
#
#        for i in xrange(0,len(w)-1):
#            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(G_history[63].GalaxyIndex[w][i].astype(int), G_history[63].CentralGalaxyIndex[w][i].astype(int),G_history[63].StellarMass[w][i], G_history[63].RadioLuminosity[w,5][i], G_history[63].BlackHoleMass[w][i], G_history[63].Temp_Gas[w][i], G_history[63].Lx_bol[w][i], G_history[63].Cooling[w][i], G_history[63].Type[w][i].astype(int)))
        DATA = np.column_stack((G_history[63].GalaxyIndex[w].astype(int), G_history[63].CentralGalaxyIndex[w].astype(int),G_history[63].Pos[w,0],G_history[63].Pos[w,1],G_history[63].Pos[w,2], G_history[63].StellarMass[w], G_history[63].RadioLuminosity[w,5], G_history[63].RadioLuminosity[w,6], G_history[63].BlackHoleMass[w], G_history[63].Temp_Gas[w], G_history[63].Lx_bol[w], G_history[63].Cooling[w], G_history[63].Type[w].astype(int),G_history[63].Vel[w,0],G_history[63].Vel[w,1],G_history[63].Vel[w,2],G_history[63].Rvir[w],G_history[63].Vvir[w],G_history[63].VelDisp[w],G_history[63].TimeSinceMajorMerger[w],G_history[63].concentration[w]))
        np.savetxt('Galaxies_MS_all.txt', DATA, delimiter="  ",header='GalaxyIndex  CentralGalaxyIndex   x  y  z  StellarMass RadioLuminosity5 RadioLuminosity6 BlackHoleMass Temp_Gas Lx_bol Cooling Type  Velx Vely Velz Rvir Vvir VelDisp TimeSinceMajorMerger concentration')
       
#       find the central galaxies or halo index at z = 0
        central_63, inv_63, inx_63= np.unique(G_history[63].CentralGalaxyIndex,return_inverse=True,return_index=True)
        
        Mass_halo_63 = np.zeros(len(central_63))
#        stellarmass_63 = np.zeros(len(central_63))
#        Lradio_63 = np.zeros(len(central_63))
#        Mass_BH = np.zeros(len(central_63))
#        Galaxytype_63    = np.zeros(len(central_63))
#        temp_63 =np.zeros(len(central_63))

#        preparing for trace into the redshift
#        for snap in xrange(41,63,23):

        snap = 41   # redshift ~ 1
    
    
#       find the central galaxies or halo index at z = 1
        central_snap, inv_snap, inx_snap= np.unique(G_history[snap].CentralGalaxyIndex,return_inverse=True,return_index=True)
        
        Mass_halo_snap = np.zeros(len(central_snap))
        
        f = open('Group_age_MS_all.txt', 'w+')
        f.write('# CentralGalaxyIndex  CentralMvir age'+"\n")
        
        for i in xrange(0,len(central_63)-1):
            Mass_halo_63[i] = G_history[63].CentralMvir[inv_63[i]]
            
            if ( (Mass_halo_63[i] > 1e3)):
               for j in xrange(0,len(central_snap)-1):
                    Mass_halo_snap[j] = G_history[snap].CentralMvir[inv_snap[j]]
                    
                    if ((central_snap[j] == central_63[i])):
                            age_alfa = Mass_halo_snap[j]/Mass_halo_63[i]
                            Central_id = central_63[i]
                            M_200 = Mass_halo_63[i] * 1e10 /self.Hubble_h
                            f.write("{0}\t{1}\t{2}\n".format(Central_id, M_200, age_alfa))

                            if (age_alfa > 0.5 ):
                               Central_id_old = central_63[i]
                               M_200_old = Mass_halo_63[i] * 1e10 /self.Hubble_h
                               print Central_id_old, age_alfa
                            elif ((age_alfa < 0.3)):
                               Central_id_young = central_63[i]
                               M_200_young = Mass_halo_63[i] * 1e10 /self.Hubble_h
                               print Central_id_young, L_radio_young, xray_young,temp_young,M_BH_young, age_alfa


        M_200 = [1,1,1,1]
        age_alfa = [1,1,1,1]
        plt.scatter(M_200, age_alfa, marker='o', s=1, color='blue', alpha=0.15, label='all')

        plt.ylabel(r'Log $L_{1.4 GHz}$ [W/Hz]')  # Set the y...
        plt.xlabel(r'Log $M_{*}$')  # and the x-axis labels
        
        # Set the x and y axis minor ticks
#        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
#        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        plt.axis([10.7, 11.7, 22, 26.6])
        

        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')



        outputFile = OutputDir + 'F_Lradio_AGE_mMS' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
#        f.close()
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def PlotHistory_Radiodensity(self, G_history):
    
        print 'Plotting Radio density evolution for all galaxies'
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        Radio_density = np.zeros((LastSnap+1-FirstSnap))
        
        for snap in xrange(FirstSnap,LastSnap+1):
            w = np.where(G_history[snap].RadioLuminosity[:,5] < 1e40)[0]
            Radio_density[snap-FirstSnap] = sum(G_history[snap].RadioLuminosity[w,5]) / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h
                                  
        z = np.array(self.redshift)
        nonzero = np.where(Radio_density > 0.0)[0]
        plt.plot(z[nonzero], np.log10(Radio_density[nonzero]), c= 'b', lw=3.0, label='Jet-model ')

        plt.ylabel(r'$\log_{10} \mathrm{Radio Luminosity\ density}\ (M_{\odot}\ \mathrm{yr}^{-1}\ \mathrm{Mpc}^{-3})$')  # Set the y...
        plt.xlabel(r'$\mathrm{redshift}$')  # and the x-axis labels
        
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
        
        plt.axis([0.0, 8.0, 19, 22])
        
        leg = plt.legend(loc='lower left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + 'G_History-Radio-density' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)
    # ---------------------------------------------------------
    
    def QjetF(self, G_history):
        
        print 'Plotting the Qjet function'
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        binwidth = 0.6  # Radio Luminosity function histogram bin width
        
        snap = 63
        w1 = np.where((G_history[snap].Qjet > 1e0)& (G_history[snap].Qjet < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0)& (G_history[snap].fcool>0))[0]
#        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        Q_jet = np.log10 (G_history[snap].Qjet[w1])
        mi_5 = np.floor(min(Q_jet)) - 3.0
        ma_5 = np.floor(max(Q_jet)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        (counts_5, binedges_5) = np.histogram(Q_jet, range=(mi_5, ma_5), bins=NB_5 )
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        #        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 2, alpha=0.55, label='Jet model, 63')
        
        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5,'c', lw=3.0, alpha=0.8, marker='s', markersize=10, label='Jet model (z = 0)')
        plt.plot(xaxeshisto_5, counts_5  / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'b--',lw=3,alpha=.8)
        
        
        snap = 48
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0))[0]
        #        if(len(w1) > dilute): w1 = sample(w1, dilute)
        
#        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        Q_jet = np.log10 (G_history[snap].Qjet[w1])
        mi_5 = np.floor(min(Q_jet)) - 3.0
        ma_5 = np.floor(max(Q_jet)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        (counts_5, binedges_5) = np.histogram(Q_jet, range=(mi_5, ma_5), bins=NB_5 )
        # Set the x-axis values to be the centre of the bins
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'm--',lw = 4, alpha=0.95, label='Jet model (z = 0.5)')
        
        snap = 41
        w1 = np.where((G_history[snap].RadioLuminosity[:,5] > 0)& (G_history[snap].RadioLuminosity[:,5] < 1e50)&(np.log10(G_history[snap].CentralMvir * 1e10/self.Hubble_h) > 0))[0]
        #        if(len(w1) > dilute): w1 = sample(w1, dilute)
        
#        Lradio1400_5 = np.log10(G_history[snap].RadioLuminosity[w1,5])
        Q_jet = np.log10 (G_history[snap].Qjet[w1])
        mi_5 = np.floor(min(Q_jet)) - 3.0
        ma_5 = np.floor(max(Q_jet)) + 3.0
        NB_5 = (ma_5 - mi_5) / binwidth
        (counts_5, binedges_5) = np.histogram(Q_jet, range=(mi_5, ma_5), bins=NB_5 )
        # Set the x-axis values to be the centre of the bins
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto_5, counts_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'r:',lw = 4, alpha=0.95, label='Jet model (z = 1.0)')
        
        
        plt.yscale('log', nonposy='clip')
        plt.axis([30, 38.5, 2.0e-7, 2.0e-4])
        
        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
        
        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} (Q_{\mathrm{jet}}\ [W])$')  # Set the y...
        
        leg = plt.legend(loc='lower left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + 'H_QjetF_allz' + OutputFormat
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
        default='model',
        help='filename base (default: model)',
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
    parser.add_option(
        '-s',
        '--snap_range',
        type='int',
        nargs=2,
        dest='SnapRange',
        default=(0, 63),
        help='first and last snapshots (default: 0 63)',
        metavar='FIRST LAST',
        )


    (opt, args) = parser.parse_args()

    if opt.DirName[-1] != '/':
        opt.DirName += '/'

#    OutputDir = opt.DirName + '/plots/'

    if not os.path.exists(OutputDir):
      os.makedirs(OutputDir)

    res = Results()

    print 'Running history...'

    FirstFile = opt.FileRange[0]
    LastFile = opt.FileRange[1]
    FirstSnap = opt.SnapRange[0]
    LastSnap = opt.SnapRange[1]

    # read in all files and put in G_history
    G_history = [0]*(LastSnap-FirstSnap+1)
    for snap in xrange(FirstSnap,LastSnap+1):

      print
      print 'SNAPSHOT NUMBER:  ', snap
      
      fin_base = opt.DirName + opt.FileName + res.redshift_file[snap]
      G_history[snap] = res.read_gals(fin_base, FirstFile, LastFile, snap)
      
    print

    res.StellarMassFunction(G_history)
    res.PlotHistory_SFRdensity(G_history)
    res.StellarMassDensityEvolution(G_history)
    res.RadioLF(G_history)
    res.Lradio_Qjet(G_history)
#    res.Lradio_AGE(G_history)
#    res.PlotHistory_Radiodensity(G_history)
#    res.QjetF(G_history)


