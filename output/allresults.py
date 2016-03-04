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
            facecolor='purple', alpha=0.25, label='Baldry et al. 2008 (z=0.1)')

        # This next line is just to get the shaded region to appear correctly in the legend
        plt.plot(xaxeshisto, counts / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, label='Baldry et al. 2008', color='purple', alpha=0.3)

        # # Cole et al. 2001 SMF (h=1.0 converted to h=0.73)
        # M = np.arange(7.0, 13.0, 0.01)
        # Mstar = np.log10(7.07*1.0e10 /self.Hubble_h/self.Hubble_h)
        # alpha = -1.18
        # phistar = 0.009 *self.Hubble_h*self.Hubble_h*self.Hubble_h
        # xval = 10.0 ** (M-Mstar)
        # yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)      
        # plt.plot(M, yval, 'g--', lw=1.5, label='Cole et al. 2001')  # Plot the SMF
        
        # Overplot the model histograms
        plt.plot(xaxeshisto, counts    / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'k-', label='Model - All')
        plt.plot(xaxeshisto, countsRED / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'r:', lw=2, label='Model - Red')
        plt.plot(xaxeshisto, countsBLU / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'b:', lw=2, label='Model - Blue')

        plt.yscale('log', nonposy='clip')
        plt.axis([8.0, 12.5, 1.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$')  # and the x-axis labels

        plt.text(12.2, 0.03, whichsimulation, size = 'large')

        leg = plt.legend(loc='lower left', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = OutputDir + '1.StellarMassFunction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def BaryonicMassFunction(self, G):

        print 'Plotting the baryonic mass function'

        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        binwidth = 0.15  # mass function histogram bin width
      
        # calculate BMF
        w = np.where(G.StellarMass + G.ColdGas > 0.0)[0]
        mass = np.log10((G.StellarMass[w] + G.ColdGas[w]) * 1.0e10 / self.Hubble_h)

        mi = np.floor(min(mass)) - 2
        ma = np.floor(max(mass)) + 2
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
       
        # Bell et al. 2003 BMF (h=1.0 converted to h=0.73)
        M = np.arange(7.0, 13.0, 0.01)
        Mstar = np.log10(5.3*1.0e10 /self.Hubble_h/self.Hubble_h)
        alpha = -1.21
        phistar = 0.0108 *self.Hubble_h*self.Hubble_h*self.Hubble_h
        xval = 10.0 ** (M-Mstar)
        yval = np.log(10.) * phistar * xval ** (alpha+1) * np.exp(-xval)
        
        if(whichimf == 0):
            # converted diet Salpeter IMF to Salpeter IMF
            plt.plot(np.log10(10.0**M /0.7), yval, 'b-', lw=2.0, label='Bell et al. 2003')  # Plot the SMF
        elif(whichimf == 1):
            # converted diet Salpeter IMF to Salpeter IMF, then to Chabrier IMF
            plt.plot(np.log10(10.0**M /0.7 /1.8), yval, 'g--', lw=1.5, label='Bell et al. 2003')  # Plot the SMF

        # Overplot the model histograms
        plt.plot(xaxeshisto, counts / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'k-', label='Model')

        plt.yscale('log', nonposy='clip')
        plt.axis([8.0, 12.5, 1.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10}\ M_{\mathrm{bar}}\ (M_{\odot})$')  # and the x-axis labels

        leg = plt.legend(loc='lower left', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = OutputDir + '2.BaryonicMassFunction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------
   
    def GasMassFunction(self, G):

        print 'Plotting the cold gas mass function'

        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        binwidth = 0.1  # mass function histogram bin width

        # calculate all
        w = np.where(G.ColdGas > 0.0)[0]
        mass = np.log10(G.ColdGas[w] * 1.0e10 / self.Hubble_h)
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
            [10.492,  -5.083]        ], dtype=np.float32)
        
        ObrRaw = np.array([
            [7.300,   -1.104],
            [7.576,   -1.302],
            [7.847,   -1.250],
            [8.133,   -1.240],
            [8.409,   -1.344],
            [8.691,   -1.479],
            [8.956,   -1.792],
            [9.231,   -2.271],
            [9.507,   -3.198],
            [9.788,   -5.062 ]        ], dtype=np.float32)

        ObrCold = np.array([
            [8.009,   -1.042],
            [8.215,   -1.156],
            [8.409,   -0.990],
            [8.604,   -1.156],
            [8.799,   -1.208],
            [9.020,   -1.333],
            [9.194,   -1.385],
            [9.404,   -1.552],
            [9.599,   -1.677],
            [9.788,   -1.812],
            [9.999,   -2.312],
            [10.172,  -2.656],
            [10.362,  -3.500],
            [10.551,  -3.635],
            [10.740,  -5.010]        ], dtype=np.float32)

        ObrCold_xval = np.log10(10**(ObrCold[:, 0])  /self.Hubble_h/self.Hubble_h)
        ObrCold_yval = (10**(ObrCold[:, 1]) * self.Hubble_h*self.Hubble_h*self.Hubble_h)
        Zwaan_xval = np.log10(10**(Zwaan[:, 0]) /self.Hubble_h/self.Hubble_h)
        Zwaan_yval = (10**(Zwaan[:, 1]) * self.Hubble_h*self.Hubble_h*self.Hubble_h)
        ObrRaw_xval = np.log10(10**(ObrRaw[:, 0])  /self.Hubble_h/self.Hubble_h)
        ObrRaw_yval = (10**(ObrRaw[:, 1]) * self.Hubble_h*self.Hubble_h*self.Hubble_h)

        plt.plot(ObrCold_xval, ObrCold_yval, color='black', lw = 7, alpha=0.25, label='Obr. \& Raw. 2009 (Cold Gas)')
        plt.plot(Zwaan_xval, Zwaan_yval, color='cyan', lw = 7, alpha=0.25, label='Zwaan et al. 2005 (HI)')
        plt.plot(ObrRaw_xval, ObrRaw_yval, color='magenta', lw = 7, alpha=0.25, label='Obr. \& Raw. 2009 (H2)')

        
        # Overplot the model histograms
        plt.plot(xaxeshisto, counts    / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'k-', label='Model - Cold Gas')

        plt.yscale('log', nonposy='clip')
        plt.axis([8.0, 11.5, 1.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{X}}\ (M_{\odot})$')  # and the x-axis labels

        leg = plt.legend(loc='lower left', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = OutputDir + '3.GasMassFunction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------
    
    def BaryonicTullyFisher(self, G):
    
        print 'Plotting the baryonic TF relationship'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
    
        # w = np.where((G.Type == 0) & (G.StellarMass + G.ColdGas > 0.0) & (G.Vmax > 0.0))[0]
        w = np.where((G.Type == 0) & (G.StellarMass + G.ColdGas > 0.0) & 
          (G.BulgeMass / G.StellarMass > 0.1) & (G.BulgeMass / G.StellarMass < 0.5))[0]
        if(len(w) > dilute): w = sample(w, dilute)
    
        mass = np.log10((G.StellarMass[w] + G.ColdGas[w]) * 1.0e10 / self.Hubble_h)
        vel = np.log10(G.Vmax[w])
                    
        plt.scatter(vel, mass, marker='o', s=1, c='k', alpha=0.5, label='Model Sb/c galaxies')
                
        # overplot Stark, McGaugh & Swatters 2009 (assumes h=0.75? ... what IMF?)
        w = np.arange(0.5, 10.0, 0.5)
        TF = 3.94*w + 1.79
        plt.plot(w, TF, 'b-', lw=2.0, label='Stark, McGaugh \& Swatters 2009')
            
        plt.ylabel(r'$\log_{10}\ M_{\mathrm{bar}}\ (M_{\odot})$')  # Set the y...
        plt.xlabel(r'$\log_{10}V_{max}\ (km/s)$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([1.4, 2.6, 8.0, 12.0])
            
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '4.BaryonicTullyFisher' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------
    
    def SpecificStarFormationRate(self, G):
    
        print 'Plotting the specific SFR'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        w = np.where(G.StellarMass > 0.01)[0]
        if(len(w) > dilute): w = sample(w, dilute)
        
        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        sSFR = np.log10( (G.SfrDisk[w] + G.SfrBulge[w]) / (G.StellarMass[w] * 1.0e10 / self.Hubble_h) )
        plt.scatter(mass, sSFR, marker='o', s=1, c='k', alpha=0.5, label='Model galaxies')
                
        # overplot dividing line between SF and passive
        w = np.arange(7.0, 13.0, 1.0)
        plt.plot(w, w/w*sSFRcut, 'b:', lw=2.0)
            
        plt.ylabel(r'$\log_{10}\ s\mathrm{SFR}\ (\mathrm{yr^{-1}})$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([8.0, 12.0, -16.0, -8.0])
            
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '5.SpecificStarFormationRate' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def GasFraction(self, G):
    
        print 'Plotting the gas fractions'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        w = np.where((G.Type == 0) & (G.StellarMass + G.ColdGas > 0.0) & 
          (G.BulgeMass / G.StellarMass > 0.1) & (G.BulgeMass / G.StellarMass < 0.5))[0]
        if(len(w) > dilute): w = sample(w, dilute)
        
        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        fraction = G.ColdGas[w] / (G.StellarMass[w] + G.ColdGas[w])
                    
        plt.scatter(mass, fraction, marker='o', s=1, c='k', alpha=0.5, label='Model Sb/c galaxies')
            
        plt.ylabel(r'$\mathrm{Cold\ Mass\ /\ (Cold+Stellar\ Mass)}$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([8.0, 12.0, 0.0, 1.0])
            
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '6.GasFraction' + OutputFormat
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
        if(len(w) > dilute): w = sample(w, dilute)
        
        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        Z = np.log10((G.MetalsColdGas[w] / G.ColdGas[w]) / 0.02) + 9.0
                    
        plt.scatter(mass, Z, marker='o', s=1, c='k', alpha=0.5, label='Model galaxies')
            
        # overplot Tremonti et al. 2003 (h=0.7)
        w = np.arange(7.0, 13.0, 0.1)
        Zobs = -1.492 + 1.847*w - 0.08026*w*w
        if(whichimf == 0):
            # Conversion from Kroupa IMF to Slapeter IMF
            plt.plot(np.log10((10**w *1.5)), Zobs, 'b-', lw=2.0, label='Tremonti et al. 2003')
        elif(whichimf == 1):
            # Conversion from Kroupa IMF to Slapeter IMF to Chabrier IMF
            plt.plot(np.log10((10**w *1.5 /1.8)), Zobs, 'b-', lw=2.0, label='Tremonti et al. 2003')
            
        plt.ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([8.0, 12.0, 8.0, 9.5])
            
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '7.Metallicity' + OutputFormat
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
                    
        plt.scatter(bulge, bh, marker='o', s=1, c='k', alpha=0.5, label='Model galaxies')
                
        # overplot Haring & Rix 2004
        w = 10. ** np.arange(20)
        BHdata = 10. ** (8.2 + 1.12 * np.log10(w / 1.0e11))
        plt.plot(np.log10(w), np.log10(BHdata), 'b-', label="Haring \& Rix 2004")

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
            
        outputFile = OutputDir + '8.BlackHoleBulgeRelationship' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------
    
    def QuiescentFraction(self, G):
    
        print 'Plotting the quiescent fraction vs stellar mass'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        groupscale = 12.5
        
        w = np.where(G.StellarMass > 0.0)[0]
        StellarMass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        CentralMvir = np.log10(G.CentralMvir[w] * 1.0e10 / self.Hubble_h)
        Type = G.Type[w]
        sSFR = (G.SfrDisk[w] + G.SfrBulge[w]) / (G.StellarMass[w] * 1.0e10 / self.Hubble_h)

        MinRange = 9.5
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

        for i in xrange(Nbins-1):
            
            w = np.where((StellarMass >= Range[i]) & (StellarMass < Range[i+1]))[0]
            if len(w) > 0:
                wQ = np.where((StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (sSFR < 10.0**sSFRcut))[0]
                Fraction.append(1.0*len(wQ) / len(w))
            else:
                Fraction.append(0.0)

            w = np.where((Type == 0) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]))[0]
            if len(w) > 0:
                wQ = np.where((Type == 0) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (sSFR < 10.0**sSFRcut))[0]
                CentralFraction.append(1.0*len(wQ) / len(w))
            else:
                CentralFraction.append(0.0)

            w = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]))[0]
            if len(w) > 0:
                wQ = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (sSFR < 10.0**sSFRcut))[0]
                SatelliteFraction.append(1.0*len(wQ) / len(w))
                wQ = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (sSFR < 10.0**sSFRcut) & (CentralMvir < groupscale))[0]
                SatelliteFractionLo.append(1.0*len(wQ) / len(w))
                wQ = np.where((Type == 1) & (StellarMass >= Range[i]) & (StellarMass < Range[i+1]) & (sSFR < 10.0**sSFRcut) & (CentralMvir > groupscale))[0]
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
        
        w = np.where(Fraction > 0)[0]
        plt.plot(Mass[w], Fraction[w], label='All')

        w = np.where(CentralFraction > 0)[0]
        plt.plot(Mass[w], CentralFraction[w], color='Blue', label='Centrals')

        w = np.where(SatelliteFraction > 0)[0]
        plt.plot(Mass[w], SatelliteFraction[w], color='Red', label='Satellites')

        w = np.where(SatelliteFractionLo > 0)[0]
        plt.plot(Mass[w], SatelliteFractionLo[w], 'r--', label='Satellites-Lo')

        w = np.where(SatelliteFractionHi > 0)[0]
        plt.plot(Mass[w], SatelliteFractionHi[w], 'r-.', label='Satellites-Hi')
        
        plt.xlabel(r'$\log_{10} M_{\mathrm{stellar}}\ (M_{\odot})$')  # Set the x-axis label
        plt.ylabel(r'$\mathrm{Quescient\ Fraction}$')  # Set the y-axis label
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([9.5, 12.0, 0.0, 1.05])
            
        leg = plt.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '9.QuiescentFraction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def BulgeMassFraction(self, G):
    
        print 'Plotting the mass fraction of galaxies'
    
        seed(2222)

        fBulge = G.BulgeMass / G.StellarMass
        fDisk = 1.0 - (G.BulgeMass) / G.StellarMass
        mass = np.log10(G.StellarMass * 1.0e10 / self.Hubble_h)
        sSFR = np.log10((G.SfrDisk + G.SfrBulge) / (G.StellarMass * 1.0e10 / self.Hubble_h))
        
        binwidth = 0.2
        shift = binwidth/2.0
        mass_range = np.arange(8.5-shift, 12.0+shift, binwidth)
        bins = len(mass_range)
        
        fBulge_ave = np.zeros(bins)
        fBulge_var = np.zeros(bins)
        fDisk_ave = np.zeros(bins)
        fDisk_var = np.zeros(bins)
        
        for i in xrange(bins-1):
            w = np.where( (mass >= mass_range[i]) & (mass < mass_range[i+1]))[0]
            # w = np.where( (mass >= mass_range[i]) & (mass < mass_range[i+1]) & (sSFR < sSFRcut))[0]
            if(len(w) > 0):
                fBulge_ave[i] = np.mean(fBulge[w])
                fBulge_var[i] = np.var(fBulge[w])
                fDisk_ave[i] = np.mean(fDisk[w])
                fDisk_var[i] = np.var(fDisk[w])

        w = np.where(fBulge_ave > 0.0)[0]
        plt.plot(mass_range[w]+shift, fBulge_ave[w], 'r-', label='bulge')
        plt.fill_between(mass_range[w]+shift, 
            fBulge_ave[w]+fBulge_var[w], 
            fBulge_ave[w]-fBulge_var[w], 
            facecolor='red', alpha=0.25)

        w = np.where(fDisk_ave > 0.0)[0]
        plt.plot(mass_range[w]+shift, fDisk_ave[w], 'k-', label='disk stars')
        plt.fill_between(mass_range[w]+shift, 
            fDisk_ave[w]+fDisk_var[w], 
            fDisk_ave[w]-fDisk_var[w], 
            facecolor='black', alpha=0.25)

        plt.axis([mass_range[0], mass_range[bins-1], 0.0, 1.05])

        plt.ylabel(r'$\mathrm{Stellar\ Mass\ Fraction}$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$')  # and the x-axis labels

        leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')

        outputFile = OutputDir + '10.BulgeMassFraction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------
    
    def BaryonFraction(self, G):
    
        print 'Plotting the average baryon fraction vs halo mass'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        HaloMass = np.log10(G.Mvir * 1.0e10 / self.Hubble_h)
        Baryons = G.StellarMass + G.ColdGas + G.HotGas + G.EjectedMass + G.IntraClusterStars + G.BlackHoleMass
        fileNr = np.floor(G.GalaxyIndex / 1.0e12)
        
        MinHalo = 11.0
        MaxHalo = 16.0
        Interval = 0.1
        Nbins = int((MaxHalo-MinHalo)/Interval)
        HaloRange = np.arange(MinHalo, MaxHalo, Interval)
        
        MeanCentralHaloMass = []
        MeanBaryonFraction = []
        MeanBaryonFractionU = []
        MeanBaryonFractionL = []

        MeanStars = []
        MeanCold = []
        MeanHot = []
        MeanEjected = []
        MeanICS = []
        MeanBH = []

        for i in xrange(Nbins-1):
            
            w1 = np.where((G.Type == 0) & (HaloMass >= HaloRange[i]) & (HaloMass < HaloRange[i+1]))[0]
            HalosFound = len(w1)
            
            if HalosFound > 2:  
                
                BaryonFraction = []
                CentralHaloMass = []
                
                Stars = []
                Cold = []
                Hot = []
                Ejected = []
                ICS = []
                BH = []
                
                for j in xrange(HalosFound):
                    
                    w2 = np.where((G.SAGEHaloIndex == G.SAGEHaloIndex[w1[j]]) & (G.SAGETreeIndex == G.SAGETreeIndex[w1[j]]) & (fileNr == fileNr[w1[j]]))[0]
                    CentralAndSatellitesFound = len(w2)
                    
                    if CentralAndSatellitesFound > 0:
                        BaryonFraction.append(sum(Baryons[w2]) / G.Mvir[w1[j]])
                        CentralHaloMass.append(np.log10(G.Mvir[w1[j]] * 1.0e10 / self.Hubble_h))

                        Stars.append(sum(G.StellarMass[w2]) / G.Mvir[w1[j]])
                        Cold.append(sum(G.ColdGas[w2]) / G.Mvir[w1[j]])
                        Hot.append(sum(G.HotGas[w2]) / G.Mvir[w1[j]])
                        Ejected.append(sum(G.EjectedMass[w2]) / G.Mvir[w1[j]])
                        ICS.append(sum(G.IntraClusterStars[w2]) / G.Mvir[w1[j]])
                        BH.append(sum(G.BlackHoleMass[w2]) / G.Mvir[w1[j]])                        
                                
                MeanCentralHaloMass.append(np.mean(CentralHaloMass))
                MeanBaryonFraction.append(np.mean(BaryonFraction))
                MeanBaryonFractionU.append(np.mean(BaryonFraction) + np.var(BaryonFraction))
                MeanBaryonFractionL.append(np.mean(BaryonFraction) - np.var(BaryonFraction))
                
                MeanStars.append(np.mean(Stars))
                MeanCold.append(np.mean(Cold))
                MeanHot.append(np.mean(Hot))
                MeanEjected.append(np.mean(Ejected))
                MeanICS.append(np.mean(ICS))
                MeanBH.append(np.mean(BH))
                
                print '  ', i, HaloRange[i], HalosFound, np.mean(BaryonFraction)
        
        plt.plot(MeanCentralHaloMass, MeanBaryonFraction, 'k-', label='TOTAL')#, color='purple', alpha=0.3)
        plt.fill_between(MeanCentralHaloMass, MeanBaryonFractionU, MeanBaryonFractionL, 
            facecolor='purple', alpha=0.25, label='TOTAL')
        
        plt.plot(MeanCentralHaloMass, MeanStars, 'k--', label='Stars')
        plt.plot(MeanCentralHaloMass, MeanCold, label='Cold', color='blue')
        plt.plot(MeanCentralHaloMass, MeanHot, label='Hot', color='red')
        plt.plot(MeanCentralHaloMass, MeanEjected, label='Ejected', color='green')
        plt.plot(MeanCentralHaloMass, MeanICS, label='ICS', color='yellow')
        # plt.plot(MeanCentralHaloMass, MeanBH, 'k:', label='BH')
        
        plt.xlabel(r'$\mathrm{Central}\ \log_{10} M_{\mathrm{vir}}\ (M_{\odot})$')  # Set the x-axis label
        plt.ylabel(r'$\mathrm{Baryon\ Fraction}$')  # Set the y-axis label
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([10.8, 15.0, 0.0, 0.23])
            
        leg = plt.legend(bbox_to_anchor=[0.99, 0.6])
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = OutputDir + '11.BaryonFraction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def SpinDistribution(self, G):
    
        print 'Plotting the spin distribution of all galaxies'

        # set up figure
        plt.figure()
        ax = plt.subplot(111)
    
        SpinParameter = np.sqrt(G.Spin[:,0]*G.Spin[:,0] + G.Spin[:,1]*G.Spin[:,1] + G.Spin[:,2]*G.Spin[:,2]) / (np.sqrt(2) * G.Vvir * G.Rvir);
        
        mi = -0.02
        ma = 0.5
        binwidth = 0.01
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(SpinParameter, range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto, counts, 'k-', label='simulation')

        plt.axis([mi, ma, 0.0, max(counts)*1.15])

        plt.ylabel(r'$\mathrm{Number}$')  # Set the y...
        plt.xlabel(r'$\mathrm{Spin\ Parameter}$')  # and the x-axis labels

        leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')

        outputFile = OutputDir + '12.SpinDistribution' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def VelocityDistribution(self, G):
    
        print 'Plotting the velocity distribution of all galaxies'
    
        seed(2222)
    
        mi = -40.0
        ma = 40.0
        binwidth = 0.5
        NB = (ma - mi) / binwidth

        # set up figure
        plt.figure()
        ax = plt.subplot(111)

        pos_x = G.Pos[:,0] / self.Hubble_h
        pos_y = G.Pos[:,1] / self.Hubble_h
        pos_z = G.Pos[:,2] / self.Hubble_h

        vel_x = G.Vel[:,0]
        vel_y = G.Vel[:,1]
        vel_z = G.Vel[:,2]

        dist_los = np.sqrt(pos_x*pos_x + pos_y*pos_y + pos_z*pos_z)
        vel_los = (pos_x/dist_los)*vel_x + (pos_y/dist_los)*vel_y + (pos_z/dist_los)*vel_z
        dist_red = dist_los + vel_los/(self.Hubble_h*100.0)

        tot_gals = len(pos_x)


        (counts, binedges) = np.histogram(vel_los/(self.Hubble_h*100.0), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto, counts / binwidth / tot_gals, 'k-', label='los-velocity')

        (counts, binedges) = np.histogram(vel_x/(self.Hubble_h*100.0), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto, counts / binwidth / tot_gals, 'r-', label='x-velocity')

        (counts, binedges) = np.histogram(vel_y/(self.Hubble_h*100.0), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto, counts / binwidth / tot_gals, 'g-', label='y-velocity')

        (counts, binedges) = np.histogram(vel_z/(self.Hubble_h*100.0), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        plt.plot(xaxeshisto, counts / binwidth / tot_gals, 'b-', label='z-velocity')


        plt.yscale('log', nonposy='clip')
        plt.axis([mi, ma, 1e-5, 0.5])
        # plt.axis([mi, ma, 0, 0.13])

        plt.ylabel(r'$\mathrm{Box\ Normalised\ Count}$')  # Set the y...
        plt.xlabel(r'$\mathrm{Velocity / H}_{0}$')  # and the x-axis labels

        leg = plt.legend(loc='upper left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')

        outputFile = OutputDir + '13.VelocityDistribution' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def MassReservoirScatter(self, G):
    
        print 'Plotting the mass in stellar, cold, hot, ejected, ICS reservoirs'
    
        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
    
        w = np.where((G.Type == 0) & (G.Mvir > 1.0) & (G.StellarMass > 0.0))[0]
        if(len(w) > dilute): w = sample(w, dilute)

        mvir = np.log10(G.Mvir[w] * 1.0e10)
        plt.scatter(mvir, np.log10(G.StellarMass[w] * 1.0e10), marker='o', s=0.3, c='k', alpha=0.5, label='Stars')
        plt.scatter(mvir, np.log10(G.ColdGas[w] * 1.0e10), marker='o', s=0.3, color='blue', alpha=0.5, label='Cold gas')
        plt.scatter(mvir, np.log10(G.HotGas[w] * 1.0e10), marker='o', s=0.3, color='red', alpha=0.5, label='Hot gas')
        plt.scatter(mvir, np.log10(G.EjectedMass[w] * 1.0e10), marker='o', s=0.3, color='green', alpha=0.5, label='Ejected gas')
        plt.scatter(mvir, np.log10(G.IntraClusterStars[w] * 1.0e10), marker='o', s=10, color='yellow', alpha=0.5, label='Intracluster stars')    

        plt.ylabel(r'$\mathrm{stellar,\ cold,\ hot,\ ejected,\ ICS\ mass}$')  # Set the y...
        plt.xlabel(r'$\log\ M_{\mathrm{vir}}\ (h^{-1}\ M_{\odot})$')  # and the x-axis labels
        
        plt.axis([10.0, 14.0, 7.5, 12.5])

        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        plt.text(13.5, 8.0, r'$\mathrm{All}')
            
        outputFile = OutputDir + '14.MassReservoirScatter' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
            
        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def SpatialDistribution(self, G):
    
        print 'Plotting the spatial distribution of all galaxies'
    
        seed(2222)
    
        plt.figure()  # New figure
    
        w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.1))[0]
        if(len(w) > dilute): w = sample(w, dilute)

        xx = G.Pos[w,0]
        yy = G.Pos[w,1]
        zz = G.Pos[w,2]

        buff = self.BoxSize*0.1

        ax = plt.subplot(221)  # 1 plot on the figure
        plt.scatter(xx, yy, marker='o', s=0.3, c='k', alpha=0.5)
        plt.axis([0.0-buff, self.BoxSize+buff, 0.0-buff, self.BoxSize+buff])

        plt.ylabel(r'$\mathrm{x}$')  # Set the y...
        plt.xlabel(r'$\mathrm{y}$')  # and the x-axis labels
        
        ax = plt.subplot(222)  # 1 plot on the figure
        plt.scatter(xx, zz, marker='o', s=0.3, c='k', alpha=0.5)
        plt.axis([0.0-buff, self.BoxSize+buff, 0.0-buff, self.BoxSize+buff])

        plt.ylabel(r'$\mathrm{x}$')  # Set the y...
        plt.xlabel(r'$\mathrm{z}$')  # and the x-axis labels
        
        ax = plt.subplot(223)  # 1 plot on the figure
        plt.scatter(yy, zz, marker='o', s=0.3, c='k', alpha=0.5)
        plt.axis([0.0-buff, self.BoxSize+buff, 0.0-buff, self.BoxSize+buff])
        plt.ylabel(r'$\mathrm{y}$')  # Set the y...
        plt.xlabel(r'$\mathrm{z}$')  # and the x-axis labels
            
        outputFile = OutputDir + '15.SpatialDistribution' + OutputFormat
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
        w = np.where((np.log10(G.Mvir * 1e10 / self.Hubble_h) >11.0))[0]
        if(len(w) > dilute): w = sample(w, dilute)
        w1 = np.where((G.Type == 0) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
        if(len(w1) > dilute): w1 = sample(w1, dilute)

        w2 = np.where((G.Type == 1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
        if(len(w2) > dilute): w2 = sample(w2, dilute)

        Lradio1400_5 = np.log10(G.RadioLuminosity[w,5])
        Q_jet_1      = np.log10(G.Qjet[w])

        # for p_value = 2.6
        Lradio1400_5_S = np.log10(G.RadioLuminosity[w2,5])
        Q_jet_S      = np.log10(G.Qjet[w2])
        Lradio1400_5_E = np.log10(G.RadioLuminosity[w1,5])
        Q_jet_E      = np.log10(G.Qjet[w1])
        #Shabala 2008
        Shabala2008 = np.array([
                               [22.087116069356547,  33.453373995710514],
                               [23.087366638467365,  34.423304206875166],
                               [24.07930737023842,   35.38708502615127],
                               [24.739908910551236,  36.02806883243122],
                               [25.39116203803964,   36.6613663549015]        ], dtype=np.float32)
                                
                                
        #Heckman & Best 2014
        Heckman2014_1 = np.array([
                                 [20.280693069306928, 33.049886621315196],
                                 [21.195544554455445, 33.83446712018141],
                                 [22.370297029702968, 34.84126984126984],
                                 [24.737128712871286, 36.87755102040816],
                                 [26.81980198019802, 38.668934240362816]        ], dtype=np.float32)
                                
        Heckman2014_2 = np.array([
                                 [20.01386138613861, 33.269085411942555],
                                 [21.275247524752473, 34.354497354497354],
                                 [22.728382838283828, 35.600151171579746],
                                 [24.768316831683165, 37.3567649281935],
                                 [26.547194719471946, 38.88662131519274]        ], dtype=np.float32)
                                
        # For fw = 10,20 convert to f= 5
        Heckman2014_xval = np.log10(10**(Heckman2014_1[:, 0])  /self.Hubble_h/self.Hubble_h)
        Heckman2014_yval = np.log10(10**(Heckman2014_1[:, 1])  / 10**(3.0/2.0)* 5**(3.0/2.0) / self.Hubble_h/self.Hubble_h)
        Heckman2014_xval_2 = np.log10(10**(Heckman2014_2[:, 0])/self.Hubble_h/self.Hubble_h)
        Heckman2014_yval_2 = np.log10(10**(Heckman2014_2[:, 1])  / 20**(3.0/2.0)* 5**(3.0/2.0) )
                                
        Shabala2008_xval = np.log10(10**(Shabala2008[:, 0]) /self.Hubble_h/self.Hubble_h)
        Shabala2008_yval = np.log10(10**(Shabala2008[:, 1]) /self.Hubble_h/self.Hubble_h)

        plt.scatter(Lradio1400_5_E, Q_jet_E, marker='o', s=10, color='blue', alpha=0.15,label='Central Galaxies (Jet-model)')
        plt.scatter(Lradio1400_5_S, Q_jet_S, marker='o', s=10, color='red', alpha=0.15,label='Satellite Galaxies (Jet-model)')

        plt.plot(Heckman2014_xval_2, Heckman2014_yval_2, 'r--', lw = 5, alpha=0.95, label='Heckman--Best (2014) ,fw = 5')
        plt.plot(Shabala2008_xval, Shabala2008_yval, 'g--', lw = 6, alpha=0.85, label='Shabala et al. (2008)')


        plt.ylabel(r'Log $Q_{jet}$ [W]')  # Set the y...
        plt.xlabel(r'Log $L_{1.4 ~GHz}$ [W/Hz]')  # and the x-axis labels
                                
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
                                
        plt.axis([22., 26.0, 33.5, 37.0])
        #        plt.axis([14.0, 25.0, 25, 37.0])
        leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
                                
        outputFile = OutputDir + '16.Lradio_Qjet' + OutputFormat
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
        
        w = np.where((G.Type == 0) & (np.log10(G.Mvir * 1e10 / self.Hubble_h) >11.0))[0]
        if(len(w) > dilute): w = sample(w, dilute)

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


        plt.scatter(R_vir, R_shocked, marker='o', s=30, color='grey', alpha=0.15, label='Central-Galaxies')
        
        plt.scatter(R_vir_Hi, R_shocked_Hi, marker='o', s=30, c='r', alpha=0.25, label='Satellite-Galaxies-Hi')
        plt.scatter(R_vir_Lo, R_shocked_Lo, marker='o', s=30, c='b', alpha=0.45, label='Satellite-Galaxies-Lo')
        
        plt.ylabel(r'$R_{shocked}$ [kpc]')  # Set the y...
        plt.xlabel(r'$R_{vir}$ [kpc]')  # and the x-axis labels
        
        plt.axis([50.0, 800, 0.0, 200.0])
        
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '17.Rshocked_Rvir' + OutputFormat
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
        if(len(w) > dilute): w = sample(w, dilute)
        mass = np.log10(G.StellarMass[w] * 1e10 / self.Hubble_h)
        Lradio1400 = np.log10(G.RadioLuminosity[w,5])
        R_shocked      = np.log10(G.Rshocked[w]* 1000.0/self.Hubble_h)
        
        w1 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>13)& (G.Rshocked >0.0001))[0]
        if(len(w1) > dilute): w1 = sample(w1, dilute)

        Lradio1400_Hi = np.log10(G.RadioLuminosity[w1,5])
        R_shocked_Hi      = np.log10(G.Rshocked[w1]* 1000.0/self.Hubble_h)
        w2 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)<13)& (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>11)& (G.Rshocked >0.001))[0]
        if(len(w2) > dilute): w2 = sample(w2, dilute)

        Lradio1400_Lo = np.log10(G.RadioLuminosity[w2,5])
        R_shocked_Lo      = np.log10(G.Rshocked[w2]* 1000.0/self.Hubble_h)
        
        
        
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
        Gendre2010_FRI= np.array([
                                    [5.1323908361550235, 2.005736894159895E22, 11.0],
                                    [5.248385471540101, 8.793977836134125E22 , 11.5],
                                    [83.03553296858263, 2.3923899116825847E25, 11.5],
                                    ], dtype=np.float32)
        Gendre2010_FRII= np.array([
                                    [4.2526476719508794, 3.128598004878272E22, 10.5],
                                    [5.065412382783366, 1.2961846939896264E23, 11.5],
                                    [17.104138564850306, 1.691126589505699E23, 11.0],
                                    [41.96295131773274, 1.583242901419614E23, 10.5],
                                    [128.55888635760584, 1.984215220139971E25, 11.5],
                                    [275.689369297745, 1.2034875257263966E25, 11.5],
                                    ], dtype=np.float32)

                                    
        Shabala2008_LS_xval = Shabala2008_LS[:,0]/2.0
        Shabala2008_LS_yval = np.log10(Shabala2008_LS[:,1])
        
        Shabala2008_LS_2_xval = Shabala2008_LS_2[:,0]/2.0
        Shabala2008_LS_2_yval = np.log10(Shabala2008_LS_2[:,1])

        Gendre2010_FRI_xval = np.log10(Gendre2010_FRI[:,0]/2.0)
        Gendre2010_FRI_yval = np.log10(Gendre2010_FRI[:,1])

        Gendre2010_FRII_xval = np.log10(Gendre2010_FRII[:,0]/2.0)
        Gendre2010_FRII_yval = np.log10(Gendre2010_FRII[:,1])

        plt.scatter(R_shocked, Lradio1400, marker='o', s=30, c='b', alpha=0.15, label='Central-Galaxies')
#        plt.scatter(R_shocked_Hi, Lradio1400_Hi, marker='o', s=30, c='r', alpha=0.25, label='Satellite-Galaxies-Hi')
#        plt.scatter(R_shocked_Lo, Lradio1400_Lo, marker='o', s=30, c='b', alpha=0.45, label='Satellite-Galaxies-Lo')
#        plt.plot(Shabala2008_LS_xval, Shabala2008_LS_yval, 'g--', lw = 7, alpha=0.95, label='Shabala et al. (2008)')
        plt.scatter(Gendre2010_FRI_xval, Gendre2010_FRI_yval, marker='^', s=199, color='brown', alpha=0.95, label='Gendre+10 (FRI)')
        plt.scatter(Gendre2010_FRII_xval, Gendre2010_FRII_yval, marker='*', s=199, color='red', alpha=0.95, label='Gendre+10 (FRII)')

        plt.ylabel(r'Log $L_{1.4~GHz}$ [W/Hz]')  # Set the y...
        plt.xlabel(r'Log $R_{Shocked}$ [kpc]')  # and the x-axis labels
        
        # Set the x and y axis minor ticks
        #        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        
        plt.axis([0, 2.5, 22, 28.0])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '18.Lradio_Rshock' + OutputFormat
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
        
        binwidth = 0.3  # Radio Luminosity function histogram bin width
        

        w1 = np.where((G.RadioLuminosity[:,5] > 0)& (G.RadioLuminosity[:,5] < 1e50)&(np.log10(G.CentralMvir * 1e10/self.Hubble_h) > 11))[0]
#        if(len(w1) > dilute): w1 = sample(w1, dilute)

        Lradio1400_5 = np.log10(G.RadioLuminosity[w1,5])
        mi_5 = np.floor(min(Lradio1400_5)) - 2.0
        ma_5 = np.floor(max(Lradio1400_5)) + 2.0
        NB_5 = (ma_5 - mi_5) / binwidth
        
        (counts_5, binedges_5) = np.histogram(Lradio1400_5, range=(mi_5, ma_5), bins=NB_5 )
        
        # Set the x-axis values to be the centre of the bins
        xaxeshisto_5 = binedges_5[:-1] + 0.5 * binwidth
        
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 10.5)
        Lradio1400_5_1 = Lradio1400_5[w]
        (counts_5_1, binedges_5) = np.histogram(Lradio1400_5_1, range=(mi_5, ma_5), bins=NB_5 )
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 10.6)
        Lradio1400_5_2 = Lradio1400_5[w]
        (counts_5_2, binedges_5) = np.histogram(Lradio1400_5_2, range=(mi_5, ma_5), bins=NB_5 )
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 10.7)
        Lradio1400_5_3 = Lradio1400_5[w]
        (counts_5_3, binedges_5) = np.histogram(Lradio1400_5_3, range=(mi_5, ma_5), bins=NB_5 )
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 10.8)
        Lradio1400_5_4 = Lradio1400_5[w]
        (counts_5_4, binedges_5) = np.histogram(Lradio1400_5_4, range=(mi_5, ma_5), bins=NB_5 )

        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 10.9)
        Lradio1400_5_5 = Lradio1400_5[w]
        (counts_5_5, binedges_5) = np.histogram(Lradio1400_5_5, range=(mi_5, ma_5), bins=NB_5 )

        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.0)
        Lradio1400_5_6 = Lradio1400_5[w]
        (counts_5_6, binedges_5) = np.histogram(Lradio1400_5_6, range=(mi_5, ma_5), bins=NB_5 )
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.1)
        Lradio1400_5_7 = Lradio1400_5[w]
        (counts_5_7, binedges_5) = np.histogram(Lradio1400_5_7, range=(mi_5, ma_5), bins=NB_5 )
        
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.2)
        Lradio1400_5_8 = Lradio1400_5[w]
        (counts_5_8, binedges_5) = np.histogram(Lradio1400_5_8, range=(mi_5, ma_5), bins=NB_5 )

        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.3)
        Lradio1400_5_9 = Lradio1400_5[w]
        (counts_5_9, binedges_5) = np.histogram(Lradio1400_5_9, range=(mi_5, ma_5), bins=NB_5 )
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.4)
        Lradio1400_5_10 = Lradio1400_5[w]
        (counts_5_10, binedges_5) = np.histogram(Lradio1400_5_10, range=(mi_5, ma_5), bins=NB_5 )
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.5)
        Lradio1400_5_11 = Lradio1400_5[w]
        (counts_5_11, binedges_5) = np.histogram(Lradio1400_5_11, range=(mi_5, ma_5), bins=NB_5 )
        w = np.where(np.log10(G.StellarMass[w1] * 1e10/self.Hubble_h)< 11.6)
        Lradio1400_5_12 = Lradio1400_5[w]
        (counts_5_12, binedges_5) = np.histogram(Lradio1400_5_12, range=(mi_5, ma_5), bins=NB_5 )

        # Shabala+ 2008 modified data used for the MCMC fitting
        shabala08 = np.array([[21.8950605434631,    -3.917819181692233],
                              [22.2948526432552,   -3.9842210188804437],
                              [22.696107350899684,  -4.486268526719551],
                              [23.09612760954614,   -4.966625533828224],
                              [23.494927367770817,  -5.450599343171034],
                              [23.898406873327,     -5.553150050591033],
                              [24.297173752764806,  -5.977470931517359],
                              [24.700680159146614,  -6.128828580369028],
                              [25.098086056073278,  -6.301885569386324],
                              [25.49358498336134,   -7.451084333515425],
                              [25.49876986841843,   -7.986147131809531] ], dtype=np.float32)
                              # Max thesis
                              
        HeckmanBest2014 = np.array([
                                   [21.997862863932504,  6.20527276386501E-5],
                                   [23.02105429376662,  2.4520014625693308E-5],
                                   [23.78657693364518,  1.190892351578673E-5],
                                   [24.14683701359509,  7.88136906640192E-6],
                                   [24.33448081357317,  6.088865023750937E-6],
                                   [24.477095857309383,  4.858107202448261E-6],
                                   [24.70229219104495,  3.1528819497691173E-6],
                                   [25.02510233978862,  1.4627878588513268E-6],
                                   [25.408018814805402,  4.577502526265837E-7],
                                   [25.68333399232533,  1.8305328229999118E-7],
                                   [26.306581928637677,  1.924150093649187E-8],
                                   [26.992412416210005,  1.552246355517814E-9],
                                    ], dtype=np.float32)
        HeckmanBest2014_radio = np.array([
                                         [22.057815903066597, 1.942144902262784E-6],
                                         [22.971383991130573, 9.281551561066824E-7],
                                         [24.449600144071947, 2.818090433815472E-7],
                                         [25.359421517370237, 1.3273584076302237E-7],
                                         [25.70272498578978, 9.652085388261188E-8],
                                         [25.90158507279573, 7.842455794744196E-8],
                                         [26.070436267861265, 6.402633276676175E-8],
                                         [26.19989560411731, 5.3547140119972934E-8],
                                         [26.35000478363883, 4.1852316941308487E-8],
                                         [26.4776055355143, 3.302824509123081E-8],
                                         [26.667153003281015, 2.148176133550013E-8],
                                         [26.873616262120898, 1.2321018544427695E-8]], dtype=np.float32)
                                                                

        shabala08_xval = np.log10(10**(shabala08[:, 0]) )
        shabala08_yval = (10**(shabala08[:, 1]))
        
                              
        Heckmanbest14_xval = np.log10(10**(HeckmanBest2014[:, 0]) )
        Heckmanbest14_yval = ((HeckmanBest2014[:, 1]))
                              
        Heckmanbest14_radio_xval = np.log10(10**(HeckmanBest2014_radio[:, 0])  )
        Heckmanbest14_radio_yval = ((HeckmanBest2014_radio[:, 1]))

        plt.plot(xaxeshisto_5, counts_5_1   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 5, alpha=0.15, label='$M_* < 10^{10.5}$')
        plt.plot(xaxeshisto_5, counts_5_2   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 5, alpha=0.35, label='$M_* < 10^{10.6}$')
        plt.plot(xaxeshisto_5, counts_5_3   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 5, alpha=0.15, label='$M_* < 10^{10.7}$')
        plt.plot(xaxeshisto_5, counts_5_4   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 5, alpha=0.35, label='$M_* < 10^{10.8}$')
        plt.plot(xaxeshisto_5, counts_5_5   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 5, alpha=0.15, label='$M_* < 10^{10.9}$')
        plt.plot(xaxeshisto_5, counts_5_6   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'p-',lw = 5, alpha=0.15, label='$M_* < 10^{11}$')
        plt.plot(xaxeshisto_5, counts_5_9   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'b-',lw = 5, alpha=0.95, label='$M_* < 10^{12}$')

        plt.plot(shabala08_xval, shabala08_yval, 'g-.', lw = 7, alpha=0.85, label='Sh08')
        plt.plot(Heckmanbest14_xval, Heckmanbest14_yval, 'r-.', lw = 7, alpha=0.85, label='HB14-JM')
#        plt.plot(Heckmanbest14_radio_xval, Heckmanbest14_radio_yval, 'r--', lw = 7, alpha=0.95, label='HB14-RL')
                              
                              
        plt.yscale('log', nonposy='clip')
        plt.axis([22, 28, 1.0e-8, 3.0e-4])
                              
        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
                              
        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'Log $L_{1.4 ~GHz}$ [W/Hz]')  # Set the y...
                              
        leg = plt.legend(loc='lower left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
                              
        outputFile = OutputDir + '19.RadioLF' + OutputFormat
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
        
        w = np.where((G.Type == 0) & (G.CentralMvir > 3000) & (G.rho_zero_Capelo > 0))[0]
        if(len(w) > dilute): w = sample(w, dilute)
        # Makino density profile 
        r_plot_i  = np.arange(-3,1,0.1)
        r_plot_Makino  = np.zeros(len(r_plot_i))
        rho_gas_Makino = np.zeros(len(r_plot_i))
        j=0
        for i in range(len(r_plot_i)):
            r_plot_Makino[i]=  G.Rvir[w][j] * (10 ** r_plot_i[i])
            rho_gas_Makino[i] = G.rho_zero_Makino[w][j]  * np.exp(-13.5 * G.b_gas[w][j]) * ((1+r_plot_Makino[i] / (G.Rs[w][j])) ** (13.5 * G.b_gas[w][j] /(r_plot_Makino[i]/(G.Rs[w][j]))))
        # Capelo density profile
        w1 = np.where((G.Type == 0) & (G.CentralMvir > 3000) & (G.rho_zero_Capelo > 0))[0]
        if(len(w1) > dilute): w1 = sample(w1, dilute)
        
        r_plot_Capelo  = np.zeros(len(r_plot_i))
        rho_gas_Capelo = np.zeros(len(r_plot_i))
        for i in range(len(r_plot_i)):
            r_plot_Capelo[i]=  G.Rvir[w1][j] * (10 ** r_plot_i[i])
            rho_gas_Capelo[i] = 1e3*G.rho_zero_Capelo[w1][j]* np.exp(-13.5 * G.b_gas[w1][j]) * ((1+r_plot_Capelo[i] / (G.Rs[w][j])) ** (13.5 * G.b_gas[w1][j] /(r_plot_Capelo[i]/(G.Rs[w1][j]))))

        r_plot_beta  = np.zeros(len(r_plot_i))
        rho_gas_beta = np.zeros(len(r_plot_i))
        for i in range(len(r_plot_i)):
            
            r_plot_beta[i]=  G.Rvir[w1][j] * (10 ** r_plot_i[i])
            rho_gas_beta[i] = G.rho_zero_Makino[w1][j]  /(1+(r_plot_beta[i]/(0.22*G.Rs[w1][j]))** 2.0)**(3*0.9*G.b_gas[w1][j]/2)
            print   np.log10(rho_gas_Makino[i]), np.log10(rho_gas_Capelo[i]),rho_gas_beta[i]


        plt.plot(r_plot_Makino/G.Rvir[w1][j], rho_gas_Makino,  'r-', lw = 4, alpha=0.5, label='Makino+98')
        plt.plot(r_plot_beta/G.Rvir[w1][j], rho_gas_beta,  'b--', lw = 4, alpha=0.65, label=(r'$\beta$ - model ( $\beta_{eff}$ = 0.9 b , $r_0$ = 0.22 $r_s$)'))
        
        plt.xscale('log', nonposy='clip')
        plt.yscale('log', nonposy='clip')
        
        #        plt.ylabel(r'$\rho  [kg / m^{3}]$ ')  # Set the y...
        plt.ylabel(r'$\rho$ ~ $[kg /m^{3}]$ ')  # Set the y...
        plt.xlabel(r'$r/r_{vir}$')  # and the x-axis labels
        
        # Set the x and y axis minor ticks
#        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
        #        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        
        plt.axis([1e-3, 9,1e-27,1e-20])
        
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')
        
        outputFile = OutputDir + '20.Density_profile' + OutputFormat
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
        
        w = np.where((G.Type == 0)&(np.log10(G.CentralMvir * 1e10 /self.Hubble_h) > 10) & (G.Cooling>39) & (G.Temp_Gas>1e4))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        E_cooling = G.Cooling[w]-40.0 #- 2*np.log10(self.Hubble_h)
        temp_x      = G.Temp_Gas[w] * 8.617328149741e-8  # [K_b T] in [kev]
#        temp_x =35.9*(G.Vvir[w]*G.Vvir[w]) / 11604.5 / 1.0e3

        plt.scatter(np.log10(temp_x), E_cooling, marker='o', s=20, c='b', alpha=0.1, label='Jet-Modle(+Heat, +Uplift) ')

        m_0, b_0  = np.polyfit(np.log10(temp_x), E_cooling, 1)
        plt.plot(np.log10(temp_x), m_0*(np.log10(temp_x))+b_0, 'b-', lw = 4,label='Linear fit', alpha=0.55)

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
        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='c', lw=2.0, alpha=0.6, marker='s', markersize=8, ls='none', label='P98 HRI', mew=1)

        w = np.where((OlumPerrU < cut*OtempP) & ((OlumPerrD < cut*OtempP)))[0]
        xplot = np.log10(OtempP[w])
        yplot = np.log10(OlumP[w])
        yerr2 = np.log10(OlumP[w]+OlumPerrU[w])-yplot
        yerr1 = yplot-np.log10(OlumP[w]-OlumPerrD[w])
        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='c', lw=2.0, alpha=0.6, marker='*', markersize=12, ls='none', label='P98 PSPC', mew=1)

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


        plt.xlabel(r' $log_{10} \ T_{gas}$ [Kev]')  # Set the y...
        plt.ylabel(r' $log_{10}$ (Net cooling $[10^{40}erg~ s^{-1}])$')  # and the x-axis labels
        
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #        plt.xscale('log', nonposy='clip')
        #        plt.yscale('log', nonposy='clip')
        plt.axis([-0.5, 1.2, -2, 7])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('small')
        
        outputFile = OutputDir + '21.Cooling_Temp' + OutputFormat
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
    res.BaryonicMassFunction(G)
    res.GasMassFunction(G)
    res.BaryonicTullyFisher(G)
    res.SpecificStarFormationRate(G)
    res.GasFraction(G)
    res.Metallicity(G)
    res.BlackHoleBulgeRelationship(G)
    res.QuiescentFraction(G)
    res.BulgeMassFraction(G)
    res.BaryonFraction(G)
    res.SpinDistribution(G)
    res.VelocityDistribution(G)
    res.MassReservoirScatter(G)
    res.SpatialDistribution(G)
    res.Lradio_Qjet(G)
    res.Rshocked_Rvir(G)
    res.Lradio_Rshock(G)
    res.RadioLF(G)
    res.Density_profile(G)
    res.cooling_Temp(G)
