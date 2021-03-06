#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

#import ROOT
import h5py as h5
import numpy as np
import pylab as plt
from random import sample, seed
from os.path import getsize as getFileSize
from mpl_toolkits.mplot3d import Axes3D
#from scipy import signal
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from scipy import signal as ss
from matplotlib.colors import LogNorm
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
            ('RadioLuminosity_lifetime'     , (np.float32, 600)),
            ('Rshocked_lifetime'            , (np.float32, 600)),
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
            ('delta'                        , np.float32),
            ('t_cool_Makino'                , np.float32)
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
            
            
            
#        coolFrac_mMS_xval      = coolFrac_mMS[:,0]
#        coolFrac_mMS_yval      = coolFrac_mMS[:,1]



        plt.plot(xaxeshisto, counts    / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h / binwidth, 'b-', lw = 2, label='Jet-model (Heat, Uplift)')
        

        plt.yscale('log', nonposy='clip')
        plt.axis([8.4, 12.7, 3.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ [\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1}]$')  # Set the y...
        plt.xlabel(r'$\log_{10} (m_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

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



#-----------------------------------------
    def Lradio_Qjet(self, G):
    
        print 'Plotting the Radio Luminosity -- Qjet relation'
        
        seed(2222)
        
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        w = np.where((np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)
        #w1 = np.where((G.Type == 0) &(G.RadioLuminosity[:,5]>0)&(G.RadioLuminosity[:,5]<1e40) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >11.0))[0]
        w1 = np.where((G.Type >= 0) &(G.RadioLuminosity[:,5]>0))[0]
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
        Q_jet_1      = G.Qjet[w]


        Lradio1400_5_S = np.log10(G.RadioLuminosity[w2,5])
        #Q_jet_S      = np.log10(G.Qjet[w2])
        Q_jet_S      = G.Qjet[w2]
        Lradio1400_5_E = np.log10(2.0* G.RadioLuminosity[w1,5])
        Q_jet_E      = G.Qjet[w1]- np.log10(2)
        
        
        Heckman2014_1 = np.array([
                                 [20.01386138613861, 33.269085411942555],
                                 [21.275247524752473, 34.354497354497354],
                                 [22.728382838283828, 35.600151171579746],
                                 [24.768316831683165, 37.3567649281935],
                                 [26.547194719471946, 38.88662131519274],
                                 ], dtype=np.float64)

        Heckman2014_2 = np.array([
                                    [20.01735733221337,  34.033570935094446],
                                    [22.775669620493755,  35.910935356386275],
                                    [24.769305876740944,  37.25960257094492],
                                    [26.916561848624763,  38.71132532014027],
                                       ], dtype=np.float64)


        # For Eq.1 and 2 in Heckman & Best (2014)
        Heckman2014_xval_1 = np.log10(10**Heckman2014_1[:, 0])
        Heckman2014_yval_1 = np.log10(10**Heckman2014_1[:, 1])
        Heckman2014_xval_2 = np.log10(10**Heckman2014_2[:, 0])
        Heckman2014_yval_2 = np.log10(10**Heckman2014_2[:, 1])
                                
 
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
        plt.plot(bins-delta/2,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (Median$~\pm~\sigma$)')
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)
        
        plt.plot(Heckman2014_xval_2, Heckman2014_yval_2, 'k-',alpha=.9, marker='o', markersize=3, label=(r'Heckman \& Best~(2014)'))
        plt.plot(Heckman2014_xval_1, Heckman2014_yval_1, 'k--',alpha=.9, marker='o', markersize=3, label=(r'Cavagnolo et al. (2010)'))


        plt.ylabel(r'$\log_{10} (Q_{jet} ~ [W])$')  # Set the y...
        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
                                
        plt.axis([20., 26.0, 33, 39.0])
        #        plt.axis([14.0, 25.0, 25, 37.0])
        leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
                                
        outputFile = OutputDir + '16_Lradio_Qjet' + OutputFormat
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
        
        w = np.where((G.Type == 0) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) >10.0)&(G.StellarMass>0)& (G.Rshocked>0.0)& (G.Rshocked<3.0) )[0]
        #if(len(w) > dilute): w = sample(w, dilute)

        R_shocked = G.Rshocked[w]* 1000.0
        R_vir      = G.Rvir[w]* 1000.0
        
        w1 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>12.5))[0]
        #if(len(w1) > dilute): w1 = sample(w1, dilute)

        R_vir_Hi      = G.Rvir[w1]* 1000.0
        R_shocked_Hi      = G.Rshocked[w1]* 1000.0
        w2 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h) <12.5))[0]
        #if(len(w2) > dilute): w2 = sample(w2, dilute)

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
        plt.plot(bins-delta/2.,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (Median$~\pm~\sigma$)')
        #running_std    = [np.sqrt(np.sum(Y[idx==k]**2)/len(Y[idx==k]) - (np.sum(Y[idx==k])*np.sum(Y[idx==k])/len(Y[idx==k])**2)) for k in range(total_bins)]
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2.,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)
        print "ID = ",bins-delta/2.,running_median,running_std

        m_0 = 1.0
        b_0 = 1.0
        x = [10,100,300,500, 1000]
        y = [10,100,300,500, 1000]
        plt.plot(x, y, 'k--', lw = 2,label='1:1')

        plt.ylabel(r'$r_{\mathrm{shock}}~ [kpc]$')  # Set the y...
        plt.xlabel(r'$R_{\mathrm{vir}}~ [kpc]$')  # and the x-axis labels

        plt.axis([100.0, 600, 0.0, 900.0])
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '17_Rshocked_Rvir' + OutputFormat
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
        
        w = np.where((G.Type == 0) & (G.Rshocked > 0.001))[0]
#        if(len(w) > dilute): w = sample(w, dilute)
        mass = np.log10(G.StellarMass[w] * 1e10 / self.Hubble_h)
        Lradio1400 = np.log10(2.0* G.RadioLuminosity[w,5])
        R_shocked      = np.log10(G.Rshocked[w]* 1000.0)
        
        w1 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>11)& (G.Rshocked >0.0001))[0]
#        if(len(w1) > dilute): w1 = sample(w1, dilute)

        Lradio1400_Hi = np.log10(G.RadioLuminosity[w1,5])
        R_shocked_Hi      = np.log10(G.Rshocked[w1]* 1000.0/self.Hubble_h)
        w2 = np.where((G.Type ==1) & (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)<13)& (np.log10(G.CentralMvir * 1e10 / self.Hubble_h)>11)& (G.Rshocked >0.001))[0]
#        if(len(w2) > dilute): w2 = sample(w2, dilute)

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
#                                    [88.6209924089638, 1.107490389102193E25, 12],
                                    ], dtype=np.float32)
        Gendre2010_FRII= np.array([
                                    [4.2526476719508794, 3.128598004878272E22, 10.5],
                                    [5.065412382783366, 1.2961846939896264E23, 11.5],
#                                    [5.114602036158313, 2.0171985260377704E24, 12.0],
                                    [17.104138564850306, 1.691126589505699E23, 11.0],
                                    [41.96295131773274, 1.583242901419614E23, 10.5],
                                    [128.55888635760584, 1.984215220139971E25, 11.5],
                                    [275.689369297745, 1.2034875257263966E25, 11.5],
                                    ], dtype=np.float32)

                                    
        Shabala2008_LS_xval = np.log10(Shabala2008_LS[:,0]/2.0)
        Shabala2008_LS_yval = np.log10(Shabala2008_LS[:,1])
        
        Shabala2008_LS_2_xval = np.log10(Shabala2008_LS_2[:,0]/2.0)
        Shabala2008_LS_2_yval = np.log10(Shabala2008_LS_2[:,1])

        Gendre2010_FRI_xval = np.log10(Gendre2010_FRI[:,0]/2.0/self.Hubble_h)
        Gendre2010_FRI_yval = np.log10(Gendre2010_FRI[:,1])

        Gendre2010_FRII_xval = np.log10(Gendre2010_FRII[:,0]/2.0/self.Hubble_h)
        Gendre2010_FRII_yval = np.log10(Gendre2010_FRII[:,1])

        x = []
        y = []
        z = []
        w = np.where((G.Type == 0) & (np.log10(G.RadioLuminosity[:,5]) >20.0) & (G.RadioLuminosity[:,5] <1e30)& (G.Rshocked >0))[0]
        for j in xrange(len(G.Rshocked[w])):
           for i in xrange(100):
              if ((G.RadioLuminosity_lifetime[j,i] > 0) & (G.Rshocked_lifetime[j,i]>0)):
                y.append((np.log10(2.0*G.RadioLuminosity_lifetime[j,i])))
                x.append(np.log10(G.Rshocked_lifetime[j,i]*1000.0))
                z.append(i)

        counts,ybins,xbins,image = plt.hist2d(x,y,weights=z,bins=300,norm=LogNorm())
        plt.colorbar()
        plt.contour(counts,linewidths=5) #,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],linewidths=3

        total_bins = 50
        X = R_shocked
        Y = Lradio1400
        bins = np.linspace(X.min(),3, total_bins)
        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (z = 0)')
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)


        plt.ylabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # Set the y...
        plt.xlabel(r'$\log_{10} (R_{\mathrm{Shock}}~ [kpc])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        #        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

        plt.axis([0, 3.2, 22, 27.0])
        
        
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
    
        outputFile = OutputDir + '19_Lradio_Rshock' + OutputFormat
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
        
        w = np.where((np.log10(G.CentralMvir* 1e10/self.Hubble_h) > 11)  & (G.RadioLuminosity[:,5]>0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)
        AGN_frac_on = 1- G.fcool[w]
        Lradio1400 = np.log10( 2* G.RadioLuminosity[w,5])
        mass      = np.log10(G.StellarMass[w] * 1e10/self.Hubble_h)
        plt.scatter( Lradio1400, mass, marker='o', s=10, color = 'grey', alpha=0.15)

        total_bins = 40
        X = Lradio1400
        Y = mass
        bins = np.linspace(X.min(),X.max(), total_bins)
#        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (Median$~\pm~\sigma$)')
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)


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

        x11 = [24, 24 ,24]
        y11 = [10.5, 11.5, 12.5]
        x22 = [22, 25 , 27]
        y22 = [11,11,11]
        
        plt.plot(x11,y11,'k--',lw=1,alpha=.5)
        plt.plot(x22,y22,'k--',lw=1,alpha=.5)
        
        plt.xlabel(r'$\log_{10} (L_{\mathrm{1.4 ~GHz}}\ [W~ Hz^{-1}])$')  # Set the y...
        plt.ylabel(r'$\log_{10} (m_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        plt.axis([ 23.0, 26.5,10.5, 12.0])
        
        
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
    
        outputFile = OutputDir + '20_Lradio_Mass' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)

# ---------------------------------------------------------

    def Temp_hist(self, G):
    
        print 'Plotting the Temp_hist relation'
        
        seed(2222)
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        w = np.where((G.Type == 0)&(G.Temp_Gas>0) & (np.log10(G.CentralMvir*1e10/self.Hubble_h) > 11)&(G.StellarMass>0))[0]
        #if(len(w) > dilute): w = sample(w, dilute)

        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
        Temp = np.log10(G.Temp_Gas[w])
        Tvir = np.log10(35.9 * G.Vvir[w] * G.Vvir[w])
        
        plt.scatter(mass,Temp/Tvir, marker='o', s=10, color = 'grey', alpha=0.05)

        total_bins = 70
        X = mass
        Y = Temp/Tvir
        bins = np.linspace(X.min(),11.6, total_bins)
        #        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (Median$~\pm~\sigma$)')
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)

        x=[9,11,12,13]
        y=[1,1,1,1]
        plt.plot(x, y, 'k--', lw=2)


        plt.ylabel(r'$T_{\mathrm{new-hot}}/T_{\mathrm{vir}} $')  # Set the y...
        plt.xlabel(r'$\log_{10} (m_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        #plt.xscale('log', nonposy='clip')
        plt.axis([10, 12.0, 0.95, 1.4])
        
        
        leg = plt.legend(loc='higher left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')




        outputFile = OutputDir + '23_Temp_hist' + OutputFormat
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
        
        binwidth = 0.6  # Radio Luminosity function histogram bin width
        w1 = np.where((G.RadioLuminosity[:,5] > 0)&(G.StellarMass >0)& (G.fcool > 0))[0]

#        if(len(w1) > dilute): w1 = sample(w1, dilute)
        Lradio1400_5 = np.log10(2.0* G.RadioLuminosity[w1,5])
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
                             
                             
#        0.15 is base on the bin center in Best 2012
        xplot= np.log10(10**(Best2012[:,0]))+0.15
        yplot = (10**Best2012[:,4])
        yerr2 = 10**(Best2012[:,4]+Best2012[:,5]) - yplot
        yerr1 = yplot - 10**(Best2012[:,4]-Best2012[:,6])
        plt.fill_between(xplot, yplot+yerr2, yplot-yerr1,facecolor='red', alpha=0.35,label= 'Best \& Heckman (2012)')
        plt.plot(xaxeshisto_5, counts_5_12   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, color='red',lw = 4, alpha=0.35, label='Best \& Heckman (2012)')

        plt.plot(xaxeshisto_5, counts_5_12   / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5,'c', lw=5.0, alpha=0.8, marker='s', markersize=13, label='Jet model')
        plt.plot(xaxeshisto_5, counts_5_12  / self.volume * self.Hubble_h*self.Hubble_h*self.Hubble_h/xaxeshisto_5, 'b--',lw=5,alpha=.8)



        AGNseen = np.zeros(len(w1))
    
        for i in xrange(len(w1)-1):
            random = np.random.randint(0, 10000)/10000.0
#            print random
            if (random < delta_duty[i]):
                AGNseen[i] = 1.0
            else:
                AGNseen[i] = 0.0
    
        w = np.where(AGNseen[:] > 0)[0]
        print len(w)
        Lradio1400_5_12 = Lradio1400_5[w]
        (counts_5_12, binedges_5) = np.histogram(Lradio1400_5_12, range=(mi_5, ma_5), bins=NB_5 )
        
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
                              
        outputFile = OutputDir + '24_RadioLF' + OutputFormat
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
        
        #w = np.where((G.Type == 0)&(np.log10(G.CentralMvir * 1e10 /self.Hubble_h) > 11) & (G.Cooling>39) & (G.Temp_Gas>1e4))[0]
        w = np.where((G.Type == 0) & (G.Cooling>0) & (G.Temp_Gas>0))[0]
#        if(len(w) > dilute): w = sample(w, dilute)

        E_cooling =np.log10(10**(G.Cooling[w]-40.0))
        temp_x      = G.Temp_Gas[w] * 8.617328149741e-8  # [K_b T] in [kev]

        plt.scatter(np.log10(temp_x), E_cooling, marker='o', s=10, color = 'grey', alpha=0.1)


        total_bins = 40
        X = np.log10(temp_x)
        Y = E_cooling
        bins = np.linspace(X.min(),1.2, total_bins)
        print bins
        delta = bins[1]-bins[0]
        idx  = np.digitize(X,bins)
        running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
        plt.plot(bins-delta/2,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (Median$~\pm~\sigma$)')
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)


        w = np.where((G.Type == 1)&(np.log10(G.CentralMvir* 1e10 /self.Hubble_h) > 11))[0]
        if(len(w) > dilute): w = sample(w, dilute)

        E_cooling = G.Cooling[w]-40.0 #- 3*np.log10(self.Hubble_h)
        temp_x      = G.Temp_Gas[w] * 8.617328149741e-8  # [K_b T] in [kev]


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
        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='c', lw=2.0, alpha=0.3, marker='o', markersize=8, ls='none', label='P98 HRI', mew=1)

        w = np.where((OlumPerrU < cut*OtempP) & ((OlumPerrD < cut*OtempP)))[0]
        xplot = np.log10(OtempP[w])
        yplot = np.log10(OlumP[w])
        yerr2 = np.log10(OlumP[w]+OlumPerrU[w])-yplot
        yerr1 = yplot-np.log10(OlumP[w]-OlumPerrD[w])
        plt.errorbar(xplot, yplot, yerr=[yerr1,yerr2], color='c', lw=2.0, alpha=0.3, marker='*', markersize=12, ls='none', label='P98 PSPC', mew=1)

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


        plt.xlabel(r' $log_{10} \ (T_{\mathrm{new-hot}}~ [Kev])$')  # Set the y...
        plt.ylabel(r' $log_{10}~(\mathrm{Net\ cooling}\ [10^{40}\mathrm{erg}~ \mathrm{s}^{-1}])$')  # and the x-axis labels

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
        
        outputFile = OutputDir + '26_Cooling_Temp' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
        
        # Add this plot to our output list
        OutputList.append(outputFile)


# ---------------------------------------------------------

    def Rshock_hist(self, G):
    
        print 'Plotting the Rshock_Rvir relation'
        
        seed(2222)
#        plt.figure(figsize=(16,6))  # New figure
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        w = np.where((G.Type == 0) & (G.Rshocked>0.0)& (G.Rshocked<3.0))[0]

        mass = np.log10(G.StellarMass[w] * 1.0e10 / self.Hubble_h)
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
        plt.plot(bins-delta/2,running_median,'b-',lw=2,alpha=.9, marker='o', markersize=3, label='Jet-model (Median$~\pm~\sigma$)')
        running_std    = [Y[idx==k].std() for k in range(total_bins)]
        plt.errorbar(bins-delta/2,running_median,running_std,color='b', lw=2.0, alpha=0.4, marker='o', markersize=3, ls='none', mew=1)
        
        
        x=[9,11,12,13]
        y=[1,1,1,1]
        plt.plot(x, y, 'k--', lw=2)
        
        plt.ylabel(r'$r_{Shock}/R_{vir}$')  # Set the y...
        plt.xlabel(r'$\log_{10} (m_{\mathrm{*}}\ [M_{\odot}])$')  # and the x-axis labels

        # Set the x and y axis minor ticks
       #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
       #ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #plt.xscale('log', nonposy='clip')
        
        plt.axis([10, 12.0, 0, 2.5])
        
        
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
        
        outputFile = OutputDir + '28_Rshock_hist' + OutputFormat
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
    res.Lradio_Qjet(G)
    res.Rshocked_Rvir(G)
    res.Lradio_Rshock(G)
    res.Lradio_Mass(G)
    res.Temp_hist(G)
    res.RadioLF(G)
    res.cooling_Temp(G)
    res.Rshock_hist(G)
