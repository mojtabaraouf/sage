"""
Ensure the data on TAO (based on an CSV download) match the actual output files of Dark Sage.  Currently set up specifically for the 2016 version of the model.
"""

import numpy as numpy


def galdtype():
    floattype = np.float32
    Nannuli = 30
    # Names of fields should be consistent with TAO
    Galdesc_full = [
                    ('Galaxy_Classification'        , np.int32),
                    ('Galaxy_ID'                    , np.int64),
                    ('Halo_Index'                   , np.int32),
                    ('Simulation_Halo_ID'           , np.int32),
                    ('Tree_Index'                   , np.int32),
                    ('Snapshot_Number'              , np.int32),
                    ('Central_Galaxy_ID'            , np.int64),
                    ('Central_Galaxy_Mvir'          , floattype),
                    ('Merger_Type'                  , np.int32), #internal
                    ('Descendant_Galaxy_Index'      , np.int32), #
                    ('Descendant_Snapshot'          , np.int32), #
                    ('Galaxy_Age'                   , floattype), #
                    ('X'                            , floattype),
                    ('Y'                            , floattype),
                    ('Z'                            , floattype),
                    ('X_Velocity'                   , floattype),
                    ('Y_Velocity'                   , floattype),
                    ('Z_Velocity'                   , floattype),
                    ('jX_Halo'                      , floattype), 
                    ('jY_Halo'                      , floattype), 
                    ('jZ_Halo'                      , floattype), 
                    ('Total_Particles'              , np.int32),
                    ('Maximum_Number_of_Particles_over_History'                  , np.int32),
                    ('Mvir'                         , floattype),
                    ('Rvir'                         , floattype),
                    ('Vvir'                         , floattype),
                    ('Vmax'                         , floattype),
                    ('Velocity_Dispersion'          , floattype),
                    ('DiscRadii'                    , (floattype, Nannuli+1)), #
                    ('Cold_Gas_Mass'                , floattype),
                    ('Total_Stellar_Mass'           , floattype),
                    ('Merger-driven_Bulge_Mass'     , floattype),
                    ('Instability-drive_Bulge_Mass' , floattype),
                    ('Hot_Gas_Mass'                 , floattype),
                    ('Ejected_Gas_Mass'             , floattype),
                    ('Black_Hole_Mass'              , floattype),
                    ('Intracluster_Stars_Mass'      , floattype),
                    ('DiscGas'                      , (floattype, Nannuli)), #
                    ('DiscStars'                    , (floattype, Nannuli)), #
                    ('X_Spin_of_Stellar_Disk'       , floattype),
                    ('Y_Spin_of_Stellar_Disk'       , floattype),
                    ('Z_Spin_of_Stellar_Disk'       , floattype),
                    ('X_Spin_of_Gas_Disk'           , floattype),
                    ('Y_Spin_of_Gas_Disk'           , floattype),
                    ('Z_Spin_of_Gas_Disk'           , floattype),
                    ('SpinMergerBulge'              , (floattype, 3)), #
                    ('StarsInSitu'                  , floattype), #
                    ('StarsInstability'             , floattype), #
                    ('StarsMergeBurst'              , floattype), #
                    ('DiscHI'                       , (floattype, Nannuli)), #
                    ('DiscH2'                       , (floattype, Nannuli)), #
                    ('DiscSFR'                      , (floattype, Nannuli)), #
                    ('AccretedGasMass'              , floattype), #
                    ('EjectedSNGasMass'             , floattype), #
                    ('EjectedQuasarGasMass'         , floattype), #
                    ('TotInstabEvents'              , np.int32), #
                    ('TotInstabEventsGas'           , np.int32), #
                    ('TotInstabEventsStar'          , np.int32), #
                    ('TotInstabAnnuliGas'           , np.int32), #
                    ('TotInstabAnnuliStar'          , np.int32), #
                    ('FirstUnstableAvGas'           , floattype), #
                    ('FirstUnstableAvStar'          , floattype), #
                    ('Metals_Cold_Gas_Mass'         , floattype),
                    ('Metals_Total_Stellar_Mass'    , floattype),
                    ('Metals_Merger_Bulge_Mass'     , floattype),
                    ('Metals_Instability_Bulge_Mass', floattype),
                    ('Metals_Hot_Gas_Mass'          , floattype),
                    ('Metals_Ejected_Gas_Mass'      , floattype),
                    ('Metals_IntraCluster_Stars_Mass', floattype),
                    ('DiscGasMetals'                , (floattype, Nannuli)), #
                    ('DiscStarsMetals'              , (floattype, Nannuli)), #
                    ('SfrDisk'                      , floattype), #
                    ('SfrBulge'                     , floattype), #
                    ('SfrDiskZ'                     , floattype), #
                    ('SfrBulgeZ'                    , floattype), #
                    ('Cooling_Scale_Radius'         , floattype),
                    ('BulgeRadius'                  , floattype), #
                    ('Hot_Gas_Cooling_Rate'         , floattype),
                    ('AGN_Heating_Rate'             , floattype),
                    ('Time_since_Last_Major_Merger' , floattype),
                    ('Supernova_Cold_Gas_Outflow_Rate', floattype),
                    ('Subhalo_Mvir_at_Infall'        , floattype),
                    ('Subhalo_Vvir_at_Infall'        , floattype),
                    ('Subhalo_Vmax_at_Infall'        , floattype)
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



# Reading of TAO data
def tao_csv(fname, keylist=None):
    with open(fname, 'r') as f: line = f.readline()
    keys = line.split(',')
    keys[-1] = keys[-1][:-1] # gets rid of \n at the end
    if keylist==None: keylist = keys
    print 'Number of properties =', len(keys)
    dict = {}
    for i, key in enumerate(keys):
        if key in ['Total_Particles', 'Maximum_Number_of_Particles_over_History', 'Snapshot_Number', 'Galaxy_Classification']:
            datatype = np.int32
        elif key in ['Galaxy_ID', 'Central_Galaxy_ID', 'Simulation_Halo_ID']:
            datatype = np.int64
        else:
            datatype = np.float32
        if key in keylist:
            print 'Reading', i, key
            dict[key] = np.loadtxt(fname, skiprows=1, usecols=(i,), dtype=datatype, delimiter=', ')
    return dict






# Still need code to actually read multiple Dark Sage files and point to directory/file
G = read_darksage()

# Read specified TAO output
T = tao_csv()


fG = np.in1d(G['Galaxy_ID'], T['Galaxy_ID'])
aG = np.argsort(G['Galaxy_ID'][fG])

fT = np.in1d(T['Galaxy_ID'], G['Galaxy_ID'])
aT = np.argsort(T['Galaxy_ID'][fT])

common_fields = np.intersect1d(G.keys(), T.keys())

print 'The following fields match?'
for field in common_fields:
    if G[field].dtype==float:
        print field, np.allclose(G[field][fG][aG], T[field][fT][aT])
    elif G[field].dtype==int:
        print field, np.all(G[field][fG][aG] == T[field][fT][aT])
