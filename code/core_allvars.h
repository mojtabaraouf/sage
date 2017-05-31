#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "core_simulation.h"

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

#define  STEPS 10         // Number of integration intervals between two snapshots 
#define  MAXGALFAC 1
#define  ALLOCPARAMETER 10.0
#define  MAX_NODE_NAME_LEN 50
#define  ABSOLUTEMAXSNAPS 1000

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define N_BINS 30
#define MIN_STARS_FOR_SN 1e-8
#define MIN_STARFORMATION 1e-10

struct GALAXY_OUTPUT  
{
  int   Type;
  long long   GalaxyIndex;
  int   HaloIndex;
  int SimulationHaloIndex;
  int   TreeIndex;
    
  int   SnapNum;
  long long CentralGalaxyIndex;
  float CentralMvir;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  float   dT;

  // properties of subhalo at the last time this galaxy was a central galaaxy 
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;
    int LenMax;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;
    
  // Radius of each annulus boundary
  float DiscRadii[N_BINS+1];

  // baryonic reservoirs 
  float ColdGas;
  float StellarMass;
  float ClassicalBulgeMass;
  float SecularBulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;
  float DiscGas[N_BINS];
  float DiscStars[N_BINS];
  float SpinStars[3];
  float SpinGas[3];
//  float SpinSecularBulge[3];
  float SpinClassicalBulge[3];
  float StarsInSitu;
  float StarsInstability;
  float StarsMergeBurst;
  float DiscHI[N_BINS];
  float DiscH2[N_BINS];
  float DiscSFR[N_BINS];
    
    // inflow/outflow tracking
//    float AccretedGasMass;
//    float EjectedSNGasMass;
//    float EjectedQuasarGasMass;
    
    // Instability tracking
//  int TotInstabEvents;
//  int TotInstabEventsGas;
//  int TotInstabEventsStar;
//  int TotInstabAnnuliGas;
//  int TotInstabAnnuliStar;
//  float FirstUnstableAvGas;
//  float FirstUnstableAvStar;
//  float TotSinkGas[N_BINS];
//  float TotSinkStar[N_BINS];
    
  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float ClassicalMetalsBulgeMass;
  float SecularMetalsBulgeMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICS;
  float DiscGasMetals[N_BINS];
  float DiscStarsMetals[N_BINS];

  // to calculate magnitudes
  float SfrDisk;
  float SfrBulge;
  float SfrDiskZ;
  float SfrBulgeZ;
  
  // misc 
  float DiskScaleRadius;
  float Cooling;
  float Heating;
  float LastMajorMerger;
  float LastMinorMerger;
  float OutflowRate;

  //infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;
};


struct GALAXY
{
  int   Type;
  int   GalaxyNr;
  int   HaloNr;
  long long  MostBoundID;
  int   SnapNum;
  int   CentralGal;
  double CentralMvir;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  double   dT;

  // properties of subhalo at the last time this galaxy was a central galaxy 
  double Pos[3];
  double Vel[3];
  int   Len;
    int LenMax;
  double Mvir;
  double deltaMvir;
  double Rvir;
  double Vvir;
  double Vmax;
    
  // Radius of each annulus boundary
  double DiscRadii[N_BINS+1];

  // baryonic reservoirs 
  double ColdGas;
  double StellarMass;
  double ClassicalBulgeMass;
  double SecularBulgeMass;
  double HotGas;
  double EjectedMass;
  double BlackHoleMass;
  double ICS;
  double DiscGas[N_BINS];
  double DiscStars[N_BINS];
  double SpinStars[3];
  double SpinGas[3];
  double SpinSecularBulge[3];
  double SpinClassicalBulge[3];
  double SpinHot[3];
  double StarsInSitu;
  double StarsInstability;
  double StarsMergeBurst;
  double DiscHI[N_BINS];
  double DiscH2[N_BINS];
  double DiscSFR[N_BINS];
    
    // inflow/outflow tracking
    double AccretedGasMass;
    double EjectedSNGasMass;
    double EjectedQuasarGasMass;
    
  // Instability tracking
  int TotInstabEvents;
  int TotInstabEventsGas;
  int TotInstabEventsStar;
  int TotInstabAnnuliGas;
  int TotInstabAnnuliStar;
  int FirstUnstableGas;
  int FirstUnstableStar;
  double TotSinkGas[N_BINS];
  double TotSinkStar[N_BINS];

  // metals
  double MetalsColdGas;
  double MetalsStellarMass;
  double ClassicalMetalsBulgeMass;
  double SecularMetalsBulgeMass;
  double MetalsHotGas;
  double MetalsEjectedMass;
  double MetalsICS;
  double DiscGasMetals[N_BINS];
  double DiscStarsMetals[N_BINS];

  // to calculate magnitudes
  double SfrDisk[STEPS];
  double SfrBulge[STEPS];
  double SfrDiskColdGas[STEPS];
  double SfrDiskColdGasMetals[STEPS];
  double SfrBulgeColdGas[STEPS];
  double SfrBulgeColdGasMetals[STEPS];

  // misc 
  double DiskScaleRadius;
    double CoolScaleRadius;
  double MergTime;
  double Cooling;
  double Heating;
  double r_heat;
  double LastMajorMerger;
  double LastMinorMerger;
  double OutflowRate;
  double TotalSatelliteBaryons;

  //infall properties
  double infallMvir;
  double infallVvir;
  double infallVmax;
}
*Gal, *HaloGal;


struct halo_aux_data   // auxiliary halo data 
{
  int DoneFlag;
  int HaloFlag;
  int NGalaxies;
  int FirstGalaxy;
}
*HaloAux;


extern int    FirstFile;    // first and last file for processing 
extern int    LastFile;

extern int    Ntrees;      // number of trees in current file 
extern int    NumGals;     // Total number of galaxies stored for current tree 
extern int    MaxGals;     // Maximum number of galaxies allowed for current tree  
extern int    FoF_MaxGals;

extern int    GalaxyCounter;     // unique galaxy ID for main progenitor line in tree

extern int    LastSnapShotNr;

extern char   OutputDir[512];
extern char   FileNameGalaxies[512];
extern char   TreeName[512];
extern char   SimulationDir[512];
extern char   FileWithSnapList[512];

extern int    TotHalos;
extern int    TotGalaxies[ABSOLUTEMAXSNAPS];
extern int    *TreeNgals[ABSOLUTEMAXSNAPS];

extern int    *FirstHaloInSnap;

extern int    *TreeNHalos;
extern int    *TreeFirstHalo;

#ifdef MPI
extern int ThisTask, NTask, nodeNameLen;
extern char *ThisNode;
#endif

extern double Omega;
extern double OmegaLambda;
extern double PartMass;
extern double Hubble_h;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

// binning information
extern double   FirstBin;
extern double   ExponentBin;

// recipe flags 
extern int    ReionizationOn;
extern int    SupernovaRecipeOn;
extern int    DiskInstabilityOn;
extern int    AGNrecipeOn;
extern int    SFprescription;
extern int    H2prescription;
extern int    GasPrecessionOn;
extern int    RamPressureOn;
extern int    HotStripOn;
extern int    HeatedToCentral;
extern int    ReincorpotationModel;

// recipe parameters 
extern double RecycleFraction;
extern double Yield;
extern double FracZleaveDisk;
extern double ReIncorporationFactor;
extern double ThreshMajorMerger;
extern double BaryonFrac;
extern double SfrEfficiency;
extern double FeedbackReheatingEpsilon;
extern double FeedbackGasSigma;
extern double FeedbackExponent;
extern double FeedbackEjectionEfficiency;
extern double RadioModeEfficiency;
extern double QuasarModeEfficiency;
extern double BlackHoleGrowthRate;
extern double H2FractionFactor;
extern double H2FractionExponent;
extern double ClumpFactor;
extern double ClumpExponent;
extern double QTotMin;
extern double GasSinkRate;
extern double ThetaThresh;
extern double DegPerTdyn;
extern double Reionization_z0;
extern double Reionization_zr;
extern double ThresholdSatDisruption;
extern double AlphaBurst;
extern double BetaBurst;

extern double UnitLength_in_cm,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s,
  UnitMass_in_g,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears, 
  G,
  Hubble,
  a0, ar;

extern int    ListOutputSnaps[ABSOLUTEMAXSNAPS];

extern double ZZ[ABSOLUTEMAXSNAPS];
extern double AA[ABSOLUTEMAXSNAPS];
extern double Age[ABSOLUTEMAXSNAPS];

extern int    MAXSNAPS;
extern int    NOUT;
extern int    Snaplistlen;

extern gsl_rng *random_generator;

extern int TreeID;
extern int FileNum;

double DiscBinEdge[N_BINS+1];
int RetroCount, ProCount;

#ifdef MINIMIZE_IO
extern char *ptr_treedata, *ptr_galaxydata, *ptr_galsnapdata[ABSOLUTEMAXSNAPS];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[ABSOLUTEMAXSNAPS], maxstorage_galsnapdata[ABSOLUTEMAXSNAPS], filled_galsnapdata[ABSOLUTEMAXSNAPS];
#endif


#endif  // #ifndef ALLVARS_H
