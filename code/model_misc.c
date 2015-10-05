#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "core_allvars.h"
#include "core_proto.h"



void init_galaxy(int p, int halonr)
{
  int j, step;

  if(halonr != Halo[halonr].FirstHaloInFOFgroup)
  {
    printf("Hah?\n");
    ABORT(1);
  }

  Gal[p].Type = 0;

  Gal[p].GalaxyNr = GalaxyCounter;
  GalaxyCounter++;
  
  Gal[p].HaloNr = halonr;
  Gal[p].MostBoundID = Halo[halonr].MostBoundID;
  Gal[p].SnapNum = Halo[halonr].SnapNum - 1;

  Gal[p].mergeType = 0;
  Gal[p].mergeIntoID = -1;
  Gal[p].mergeIntoSnapNum = -1;
  Gal[p].dT = -1.0;

  for(j = 0; j < 3; j++)
  {
    Gal[p].Pos[j] = Halo[halonr].Pos[j];
    Gal[p].Vel[j] = Halo[halonr].Vel[j];
	Gal[p].SpinStars[j] = Halo[halonr].Spin[j] / pow(pow(Halo[halonr].Spin[0], 2.0) + pow(Halo[halonr].Spin[1], 2.0) + pow(Halo[halonr].Spin[2], 2.0), 0.5);
	Gal[p].SpinGas[j] = Halo[halonr].Spin[j] / pow(pow(Halo[halonr].Spin[0], 2.0) + pow(Halo[halonr].Spin[1], 2.0) + pow(Halo[halonr].Spin[2], 2.0), 0.5);
  }

  Gal[p].Len = Halo[halonr].Len;
  Gal[p].Vmax = Halo[halonr].Vmax;
  Gal[p].Vvir = get_virial_velocity(halonr);
  Gal[p].Mvir = get_virial_mass(halonr);
  Gal[p].Rvir = get_virial_radius(halonr);

  Gal[p].deltaMvir = 0.0;

  Gal[p].ColdGas = 0.0;
  Gal[p].StellarMass = 0.0;
  Gal[p].ClassicalBulgeMass = 0.0;
  Gal[p].SecularBulgeMass = 0.0;
  Gal[p].HotGas = 0.0;
  Gal[p].EjectedMass = 0.0;
  Gal[p].BlackHoleMass = 0.0;
  Gal[p].ICS = 0.0;
    
  Gal[p].StarsInSitu = 0.0;
  Gal[p].StarsInstability = 0.0;
  Gal[p].StarsMergeBurst = 0.0;

  Gal[p].MetalsColdGas = 0.0;
  Gal[p].MetalsStellarMass = 0.0;
  Gal[p].ClassicalMetalsBulgeMass = 0.0;
  Gal[p].SecularMetalsBulgeMass = 0.0;
  Gal[p].MetalsHotGas = 0.0;
  Gal[p].MetalsEjectedMass = 0.0;
  Gal[p].MetalsICS = 0.0;
  
  for(j=0; j<30; j++)
  {
	Gal[p].DiscGas[j] = 0.0;
	Gal[p].DiscStars[j] = 0.0;
	Gal[p].DiscGasMetals[j] = 0.0;
	Gal[p].DiscStarsMetals[j] = 0.0;
  }

  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrBulge[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
    Gal[p].SfrBulgeColdGas[step] = 0.0;
    Gal[p].SfrBulgeColdGasMetals[step] = 0.0;
  }

  Gal[p].DiskScaleRadius = get_disk_radius(halonr, p);
  Gal[p].ClassicalBulgeRadius = 0.0;
  Gal[p].MergTime = 999.9;
  Gal[p].Cooling = 0.0;
  Gal[p].Heating = 0.0;
  Gal[p].r_heat = 0.0;
  Gal[p].LastMajorMerger = -1.0;
  Gal[p].OutflowRate = 0.0;
  Gal[p].TotalSatelliteBaryons = 0.0;

  Gal[p].infallMvir = -1.0;  //infall properties
  Gal[p].infallVvir = -1.0;
  Gal[p].infallVmax = -1.0;
  
}



double get_disk_radius(int halonr, int p)
{
  // See Mo, Shude & White (1998) eq12, and using a Bullock style lambda.
  double SpinMagnitude, SpinParameter;
  
	SpinMagnitude = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] + 
		Halo[halonr].Spin[1] * Halo[halonr].Spin[1] + Halo[halonr].Spin[2] * Halo[halonr].Spin[2]);
  
  // trim the extreme tail of the spin distribution for more a realistic r_s
  //if(SpinMagnitude > 1.5) SpinMagnitude = 1.5;
  
  SpinParameter = SpinMagnitude / (1.414 * Gal[p].Vvir * Gal[p].Rvir);
    
  return (SpinParameter / 1.414) * Gal[p].Rvir;

}



double get_metallicity(double gas, double metals)
{
  double metallicity;

  if(metals>gas)
	printf("get_metallicity report: metals, gas/stars = %e, %e\n", metals, gas);
  //assert(gas>metals);

  if(gas > 0.0 && metals > 0.0)
  {
    metallicity = metals / gas;
    if(metallicity < 1.0)
      return metallicity;
    else
      return 1.0;
  }
  else
    return 0.0;

}



double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}



double get_virial_mass(int halonr)
{
  if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir >= 0.0)
    return Halo[halonr].Mvir;   /* take spherical overdensity mass estimate */ 
  else
    return Halo[halonr].Len * PartMass;
}



double get_virial_velocity(int halonr)
{
	double Rvir;
	
	Rvir = get_virial_radius(halonr);
	
  if(Rvir > 0.0)
		return sqrt(G * get_virial_mass(halonr) / Rvir);
	else
		return 0.0;
}



double get_virial_radius(int halonr)
{
  // return Halo[halonr].Rvir;  // Used for Bolshoi

  double zplus1, hubble_of_z_sq, rhocrit, fac;
  
  zplus1 = 1 + ZZ[Halo[halonr].SnapNum];
  hubble_of_z_sq =
    Hubble * Hubble *(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
    OmegaLambda);
  
  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
  
  return cbrt(get_virial_mass(halonr) * fac);
}


double get_disc_gas(int p)
{
	double DiscGasSum, DiscMetalsSum;
	int l;
	
	DiscGasSum = DiscMetalsSum = 0.0;
	for(l=0; l<30; l++)
	{
		DiscGasSum += Gal[p].DiscGas[l];
		DiscMetalsSum += Gal[p].DiscGasMetals[l];
	}
    
    if((Gal[p].ColdGas < 1e-15 && Gal[p].ColdGas!=0.0) || (DiscGasSum < 1e-15 && DiscGasSum!=0.0))
    {
        printf("get_disc_gas initial DiscGasSum, ColdGas = %e, %e\n", DiscGasSum, Gal[p].ColdGas);
        Gal[p].ColdGas = 0.0;
        Gal[p].MetalsColdGas = 0.0;
        for(l=0; l<30; l++)
        {
            Gal[p].DiscGas[l] = 0.0;
            Gal[p].DiscGasMetals[l] = 0.0;
        }
        DiscGasSum = 0.0;
    }

	if(DiscGasSum>1.001*Gal[p].ColdGas || DiscGasSum<Gal[p].ColdGas/1.001)
	{
		printf("get_disc_gas report ... DiscSum, ColdGas =  %e, %e\n", DiscGasSum, Gal[p].ColdGas);
		printf("get_disc_gas report ... MetalsSum, ColdMetals =  %e, %e\n", DiscMetalsSum, Gal[p].MetalsColdGas);
		
		if(Gal[p].ColdGas<=1e-15)
		{
			for(l=0; l<30; l++)
				Gal[p].DiscGas[l] = 0.0; // Sometimes a tiny non-zero difference can creep in (probably due to projecting discs).  This just takes care of that.
			DiscGasSum = 0.0;
			Gal[p].ColdGas = 0.0;
		}
        
        if(Gal[p].ColdGas<=1e-15 || Gal[p].MetalsColdGas<=0.0)
        {
            for(l=0; l<30; l++)
                Gal[p].DiscGasMetals[l] = 0.0;
            DiscMetalsSum = 0.0;
            Gal[p].MetalsColdGas = 0.0;
        }
		
		if((DiscGasSum<1.01*Gal[p].ColdGas && DiscGasSum>Gal[p].ColdGas/1.01) || (DiscGasSum==0.0 && Gal[p].ColdGas<1e-10) || DiscGasSum < 1e-10)
		{
			Gal[p].ColdGas = DiscGasSum; // If difference is small, just set the numbers to be the same to prevent small errors from blowing up
			Gal[p].MetalsColdGas = DiscMetalsSum;
		}
		
  	}
	return DiscGasSum;
}

double get_disc_stars(int p)
{
    double DiscStarSum, DiscAndBulge;
    int l;
    
    DiscStarSum = 0.0;
    for(l=0; l<30; l++)
        DiscStarSum += Gal[p].DiscStars[l];
    
    DiscAndBulge = DiscStarSum + Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
    
    if(DiscAndBulge>1.001*Gal[p].StellarMass || DiscAndBulge<Gal[p].StellarMass/1.001)
    {
        printf("get_disc_stars report %e\t%e\n", DiscAndBulge, Gal[p].StellarMass);
        if(DiscAndBulge<1.01*Gal[p].StellarMass && DiscAndBulge>Gal[p].StellarMass/1.01)
            Gal[p].StellarMass = DiscAndBulge; // If difference is small, just set the numbers to be the same to prevent small errors from blowing up
        if(Gal[p].StellarMass <= Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass)
        {
            for(l=0; l<30; l++)
                Gal[p].DiscStars[l] = 0.0;
            DiscStarSum = 0.0;
        }
    }
    return DiscStarSum;
}

double get_disc_ang_mom(int p, int type)
{
	// type=0 for gas, type=1 for stars
	double J_sum;
	int l;
	
	J_sum = 0.0;
	if(type==0)
	{
		for(l=0; l<30; l++)
			J_sum += Gal[p].DiscGas[l] * pow((pow(DiscBinEdge[l],2.0) + pow(DiscBinEdge[l+1],2.0))/2.0, 0.5);
	}
	else if(type==1)
	{
		for(l=0; l<30; l++)
			J_sum += Gal[p].DiscStars[l] * pow((pow(DiscBinEdge[l],2.0) + pow(DiscBinEdge[l+1],2.0))/2.0, 0.5);
	}
	
	return J_sum;
}