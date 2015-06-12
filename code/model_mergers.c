#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double estimate_merging_time(int sat_halo, int mother_halo, int ngal)
{
  double coulomb, mergtime, SatelliteMass, SatelliteRadius;

  if(sat_halo == mother_halo) 
  {
    printf("\t\tSnapNum, Type, IDs, sat radius:\t%i\t%i\t%i\t%i\t--- sat/cent have the same ID\n", 
      Gal[ngal].SnapNum, Gal[ngal].Type, sat_halo, mother_halo);
    return -1.0;
  }
  
  coulomb = log(Halo[mother_halo].Len / ((double) Halo[sat_halo].Len) + 1);

  SatelliteMass = get_virial_mass(sat_halo) + Gal[ngal].StellarMass + Gal[ngal].ColdGas;
  SatelliteRadius = get_virial_radius(mother_halo);

  if(SatelliteMass > 0.0 && coulomb > 0.0)
    mergtime = 2.0 *
    1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halo) / (coulomb * G * SatelliteMass);
  else
    mergtime = -1.0;
  
  return mergtime;

}



void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step)
{
  double mi, ma, mass_ratio, central_bulge_fraction;
  double R1, R2, Eini1, Eini2, Eorb, Erad;
  double disc_mass_ratio[30];
  //int i;

  // calculate mass ratio of merging galaxies 
  if(Gal[p].StellarMass + Gal[p].ColdGas <
    Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas)
  {
    mi = Gal[p].StellarMass + Gal[p].ColdGas;
    ma = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
  }
  else
  {
    mi = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
    ma = Gal[p].StellarMass + Gal[p].ColdGas;
  }

  if(ma > 0)
    mass_ratio = mi / ma;
  else
    mass_ratio = 1.0;
	
  if(mass_ratio<0.0)
  {
	mass_ratio = 0.0;
	printf("Had to correct mass_ratio < 0.0");
  }
	
	if(Gal[merger_centralgal].StellarMass > 0.0)
		central_bulge_fraction = Gal[merger_centralgal].ClassicalBulgeMass / Gal[merger_centralgal].StellarMass;
	else
		central_bulge_fraction = 0.0;

  // pre-merger information needed to calculate the final classical bulge radius below
  if( mass_ratio > ThreshMajorMerger || central_bulge_fraction > 0.5)
  {
    if( central_bulge_fraction > 0.5)
      R1 = Gal[merger_centralgal].ClassicalBulgeRadius;
    else
      R1 = dmax(3.0 * Gal[merger_centralgal].DiskScaleRadius, Gal[merger_centralgal].ClassicalBulgeRadius);

    if(Gal[p].ClassicalBulgeMass / Gal[p].StellarMass > 0.5)
      R2 = Gal[p].ClassicalBulgeRadius;
    else
      R2 = dmax(3.0 * Gal[p].DiskScaleRadius, Gal[p].ClassicalBulgeRadius);
    
    if(R1 > 0.0)
      Eini1 = G * pow(Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas, 2.0) / R1;
    else 
      Eini1 = 0.0;
    
    if(R2 > 0.0)
      Eini2 = G * pow(Gal[p].Mvir + Gal[p].StellarMass + Gal[p].ColdGas, 2.0) / R2;
    else 
      Eini2 = 0.0;
    
    if(R1 + R2 > 0.0)
      Eorb = G * ma * mi / (R1 + R2);
    else
      Eorb = 0.0;
    
    if(ma + mi > 0.0)
      Erad = 2.75 * (Eini1 + Eini2) * (Gal[merger_centralgal].ColdGas + Gal[p].ColdGas) / (ma + mi);
    else
      Erad = 0.0;
  }

  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0 && Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[merger_centralgal].HotGas == Gal[merger_centralgal].HotGas);

  add_galaxies_together(merger_centralgal, p, disc_mass_ratio);

  // Initiate quasar accretion and feedback
  if(AGNrecipeOn>0)
    grow_black_hole(merger_centralgal, mass_ratio);
  
  collisional_starburst_recipe(disc_mass_ratio, merger_centralgal, centralgal, time, dt, halonr, 0, step, mass_ratio);

  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0 && Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[merger_centralgal].HotGas == Gal[merger_centralgal].HotGas);

  if(mass_ratio > ThreshMajorMerger)
  {
    make_bulge_from_burst(merger_centralgal);
    Gal[merger_centralgal].LastMajorMerger = time;
    Gal[p].mergeType = 2;  // Mark as major merger
  }
  else
    Gal[p].mergeType = 1;  // Mark as minor merger
	
  if(mass_ratio > ThreshMajorMerger || central_bulge_fraction > 0.5)
  {
	// Calculate the post-merger bulge radius
    Gal[merger_centralgal].ClassicalBulgeRadius = 
    G * pow(Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas, 2.0) / (Eini1 + Eini2 + Eorb + Erad);
  }
}



void grow_black_hole(int merger_centralgal, double mass_ratio)
{
  double BHaccrete, metallicity, accrete_ratio, DiscGasSum;
  int k;

  if(Gal[merger_centralgal].ColdGas > 0.0)
  {
    BHaccrete = BlackHoleGrowthRate * mass_ratio / 
      (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;

    // Cannot accrete more gas than is available
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;
  
	DiscGasSum = get_disc_gas(merger_centralgal);
	assert(DiscGasSum <= 1.01*Gal[merger_centralgal].ColdGas && DiscGasSum >= Gal[merger_centralgal].ColdGas/1.01);

	for(k=0; k<30; k++)
    {
	  accrete_ratio = Gal[merger_centralgal].DiscGas[k] / DiscGasSum;
	  metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	  Gal[merger_centralgal].DiscGas[k] -= BHaccrete * accrete_ratio;
	  Gal[merger_centralgal].DiscGasMetals[k] -= BHaccrete * accrete_ratio * metallicity;
	  Gal[merger_centralgal].ColdGas -= BHaccrete * accrete_ratio;
	  Gal[merger_centralgal].MetalsColdGas -= BHaccrete * accrete_ratio * metallicity;
	  Gal[merger_centralgal].BlackHoleMass += BHaccrete * accrete_ratio;
    }

    quasar_mode_wind(merger_centralgal, BHaccrete);
  }
}



void quasar_mode_wind(int gal, float BHaccrete)
{
  float quasar_energy, cold_gas_energy, hot_gas_energy;
  int k;
  
  // work out total energies in quasar wind (eta*m*c^2), cold and hot gas (1/2*m*Vvir^2)
  quasar_energy = QuasarModeEfficiency * 0.1 * BHaccrete * (C / UnitVelocity_in_cm_per_s) * (C / UnitVelocity_in_cm_per_s);
  cold_gas_energy = 0.5 * Gal[gal].ColdGas * Gal[gal].Vvir * Gal[gal].Vvir;
  hot_gas_energy = 0.5 * Gal[gal].HotGas * Gal[gal].Vvir * Gal[gal].Vvir;
   
  // compare quasar wind and cold gas energies and eject cold
  if(quasar_energy > cold_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].ColdGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsColdGas;
   
    Gal[gal].ColdGas = 0.0;
    Gal[gal].MetalsColdGas = 0.0;

	for(k=0; k<30; k++)
    {
	  Gal[gal].DiscGas[k] = 0.0;
	  Gal[gal].DiscGasMetals[k] = 0.0;
    }
  }
  
  // compare quasar wind and cold+hot gas energies and eject hot
  if(quasar_energy > cold_gas_energy + hot_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].HotGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsHotGas;
   
    Gal[gal].HotGas = 0.0;
    Gal[gal].MetalsHotGas = 0.0;
  }
}



void add_galaxies_together(int t, int p, double *disc_mass_ratio)
{
  int step, i;
  double DiscGasSum;
  
  for(i=0; i<30; i++)
  {
	if(Gal[t].DiscGas[i]>0.0 && Gal[p].DiscGas[i]>0.0)
	  disc_mass_ratio[i] = Gal[p].DiscGas[i] / Gal[t].DiscGas[i];
	else
	  disc_mass_ratio[i] = 0.0;
	
	Gal[t].DiscGas[i] += Gal[p].DiscGas[i];
	Gal[t].DiscGasMetals[i] += Gal[p].DiscGasMetals[i];
	Gal[t].ColdGas += Gal[p].DiscGas[i];
	Gal[t].MetalsColdGas += Gal[p].DiscGasMetals[i];
  }

  Gal[t].StellarMass += Gal[p].StellarMass;
  Gal[t].MetalsStellarMass += Gal[p].MetalsStellarMass;

  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;
  
  Gal[t].EjectedMass += Gal[p].EjectedMass;
  Gal[t].MetalsEjectedMass += Gal[p].MetalsEjectedMass;
  
  Gal[t].ICS += Gal[p].ICS;
  Gal[t].MetalsICS += Gal[p].MetalsICS;

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;

  // Add merger to bulge
  if(Gal[t].StellarMass > 0.0 && Gal[t].ClassicalBulgeMass / Gal[t].StellarMass > 0.5)
  {
	Gal[t].ClassicalBulgeMass += Gal[p].StellarMass;
	Gal[t].ClassicalMetalsBulgeMass += Gal[p].MetalsStellarMass;		
  }
  else
  {
	Gal[t].SecularBulgeMass += Gal[p].StellarMass;
	Gal[t].SecularMetalsBulgeMass += Gal[p].MetalsStellarMass;				
  }

  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulge[step] += Gal[p].SfrDisk[step] + Gal[p].SfrBulge[step];
    Gal[t].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step] + Gal[p].SfrBulgeColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step] + Gal[p].SfrBulgeColdGasMetals[step];
  }

  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
}



void make_bulge_from_burst(int p)
{
  int step, i;
  
  // generate bulge 
  Gal[p].ClassicalBulgeMass = Gal[p].StellarMass;
  Gal[p].ClassicalMetalsBulgeMass = Gal[p].MetalsStellarMass;
  
  Gal[p].SecularBulgeMass = 0.0;
  Gal[p].SecularMetalsBulgeMass = 0.0;

  // Remove stars from the disc shells
  for(i=0; i<30; i++)
  {
	Gal[p].DiscStars[i] = 0.0;
	Gal[p].DiscStarsMetals[i] = 0.0;
  }  

  // update the star formation rate 
  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrBulge[step] += Gal[p].SfrDisk[step];
    Gal[p].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step];
    Gal[p].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step];
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
  }
}



void collisional_starburst_recipe(double disc_mass_ratio[30], int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step, double mass_ratio)
{
 double stars, reheated_mass, ejected_mass, fac, metallicity, CentralVvir, eburst, Sigma_0gas, area, stars_sum;
 double NewStars[30], NewStarsMetals[30];
 int k;

 // This is the major and minor merger starburst recipe of Somerville et al. 2001. 
 // The coefficients in eburst are taken from TJ Cox's PhD thesis and should be more 
 // accurate then previous. 

 stars_sum = 0.0;

 if(Gal[merger_centralgal].ColdGas>0)
 {
  CentralVvir = Gal[centralgal].Vvir;

  // update the star formation rate 
  Gal[merger_centralgal].SfrBulgeColdGas[step] += Gal[merger_centralgal].ColdGas;
  Gal[merger_centralgal].SfrBulgeColdGasMetals[step] += Gal[merger_centralgal].MetalsColdGas;

  for(k=0; k<30; k++)
  {
	// the bursting fraction 
    if(mode == 1)
      eburst = disc_mass_ratio[k];
    else
      eburst = 0.56 * pow(disc_mass_ratio[k], 0.7);

    stars = eburst * Gal[merger_centralgal].DiscGas[k];
    if(stars < 0.0)
      stars = 0.0;

	if(stars > Gal[merger_centralgal].DiscGas[k])
      stars = Gal[merger_centralgal].DiscGas[k];
	
	// this bursting results in SN feedback on the cold/hot gas 
    if(SupernovaRecipeOn == 1 && Gal[merger_centralgal].DiscGas[k] > 0.0 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
		if(Gal[merger_centralgal].Vvir > 0.0)
		  area = M_PI * (pow(DiscBinEdge[k+1]/Gal[merger_centralgal].Vvir, 2.0) - pow(DiscBinEdge[k]/Gal[merger_centralgal].Vvir, 2.0));
		else
		  area = M_PI * (pow(DiscBinEdge[k+1]/Gal[merger_centralgal].Vmax, 2.0) - pow(DiscBinEdge[k]/Gal[merger_centralgal].Vmax, 2.0));
		Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
        reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[merger_centralgal].DiscGas[k]/area/1.3);
		
		// can't use more cold gas than is available! so balance SF and feedback 
	    if((stars + reheated_mass) > Gal[merger_centralgal].DiscGas[k] && (stars + reheated_mass) > 0.0)
	    {
	      fac = Gal[merger_centralgal].DiscGas[k] / (stars + reheated_mass);
	      stars *= fac;
	      reheated_mass *= fac;
	    }
	
	    if(stars<1e-8)
	    {
		  stars = 1e-8;
		  reheated_mass = Gal[merger_centralgal].DiscGas[k] - (1-RecycleFraction)*stars;
	    }
	
	    ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (CentralVvir * CentralVvir) - FeedbackReheatingEpsilon) * stars;
	    if(ejected_mass < 0.0)
	        ejected_mass = 0.0;
	
		assert(RecycleFraction*stars+reheated_mass <= 1.001*Gal[merger_centralgal].DiscGas[k]);
	  }

	  else
	  {
		reheated_mass = RecycleFraction * stars;
		ejected_mass = 0.0;
	  }
	}  
    else
	{
	  stars=0.0;
      reheated_mass = 0.0;
	  ejected_mass = 0.0;
	}
	if(reheated_mass!=reheated_mass || reheated_mass<0.0)
		printf("reheated_mass, stars, fac, DiscGas -- %e\t%e\t%e\t%e\n", reheated_mass, stars, fac, Gal[merger_centralgal].DiscGas[k]);	
	assert(reheated_mass >= 0.0);

	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	NewStars[k] = (1 - RecycleFraction) * stars;
	NewStarsMetals[k] = (1 - RecycleFraction) * metallicity * stars;
    update_from_star_formation(merger_centralgal, stars, metallicity, k);

    if(reheated_mass > Gal[merger_centralgal].DiscGas[k] && reheated_mass < 1.01*Gal[merger_centralgal].DiscGas[k])
	  reheated_mass = Gal[merger_centralgal].DiscGas[k];

	// update from feedback
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	update_from_feedback(merger_centralgal, centralgal, reheated_mass, ejected_mass, metallicity, k);
 
    // Inject new metals from SN II
	if(SupernovaRecipeOn == 1 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
	    Gal[merger_centralgal].DiscGasMetals[k] += Yield * stars;
	    Gal[merger_centralgal].MetalsColdGas += Yield * stars;
  	  }
	  else
		Gal[merger_centralgal].MetalsHotGas += Yield * stars;
	}
	assert(Gal[merger_centralgal].DiscGasMetals[k]<=Gal[merger_centralgal].DiscGas[k]);
	
	stars_sum += stars;
  }

  // Sum stellar discs together
  combine_stellar_discs(merger_centralgal, NewStars, NewStarsMetals);

  Gal[merger_centralgal].SfrBulge[step] += stars_sum / dt;

  // check for disk instability
  if(DiskInstabilityOn && mode == 0)
    if(mass_ratio < ThreshMajorMerger)
      check_disk_instability(merger_centralgal, centralgal, halonr, time, dt, step);
 }
}



void disrupt_satellite_to_ICS(int centralgal, int gal)
{  
  Gal[centralgal].HotGas += Gal[gal].ColdGas + Gal[gal].HotGas;
  Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas + Gal[gal].MetalsHotGas;
  
  Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
  Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
  
  Gal[centralgal].ICS += Gal[gal].ICS;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsICS;

  Gal[centralgal].ICS += Gal[gal].StellarMass;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsStellarMass;
  
  // what should we do with the disrupted satellite BH?
  Gal[gal].mergeType = 4;  // mark as disruption to the ICS
}