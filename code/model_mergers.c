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

  SatelliteMass = get_virial_mass(sat_halo, ngal) + Gal[ngal].StellarMass + Gal[ngal].ColdGas;
  SatelliteRadius = get_virial_radius(mother_halo, ngal); // Feeding ngal here is right, but it's better than nothing.  Frankly, why doesn't this use coords to get the satellite radius anyway?

  if(SatelliteMass > 0.0 && coulomb > 0.0)
    mergtime = 2.0 *
    1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halo, ngal) / (coulomb * G * SatelliteMass);
  else
    mergtime = -1.0;
  
  return mergtime;

}



void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step)
{
  double mi, ma, mass_ratio, central_bulge_fraction;
  double R1, R2, Eini1, Eini2, Eorb, Erad;
  double disc_mass_ratio[N_BINS], PostRetroGas[N_BINS];
  double DiscGasSum, DiscStarsSum, StarChannelSum;
  int i, s;

    // These lines are really just to ensure everything is in order before merging
    DiscStarsSum = get_disc_stars(merger_centralgal);
    DiscStarsSum = get_disc_stars(p);
    
	for(i=0; i<N_BINS; i++) 
	{
//		if (Gal[merger_centralgal].DiscStarsMetals[i] > Gal[merger_centralgal].DiscStars[i])
//			printf("DiscStars, Metals = %e, %e\n", Gal[merger_centralgal].DiscStars[i], Gal[merger_centralgal].DiscStarsMetals[i]);
//		assert(Gal[merger_centralgal].DiscStarsMetals[i] <= Gal[merger_centralgal].DiscStars[i]);
//		
//		if (Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i])
//			printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
//		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
        
        disc_mass_ratio[i] = 0.0;
	}

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

 // GIVEN THIS ISN'T ACTUALLY USED IN THE MODEL, I CAN PROBABLY DELETE IT
  // pre-merger information needed to calculate the final classical bulge radius below
  if( mass_ratio > ThreshMajorMerger || central_bulge_fraction > 0.5)
  {
    if( central_bulge_fraction > 0.5 && Gal[merger_centralgal].ClassicalBulgeRadius > 0.0)
      R1 = Gal[merger_centralgal].ClassicalBulgeRadius;
    else
      R1 = dmax(3.0 * Gal[merger_centralgal].DiskScaleRadius, Gal[merger_centralgal].ClassicalBulgeRadius);

    if(Gal[p].ClassicalBulgeMass / Gal[p].StellarMass > 0.5 && Gal[p].ClassicalBulgeRadius > 0.0)
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

  DiscGasSum = get_disc_gas(merger_centralgal);
  assert(DiscGasSum <= 1.01*Gal[merger_centralgal].ColdGas && DiscGasSum >= Gal[merger_centralgal].ColdGas/1.01);
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);

  add_galaxies_together(merger_centralgal, p, mass_ratio, disc_mass_ratio, centralgal, dt, PostRetroGas);
  
    for(i=0; i<N_BINS; i++) assert(disc_mass_ratio[i] <= 1.0 && disc_mass_ratio[i]>=0.0);
    
  DiscGasSum = get_disc_gas(merger_centralgal);
  assert(DiscGasSum <= 1.01*Gal[merger_centralgal].ColdGas && DiscGasSum >= Gal[merger_centralgal].ColdGas/1.01);
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
  for(i=N_BINS-1; i>=0; i--) assert(Gal[merger_centralgal].DiscGasMetals[i] <= Gal[merger_centralgal].DiscGas[i]);

    for(i=0; i<N_BINS; i++) assert(disc_mass_ratio[i] <= 1.0);

  collisional_starburst_recipe(disc_mass_ratio, merger_centralgal, centralgal, time, dt, halonr, 0, step, mass_ratio);
  for(i=N_BINS-1; i>=0; i--) assert(Gal[merger_centralgal].DiscGasMetals[i] <= Gal[merger_centralgal].DiscGas[i]);

  double BHaccrete = grow_black_hole(merger_centralgal, disc_mass_ratio);
  for(i=N_BINS-1; i>=0; i--) assert(Gal[merger_centralgal].DiscGasMetals[i] <= Gal[merger_centralgal].DiscGas[i]);

  if(AGNrecipeOn>0)
	quasar_mode_wind(p, BHaccrete, centralgal);

  DiscGasSum = get_disc_gas(merger_centralgal);
  assert(DiscGasSum <= 1.01*Gal[merger_centralgal].ColdGas && DiscGasSum >= Gal[merger_centralgal].ColdGas/1.01);

  // Check whether any retrograde gas is left over
  double unstable_gas, metallicity, stars, net_stars, j_bin;
  for(i=N_BINS-1; i>=0; i--)
  {
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
	assert(Gal[merger_centralgal].DiscGasMetals[i] <= Gal[merger_centralgal].DiscGas[i]);
	
	if(PostRetroGas[i] < 0.99*Gal[merger_centralgal].DiscGas[i] && Gal[merger_centralgal].DiscGas[i]-PostRetroGas[i] > 1e-10)
	{
		unstable_gas = Gal[merger_centralgal].DiscGas[i] - PostRetroGas[i];
        
        if(!(Gal[merger_centralgal].DiscRadii[i]!=INFINITY && Gal[merger_centralgal].DiscRadii[i]==Gal[merger_centralgal].DiscRadii[i]))
        {
            printf("i = %d\n", i);
            printf("Mvir = %e\n", Gal[merger_centralgal].Mvir);
            printf("radius before = %e\n", Gal[merger_centralgal].DiscRadii[i]);
            update_disc_radii(merger_centralgal);
            printf("radius after = %e\n", Gal[merger_centralgal].DiscRadii[i]);
        }
        
        assert(Gal[merger_centralgal].DiscRadii[i]!=INFINITY && Gal[merger_centralgal].DiscRadii[i]==Gal[merger_centralgal].DiscRadii[i]);
        stars = deal_with_unstable_gas(unstable_gas, merger_centralgal, i, Gal[merger_centralgal].Vvir, metallicity, centralgal, 0, Gal[merger_centralgal].DiscRadii[i], Gal[merger_centralgal].DiscRadii[i+1]);
        
        if(stars>=MIN_STARS_FOR_SN)
            net_stars = (1 - RecycleFraction) * stars;
        else
            net_stars = stars;
        
        Gal[merger_centralgal].StellarMass += net_stars;
        Gal[merger_centralgal].MetalsStellarMass += metallicity * net_stars;
        Gal[merger_centralgal].StarsMergeBurst += net_stars;
        Gal[merger_centralgal].SfrBulge[step] += stars / dt;
        StarChannelSum = get_channel_stars(merger_centralgal);
        assert(Gal[merger_centralgal].StellarMass >= (Gal[merger_centralgal].StarsInSitu+Gal[merger_centralgal].StarsInstability+Gal[merger_centralgal].StarsMergeBurst)/1.01 && Gal[merger_centralgal].StellarMass <= (Gal[merger_centralgal].StarsInSitu+Gal[merger_centralgal].StarsInstability+Gal[merger_centralgal].StarsMergeBurst)*1.01);
        
        // Add the new stars from the retrograde starburst to the classical bulge
        j_bin = (DiscBinEdge[i]+DiscBinEdge[i+1])/2;
        if(net_stars>0)
        {
            double spinmagpre = pow(pow(Gal[merger_centralgal].SpinClassicalBulge[0],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[1],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[2],2.0), 0.5);
            for(s=0; s<3; s++)
            {
                assert(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY);
                // Need to know magnitude of j of newstars, not just direction
                Gal[merger_centralgal].SpinClassicalBulge[s] = (Gal[merger_centralgal].SpinClassicalBulge[s]*Gal[merger_centralgal].ClassicalBulgeMass + Gal[merger_centralgal].SpinGas[s]*net_stars*j_bin)/(Gal[merger_centralgal].ClassicalBulgeMass+net_stars);
                if(!(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY))
                {
                    printf("205: Reset classical bulge spin to zero\n");
                    Gal[merger_centralgal].SpinClassicalBulge[s] = 0.0;
                }
                assert(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY);
            }
            double spinmag = pow(pow(Gal[merger_centralgal].SpinClassicalBulge[0],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[1],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[2],2.0), 0.5);
//            if(spinmag<1.01 && spinmag>0.99 && !(spinmagpre<1.01 && spinmagpre>0.99))
//                printf("\nspinmag %e at 211 with i=%d\n", spinmag, i);
            
            Gal[merger_centralgal].ClassicalBulgeMass += net_stars;
            Gal[merger_centralgal].ClassicalMetalsBulgeMass += metallicity * net_stars;
        }

	}
  }

  if(mass_ratio > ThreshMajorMerger)
  {
    stars_to_bulge(merger_centralgal, p);
    Gal[merger_centralgal].LastMajorMerger = time;
    Gal[p].mergeType = 2;  // Mark as major merger
  }
  else
    Gal[p].mergeType = 1;  // Mark as minor merger


  if(DiskInstabilityOn>0)
  	check_disk_instability(merger_centralgal, centralgal, time, dt, step);

  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0 && Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[merger_centralgal].HotGas == Gal[merger_centralgal].HotGas);
  DiscGasSum = get_disc_gas(merger_centralgal);
  assert(DiscGasSum <= 1.01*Gal[merger_centralgal].ColdGas && DiscGasSum >= Gal[merger_centralgal].ColdGas/1.01);
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
  if(mass_ratio > ThreshMajorMerger || central_bulge_fraction > 0.5 && ma>0.0)
  {
	// Calculate the post-merger bulge radius
    double Efinal = Eini1 + Eini2 + Eorb + Erad;
    if(Efinal>0.0)
        Gal[merger_centralgal].ClassicalBulgeRadius = G * pow(Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas, 2.0) / Efinal;
    else
        Gal[merger_centralgal].ClassicalBulgeRadius = 0.0;
      
      
    double BR = Gal[merger_centralgal].ClassicalBulgeRadius;
      
    if(Gal[merger_centralgal].ClassicalBulgeMass>0.0 && BR <= 0.0)
    {
        printf("BulgeMass, BulgeRadius = %e, %e\n", Gal[merger_centralgal].ClassicalBulgeMass, BR);
        printf("Efinal, centralfraction = %e, %e\n", Efinal, central_bulge_fraction);
        printf("Disk radii central, sat = %e, %e\n", Gal[merger_centralgal].DiskScaleRadius, Gal[p].DiskScaleRadius);
        //ABORT(0);
    }
      
      if(BR!=BR || BR==INFINITY)
      {
          printf("BulgeRadius = %e\n", BR);
          printf("SM, CGM = %e, %e\n", Gal[merger_centralgal].StellarMass, Gal[merger_centralgal].ColdGas);
          printf("Eini1, Eini2, Eorb, Erad = %e, %e, %e, %e\n", Eini1, Eini2, Eorb, Erad);
          //ABORT(0);
      }
  }
}



double grow_black_hole(int merger_centralgal, double* disc_mass_ratio)
{
  double BHaccrete, BHaccrete_tot, metallicity;//, accrete_ratio, DiscGasSum;
  int i;

  BHaccrete_tot = 0.0;

  for(i=0; i<N_BINS; i++)
  {
	if(Gal[merger_centralgal].DiscGas[i] > 0.0)
	{
		BHaccrete = BlackHoleGrowthRate * disc_mass_ratio[i] / (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].DiscGas[i];
        if(disc_mass_ratio[i] > 1.0 || disc_mass_ratio[i]!=disc_mass_ratio[i] || disc_mass_ratio[i]<0.0) printf("i, disc_mass_ratio[i] = %d, %e\n", i, disc_mass_ratio[i]);
		assert(disc_mass_ratio[i]<=1.0);
        
        if(BHaccrete<0.0)
        {
            printf("BHaccrete, disc_mass_ratio, DiscGas, Vvir, BlackHoleGrowthRate = %e, %e, %e, %e, %e\n", BHaccrete, disc_mass_ratio[i], Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].Vvir, BlackHoleGrowthRate);
        }
        
        assert(BHaccrete>=0.0);
		if(BHaccrete > Gal[merger_centralgal].DiscGas[i]) // This could only be possible if BlackHoleGrowthRate is set to >1.0, which shouldn't happen...
		{
			BHaccrete_tot += Gal[merger_centralgal].DiscGas[i];
			Gal[merger_centralgal].ColdGas -= Gal[merger_centralgal].DiscGas[i];
			Gal[merger_centralgal].MetalsColdGas -= Gal[merger_centralgal].DiscGasMetals[i];
			Gal[merger_centralgal].DiscGas[i] = 0.0;
			Gal[merger_centralgal].DiscGasMetals[i] = 0.0;
		}
		else
		{
			BHaccrete_tot += BHaccrete;
			metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
			Gal[merger_centralgal].DiscGas[i] -= BHaccrete;
			Gal[merger_centralgal].DiscGasMetals[i] -= BHaccrete * metallicity;
			if(Gal[merger_centralgal].DiscGasMetals[i]<0.0) Gal[merger_centralgal].DiscGasMetals[i] = 0.0;
			Gal[merger_centralgal].ColdGas -= BHaccrete;
			Gal[merger_centralgal].MetalsColdGas -= BHaccrete * metallicity;
			
		}
	}
  }

  Gal[merger_centralgal].BlackHoleMass += BHaccrete_tot;
    assert(Gal[merger_centralgal].BlackHoleMass>=0.0);
  return BHaccrete_tot;
}



void quasar_mode_wind(int p, float BHaccrete, int centralgal)
{ // I should probably out through the centralgal ID here
  double quasar_energy, cold_gas_energy, hot_gas_energy, DiscGasSum, cold_gas_energy_tot;
  int k;

  check_ejected(p);
  assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass); // Really should be centralgal I'm considering...

  // work out total energies in quasar wind (eta*m*c^2), cold and hot gas (1/2*m*Vvir^2)
  quasar_energy = QuasarModeEfficiency * 0.1 * BHaccrete * (C / UnitVelocity_in_cm_per_s) * (C / UnitVelocity_in_cm_per_s);
  cold_gas_energy_tot = 0.5 * Gal[p].ColdGas * Gal[p].Vvir * Gal[p].Vvir;

  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);

  if(quasar_energy > cold_gas_energy_tot)
  {
    Gal[p].EjectedQuasarGasMass += Gal[p].ColdGas;
      
    if(HeatedToCentral)
    {
        Gal[centralgal].EjectedMass += Gal[p].ColdGas;
        Gal[centralgal].MetalsEjectedMass += Gal[p].MetalsColdGas;
    }
    else
    {
        Gal[p].EjectedMass += Gal[p].ColdGas;
        Gal[p].MetalsEjectedMass += Gal[p].MetalsColdGas;
    }
      
	Gal[p].ColdGas = 0.0;
	Gal[p].MetalsColdGas = 0.0;
	
	for(k=0; k<N_BINS; k++)
	{
		Gal[p].DiscGas[k] = 0.0;
		Gal[p].DiscGasMetals[k] = 0.0;
	}
	
	quasar_energy -= cold_gas_energy_tot;
  }
  else
  {
	for(k=0; k<N_BINS; k++)
	{
		cold_gas_energy = 0.5 * Gal[p].DiscGas[k] * Gal[p].Vvir * Gal[p].Vvir;
		if(quasar_energy >= cold_gas_energy && cold_gas_energy > 0.0)
		{
            Gal[p].EjectedQuasarGasMass += Gal[p].DiscGas[k];
            if(HeatedToCentral)
            {
                Gal[centralgal].EjectedMass += Gal[p].DiscGas[k];
                Gal[centralgal].MetalsEjectedMass += Gal[p].DiscGasMetals[k];
            }
            else
            {
                Gal[p].EjectedMass += Gal[p].DiscGas[k];
                Gal[p].MetalsEjectedMass += Gal[p].DiscGasMetals[k];
            }
			Gal[p].ColdGas -= Gal[p].DiscGas[k];
			Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[k];
			Gal[p].DiscGas[k] = 0.0;
			Gal[p].DiscGasMetals[k] = 0.0;
			quasar_energy -= cold_gas_energy;
		}
		else if(quasar_energy > 0.0 && cold_gas_energy > 0.0)
		{
            Gal[p].EjectedQuasarGasMass += Gal[p].DiscGas[k] * quasar_energy/cold_gas_energy;
            if(HeatedToCentral)
            {
                Gal[centralgal].EjectedMass += Gal[p].DiscGas[k] * quasar_energy/cold_gas_energy;
                Gal[centralgal].MetalsEjectedMass += Gal[p].DiscGasMetals[k] * quasar_energy/cold_gas_energy;
            }
            else
            {
                Gal[p].EjectedMass += Gal[p].DiscGas[k] * quasar_energy/cold_gas_energy;
                Gal[p].MetalsEjectedMass += Gal[p].DiscGasMetals[k] * quasar_energy/cold_gas_energy;
            }
			Gal[p].ColdGas -= Gal[p].DiscGas[k] * quasar_energy/cold_gas_energy;
			Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[k] * quasar_energy/cold_gas_energy;
			Gal[p].DiscGas[k] *= (1 - quasar_energy/cold_gas_energy);
			Gal[p].DiscGasMetals[k] *= (1 - quasar_energy/cold_gas_energy);
			assert(Gal[p].DiscGasMetals[k] <= Gal[p].DiscGas[k]);
			quasar_energy = 0.0;
			break;
		}
		else if(cold_gas_energy > 0.0 || quasar_energy < 0.0)
		{
			quasar_energy = 0.0;
			break;
		}
	}

	DiscGasSum = get_disc_gas(p);
	assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
  }

   
  hot_gas_energy = 0.5 * Gal[p].HotGas * Gal[p].Vvir * Gal[p].Vvir;
  
  // compare quasar wind and cold+hot gas energies and eject hot
  if(quasar_energy > hot_gas_energy)
  {
      if(HeatedToCentral)
      {
          Gal[centralgal].EjectedMass += Gal[p].HotGas;
          Gal[centralgal].MetalsEjectedMass += Gal[p].MetalsHotGas;
      }
      else
      {
          Gal[p].EjectedMass += Gal[p].HotGas;
          Gal[p].MetalsEjectedMass += Gal[p].MetalsHotGas;
      }
   
    Gal[p].HotGas = 0.0;
    Gal[p].MetalsHotGas = 0.0;
  }
  check_ejected(p);
  assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
}



void add_galaxies_together(int t, int p, double mass_ratio, double *disc_mass_ratio, int centralgal, double dt, double *PostRetroGas)
{
  int step, i, s;
  double DiscGasSum, CentralGasOrig, ExpFac, StarChannelSum;
    ExpFac = AA[Gal[t].SnapNum]; // Expansion factor needed for determining physical distances for calculating j
    
	CentralGasOrig = get_disc_gas(t);
	assert(CentralGasOrig <= 1.01*Gal[t].ColdGas && CentralGasOrig >= Gal[t].ColdGas/1.01);

  	DiscGasSum = get_disc_gas(p);
	assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
	assert(Gal[t].ColdGas >= Gal[t].MetalsColdGas);
	assert(Gal[t].StellarMass >= Gal[t].MetalsStellarMass);
	
	for(i=0; i<N_BINS; i++) 
	{
		if(Gal[t].DiscStarsMetals[i] > Gal[t].DiscStars[i])
			printf("DiscStars, Metals = %e, %e\n", Gal[t].DiscStars[i], Gal[t].DiscStarsMetals[i]);
		assert(Gal[t].DiscStarsMetals[i] <= Gal[t].DiscStars[i]);
        assert(disc_mass_ratio[i]==0.0);
	}
	
  Gal[t].ColdGas += DiscGasSum;
  Gal[t].MetalsColdGas += Gal[p].MetalsColdGas;

  if(mass_ratio<ThreshMajorMerger) // Minor mergers, combine discs by conserving angular momentum
  {
	// Satellite's specific angular momentum
	double sat_sam[3];
	sat_sam[0] = (Gal[p].Pos[1]-Gal[t].Pos[1])*(Gal[p].Vel[2]-Gal[t].Vel[2])*ExpFac - (Gal[p].Pos[2]-Gal[t].Pos[2])*(Gal[p].Vel[1]-Gal[t].Vel[1])*ExpFac;
	sat_sam[1] = (Gal[p].Pos[2]-Gal[t].Pos[2])*(Gal[p].Vel[0]-Gal[t].Vel[0])*ExpFac - (Gal[p].Pos[0]-Gal[t].Pos[0])*(Gal[p].Vel[2]-Gal[t].Vel[2])*ExpFac;
	sat_sam[2] = (Gal[p].Pos[0]-Gal[t].Pos[0])*(Gal[p].Vel[1]-Gal[t].Vel[1])*ExpFac - (Gal[p].Pos[1]-Gal[t].Pos[1])*(Gal[p].Vel[0]-Gal[t].Vel[0])*ExpFac;
	
	double sat_sam_mag, cos_angle_sat_disc, sat_sam_max, sat_sam_min;
	int i_min, i_max, bin_num;
	
    sat_sam_mag = pow(sat_sam[0]*sat_sam[0] + sat_sam[1]*sat_sam[1] + sat_sam[2]*sat_sam[2], 0.5);
      
	//if(CentralGasOrig > 0.0 && Gal[p].ColdGas > 0.0)
	if(Gal[p].ColdGas > 0.0)
	{
		cos_angle_sat_disc = (Gal[t].SpinGas[0]*sat_sam[0] + Gal[t].SpinGas[1]*sat_sam[1] + Gal[t].SpinGas[2]*sat_sam[2]) / sat_sam_mag; // Angle between ang mom of satellite and central's disc
		sat_sam_mag *= fabs(cos_angle_sat_disc); // Project satellite's (gas) angular momentum onto central's disc
	
		// Consider that the satellite will have rotation and hence it will have a distribution of angular momentum to contribute
		sat_sam_max =  sat_sam_mag  +  Gal[p].Vvir * fabs(cos_angle_sat_disc) * pow(pow(Gal[p].Pos[0]-Gal[t].Pos[0], 2.0) + pow(Gal[p].Pos[1]-Gal[t].Pos[1], 2.0) + pow(Gal[p].Pos[2]-Gal[t].Pos[2], 2.0), 0.5);
		sat_sam_min = 2.0*sat_sam_mag - sat_sam_max;
		if(sat_sam_min<0.0)
			sat_sam_min = 0.0;
		
        if(cos_angle_sat_disc < 0.0)
            RetroCount += 1;
        else
            ProCount += 1;
        
		i_min=0;
		while(DiscBinEdge[i_min]<=sat_sam_min)
		{
			i_min++;
			if(i_min==N_BINS) break;
		}
		i_min -= 1;
        
		i_max=i_min;
		while(DiscBinEdge[i_max]<=sat_sam_max)
		{
			i_max++;
			if(i_max==N_BINS) break;
		}
	
		bin_num = i_max - i_min; // How many bins the satellite's gas will be added to in the main disc
	
		double gas_added = 0.0;
	
		for(i=0; i<N_BINS; i++)
		{
			if(i<i_min || i>=i_max)
			{
				disc_mass_ratio[i] = 0.0; // Probably redundant line, given all initialised at 0.
				PostRetroGas[i] = Gal[t].DiscGas[i];
			}
			else
			{
                if(Gal[t].DiscGas[i] > 0.0 && bin_num > 0)
                    disc_mass_ratio[i] = Gal[p].ColdGas / bin_num / Gal[t].DiscGas[i];
                else
                    disc_mass_ratio[i] = 0.0;
				if(disc_mass_ratio[i] > 1.0) disc_mass_ratio[i] = 1.0/disc_mass_ratio[i];
                assert(disc_mass_ratio[i] <= 1.0 && disc_mass_ratio[i]>=0.0);
				Gal[t].DiscGas[i] += Gal[p].ColdGas / bin_num;
				gas_added += Gal[p].ColdGas / bin_num;
				Gal[t].DiscGasMetals[i] += Gal[p].MetalsColdGas / bin_num;
				assert(Gal[t].DiscGasMetals[i] <= Gal[t].DiscGas[i]);
				
				if(cos_angle_sat_disc < 0.0)
					PostRetroGas[i] = Gal[t].DiscGas[i] - 2.0*Gal[p].ColdGas / bin_num;
				else
					PostRetroGas[i] = Gal[t].DiscGas[i];
					
				if(PostRetroGas[i] < 0.0) 
					PostRetroGas[i] = 0.0;
			}
		}
	
		assert(gas_added <= 1.01*Gal[p].ColdGas && gas_added >= Gal[p].ColdGas/1.01);
		DiscGasSum = get_disc_gas(p);
		assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
		DiscGasSum = get_disc_gas(t);
		assert(DiscGasSum <= 1.01*Gal[t].ColdGas && DiscGasSum >= Gal[t].ColdGas/1.01);

    }
    else
        for(i=0; i<N_BINS; i++) PostRetroGas[i] = Gal[t].DiscGas[i];
        

    // Now just going to throw stars into the bulge.
    if(Gal[p].StellarMass>0.0)
    {
//        for(s=0; s<3; s++)
//        {
//            assert(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY);
//            Gal[t].SpinClassicalBulge[s] = (Gal[t].SpinClassicalBulge[s]*Gal[t].ClassicalBulgeMass + sat_sam[s]*Gal[p].StellarMass)/(Gal[t].ClassicalBulgeMass+Gal[p].StellarMass);//sat_sam_mag;
//            if(!(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY))
//            {
//                printf("548: Reset classical bulge spin to zero\n");
//                Gal[t].SpinClassicalBulge[s] = 0.0;
//            }
//            assert(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY);
//        }
//        double spinmag = pow(pow(Gal[t].SpinClassicalBulge[0],2.0)+pow(Gal[t].SpinClassicalBulge[1],2.0)+pow(Gal[t].SpinClassicalBulge[2],2.0), 0.5);
//        if(spinmag<1.01 && spinmag>0.99)
//        {
//            printf("spinmag %e at 537\n, spinmag");
//            printf("sat_sam = %e, %e, %e\n", sat_sam[0], sat_sam[1], sat_sam[2]);
//        }
        Gal[t].ClassicalBulgeMass += Gal[p].StellarMass;
        Gal[t].ClassicalMetalsBulgeMass += Gal[p].MetalsStellarMass;
    }
      
  }
  else // Major mergers -- a more complex treatment of the gas could be done in future versions
  {
	
	DiscGasSum = get_disc_gas(p);
	assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
	if(Gal[p].ColdGas > 0.0)
	{
		double new_spin_mag, cos_angle_t, cos_angle_p;
		double NewSpin[3];
		double NewDiscT[N_BINS], NewDiscP[N_BINS], NewDiscMetalsT[N_BINS], NewDiscMetalsP[N_BINS];
	
		// Determine spin of new gaseous disc
		for(i=0; i<3; i++) NewSpin[i] = Gal[t].SpinGas[i]*CentralGasOrig + Gal[p].SpinGas[i]*Gal[p].ColdGas;
		new_spin_mag = pow(NewSpin[0]*NewSpin[0] + NewSpin[1]*NewSpin[1] + NewSpin[2]*NewSpin[2], 0.5);
		for(i=0; i<3; i++) NewSpin[i] /= new_spin_mag;
	
		DiscGasSum = get_disc_gas(p);
		assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
		cos_angle_t = Gal[t].SpinGas[0]*NewSpin[0] + Gal[t].SpinGas[1]*NewSpin[1] + Gal[t].SpinGas[2]*NewSpin[2];
		cos_angle_p = Gal[p].SpinGas[0]*NewSpin[0] + Gal[p].SpinGas[1]*NewSpin[1] + Gal[p].SpinGas[2]*NewSpin[2];
		
		DiscGasSum = get_disc_gas(p);
		assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
		
		project_disc(Gal[t].DiscGas, cos_angle_t, t, NewDiscT);		
		project_disc(Gal[p].DiscGas, cos_angle_p, p, NewDiscP);
		project_disc(Gal[t].DiscGasMetals, cos_angle_t, t, NewDiscMetalsT);		
		project_disc(Gal[p].DiscGasMetals, cos_angle_p, p, NewDiscMetalsP);
        
		for(i=0; i<N_BINS; i++)
		{
			Gal[p].DiscGas[i] = NewDiscP[i];
			Gal[p].DiscGasMetals[i] = NewDiscMetalsP[i]; // Evidently I need these to prevent an error -- project_gas must actually change the DiscGas values.
			Gal[t].DiscGas[i] = NewDiscT[i] + NewDiscP[i];
			Gal[t].DiscGasMetals[i] = NewDiscMetalsT[i] + NewDiscMetalsP[i];
			assert(Gal[t].DiscGasMetals[i] <= Gal[t].DiscGas[i]);
			assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
			
            if(NewDiscP[i] > 0.0)
                disc_mass_ratio[i] = NewDiscT[i] / NewDiscP[i];
            else
                disc_mass_ratio[i] = 0.0;
            
			if(disc_mass_ratio[i] > 1.0)
				disc_mass_ratio[i] = 1.0 / disc_mass_ratio[i];
            
            assert(disc_mass_ratio[i]<=1.0 && disc_mass_ratio[i]>=0.0);
		}
        
		DiscGasSum = get_disc_gas(p);
		assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
		DiscGasSum = get_disc_gas(t);
		assert(DiscGasSum <= 1.01*Gal[t].ColdGas && DiscGasSum >= Gal[t].ColdGas/1.01);
		
		// Output expected mass of each annulus after retrograde gas is dealt with
		for(i=0; i<N_BINS; i++)
		{
			if(cos_angle_t < 0.0)
				PostRetroGas[i] = NewDiscP[i] - NewDiscT[i];
			else if(cos_angle_p < 0.0)
				PostRetroGas[i] = NewDiscT[i] - NewDiscP[i];
			else
				PostRetroGas[i] = Gal[t].DiscGas[i];
				
			if(PostRetroGas[i] < 0.0) 
				PostRetroGas[i] = 0.0;
		}
        
        // Set the new spin direction of the gas
        for(i=0; i<3; i++) Gal[t].SpinGas[i] = NewSpin[i];
    }
    else
        for(i=0; i<N_BINS; i++) PostRetroGas[i] = Gal[t].DiscGas[i];

	DiscGasSum = get_disc_gas(t);
	assert(DiscGasSum <= 1.01*Gal[t].ColdGas && DiscGasSum >= Gal[t].ColdGas/1.01);

      // Set spin of classical bulge.  The mass itself will be transfered there in stars_to_bulge
      // First need to know the pos and vel of the centre of momentum
      double COM[3], vCOM[3], j_t[3], j_p[3], m_t, m_p;
      m_t = dmax(Gal[t].Mvir, Gal[t].StellarMass+Gal[t].HotGas+Gal[t].ColdGas+Gal[t].BlackHoleMass);
      m_p = dmax(Gal[p].Mvir, Gal[p].StellarMass+Gal[p].HotGas+Gal[p].ColdGas+Gal[p].BlackHoleMass);
      for(i=0; i<3; i++)
      {
         COM[i] = (Gal[t].Pos[i]*m_t + Gal[p].Pos[i]*m_p)/(m_t+m_p);
        vCOM[i] = (Gal[t].Vel[i]*m_t + Gal[p].Vel[i]*m_p)/(m_t+m_p);
      }

      j_p[0] = (Gal[p].Pos[1]-COM[1])*(Gal[p].Vel[2]-vCOM[2])*ExpFac - (Gal[p].Pos[2]-COM[2])*(Gal[p].Vel[1]-vCOM[1])*ExpFac;
      j_p[1] = (Gal[p].Pos[2]-COM[2])*(Gal[p].Vel[0]-vCOM[0])*ExpFac - (Gal[p].Pos[0]-COM[0])*(Gal[p].Vel[2]-vCOM[2])*ExpFac;
      j_p[2] = (Gal[p].Pos[0]-COM[0])*(Gal[p].Vel[1]-vCOM[1])*ExpFac - (Gal[p].Pos[1]-COM[1])*(Gal[p].Vel[0]-vCOM[0])*ExpFac;
      j_t[0] = (Gal[t].Pos[1]-COM[1])*(Gal[t].Vel[2]-vCOM[2])*ExpFac - (Gal[t].Pos[2]-COM[2])*(Gal[t].Vel[1]-vCOM[1])*ExpFac;
      j_t[1] = (Gal[t].Pos[2]-COM[2])*(Gal[t].Vel[0]-vCOM[0])*ExpFac - (Gal[t].Pos[0]-COM[0])*(Gal[t].Vel[2]-vCOM[2])*ExpFac;
      j_t[2] = (Gal[t].Pos[0]-COM[0])*(Gal[t].Vel[1]-vCOM[1])*ExpFac - (Gal[t].Pos[1]-COM[1])*(Gal[t].Vel[0]-vCOM[0])*ExpFac;
      
      double spinmagpre = pow(pow(Gal[t].SpinClassicalBulge[0],2.0)+pow(Gal[t].SpinClassicalBulge[1],2.0)+pow(Gal[t].SpinClassicalBulge[2],2.0), 0.5);
      for(s=0; s<3; s++)
      {
          assert(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY);
//          Gal[t].SpinClassicalBulge[s] = (Gal[p].StellarMass*j_p[s] + Gal[t].StellarMass*j_t[s] + Gal[t].SecularBulgeMass*Gal[t].SpinSecularBulge[s] + Gal[p].SecularBulgeMass*Gal[p].SpinSecularBulge[s] + Gal[t].ClassicalBulgeMass*Gal[t].SpinClassicalBulge[s] + Gal[p].ClassicalBulgeMass*Gal[p].SpinClassicalBulge[s] + get_disc_ang_mom(p,1)*Gal[p].SpinStars[s] + get_disc_ang_mom(t,1)*Gal[t].SpinStars[s]) / (Gal[p].StellarMass+Gal[t].StellarMass);
          Gal[t].SpinClassicalBulge[s] = (Gal[t].SecularBulgeMass*Gal[t].SpinSecularBulge[s] + Gal[p].SecularBulgeMass*Gal[p].SpinSecularBulge[s] + Gal[t].ClassicalBulgeMass*Gal[t].SpinClassicalBulge[s] + Gal[p].ClassicalBulgeMass*Gal[p].SpinClassicalBulge[s] + get_disc_ang_mom(p,1)*Gal[p].SpinStars[s] + get_disc_ang_mom(t,1)*Gal[t].SpinStars[s]) / (Gal[p].StellarMass+Gal[t].StellarMass);
          if(!(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY))
          {
              //printf("SM_p, SM_t, j_p, j_t = %e, %e, %e, %e\n", Gal[p].StellarMass, Gal[t].StellarMass, j_p[s], j_t[s]);
              Gal[t].SpinClassicalBulge[s] = 0.0;
          }
          assert(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY);

      }
      double spinmag = pow(pow(Gal[t].SpinClassicalBulge[0],2.0)+pow(Gal[t].SpinClassicalBulge[1],2.0)+pow(Gal[t].SpinClassicalBulge[2],2.0), 0.5);
      if(spinmag<1.01 && spinmag>0.99 && !(spinmagpre<1.01 && spinmagpre>0.99))
      {
          printf("spinmag %e at 654\n", spinmag);
          printf("j_t = %e, %e, %e\n", j_t[0], j_t[1], j_t[2]);
          printf("j_p = %e, %e, %e\n", j_p[0], j_p[1], j_p[2]);
      }
      
    Gal[t].ClassicalBulgeMass += Gal[p].StellarMass;
	Gal[t].ClassicalMetalsBulgeMass += Gal[p].MetalsStellarMass;
      
  }


  Gal[t].StarsInSitu += Gal[p].StarsInSitu;
  Gal[t].StarsInstability += Gal[p].StarsInstability;
  Gal[t].StarsMergeBurst += Gal[p].StarsMergeBurst;
    
    Gal[t].AccretedGasMass += Gal[p].AccretedGasMass;
    Gal[t].EjectedSNGasMass += Gal[p].EjectedSNGasMass;
    Gal[t].EjectedQuasarGasMass += Gal[p].EjectedQuasarGasMass;

  Gal[t].StellarMass += Gal[p].StellarMass;
  Gal[t].MetalsStellarMass += Gal[p].MetalsStellarMass;

    StarChannelSum = get_channel_stars(t);
    assert(Gal[t].StellarMass >= (Gal[t].StarsInSitu+Gal[t].StarsInstability+Gal[t].StarsMergeBurst)/1.01 && Gal[t].StellarMass <= (Gal[t].StarsInSitu+Gal[t].StarsInstability+Gal[t].StarsMergeBurst)*1.01);

    
  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;
  
  Gal[t].EjectedMass += Gal[p].EjectedMass;
  Gal[t].MetalsEjectedMass += Gal[p].MetalsEjectedMass;
  
  Gal[t].ICS += Gal[p].ICS;
  Gal[t].MetalsICS += Gal[p].MetalsICS;

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;
    assert(Gal[t].BlackHoleMass>=0.0);

  

  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulge[step] += Gal[p].SfrDisk[step] + Gal[p].SfrBulge[step];
    Gal[t].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step] + Gal[p].SfrBulgeColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step] + Gal[p].SfrBulgeColdGasMetals[step];
  }

  DiscGasSum = get_disc_gas(t);
  assert(DiscGasSum <= 1.01*Gal[t].ColdGas && DiscGasSum >= Gal[t].ColdGas/1.01);

	for(i=0; i<N_BINS; i++) 
	{
		if (Gal[t].DiscStarsMetals[i] > Gal[t].DiscStars[i])
			printf("DiscStars, Metals = %e, %e\n", Gal[t].DiscStars[i], Gal[t].DiscStarsMetals[i]);
		assert(Gal[t].DiscStarsMetals[i] <= Gal[t].DiscStars[i]);
        
        // This should already be taken care of above, but for whatever reason I needed to add it here to actually work.
        if(disc_mass_ratio[i] > 1.0) disc_mass_ratio[i] = 1.0/disc_mass_ratio[i];
        assert(disc_mass_ratio[i]<=1.0 && disc_mass_ratio[i]>=0.0);
	}
}



void stars_to_bulge(int t, int p)
{
  int step, i;
    
  // generate bulge 
  Gal[t].ClassicalBulgeMass = Gal[t].StellarMass;
  Gal[t].ClassicalMetalsBulgeMass = Gal[t].MetalsStellarMass;
  
  Gal[t].SecularBulgeMass = 0.0;
  Gal[t].SecularMetalsBulgeMass = 0.0;
  for(i=0; i<3; i++) Gal[t].SpinSecularBulge[i] = 0.0;
    
  // Remove stars from the disc annuli
  for(i=0; i<N_BINS; i++)
  {
	Gal[t].DiscStars[i] = 0.0;
	Gal[t].DiscStarsMetals[i] = 0.0;
  }  

    
  // update the star formation rate 
  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulge[step] += Gal[t].SfrDisk[step];
    Gal[t].SfrBulgeColdGas[step] += Gal[t].SfrDiskColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[t].SfrDiskColdGasMetals[step];
    Gal[t].SfrDisk[step] = 0.0;
    Gal[t].SfrDiskColdGas[step] = 0.0;
    Gal[t].SfrDiskColdGasMetals[step] = 0.0;
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


void collisional_starburst_recipe(double disc_mass_ratio[N_BINS], int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step, double mass_ratio)
{
 double stars, reheated_mass, ejected_mass, fac, metallicity, CentralVvir, eburst, Sigma_0gas, area, stars_sum, metals_stars, stars_angmom;
 double r_inner, r_outer, StarChannelSum, j_bin;
 //double NewStars[N_BINS], NewStarsMetals[N_BINS];
 int k, s;

 // This is the major and minor merger starburst recipe of Somerville et al. 2001. 
 // The coefficients in eburst are taken from TJ Cox's PhD thesis and should be more 
 // accurate then previous. The recipe has been modified to function for each annulus.
    
    double ejected_sum = 0.0;
    double metals_stars_sum = 0.0;

 stars_sum = 0.0;
    stars_angmom = 0.0;
 assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

 if(Gal[merger_centralgal].ColdGas>0)
 {
  CentralVvir = Gal[centralgal].Vvir;

  // update the star formation rate 
  Gal[merger_centralgal].SfrBulgeColdGas[step] += Gal[merger_centralgal].ColdGas;
  Gal[merger_centralgal].SfrBulgeColdGasMetals[step] += Gal[merger_centralgal].MetalsColdGas;

  for(k=0; k<N_BINS; k++)
  {
      if(disc_mass_ratio[k] > 1.0 || disc_mass_ratio[k]!=disc_mass_ratio[k]) printf("i, disc_mass_ratio[i] = %d, %e\n", k, disc_mass_ratio[k]);
      assert(disc_mass_ratio[k] <= 1.0);
      
	// the bursting fraction 
    if(mode == 1)
      eburst = disc_mass_ratio[k];
    else
      eburst = 0.56 * pow(disc_mass_ratio[k], 0.7);

    stars = eburst * Gal[merger_centralgal].DiscGas[k];
    if(stars < MIN_STARFORMATION)
      stars = 0.0;

	if(stars > Gal[merger_centralgal].DiscGas[k])
      stars = Gal[merger_centralgal].DiscGas[k];
	
	// this bursting results in SN feedback on the cold/hot gas 
    if(SupernovaRecipeOn>0 && Gal[merger_centralgal].DiscGas[k] > 0.0 && stars>=MIN_STARS_FOR_SN)//1e-9)
	{
//	  if(stars>MIN_STARS_FOR_SN)
//	  {
        r_inner = Gal[merger_centralgal].DiscRadii[k];
        r_outer = Gal[merger_centralgal].DiscRadii[k+1];
          
        area = M_PI * (r_outer*r_outer - r_inner*r_inner);
        
        if(SupernovaRecipeOn == 1)
        {
            Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
            reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[merger_centralgal].DiscGas[k]/area/1.3);
        }
        else if(SupernovaRecipeOn == 2)
            reheated_mass = FeedbackReheatingEpsilon * stars;
		
		// can't use more cold gas than is available! so balance SF and feedback
	    if((stars + reheated_mass) > Gal[merger_centralgal].DiscGas[k] && (stars + reheated_mass) > 0.0)
	    {
	      fac = Gal[merger_centralgal].DiscGas[k] / (stars + reheated_mass);
	      stars *= fac;
	      reheated_mass *= fac;
	    }
	
	    if(stars<MIN_STARS_FOR_SN)
	    {
		  stars = MIN_STARS_FOR_SN;
		  reheated_mass = Gal[merger_centralgal].DiscGas[k] - stars; // Used to have (1-RecycleFraction)* in front of stars here, but changed philosophy
	    }
	
        ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (CentralVvir * CentralVvir) - FeedbackReheatingEpsilon) * stars;
	    if(ejected_mass < 0.0)
	        ejected_mass = 0.0;
          
//		assert(RecycleFraction*stars+reheated_mass <= 1.01*Gal[merger_centralgal].DiscGas[k]);
//	  }

//	  else
//	  {
//		reheated_mass = RecycleFraction * stars;
//		ejected_mass = 0.0;
//	  }
	}
    else
	{
	  //stars=0.0;
      reheated_mass = 0.0;
	  ejected_mass = 0.0;
	}
      
      ejected_sum += ejected_mass;
      
	if(reheated_mass!=reheated_mass || reheated_mass<0.0)
		printf("reheated_mass, stars, fac, DiscGas -- %e\t%e\t%e\t%e\n", reheated_mass, stars, fac, Gal[merger_centralgal].DiscGas[k]);	
	assert(reheated_mass >= 0.0);

//      if(stars==MIN_STARS_FOR_SN)
//          printf("stars at min for SF: annulus gas, stars, reheated = %e, %e, %e\n", Gal[merger_centralgal].DiscGas[k], stars, reheated_mass);
      
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	assert(Gal[merger_centralgal].DiscGasMetals[k] <= Gal[merger_centralgal].DiscGas[k]);
    metals_stars = metallicity * stars;
    
	//NewStars[k] = (1 - RecycleFraction) * stars;
	//NewStarsMetals[k] = (1 - RecycleFraction) * metallicity * stars;
    update_from_star_formation(merger_centralgal, stars, metallicity, k);

    if(reheated_mass > Gal[merger_centralgal].DiscGas[k] && reheated_mass < 1.01*Gal[merger_centralgal].DiscGas[k])
	  reheated_mass = Gal[merger_centralgal].DiscGas[k];

	// update from feedback
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	assert(Gal[merger_centralgal].DiscGasMetals[k] <= Gal[merger_centralgal].DiscGas[k]);
      assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
	update_from_feedback(merger_centralgal, centralgal, reheated_mass, metallicity, k);
 
    // Inject new metals from SN II
	if(SupernovaRecipeOn>0 && stars>=MIN_STARS_FOR_SN)
	{
	    Gal[merger_centralgal].DiscGasMetals[k] += Yield * stars*(1-get_metallicity(stars,metals_stars));
	    Gal[merger_centralgal].MetalsColdGas += Yield * stars*(1-get_metallicity(stars,metals_stars));
	}
      
    if(!(Gal[merger_centralgal].DiscGasMetals[k]<=Gal[merger_centralgal].DiscGas[k]))
          printf("metals, gas = %e, %e\n", Gal[merger_centralgal].DiscGasMetals[k], Gal[merger_centralgal].DiscGas[k]);
      
	assert(Gal[merger_centralgal].DiscGasMetals[k]<=Gal[merger_centralgal].DiscGas[k]);
	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
	
      j_bin = (DiscBinEdge[k]+DiscBinEdge[k+1])/2;
    if(stars>=MIN_STARS_FOR_SN)
    {
        stars_sum += (1 - RecycleFraction) * stars;
        metals_stars_sum += (1 - RecycleFraction) * metals_stars;
        stars_angmom += (1 - RecycleFraction) * stars * j_bin;
    }
    else
    {
        stars_sum += stars;
        metals_stars_sum += metals_stars;
        stars_angmom += stars * j_bin;
    }
    Gal[merger_centralgal].DiscSFR[k] += stars / dt;
  }
     
  if(ejected_sum>0.0)
      update_from_ejection(centralgal, ejected_sum);

  // Sum stellar discs together
  //combine_stellar_discs(merger_centralgal, NewStars, NewStarsMetals);
     
  if(stars_sum>0)
  {
      double spinmagpre = pow(pow(Gal[merger_centralgal].SpinClassicalBulge[0],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[1],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[2],2.0), 0.5);
     // Update bulge spin
     for(s=0; s<3; s++)
     {
         assert(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY);
         Gal[merger_centralgal].SpinClassicalBulge[s] = (Gal[merger_centralgal].SpinClassicalBulge[s]*Gal[merger_centralgal].ClassicalBulgeMass + Gal[merger_centralgal].SpinGas[s]*stars_angmom) / (Gal[merger_centralgal].ClassicalBulgeMass+stars_sum);
         if(!(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY))
         {
             printf("BulgeSpin, BulgeMass, SpinGas, stars_sum = %e, %e, %e, %e\n", Gal[merger_centralgal].SpinClassicalBulge[s], Gal[merger_centralgal].ClassicalBulgeMass, Gal[merger_centralgal].SpinGas[s], stars_sum);
             Gal[merger_centralgal].SpinClassicalBulge[s] = 0.0;
         }
         assert(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY);
     }
      double spinmag = pow(pow(Gal[merger_centralgal].SpinClassicalBulge[0],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[1],2.0)+pow(Gal[merger_centralgal].SpinClassicalBulge[2],2.0), 0.5);
      //if(spinmag<1.01 && spinmag>0.99 && !(spinmagpre<1.01 && spinmagpre>0.99)) printf("spinmag %e at 925\n", spinmag);
      
     // Now adding all new stars directly to the bulge
     Gal[merger_centralgal].StellarMass += stars_sum; // Recycling fraction already taken into account when adding to stars_sum etc above
     Gal[merger_centralgal].ClassicalBulgeMass += stars_sum;
     Gal[merger_centralgal].MetalsStellarMass += metals_stars_sum;
     Gal[merger_centralgal].ClassicalMetalsBulgeMass += metals_stars_sum;
  }

  Gal[merger_centralgal].SfrBulge[step] += stars_sum / dt; // This can probably be handled better.  It's a bit complicated though (as to whether it's Sfr in the bulge or disc)
  Gal[merger_centralgal].StarsMergeBurst += stars_sum;
     
  StarChannelSum = get_channel_stars(merger_centralgal);
  assert(Gal[merger_centralgal].StellarMass >= (Gal[merger_centralgal].StarsInSitu+Gal[merger_centralgal].StarsInstability+Gal[merger_centralgal].StarsMergeBurst)/1.01 && Gal[merger_centralgal].StellarMass <= (Gal[merger_centralgal].StarsInSitu+Gal[merger_centralgal].StarsInstability+Gal[merger_centralgal].StarsMergeBurst)*1.01);

  // check for disk instability
//   if(DiskInstabilityOn && mode == 0)
//     check_disk_instability(merger_centralgal, centralgal, time, dt, step);
  //   if(mass_ratio < ThreshMajorMerger)
 }
}



