#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double cooling_recipe(int gal, double dt)
{
  double tcool, x, logZ, lambda, rcool, rho_rcool, rho0, temp, coolingGas;

  if(Gal[gal].HotGas > 0.0 && Gal[gal].Vvir > 0.0)
  {
    tcool = Gal[gal].Rvir / Gal[gal].Vvir;
    temp = 35.9 * Gal[gal].Vvir * Gal[gal].Vvir;         // in Kelvin 

    if(Gal[gal].MetalsHotGas > 0)
      logZ = log10(Gal[gal].MetalsHotGas / Gal[gal].HotGas);
    else
      logZ = -10.0;

    lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
    x = PROTONMASS * BOLTZMANN * temp / lambda;        // now this has units sec g/cm^3  
    x /= (UnitDensity_in_cgs * UnitTime_in_s);         // now in internal units 
    rho_rcool = x / tcool * 0.885;  // 0.885 = 3/2 * mu, mu=0.59 for a fully ionized gas

    // an isothermal density profile for the hot gas is assumed here 
    rho0 = Gal[gal].HotGas / (4 * M_PI * Gal[gal].Rvir);
    rcool = sqrt(rho0 / rho_rcool);

    if(rcool > Gal[gal].Rvir)
    {
      // infall dominated regime 
      coolingGas = Gal[gal].HotGas / (Gal[gal].Rvir / Gal[gal].Vvir) * dt;
//      Gal[gal].CoolScaleRadius = pow(10, 0.23*log10(1.414*Gal[gal].DiskScaleRadius/Gal[gal].Rvir) - 0.67-0.18) * Gal[gal].Rvir; // Stevens et al. (2017)
    }
    else
    {
      // hot phase regime 
      coolingGas = (Gal[gal].HotGas / Gal[gal].Rvir) * (rcool / (2.0 * tcool)) * dt;
//      Gal[gal].CoolScaleRadius = 1.0*Gal[gal].DiskScaleRadius;
    }
      Gal[gal].CoolScaleRadius = 1.0*Gal[gal].DiskScaleRadius;
      
    if(coolingGas > Gal[gal].HotGas)
      coolingGas = Gal[gal].HotGas;
    else if(coolingGas < 0.0)
      coolingGas = 0.0;

    if(AGNrecipeOn > 0 && coolingGas > 0.0)
		coolingGas = do_AGN_heating(coolingGas, gal, dt, x, rcool);

    if (coolingGas > 0.0)
      Gal[gal].Cooling += 0.5 * coolingGas * Gal[gal].Vvir * Gal[gal].Vvir;
  }
  else
    coolingGas = 0.0;
    
  assert(coolingGas >= 0.0);
  return coolingGas;

}



double do_AGN_heating(double coolingGas, int p, double dt, double x, double rcool)
{
  double AGNrate, EDDrate, AGNaccreted, AGNcoeff, AGNheating, metallicity, r_heat_new;

  // first update the cooling rate based on the past AGN heating
  if(Gal[p].r_heat < rcool)
	coolingGas = (1.0 - Gal[p].r_heat / rcool) * coolingGas;
  else
	coolingGas = 0.0;
	
  assert(coolingGas >= 0.0);

    //
  if(Gal[p].HotGas > 0.0)
  {

    if(AGNrecipeOn == 2)
    {
      // Bondi-Hoyle accretion recipe
      AGNrate = (2.5 * M_PI * G) * (0.375 * 0.6 * x) * Gal[p].BlackHoleMass * RadioModeEfficiency;
    }
    else if(AGNrecipeOn == 3)
    {
      // Cold cloud accretion: trigger: rBH > 1.0e-4 Rsonic, and accretion rate = 0.01% cooling rate 
      if(Gal[p].BlackHoleMass > 0.0001 * Gal[p].Mvir * pow(rcool/Gal[p].Rvir, 3.0))
        AGNrate = 0.0001 * coolingGas / dt;
      else
        AGNrate = 0.0;
    }
    else
    {
      // empirical (standard) accretion recipe 
      if(Gal[p].Mvir > 0.0)
        AGNrate = RadioModeEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
          * (Gal[p].BlackHoleMass / 0.01) * pow(Gal[p].Vvir / 200.0, 3.0)
            * ((Gal[p].HotGas / Gal[p].Mvir) / 0.1);
      else
        AGNrate = RadioModeEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
          * (Gal[p].BlackHoleMass / 0.01) * pow(Gal[p].Vvir / 200.0, 3.0);
    }
    
    // Eddington rate 
    EDDrate = 1.3e48 * Gal[p].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;

    // accretion onto BH is always limited by the Eddington rate 
    if(AGNrate > EDDrate)
      AGNrate = EDDrate;

    // accreted mass onto black hole 
    AGNaccreted = AGNrate * dt;

    // cannot accrete more mass than is available! 
    if(AGNaccreted > Gal[p].HotGas)
      AGNaccreted = Gal[p].HotGas;

    // coefficient to heat the cooling gas back to the virial temperature of the halo 
    // 1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s 
    AGNcoeff = (1.34e5 / Gal[p].Vvir) * (1.34e5 / Gal[p].Vvir);

    // cooling mass that can be suppresed from AGN heating 
    AGNheating = AGNcoeff * AGNaccreted;

    // limit heating to cooling rate 
    if(AGNheating > coolingGas)
    {
      AGNaccreted = coolingGas / AGNcoeff;
      AGNheating = coolingGas;
    }

    // accreted mass onto black hole
    metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
	assert(Gal[p].MetalsHotGas <= Gal[p].HotGas);
    Gal[p].BlackHoleMass += AGNaccreted;
      assert(Gal[p].BlackHoleMass>=0.0);
    Gal[p].HotGas -= AGNaccreted;
    Gal[p].MetalsHotGas -= metallicity * AGNaccreted;

		// update the heating radius as needed
		if(Gal[p].r_heat < rcool && coolingGas > 0.0)
		{
            r_heat_new = (AGNheating / coolingGas) * rcool;
			
			if(r_heat_new > Gal[p].r_heat)
				Gal[p].r_heat = r_heat_new;
		}

		if (AGNheating > 0.0)
			Gal[p].Heating += 0.5 * AGNheating * Gal[p].Vvir * Gal[p].Vvir;
  }

  return coolingGas;

}



// This cools the gas onto the correct disc bins
void cool_gas_onto_galaxy(int p, int centralgal, double coolingGas, double dt, int step)
{
  double metallicity, coolingGasBin, coolingGasBinSum, DiscGasSum, cos_angle_disc_new, cos_angle_halo_new, ratio_last_bin, high_bound, disc_spin_mag, J_disc, J_cool, SpinMag;
//  double r_inner, r_outer;
  double DiscNewSpin[3];
  double OldDisc[N_BINS], OldDiscMetals[N_BINS], RetroGas[N_BINS];
  int i, j, k, j_old;
  double jfrac1, jfrac2;
  double rfrac1, rfrac2;

  // Check that Cold Gas has been treated properly prior to this function
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);

  disc_spin_mag = pow(pow(Gal[p].SpinGas[0], 2.0) + pow(Gal[p].SpinGas[1], 2.0) + pow(Gal[p].SpinGas[2], 2.0), 0.5);
  assert(disc_spin_mag==disc_spin_mag);
    
  if(disc_spin_mag>0.0)
  {
    for(i=0; i<3; i++)
	{
		Gal[p].SpinGas[i] /= disc_spin_mag; // Ensure the disc spin magnitude is normalised
		DiscNewSpin[i] = Gal[p].SpinGas[i]; // Default value as the reverse line occurs later
	}
  }

  // add the fraction 1/STEPS of the total cooling gas to the cold disk 
  if(coolingGas > 0.0)
  {
	coolingGasBinSum = 0.0;
	metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
	assert(Gal[p].MetalsHotGas <= Gal[p].HotGas);
	
	// Get ang mom of cooling gas in its native orientations
	J_cool = 2.0 * coolingGas * Gal[p].Vvir * Gal[p].CoolScaleRadius;
	
	if(Gal[p].ColdGas > 0.0)
	{
		// Get magnitude of ang mom of disc currently in native orientation 
		J_disc = get_disc_ang_mom(p, 0);
		
		// Determine orientation of disc after cooling
		for(i=0; i<3; i++)
			DiscNewSpin[i] = Gal[p].SpinHot[i]*J_cool + Gal[p].SpinGas[i]*J_disc; // Not normalised yet

		disc_spin_mag = pow(pow(DiscNewSpin[0], 2.0) + pow(DiscNewSpin[1], 2.0) + pow(DiscNewSpin[2], 2.0), 0.5);
		for(i=0; i<3; i++)
			DiscNewSpin[i] /= disc_spin_mag; // Normalise it now
		cos_angle_disc_new = Gal[p].SpinGas[0]*DiscNewSpin[0] + Gal[p].SpinGas[1]*DiscNewSpin[1] + Gal[p].SpinGas[2]*DiscNewSpin[2];
		cos_angle_halo_new = Gal[p].SpinHot[0]*DiscNewSpin[0] + Gal[p].SpinHot[1]*DiscNewSpin[1] + Gal[p].SpinHot[2]*DiscNewSpin[2];
	}
	else
	{
		cos_angle_disc_new = 1.0;
		cos_angle_halo_new = 1.0;
		J_disc = 0.0;
        SpinMag = pow(pow(Halo[Gal[p].HaloNr].Spin[0], 2.0) + pow(Halo[Gal[p].HaloNr].Spin[1], 2.0) + pow(Halo[Gal[p].HaloNr].Spin[2], 2.0), 0.5);
        for(i=0; i<3; i++)
        {
            if(SpinMag>0)
                DiscNewSpin[i] = Halo[Gal[p].HaloNr].Spin[i] / SpinMag;
            else
                DiscNewSpin[i] = 0.0;
        }
	}
	
    assert(cos_angle_disc_new==cos_angle_disc_new);
		
	if(cos_angle_disc_new != 1.0 && cos_angle_disc_new != 0.0)
	{
		// Project current disc to new orientation.  Could use the project_disc function here.
		for(i=0; i<N_BINS; i++)
		{
			OldDisc[i] = Gal[p].DiscGas[i];
			OldDiscMetals[i] = Gal[p].DiscGasMetals[i];
		}
		j_old = 0;
	
		for(i=0; i<N_BINS; i++)
		{
			high_bound = DiscBinEdge[i+1] / fabs(cos_angle_disc_new);
			j = j_old;
			
			while(DiscBinEdge[j]<=high_bound)
			{
				j++;
				if(j==N_BINS) break;
			} 
			j -= 1;
			
			Gal[p].DiscGas[i] = 0.0;
			Gal[p].DiscGasMetals[i] = 0.0;
			for(k=j_old; k<j; k++) 
			{
				Gal[p].DiscGas[i] += OldDisc[k];
				Gal[p].DiscGasMetals[i] += OldDiscMetals[k];
				assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
				OldDisc[k] = 0.0;
				OldDiscMetals[k] = 0.0;
			}
			if(i!=N_BINS-1)
			{
				if(j!=N_BINS-1)
                {
					ratio_last_bin = pow((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]), 2.0);
					assert(ratio_last_bin<=1.0);
                }
				else if(high_bound < Gal[p].Rvir/Gal[p].Vvir)
                {
					ratio_last_bin = pow((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir/Gal[p].Vvir-DiscBinEdge[j]), 2.0);
					assert(ratio_last_bin<=1.0);
                }
				else
					ratio_last_bin = 1.0;
				Gal[p].DiscGas[i] += ratio_last_bin * OldDisc[j];
				Gal[p].DiscGasMetals[i] += ratio_last_bin * OldDiscMetals[j];
				OldDisc[j] -= ratio_last_bin * OldDisc[j];
				OldDiscMetals[j] -= ratio_last_bin * OldDiscMetals[j];
			}
			else
			{
				Gal[p].DiscGas[i] = OldDisc[i];
				Gal[p].DiscGasMetals[i] = OldDiscMetals[i]; // Shouldn't need to set the Old stuff to zero for this last bit, as it'll just get overwritten by the next galaxy
			}
			assert(Gal[p].DiscGas[i]>=0.0);
			assert(Gal[p].DiscGasMetals[i]>=0.0);
			if(cos_angle_disc_new<0.0)
				RetroGas[i] = Gal[p].DiscGas[i];
			j_old = j;
		}
	}
		
	for(i=0; i<N_BINS; i++) assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);	
	
    if(coolingGas < Gal[p].HotGas)
    {
	  for(i=0; i<N_BINS; i++)
      {
          if(Gal[p].DiskScaleRadius==Gal[p].CoolScaleRadius)
          {
              jfrac1 = DiscBinEdge[i] / (Gal[p].Vvir * Gal[p].DiskScaleRadius);
              jfrac2 = DiscBinEdge[i+1] / (Gal[p].Vvir * Gal[p].DiskScaleRadius);
              coolingGasBin = coolingGas * ((jfrac1+1.0)*exp(-jfrac1) - (jfrac2+1.0)*exp(-jfrac2));
          }
          else
          {
              rfrac1 = Gal[p].DiscRadii[i] / Gal[p].CoolScaleRadius;
              rfrac2 = Gal[p].DiscRadii[i+1] / Gal[p].CoolScaleRadius;
              coolingGasBin = coolingGas * ((rfrac1+1.0)*exp(-rfrac1) - (rfrac2+1.0)*exp(-rfrac2));
          }
          
		if(coolingGasBin + coolingGasBinSum > coolingGas || i==N_BINS-1)
		  coolingGasBin = coolingGas - coolingGasBinSum;
          
	    Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
		if(cos_angle_halo_new<0.0)
			RetroGas[i] = coolingGasBin;
          
        
		assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
	  }
	
	  assert(coolingGasBinSum <= 1.01*coolingGas && coolingGasBinSum >= coolingGas/1.01);

      Gal[p].ColdGas += coolingGas;
      Gal[p].MetalsColdGas += metallicity * coolingGas;
      Gal[p].HotGas -= coolingGas;
      Gal[p].MetalsHotGas -= metallicity * coolingGas;
    }
    else
    {
	  for(i=0; i<N_BINS; i++)
      {
        if(Gal[p].DiskScaleRadius==Gal[p].CoolScaleRadius)
        {
            jfrac1 = DiscBinEdge[i] / (Gal[p].Vvir * Gal[p].DiskScaleRadius);
            jfrac2 = DiscBinEdge[i+1] / (Gal[p].Vvir * Gal[p].DiskScaleRadius);
            coolingGasBin = coolingGas * ((jfrac1+1.0)*exp(-jfrac1) - (jfrac2+1.0)*exp(-jfrac2));
        }
        else
        {
            rfrac1 = Gal[p].DiscRadii[i] / Gal[p].CoolScaleRadius;
            rfrac2 = Gal[p].DiscRadii[i+1] / Gal[p].CoolScaleRadius;
            coolingGasBin = Gal[p].HotGas * ((rfrac1+1.0)*exp(-rfrac1) - (rfrac2+1.0)*exp(-rfrac2));
        }
            
        assert(coolingGasBin>=0.0);
		if(coolingGasBin + coolingGasBinSum > coolingGas || i==N_BINS-1)
		  coolingGasBin = coolingGas - coolingGasBinSum;

		Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
		if(cos_angle_halo_new<0.0)
			RetroGas[i] = coolingGasBin;
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
      }

	  assert(coolingGasBinSum <= 1.001*coolingGas && coolingGasBinSum >= coolingGas/1.001);

      Gal[p].ColdGas += Gal[p].HotGas;
      Gal[p].MetalsColdGas += Gal[p].MetalsHotGas;
      Gal[p].HotGas = 0.0;
      Gal[p].MetalsHotGas = 0.0;
    }
			
    // Set spin of new disc
	for(i=0; i<3; i++) 
		Gal[p].SpinGas[i] = DiscNewSpin[i];
	
	if(cos_angle_disc_new < -1e-5 || cos_angle_halo_new < -1e-5) 
		retrograde_gas_collision(p, RetroGas, cos_angle_halo_new, cos_angle_disc_new, J_disc, J_cool);
	
  
  }
  assert(Gal[p].ColdGas == Gal[p].ColdGas);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
}


void retrograde_gas_collision(int p, double RetroGas[N_BINS], double cos_angle_halo_new, double cos_angle_disc_new, double J_disc, double J_cool)
{
	double J_sum, J_retro, bin_ratio;
	double NewDisc[N_BINS];
	double NewDiscMetals[N_BINS];
	int i;
	
	J_sum = J_disc*fabs(cos_angle_disc_new) + J_cool*fabs(cos_angle_halo_new);
	if(cos_angle_disc_new<0.0)
		J_retro = J_disc*fabs(cos_angle_disc_new);
	else if(cos_angle_halo_new<0.0)
		J_retro = J_cool*fabs(cos_angle_halo_new);
	else{
		J_retro = 0.0;
		printf("retrograde_gas_collision entered despite no retrograde gas");}
		
	assert(J_sum >= 2.0*J_retro);
	
	// Change the bin edges by the ratio of what the ang mom should be to the actual current ang mom
	bin_ratio = (J_sum - 2.0*J_retro) / J_sum;
	
	for(i=0; i<N_BINS; i++) assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	
	project_disc(Gal[p].DiscGas, bin_ratio, p, NewDisc);
	project_disc(Gal[p].DiscGasMetals, bin_ratio, p, NewDiscMetals);
	
	for(i=0; i<N_BINS; i++)
	{
		Gal[p].DiscGas[i] = NewDisc[i];
		Gal[p].DiscGasMetals[i] = NewDiscMetals[i];
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	}
}
