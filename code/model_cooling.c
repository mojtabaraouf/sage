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

  //printf("HotGas in cooling_recipe\n");
  //printf("%e", Gal[gal].HotGas);
  //printf("\n");

  if(Gal[gal].HotGas > 1.0e-6)
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
      // infall dominated regime 
      coolingGas = Gal[gal].HotGas / (Gal[gal].Rvir / Gal[gal].Vvir) * dt; 
    else
      // hot phase regime 
      coolingGas = (Gal[gal].HotGas / Gal[gal].Rvir) * (rcool / (2.0 * tcool)) * dt;

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

  return coolingGas;

}



double do_AGN_heating(double coolingGas, int centralgal, double dt, double x, double rcool)
{
  double AGNrate, EDDrate, AGNaccreted, AGNcoeff, AGNheating, metallicity, r_heat_new;


  // first update the cooling rate based on the past AGN heating
  if(Gal[centralgal].r_heat < rcool)
	coolingGas = (1.0 - Gal[centralgal].r_heat / rcool) * coolingGas;
  else
	coolingGas = 0.0;
	
  assert(coolingGas >= 0.0);

  if(Gal[centralgal].HotGas > 0.0)
  {

    if(AGNrecipeOn == 2)
    {
      // Bondi-Hoyle accretion recipe
      AGNrate = (2.5 * M_PI * G) * (0.375 * 0.6 * x) * Gal[centralgal].BlackHoleMass * RadioModeEfficiency;
    }
    else if(AGNrecipeOn == 3)
    {
      // Cold cloud accretion: trigger: rBH > 1.0e-4 Rsonic, and accretion rate = 0.01% cooling rate 
      if(Gal[centralgal].BlackHoleMass > 0.0001 * Gal[centralgal].Mvir * pow(rcool/Gal[centralgal].Rvir, 3.0))
        AGNrate = 0.0001 * coolingGas / dt;
      else
        AGNrate = 0.0;
    }
    else
    {
      // empirical (standard) accretion recipe 
      if(Gal[centralgal].Mvir > 0.0)
        AGNrate = RadioModeEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
          * (Gal[centralgal].BlackHoleMass / 0.01) * pow(Gal[centralgal].Vvir / 200.0, 3.0)
            * ((Gal[centralgal].HotGas / Gal[centralgal].Mvir) / 0.1);
      else
        AGNrate = RadioModeEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
          * (Gal[centralgal].BlackHoleMass / 0.01) * pow(Gal[centralgal].Vvir / 200.0, 3.0);
    }
    
    // Eddington rate 
    EDDrate = 1.3e48 * Gal[centralgal].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;

    // accretion onto BH is always limited by the Eddington rate 
    if(AGNrate > EDDrate)
      AGNrate = EDDrate;

    // accreted mass onto black hole 
    AGNaccreted = AGNrate * dt;

    // cannot accrete more mass than is available! 
    if(AGNaccreted > Gal[centralgal].HotGas)
      AGNaccreted = Gal[centralgal].HotGas;

    // coefficient to heat the cooling gas back to the virial temperature of the halo 
    // 1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s 
    AGNcoeff = (1.34e5 / Gal[centralgal].Vvir) * (1.34e5 / Gal[centralgal].Vvir);

    // cooling mass that can be suppresed from AGN heating 
    AGNheating = AGNcoeff * AGNaccreted;

    // limit heating to cooling rate 
    if(AGNheating > coolingGas)
    {
      AGNaccreted = coolingGas / AGNcoeff;
      AGNheating = coolingGas;
    }

    // accreted mass onto black hole 
    metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
    Gal[centralgal].BlackHoleMass += AGNaccreted;
    Gal[centralgal].HotGas -= AGNaccreted;
    Gal[centralgal].MetalsHotGas -= metallicity * AGNaccreted;

		// update the heating radius as needed
		if(Gal[centralgal].r_heat < rcool && coolingGas > 0.0)
		{
			r_heat_new = HeatingRadiusScale * (AGNheating / coolingGas) * rcool;
			if(r_heat_new > Gal[centralgal].r_heat)
				Gal[centralgal].r_heat = r_heat_new;
		}

		if (AGNheating > 0.0)
			Gal[centralgal].Heating += 0.5 * AGNheating * Gal[centralgal].Vvir * Gal[centralgal].Vvir;
  }

  return coolingGas;

}



// This cools the gas onto the correct disc bins
void cool_gas_onto_galaxy(int p, int centralgal, double coolingGas, double dt, int step)
{
  double metallicity, coolingGasBin, coolingGasBinSum, DiscGasSum, DiscGasSum_new, cos_angle_disc_new, cos_angle_halo_new, ratio_last_bin, high_bound, disc_spin_mag, J_disc, J_cool;
  double HaloSpin[3], DiscNewSpin[3];
  double OldDisc[30], OldDiscMetals[30], RetroGas[30];
  int i, j, k, j_old;

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

  assert(Gal[p].ColdGas == Gal[p].ColdGas);

  // add the fraction 1/STEPS of the total cooling gas to the cold disk 
  if(coolingGas > 0.0)
  {
	coolingGasBinSum = 0.0;
	metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
	
	// Get ang mom of cooling gas in its native orientations
	J_cool = 2.0 * coolingGas * Gal[p].Vvir * Gal[p].DiskScaleRadius;
	
	if(Gal[p].ColdGas > 0.0)
	{
		// Get magnitude of ang mom of disc currently in native orientation 
		J_disc = get_disc_ang_mom(p, 0);
		
		// Determine orientation of disc after cooling
		for(i=0; i<3; i++)
		{
			HaloSpin[i] = Halo[Gal[p].HaloNr].Spin[i] / pow(pow(Halo[Gal[p].HaloNr].Spin[0], 2.0) + pow(Halo[Gal[p].HaloNr].Spin[1], 2.0) + pow(Halo[Gal[p].HaloNr].Spin[2], 2.0), 0.5);
			DiscNewSpin[i] = HaloSpin[i]*J_cool + Gal[p].SpinGas[i]*J_disc; // Not normalised yet
		}

		disc_spin_mag = pow(pow(DiscNewSpin[0], 2.0) + pow(DiscNewSpin[1], 2.0) + pow(DiscNewSpin[2], 2.0), 0.5);
		for(i=0; i<3; i++)
			DiscNewSpin[i] /= disc_spin_mag; // Normalise it now
		cos_angle_disc_new = Gal[p].SpinGas[0]*DiscNewSpin[0] + Gal[p].SpinGas[1]*DiscNewSpin[1] + Gal[p].SpinGas[2]*DiscNewSpin[2];
		cos_angle_halo_new = HaloSpin[0]*DiscNewSpin[0] + HaloSpin[1]*DiscNewSpin[1] + HaloSpin[2]*DiscNewSpin[2];
	}
	else
	{
		cos_angle_disc_new = 1.0;
		cos_angle_halo_new = 1.0;
		J_disc = 0.0;
	}
	
    assert(Gal[p].ColdGas == Gal[p].ColdGas);
		
	if(cos_angle_disc_new != 1.0)
	{
		// Project current disc to new orientation
		for(i=0; i<30; i++)
		{
			OldDisc[i] = Gal[p].DiscGas[i];
			OldDiscMetals[i] = Gal[p].DiscGasMetals[i];
		}
		j_old = 0;
	
		for(i=0; i<30; i++)
		{
			high_bound = DiscBinEdge[i+1] / fabs(cos_angle_disc_new);
			j = j_old;
			
			while(DiscBinEdge[j]<=high_bound)
			{
				j++;
				if(j==30) break;
			} 
			j -= 1;
			
			Gal[p].DiscGas[i] = 0.0;
			Gal[p].DiscGasMetals[i] = 0.0;
			for(k=j_old; k<j; k++) 
			{
				Gal[p].DiscGas[i] += OldDisc[k];
				Gal[p].DiscGasMetals[i] += OldDiscMetals[k];
				OldDisc[k] = 0.0;
				OldDiscMetals[k] = 0.0;
			}
			if(i!=29)
			{
				if(j!=29){
					ratio_last_bin = pow((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]), 2.0);
					assert(ratio_last_bin<=1.0);}
				else if(high_bound < Gal[p].Rvir/Gal[p].Vvir){
					ratio_last_bin = pow((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir/Gal[p].Vvir-DiscBinEdge[j]), 2.0);
					assert(ratio_last_bin<=1.0);}
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
	assert(Gal[p].ColdGas == Gal[p].ColdGas);
		
	DiscGasSum_new = get_disc_gas(p);
	assert(DiscGasSum <= 1.001*DiscGasSum_new && DiscGasSum >= DiscGasSum_new/1.001);	
	
    if(coolingGas < Gal[p].HotGas)
    {
	  for(i=0; i<30; i++)
      {
		coolingGasBin = (coolingGas / Gal[p].DiskScaleRadius) * (exp(-DiscBinEdge[i]/fabs(cos_angle_halo_new)/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i]/fabs(cos_angle_halo_new)/Gal[p].Vvir + Gal[p].DiskScaleRadius) - exp(-DiscBinEdge[i+1]/fabs(cos_angle_halo_new)/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i+1]/fabs(cos_angle_halo_new)/Gal[p].Vvir + Gal[p].DiskScaleRadius));
		if(coolingGasBin + coolingGasBinSum > coolingGas)
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
	  assert(Gal[p].ColdGas == Gal[p].ColdGas);

    }
    else
    {
	  for(i=0; i<30; i++)
      {
		coolingGasBin = (Gal[p].HotGas / Gal[p].DiskScaleRadius) * (exp(-DiscBinEdge[i]/fabs(cos_angle_halo_new)/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i]/fabs(cos_angle_halo_new)/Gal[p].Vvir + Gal[p].DiskScaleRadius) - exp(-DiscBinEdge[i+1]/fabs(cos_angle_halo_new)/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i+1]/fabs(cos_angle_halo_new)/Gal[p].Vvir + Gal[p].DiskScaleRadius));
		assert(coolingGasBin>=0.0);
		if(coolingGasBin + coolingGasBinSum > coolingGas)
		  coolingGasBin = coolingGas - coolingGasBinSum;

		Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
		if(cos_angle_halo_new<0.0)
			RetroGas[i] = coolingGasBin;
		assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
      }

	  assert(coolingGasBinSum <= 1.001*coolingGas && coolingGasBinSum >= coolingGas/1.001);

      Gal[p].ColdGas += Gal[p].HotGas;
      Gal[p].MetalsColdGas += Gal[p].MetalsHotGas;
      Gal[p].HotGas = 0.0;
      Gal[p].MetalsHotGas = 0.0;
	  assert(Gal[p].ColdGas == Gal[p].ColdGas);
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


void retrograde_gas_collision(int p, double RetroGas[30], double cos_angle_halo_new, double cos_angle_disc_new, double J_disc, double J_cool)
{
	double J_sum, J_retro, bin_ratio;
	double NewDisc[30];
	double NewDiscMetals[30];
	int i;
	
	//, RetroSum,  high_bound, ratio_last_bin;
	//double OldDisc[30], OldDiscMetals[30];
	//int i, j, j_old, k;
	
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
	
	project_disc(Gal[p].DiscGas, bin_ratio, p, NewDisc);
	project_disc(Gal[p].DiscGasMetals, bin_ratio, p, NewDiscMetals);
	
	for(i=0; i<30; i++)
	{
		Gal[p].DiscGas[i] = NewDisc[i];
		Gal[p].DiscGasMetals[i] = NewDiscMetals[i];
	}
	
	// Project the smaller bins back to the original bins
	// for(i=0; i<30; i++)
	// {
	// 	OldDisc[i] = Gal[p].DiscGas[i];
	// 	OldDiscMetals[i] = Gal[p].DiscGasMetals[i];
	// }
	// j_old = 0;
	// 
	// for(i=0; i<30; i++)
	// {
	// 	high_bound = DiscBinEdge[i+1] / bin_ratio;
	// 	j = j_old;
	// 	while(DiscBinEdge[j]<=high_bound)
	// 	{
	// 		j++;
	// 		if(j==30) break;
	// 	} 
	// 	j -= 1;
	// 	
	// 	Gal[p].DiscGas[i] = 0.0;
	// 	Gal[p].DiscGasMetals[i] = 0.0;
	// 	for(k=j_old; k<j; k++) 
	// 	{
	// 		Gal[p].DiscGas[i] += OldDisc[k];
	// 		Gal[p].DiscGasMetals[i] += OldDiscMetals[k];
	// 		OldDisc[k] = 0.0;
	// 		OldDiscMetals[k] = 0.0;
	// 	}
	// 	
	// 	if(i!=29)
	// 	{
	// 		if(j!=29){
	// 			ratio_last_bin = pow((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]), 2.0);
	// 			assert(ratio_last_bin<=1.0);}
	// 		else if(high_bound < Gal[p].Rvir/Gal[p].Vvir){
	// 			ratio_last_bin = pow((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir/Gal[p].Vvir-DiscBinEdge[j]), 2.0);
	// 			assert(ratio_last_bin<=1.0);}
	// 		else
	// 			ratio_last_bin = 1.0;
	// 		Gal[p].DiscGas[i] += ratio_last_bin * OldDisc[j];
	// 		Gal[p].DiscGasMetals[i] += ratio_last_bin * OldDiscMetals[j];
	// 		OldDisc[j] -= ratio_last_bin * OldDisc[j];
	// 		OldDiscMetals[j] -= ratio_last_bin * OldDiscMetals[j];
	// 	}
	// 	else
	// 	{
	// 		Gal[p].DiscGas[i] = OldDisc[i];
	// 		Gal[p].DiscGasMetals[i] = OldDiscMetals[i];
	// 	}
	// 	assert(Gal[p].DiscGas[i]>=0.0);
	// 	assert(Gal[p].DiscGasMetals[i]>=0.0);
	// 	j_old = j;
	// }
}
