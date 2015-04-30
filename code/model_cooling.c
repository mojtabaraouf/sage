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
	{
		do_AGN_heating(coolingGas, gal, dt, x, rcool);

		// update the cooling rate based on current AGN heating
		if(Gal[gal].r_heat < rcool)
			coolingGas = (1.0 - Gal[gal].r_heat / rcool) * coolingGas;
		else
			coolingGas = 0.0;
	}
	
    if(coolingGas < 0.0)
      coolingGas = 0.0;

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
  }
  else
    AGNheating = 0.0;

  if (AGNheating > 0.0)
    Gal[centralgal].Heating += 0.5 * AGNheating * Gal[centralgal].Vvir * Gal[centralgal].Vvir;
  
  // update the heating radius as needed
  r_heat_new = (AGNheating / coolingGas) * rcool;
  if(r_heat_new > Gal[centralgal].r_heat)
	  Gal[centralgal].r_heat = r_heat_new;

  return AGNheating;

}



// This cools the gas onto the correct disc bins
void cool_gas_onto_galaxy(int p, double coolingGas)
{
  double metallicity, coolingGasBin, coolingGasBinSum, DiscGasSum;
  int i;

  // Check that Cold Gas has been treated properly prior to this function
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum < 1.001*Gal[p].ColdGas && DiscGasSum > Gal[p].ColdGas/1.001);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);

  // add the fraction 1/STEPS of the total cooling gas to the cold disk 
  if(coolingGas > 0.0)
  {
	coolingGasBinSum = 0.0;
	metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
    if(coolingGas < Gal[p].HotGas)
    {
	  for(i=0; i<30; i++)
      {
		coolingGasBin = (coolingGas / Gal[p].DiskScaleRadius) * (exp(-DiscBinEdge[i]/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i]/Gal[p].Vvir + Gal[p].DiskScaleRadius) - exp(-DiscBinEdge[i+1]/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i+1]/Gal[p].Vvir + Gal[p].DiskScaleRadius));
		if(coolingGasBin + coolingGasBinSum > coolingGas)
		  coolingGasBin = coolingGas - coolingGasBinSum;
		
	    Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
		assert(Gal[p].DiscGasMetals[i]<Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
	  }
	
	  assert(coolingGasBinSum < 1.001*coolingGas && coolingGasBinSum > coolingGas/1.001);

      Gal[p].ColdGas += coolingGas;
      Gal[p].MetalsColdGas += metallicity * coolingGas;
      Gal[p].HotGas -= coolingGas;
      Gal[p].MetalsHotGas -= metallicity * coolingGas;
    }
    else
    {
	  for(i=0; i<30; i++)
      {
		coolingGasBin = (Gal[p].HotGas / Gal[p].DiskScaleRadius) * (exp(-DiscBinEdge[i]/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i]/Gal[p].Vvir + Gal[p].DiskScaleRadius) - exp(-DiscBinEdge[i+1]/Gal[p].Vvir/Gal[p].DiskScaleRadius)*(DiscBinEdge[i+1]/Gal[p].Vvir + Gal[p].DiskScaleRadius));
		if(coolingGasBin + coolingGasBinSum > coolingGas)
		  coolingGasBin = coolingGas - coolingGasBinSum;

		Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
		assert(Gal[p].DiscGasMetals[i]<Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
      }

	  assert(coolingGasBinSum < 1.001*coolingGas && coolingGasBinSum > coolingGas/1.001);

      Gal[p].ColdGas += Gal[p].HotGas;
      Gal[p].MetalsColdGas += Gal[p].MetalsHotGas;
      Gal[p].HotGas = 0.0;
      Gal[p].MetalsHotGas = 0.0;
    }
  }

  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
}