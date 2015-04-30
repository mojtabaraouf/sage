#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void check_disk_instability(int p, int centralgal, int halonr, double time, double dt, int step)
{
  double Mcrit, gas_fraction, unstable_gas, unstable_gas_fraction, unstable_stars, diskmass, metallicity, DiscGasSum;
  double star_fraction, ring_fraction, GasBeforeBH;
  int j;

  // Here we calculate the stability of the stellar and gaseous disk as discussed in Mo, Mao & White (1998).
  // For unstable stars and gas, we transfer the required ammount to the bulge to make the disk stable again

  // Check that Cold Gas has been treated properly prior to this function
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum < 1.001*Gal[p].ColdGas || DiscGasSum > Gal[p].ColdGas/1.001);

  // Disk mass has to be > 0.0 !
  diskmass = Gal[p].ColdGas + (Gal[p].StellarMass - Gal[p].ClassicalBulgeMass - Gal[p].SecularBulgeMass);
  if(diskmass > 0.0)
  {
    // calculate critical disk mass
    Mcrit = Gal[p].Vmax * Gal[p].Vmax * (3.0 * Gal[p].DiskScaleRadius) / G;
    if(Mcrit > diskmass)
      Mcrit = diskmass;
    
    // use Disk mass here !
    gas_fraction   = Gal[p].ColdGas / diskmass;
    unstable_gas   = gas_fraction * (diskmass - Mcrit);
    star_fraction  = 1.0 - gas_fraction;
    unstable_stars = star_fraction * (diskmass - Mcrit);

    // add excess stars to the bulge
    if(unstable_stars > 0.0)
    {
      for(j=0; j<30; j++)
	  {
		metallicity = get_metallicity(Gal[p].DiscStars[j], Gal[p].DiscStarsMetals[j]);
	    ring_fraction = Gal[p].DiscStars[j] / (star_fraction*diskmass);
		Gal[p].DiscStars[j] -= unstable_stars * ring_fraction;
		Gal[p].DiscStarsMetals[j] -= metallicity * unstable_stars * ring_fraction;
		Gal[p].SecularBulgeMass += unstable_stars * ring_fraction;
	    Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars * ring_fraction;
	  }
    }  
      
    // Need to fix this. Excluded for now.
    // Gal[p].mergeType = 3;  // mark as disk instability partial mass transfer
    // Gal[p].mergeIntoID = NumGals + p - 1;      
      
	assert((Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass) / Gal[p].StellarMass < 1.0001);
	assert((Gal[p].ClassicalMetalsBulgeMass + Gal[p].SecularMetalsBulgeMass) / Gal[p].MetalsStellarMass < 1.0001);

    // burst excess gas and feed black hole (really need a dedicated model for bursts and BH growth here)
    if(unstable_gas > 0.0)
    {
	  assert(unstable_gas/Gal[p].ColdGas < 1.0001);

	  unstable_gas_fraction = unstable_gas / Gal[p].ColdGas;
	  GasBeforeBH = Gal[p].ColdGas;
		
      if(AGNrecipeOn > 0)
        grow_black_hole(p, unstable_gas_fraction);
	
      // MY NEW SIMPLE WAY TO DEAL WITH AN INSTABILITY
      unstable_gas = unstable_gas - (GasBeforeBH - Gal[p].ColdGas); // Update the unstable gas after BH accretion
      if(DiscGasSum>0 && unstable_gas>0)
	  {
	    for(j=0; j<30; j++)
	    {	
		  assert(Gal[p].DiscGasMetals[j]<Gal[p].DiscGas[j]);
		  metallicity = get_metallicity(Gal[p].DiscGas[j], Gal[p].DiscGasMetals[j]);
		  ring_fraction = Gal[p].DiscGas[j] / DiscGasSum;
		
		  Gal[p].DiscGas[j] -= ring_fraction * unstable_gas;
		  Gal[p].DiscGasMetals[j] -= metallicity * ring_fraction * unstable_gas;
		  Gal[p].MetalsColdGas -= metallicity * ring_fraction * unstable_gas;
		  Gal[p].ColdGas -= ring_fraction * unstable_gas;

		  Gal[p].StellarMass += ring_fraction * unstable_gas * (1-RecycleFraction);
		  Gal[p].MetalsStellarMass += metallicity * ring_fraction * unstable_gas * (1-RecycleFraction);
		  Gal[p].SecularBulgeMass += ring_fraction * unstable_gas * (1-RecycleFraction);
	      Gal[p].SecularMetalsBulgeMass += metallicity * ring_fraction * unstable_gas * (1-RecycleFraction);

		  Gal[p].HotGas += ring_fraction * unstable_gas * RecycleFraction;
		  Gal[p].MetalsHotGas += metallicity * ring_fraction * unstable_gas * RecycleFraction;
	  	  Gal[p].MetalsHotGas += Yield * ring_fraction * unstable_gas;

		  Gal[p].SfrBulgeColdGas[step] += ring_fraction * unstable_gas;
		  Gal[p].SfrBulgeColdGasMetals[step] += metallicity * ring_fraction * unstable_gas;
		  Gal[p].SfrBulge[step] += ring_fraction * unstable_gas / dt;
	    }
	  }
	}
  }
}


