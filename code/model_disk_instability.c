#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "core_allvars.h"
#include "core_proto.h"



void check_disk_instability(int p, int centralgal, int halonr, double time, double dt, int step)
{
  double Mcrit, gas_fraction, unstable_gas, unstable_gas_fraction, unstable_stars, diskmass, metallicity, DiscGasSum, DiscDiff;
  double star_fraction, ring_fraction, GasBeforeBH;
  int j;

  // Here we calculate the stability of the stellar and gaseous disk as discussed in Mo, Mao & White (1998).
  // For unstable stars and gas, we transfer the required ammount to the bulge to make the disk stable again

  // Check that Cold Gas has been treated properly prior to this function
  DiscGasSum = 0.0;
  for(j=0; j<30; j++)
	DiscGasSum += Gal[p].DiscGas[j];
  if(DiscGasSum > 1.001*Gal[p].ColdGas || DiscGasSum < Gal[p].ColdGas/1.001)
  {
	printf("Gas uneven at start of instability check....%e\t%e\n", DiscGasSum, Gal[p].ColdGas);
	//ABORT(1);
  }
  
  DiscDiff = DiscGasSum - Gal[p].ColdGas;


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
      // Use disk metallicity here !
      //metallicity = get_metallicity(Gal[p].StellarMass - (Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass), Gal[p].MetalsStellarMass - (Gal[p].ClassicalMetalsBulgeMass + Gal[p].SecularMetalsBulgeMass));

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
      
      if ((Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass) / Gal[p].StellarMass > 1.0001 || (Gal[p].ClassicalMetalsBulgeMass + Gal[p].SecularMetalsBulgeMass) / Gal[p].MetalsStellarMass > 1.0001)
	    {
        printf("Instability: Mbulge > Mtot (stars or metals)\t%e\t%e\t%e\t%e\t%e\n", (Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass), Gal[p].StellarMass, (Gal[p].ClassicalMetalsBulgeMass + Gal[p].SecularMetalsBulgeMass), Gal[p].MetalsStellarMass, unstable_stars);
        // ABORT(96);
      }

    // burst excess gas and feed black hole (really need a dedicated model for bursts and BH growth here)
    if(unstable_gas > 0.0)
    {
      if(unstable_gas/Gal[p].ColdGas > 1.0001)
      {
        printf("unstable_gas > Gal[p].ColdGas\t%e\t%e\n", unstable_gas, Gal[p].ColdGas);
        // ABORT(97);
      }


	  unstable_gas_fraction = unstable_gas / Gal[p].ColdGas;
	  GasBeforeBH = Gal[p].ColdGas;

      //for(j=0; j<30; j++)
	//	if(Gal[p].DiscGas[j]>0.0 && unstable_gas > 0.0)
	  //    unstable_gas_fraction_arr[j] = unstable_gas_fraction; // Make the unstable fraction the same for each disc annulus
	    //else
		  //unstable_gas_fraction_arr[j] = 0.0;
		
      if(AGNrecipeOn > 0)
        grow_black_hole(p, unstable_gas_fraction);

	  // Check that Cold Gas has been treated properly by to this function
	  DiscGasSum = 0.0;
	  for(j=0; j<30; j++)
		DiscGasSum += Gal[p].DiscGas[j];
	  if(DiscGasSum > 1.02*Gal[p].ColdGas || DiscGasSum < Gal[p].ColdGas/1.02)
	  {
		printf("Gas uneven at after growing BH during instability check....%e\t%e\n", DiscGasSum, Gal[p].ColdGas);
		//ABORT(1);
	  }
	
	
	  // MY NEW SIMPLE WAY TO DEAL WITH AN INSTABILITY
	  unstable_gas = unstable_gas - (GasBeforeBH - Gal[p].ColdGas); // Update the unstable gas after BH accretion
	  //unstable_gas_fraction = unstable_gas / Gal[p].ColdGas;
	  if(DiscGasSum>0 && unstable_gas>0)
	  {
	    for(j=0; j<30; j++)
	    {
		  metallicity = get_metallicity(Gal[p].DiscStars[j], Gal[p].DiscStarsMetals[j]);
		  ring_fraction = Gal[p].DiscGas[j] / DiscGasSum;
		
		  Gal[p].DiscGas[j] -= ring_fraction * unstable_gas;
		  Gal[p].DiscGasMetals[j] -= metallicity * ring_fraction * unstable_gas;
		  Gal[p].MetalsColdGas -= metallicity * ring_fraction * unstable_gas;
		  Gal[p].ColdGas -= ring_fraction * unstable_gas;

		  Gal[p].StellarMass += ring_fraction * unstable_gas * (1-RecycleFraction);
		  Gal[p].MetalsStellarMass += metallicity * ring_fraction * unstable_gas * (1-RecycleFraction);
		  Gal[p].SecularBulgeMass += unstable_stars * ring_fraction * (1-RecycleFraction);
	      Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars * ring_fraction * (1-RecycleFraction);

		  Gal[p].HotGas += ring_fraction * unstable_gas * RecycleFraction;
		  Gal[p].MetalsHotGas += metallicity * ring_fraction * unstable_gas * RecycleFraction;
		  Gal[p].MetalsHotGas += Yield * ring_fraction * unstable_gas;

		  Gal[p].SfrBulgeColdGas[step] += ring_fraction * unstable_gas;
		  Gal[p].SfrBulgeColdGasMetals[step] += metallicity * ring_fraction * unstable_gas;
		  Gal[p].SfrBulge[step] += ring_fraction * unstable_gas / dt;
		
		  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
		  {
		    printf("HotGas in instability...%e\n", Gal[p].HotGas);
			printf("j, ring_fraction, metallicity, unstable_gas, DiscGasSum\n");
			printf("%d\t%e\t%e\t%e\t%e\n", j, ring_fraction, metallicity, unstable_gas, DiscGasSum);
		    //ABORT(1);
		  }
		}
	  }
	
	
	  //printf("Calling starburst from instability\n");
      //collisional_starburst_recipe(unstable_gas_fraction_arr, p, centralgal, time, dt, halonr, 1, step, unstable_gas_fraction);
    }

  // Check that Cold Gas has been treated properly by to this function
  DiscGasSum = 0.0;
  for(j=0; j<30; j++)
	DiscGasSum += Gal[p].DiscGas[j];
	
  //if(abs(DiscDiff) > 1.001*abs(DiscGasSum-Gal[p].ColdGas) || abs(DiscDiff) < abs(DiscGasSum-Gal[p].ColdGas)/1.001 || DiscDiff*(DiscGasSum-Gal[p].ColdGas)<0.0)
  //{
	//printf("Instability inducing absolute differences between ColdGas and DiscGas......%e\n", DiscDiff/(DiscGasSum-Gal[p].ColdGas));
	//ABORT(1);
  //}

  if(DiscGasSum > 1.001*Gal[p].ColdGas || DiscGasSum < Gal[p].ColdGas/1.001)
  {
	printf("Gas uneven at end of instability check....%e\t%e\n", DiscGasSum, Gal[p].ColdGas);
	//ABORT(1);
  }



  }

}
