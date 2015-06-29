#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


void check_disk_instability(int p, int centralgal, double time, double dt, int step)
{
	// New treatment of instabilities based on the Toomre Q parameter
	double Q_star, Q_gas, Sigma_star, Sigma_gas, area, radius, V_vir;
	double unstable_gas, unstable_stars, gas_sink, metallicity, gas_sink_sum, gas_sf;
	double stars, reheated_mass, ejected_mass, stars_sum, Sigma_0gas, fac;
	double NewStars[30], NewStarsMetals[30];
	int i;
	
	if(Gal[p].Vvir>0.0)
		V_vir = Gal[p].Vvir;
	else
		V_vir = Gal[p].Vmax;
	
	// Deal with gaseous instabilities
	stars_sum = 0.0;
	gas_sink_sum = 0.0;
	for(i=29; i>=0; i--)
	{
		area = M_PI * (pow(DiscBinEdge[i+1]/V_vir, 2.0) - pow(DiscBinEdge[i]/V_vir, 2.0)) * pow(CM_PER_MPC, 2.0);
		Sigma_gas = Gal[p].DiscGas[i] * 1e10*SOLAR_MASS / area;
		radius = pow((DiscBinEdge[i]*DiscBinEdge[i] + DiscBinEdge[i+1]*DiscBinEdge[i+1])/2.0, 0.5)/V_vir * CM_PER_MPC;
		Q_gas = pow(V_vir*1e5, 2.0) / (M_PI * GRAVITY * radius * Sigma_gas);
		
		if(Q_gas<1.0)
		{
			unstable_gas = Gal[p].DiscGas[i] - pow(V_vir*1e5, 2.0)*area / (M_PI * GRAVITY * radius * 1e10*SOLAR_MASS);
			metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
			
			// Let gas sink
			gas_sink = BlackHoleGrowthRate * unstable_gas / (1.0 + pow(280.0 / V_vir, 2.0));
			Gal[p].DiscGas[i] -= gas_sink;
			Gal[p].DiscGasMetals[i] -= metallicity * gas_sink;
			gas_sink_sum += gas_sink;
			// if(i!=0)
			// 			{
			// 				Gal[p].DiscGas[i-1] += gas_sink;
			// 				Gal[p].DiscGasMetals[i-1] += metallicity * gas_sink;
			// 			}
			// 			else
			// 			{
				Gal[p].BlackHoleMass += gas_sink;
				Gal[p].ColdGas -= gas_sink;
				Gal[p].MetalsColdGas -= metallicity * gas_sink;
			// }
		
			// Calculate new stars formed in that annulus
			gas_sf = unstable_gas - gas_sink;
			stars = unstable_gas - gas_sink;
			if(Gal[p].DiscGas[i] > 0.0 && stars > 0.0) // Quasar feedback could blow out the unstable gas
			{
				if(SupernovaRecipeOn == 1)
				{
					Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
			        reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area/1.3*pow(CM_PER_MPC, 2.0));
			
					// Can't use more cold gas than is available, so balance SF and feedback 
				    if((stars + reheated_mass) > gas_sf && (stars + reheated_mass) > 0.0)
				    {
				    	fac = gas_sf / (stars + reheated_mass);
				    	stars *= fac;
				    	reheated_mass *= fac;
				    }
				
					if(stars<1e-8)
				    {
						if(unstable_gas - gas_sink >= 1e-8)
						{
				    		stars = 1e-8;
							reheated_mass = gas_sf - (1-RecycleFraction)*stars;
						}
						else
						{
							stars = gas_sf;
							reheated_mass = 0.0;
						}
						ejected_mass = 0.0;
				    }
					else
					{
						ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (V_vir * V_vir) - FeedbackReheatingEpsilon) * stars;
					    if(ejected_mass < 0.0)
					        ejected_mass = 0.0;
					}
				}
				else
				{
					reheated_mass = 0.0;// not a great treatment right now
					ejected_mass = 0.0;
				}
				
				stars_sum += stars;
				
				NewStars[i] = (1 - RecycleFraction) * stars;
				NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
			    update_from_star_formation(p, stars, metallicity, i);
			
				if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
				  reheated_mass = Gal[p].DiscGas[i];
				
				metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
			    update_from_feedback(p, centralgal, reheated_mass, ejected_mass, metallicity, i);
			
				// Update metals from SN II feedback
				Gal[p].DiscGasMetals[i] += Yield * stars;
			    Gal[p].MetalsColdGas += Yield * stars;
			}
						
			//if(i==0 && AGNrecipeOn > 0)  // Deal with quasar feedback
				
			
		}
	}
	if(gas_sink_sum > 0.0 && AGNrecipeOn > 0)
		quasar_mode_wind(p, gas_sink_sum);
	
	// Merge new-star disc with previous stellar disc
	if(stars_sum>0.0)
	{
		combine_stellar_discs(p, NewStars, NewStarsMetals);
		Gal[p].SfrDisk[step] += stars_sum / dt; // Some of these stars may quickly be transferred to the bulge, so simply updating SfrDisk might be crude
	}
				  
	
	// Deal with stellar instabilities
	for(i=29; i>=0; i--)
	{
		area = M_PI * (pow(DiscBinEdge[i+1]/V_vir, 2.0) - pow(DiscBinEdge[i]/V_vir, 2.0)) * pow(CM_PER_MPC, 2.0);
		Sigma_star = Gal[p].DiscStars[i] * 1e10*SOLAR_MASS / area;
		radius = pow((DiscBinEdge[i]*DiscBinEdge[i] + DiscBinEdge[i+1]*DiscBinEdge[i+1])/2.0, 0.5)/V_vir * CM_PER_MPC;
		Q_star = 0.126 * pow(V_vir*1e5, 2.0) / (GRAVITY * radius * Sigma_star);
		
		if(Q_star<1.0)
		{
			unstable_stars = Gal[p].DiscStars[i] - 0.126*pow(V_vir*1e5, 2.0)*area / (GRAVITY * 1e10*SOLAR_MASS * radius);
			metallicity = get_metallicity(Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			Gal[p].DiscStars[i] -= unstable_stars;
			Gal[p].DiscStarsMetals[i] -= metallicity * unstable_stars;
			if(i!=0)
			{
				Gal[p].DiscStars[i-1] += unstable_stars;
				Gal[p].DiscStarsMetals[i-1] += metallicity * unstable_stars;
			}
			else
			{
				Gal[p].SecularBulgeMass += unstable_stars;
				Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars;
			}
		}
	}
}


// THIS IS NO LONGER USED
// void check_disk_instability_old(int p, int centralgal, int halonr, double time, double dt, int step)
// {
//   double Mcrit, gas_fraction, unstable_gas, unstable_gas_fraction, unstable_stars, diskmass, metallicity, DiscGasSum, DiscStarSum;
//   double star_fraction, ring_fraction, GasBeforeBH;
//   int j;
// 
//   // Here we calculate the stability of the stellar and gaseous disk as discussed in Mo, Mao & White (1998).
//   // For unstable stars and gas, we transfer the required ammount to the bulge to make the disk stable again
// 
//   // Check that Cold Gas has been treated properly prior to this function
//   DiscGasSum = get_disc_gas(p);
//   assert(DiscGasSum <= 1.001*Gal[p].ColdGas || DiscGasSum >= Gal[p].ColdGas/1.001);
// 
//   // Disk mass has to be > 0.0 !
//   diskmass = Gal[p].ColdGas + (Gal[p].StellarMass - Gal[p].ClassicalBulgeMass - Gal[p].SecularBulgeMass);
//   if(diskmass > 0.0)
//   {
//     // calculate critical disk mass
//     Mcrit = Gal[p].Vmax * Gal[p].Vmax * (3.0 * Gal[p].DiskScaleRadius) / G;
//     if(Mcrit > diskmass)
//       Mcrit = diskmass;
//     
//     // use Disk mass here !
//     gas_fraction   = Gal[p].ColdGas / diskmass;
//     unstable_gas   = gas_fraction * (diskmass - Mcrit);
//     star_fraction  = 1.0 - gas_fraction;
//     unstable_stars = star_fraction * (diskmass - Mcrit);
// 
//     // add excess stars to the bulge
//     if(unstable_stars > 0.0)
//     {
//       for(j=0; j<30; j++)
// 	  {
// 		metallicity = get_metallicity(Gal[p].DiscStars[j], Gal[p].DiscStarsMetals[j]);
// 	    ring_fraction = Gal[p].DiscStars[j] / (star_fraction*diskmass);
// 		Gal[p].DiscStars[j] -= unstable_stars * ring_fraction;
// 		Gal[p].DiscStarsMetals[j] -= metallicity * unstable_stars * ring_fraction;
// 		Gal[p].SecularBulgeMass += unstable_stars * ring_fraction;
// 	    Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars * ring_fraction;
// 	  }
//     }  
//       
//     // Need to fix this. Excluded for now.
//     // Gal[p].mergeType = 3;  // mark as disk instability partial mass transfer
//     // Gal[p].mergeIntoID = NumGals + p - 1;      
//     
//     DiscStarSum = get_disc_stars(p);
//     assert((Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass) <= 1.001*Gal[p].StellarMass);
//     //if(DiscStarSum > 1.001*(Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass) || DiscStarSum < (Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass)/1.001)
//         //printf("DiscStarSum, StellarMass, BulgeSum, diff = %e, %e, %e, %e\n", DiscStarSum, Gal[p].StellarMass, Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass, Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass);
//     if(DiscStarSum>0.0) assert(DiscStarSum <= 1.001*(Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass) && DiscStarSum >= (Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass)/1.001);
// 	assert((Gal[p].ClassicalMetalsBulgeMass + Gal[p].SecularMetalsBulgeMass) <= 1.001*Gal[p].MetalsStellarMass);
// 
//     // burst excess gas and feed black hole (really need a dedicated model for bursts and BH growth here)
//     if(unstable_gas > 0.0)
//     {
// 	  assert(unstable_gas/Gal[p].ColdGas < 1.0001);
// 
// 	  unstable_gas_fraction = unstable_gas / Gal[p].ColdGas;
// 	  GasBeforeBH = Gal[p].ColdGas;
// 		
//       //if(AGNrecipeOn > 0)
//         //grow_black_hole(p, unstable_gas_fraction);
// 	
//       // MY NEW SIMPLE WAY TO DEAL WITH AN INSTABILITY
//       //unstable_gas = unstable_gas - (GasBeforeBH - Gal[p].ColdGas); // Update the unstable gas after BH accretion
//       if(DiscGasSum>0 && unstable_gas>0)
// 	  {
// 	    for(j=0; j<30; j++)
// 	    {	
// 		  assert(Gal[p].DiscGasMetals[j]<=Gal[p].DiscGas[j]);
// 		  metallicity = get_metallicity(Gal[p].DiscGas[j], Gal[p].DiscGasMetals[j]);
// 		  ring_fraction = Gal[p].DiscGas[j] / DiscGasSum;
// 		
// 		  Gal[p].DiscGas[j] -= ring_fraction * unstable_gas;
// 		  Gal[p].DiscGasMetals[j] -= metallicity * ring_fraction * unstable_gas;
// 		  Gal[p].MetalsColdGas -= metallicity * ring_fraction * unstable_gas;
// 		  Gal[p].ColdGas -= ring_fraction * unstable_gas;
// 
// 		  Gal[p].StellarMass += ring_fraction * unstable_gas * (1-RecycleFraction);
// 		  Gal[p].MetalsStellarMass += metallicity * ring_fraction * unstable_gas * (1-RecycleFraction);
// 		  Gal[p].SecularBulgeMass += ring_fraction * unstable_gas * (1-RecycleFraction);
// 	      Gal[p].SecularMetalsBulgeMass += metallicity * ring_fraction * unstable_gas * (1-RecycleFraction);
// 
// 		  Gal[p].HotGas += ring_fraction * unstable_gas * RecycleFraction;
// 		  Gal[p].MetalsHotGas += metallicity * ring_fraction * unstable_gas * RecycleFraction;
// 	  	  Gal[p].MetalsHotGas += Yield * ring_fraction * unstable_gas;
// 
// 		  Gal[p].SfrBulgeColdGas[step] += ring_fraction * unstable_gas;
// 		  Gal[p].SfrBulgeColdGasMetals[step] += metallicity * ring_fraction * unstable_gas;
// 		  Gal[p].SfrBulge[step] += ring_fraction * unstable_gas / dt;
// 	    }
// 	  }
// 	}
//   }
// }
// 
// 
