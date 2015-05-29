#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step)
{
  double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, SFE_H2, f_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscPre, ColdPre;
  int i;

  // Checks that the deconstructed disc is being treated properly and not generating NaNs
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);

  f_H2_const = 1.306e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*0.92);
  SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR;

  // Initialise variables
  strdot = 0.0;
  stars_sum = 0.0;

  Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
  Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas;

  for(i=0; i<30; i++)
  {
	if(Gal[p].Vvir>0) //There can be galaxies passed through they don't appear to be real
	{
		area = M_PI * (pow(DiscBinEdge[i+1]/Gal[p].Vvir, 2.0) - pow(DiscBinEdge[i]/Gal[p].Vvir, 2.0));
		f_H2 = f_H2_const * pow(pow(Gal[p].DiscGas[i]/area, 2.0) + 0.1*Gal[p].DiscGas[i]/area * pow(Gal[p].DiscStars[i]*Gal[p].DiscStars[0], 0.5)/area, 0.92);
		if(f_H2 > 0.0)
		{
		  assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
		  f_H2 = 0.75 * 1.0/(1.0/f_H2 + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i]); //Changes f_H2 from being H2/HI to H2/Cold Gas
		}
		else
		  f_H2 = 0.0;
	
		strdot = SfrEfficiency * SFE_H2 * f_H2 * Gal[p].DiscGas[i]/1.3;
    }

	stars = strdot * dt;
	
	if(stars < 0.0)
	  stars = 0.0;

    if(stars > Gal[p].DiscGas[i])
      stars = Gal[p].DiscGas[i];

    if(SupernovaRecipeOn == 1 && Gal[p].DiscGas[i] > 0.0 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
		if(Gal[p].Vvir > 0.0)
		  area = M_PI * (pow(DiscBinEdge[i+1]/Gal[p].Vvir, 2.0) - pow(DiscBinEdge[i]/Gal[p].Vvir, 2.0));
		else
		  area = M_PI * (pow(DiscBinEdge[i+1]/Gal[p].Vmax, 2.0) - pow(DiscBinEdge[i]/Gal[p].Vmax, 2.0));
		Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
        reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area/1.3);

		// Can't use more cold gas than is available, so balance SF and feedback 
	    if((stars + reheated_mass) > Gal[p].DiscGas[i] && (stars + reheated_mass) > 0.0)
	    {
	      fac = Gal[p].DiscGas[i] / (stars + reheated_mass);
	      stars *= fac;
	      reheated_mass *= fac;
	    }

	    if(stars<1e-8)
	    {
	      stars = 1e-8;
		  reheated_mass = Gal[p].DiscGas[i] - (1-RecycleFraction)*stars;
	    }
	
	    ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) - FeedbackReheatingEpsilon) * stars;
	    if(ejected_mass < 0.0)
	        ejected_mass = 0.0;
	
		assert(RecycleFraction*stars+reheated_mass < 1.001*Gal[p].DiscGas[i]);
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

	stars_sum += stars;

	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

    // Update for star formation
    metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    update_from_star_formation(p, stars, metallicity, i);

	if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
	  reheated_mass = Gal[p].DiscGas[i];

	// These checks ensure numerical uncertainties don't blow up	
    assert(abs(Gal[p].ColdGas-ColdPre) <= 1.001*abs(Gal[p].DiscGas[i]-DiscPre) && abs(Gal[p].ColdGas-ColdPre) >= 0.999*abs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
 
	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

    // Update from SN feedback
	metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    update_from_feedback(p, centralgal, reheated_mass, ejected_mass, metallicity, i);

	assert(abs(Gal[p].ColdGas-ColdPre) <= 1.001*abs(Gal[p].DiscGas[i]-DiscPre) && abs(Gal[p].ColdGas-ColdPre) >= 0.999*abs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);

	// Inject new metals from SN II
	if(SupernovaRecipeOn == 1 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
	    Gal[p].DiscGasMetals[i] += Yield * stars;
	    Gal[p].MetalsColdGas += Yield * stars;
  	  }
	  else
		Gal[p].MetalsHotGas += Yield * stars;
	}
	assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
  }

  // Update the star formation rate 
  Gal[p].SfrDisk[step] += stars_sum / dt;

  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);

  // Check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, halonr, time, dt, step);
}



void update_from_star_formation(int p, double stars, double metallicity, int i)
{
  //double DiscNewSpin[3], OldDisc[30], OldDiscMetals[30];
	
  // update gas and metals from star formation 
  Gal[p].DiscGas[i] -= (1 - RecycleFraction) * stars;
  Gal[p].DiscGasMetals[i] -= metallicity * (1 - RecycleFraction) * stars;

 // Here's where I need to consider the offset of the gas and star discs
  Gal[p].DiscStars[i] += (1 - RecycleFraction) * stars;
  Gal[p].DiscStarsMetals[i] += metallicity * (1 - RecycleFraction) * stars;
  
  Gal[p].ColdGas -= (1 - RecycleFraction) * stars;
  Gal[p].MetalsColdGas -= metallicity * (1 - RecycleFraction) * stars;
  Gal[p].StellarMass += (1 - RecycleFraction) * stars;
  Gal[p].MetalsStellarMass += metallicity * (1 - RecycleFraction) * stars;

  if(Gal[p].DiscGas[i] < 0.0){
	printf("DiscGas in update_SF...%e\n", Gal[p].DiscGas[i]);
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;}

  if(Gal[p].ColdGas < 0.0){
	printf("DiscGas in update_SF...%e\n", Gal[p].DiscGas[i]);
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;}
  
}



void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double metallicity, int i)
{
  double metallicityHot;

  // Check first just to be sure 
  assert(reheated_mass <= Gal[p].DiscGas[i]);

  if(SupernovaRecipeOn == 1)
  {
    Gal[p].ColdGas -= reheated_mass;
    Gal[p].DiscGas[i] -= reheated_mass;
    Gal[centralgal].HotGas += reheated_mass;

	if(Gal[p].DiscGas[i]>0.0)
	{
	  Gal[p].MetalsColdGas -= metallicity * reheated_mass;
      Gal[centralgal].MetalsHotGas += metallicity * reheated_mass;
      Gal[p].DiscGasMetals[i] -= metallicity * reheated_mass;
	}
	else
	{
	  Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[i];
      Gal[centralgal].MetalsHotGas += Gal[p].DiscGasMetals[i];
	  Gal[p].DiscGasMetals[i] = 0.0;
    }

    if(ejected_mass > Gal[centralgal].HotGas)
      ejected_mass = Gal[centralgal].HotGas;

    Gal[centralgal].HotGas -= ejected_mass;
    Gal[centralgal].EjectedMass += ejected_mass;

	if(Gal[centralgal].HotGas>0.0)
	{
	  metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
      Gal[centralgal].MetalsHotGas -= metallicityHot * ejected_mass;
      Gal[centralgal].MetalsEjectedMass += metallicityHot * ejected_mass;
	}
	else
	{
      Gal[centralgal].MetalsEjectedMass += Gal[centralgal].MetalsHotGas;
      Gal[centralgal].MetalsHotGas = 0.0;
	}

    Gal[p].OutflowRate += reheated_mass;    
  }

  if(Gal[p].DiscGas[i] < 0.0)
  {
    printf("DiscGas in update_feedback...%e\n", Gal[p].DiscGas[i]);
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;
  }

  if(Gal[p].ColdGas < 0.0)
  {
    printf("ColdGas in update_feedback...%e\n", Gal[p].ColdGas);
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;
  }

  if(Gal[p].HotGas < 0.0)
  {
	printf("HotGas in update_feedback...%e\n", Gal[p].HotGas);
	Gal[p].HotGas=0.0;
	Gal[p].MetalsHotGas=0.0;
  }
}


