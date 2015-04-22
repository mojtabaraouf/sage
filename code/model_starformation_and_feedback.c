#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "core_allvars.h"
#include "core_proto.h"



void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step)
{
  double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, SFE_H2, f_H2, f_H2_const, Sigma_0gas;
  int i;

  f_H2_const = 1.306e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*0.92);
  SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
  {
    printf("HotGas initial starformation_and_feedback...%e\n", Gal[p].HotGas);
    ABORT(1);
  }
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
		  f_H2 = 3.0 * 1.0/(1.0/f_H2 + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i]); //Changes f_H2 from being H2/HI to H2/Cold Gas
		else
	      f_H2 = 0.0;
		strdot = SFE_H2 * f_H2 * Gal[p].DiscGas[i];
	
		if(area!=area || area<0)
		{
			printf("DiscBinEdge, DiscBinEdge, Vvir, area...%e\t%e\t%e\t%e\n", DiscBinEdge[i], DiscBinEdge[i+1], Gal[p].Vvir, area);
			printf("Mvir, ColdGas...%e\t%e\n", Gal[p].Mvir, Gal[p].ColdGas);
			ABORT(1);
		}
	
		//printf("Disc properties\n");
		//printf("%e", Gal[p].DiscGas[i]);
		//printf("\n");
		//printf("%e", Gal[p].DiscStars[i]);
		//printf("\n");
		//printf("%e", Gal[p].DiscStars[0]);
		//printf("\n");
	

    }

	stars = strdot * dt;
	
	if(stars!=stars || stars<-1e-30)
	{
		printf("DiscGas, metals, area, f_H2...%e\t%e\t%e\t%e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i], area, f_H2);
		printf("stars, strdot, dt...%e\t%e\t%e\n", stars, strdot, dt);
		ABORT(1);
	}
	
	if(stars < 0.0)
	  stars = 0.0;

	if(SupernovaRecipeOn == 1 && Gal[p].DiscGas[i] > 0.0)
	{
	  Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
	  reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area);
	}
	else
	  reheated_mass = 0.0;

	if(reheated_mass < 0.0)
	{
	  printf("Something strange here (SF1)....\n");
	  ABORT(32);
	  reheated_mass = 0.0;
	}

	// cant use more cold gas than is available! so balance SF and feedback 
	if((stars + reheated_mass) > Gal[p].DiscGas[i] && (stars + reheated_mass) > 0.0)
	{
	  fac = Gal[p].DiscGas[i] / (stars + reheated_mass);
	  stars *= fac;
	  reheated_mass *= fac;
	}

	// determine ejection
	if(SupernovaRecipeOn == 1)
	{
	  if(Gal[centralgal].Vvir > 0.0)
		ejected_mass = 
			(FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) -
				FeedbackReheatingEpsilon) * stars;
	  else
		ejected_mass = 0.0;

      if(ejected_mass < 0.0)
        ejected_mass = 0.0;
	  }
    else
      ejected_mass = 0.0;

	stars_sum += stars;

    if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
    {
      printf("HotGas starformation_and_feedback (2)...%e\n", Gal[p].HotGas);
      ABORT(1);
    }

    // update for star formation
    metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    update_from_star_formation(p, stars, metallicity, i);

    if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
    {
      printf("HotGas starformation_and_feedback (3)...%e\n", Gal[p].HotGas);
      ABORT(1);
    }

    // recompute the metallicity of the cold phase
    metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);

	if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass > 1.0e-8)
	  printf("reheated mass too high in standard SF recipe....%e\n", reheated_mass/Gal[p].DiscGas[i]);

    // update from SN feedback 
    update_from_feedback(p, centralgal, reheated_mass, ejected_mass, metallicity, i);

    if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
    {
      printf("HotGas starformation_and_feedback (4)...%e\n", Gal[p].HotGas);
	  printf("reheated, ejected...%e\t%e\n", reheated_mass, ejected_mass);
      ABORT(1);
    }

	// Inject new metals from SN II
	if(SupernovaRecipeOn == 1)
	{
	  Gal[p].DiscGasMetals[i] += Yield * stars;
	  Gal[p].MetalsColdGas += Yield * stars;
	}
  }

  // update the star formation rate 
  Gal[p].SfrDisk[step] += stars_sum / dt;

  // check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, halonr, time, dt, step);

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
  {
    printf("HotGas final starformation_and_feedback...%e\n", Gal[p].HotGas);
    ABORT(1);
  }
}



void update_from_star_formation(int p, double stars, double metallicity, int i)
{
  // update gas and metals from star formation 
  Gal[p].DiscGas[i] -= (1 - RecycleFraction) * stars;
  Gal[p].DiscGasMetals[i] -= metallicity * (1 - RecycleFraction) * stars;
  Gal[p].DiscStars[i] += (1 - RecycleFraction) * stars;
  Gal[p].DiscStarsMetals[i] += metallicity * (1 - RecycleFraction) * stars;
  
  Gal[p].ColdGas -= (1 - RecycleFraction) * stars;
  Gal[p].MetalsColdGas -= metallicity * (1 - RecycleFraction) * stars;
  Gal[p].StellarMass += (1 - RecycleFraction) * stars;
  Gal[p].MetalsStellarMass += metallicity * (1 - RecycleFraction) * stars;

  if(Gal[p].DiscGas[i] < 0){
	printf("DiscGas in update_SF...%e\n", Gal[p].DiscGas[i]);
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;}

  if(Gal[p].ColdGas < 0){
	printf("DiscGas in update_SF...%e\n", Gal[p].DiscGas[i]);
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;}
  
}



void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double metallicity, int i)
{
  double metallicityHot;

  // check first just to be sure 
  if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass > 1.0e-8)
  {
    printf("Something strange here (SF2)....%d\t%e\t%e\n", i, reheated_mass, Gal[p].ColdGas);
    ABORT(19);
    reheated_mass = Gal[p].ColdGas;
  }

  if(SupernovaRecipeOn == 1)
  {
    Gal[p].ColdGas -= reheated_mass;
    Gal[p].MetalsColdGas -= metallicity * reheated_mass;

    Gal[p].DiscGas[i] -= reheated_mass;
    Gal[p].DiscGasMetals[i] -= metallicity * reheated_mass;

    Gal[centralgal].HotGas += reheated_mass;
    Gal[centralgal].MetalsHotGas += metallicity * reheated_mass;

    if(ejected_mass > Gal[centralgal].HotGas)
      ejected_mass = Gal[centralgal].HotGas;
    metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);

    Gal[centralgal].HotGas -= ejected_mass;
    Gal[centralgal].MetalsHotGas -= metallicityHot * ejected_mass;
    Gal[centralgal].EjectedMass += ejected_mass;
    Gal[centralgal].MetalsEjectedMass += metallicityHot * ejected_mass;

    Gal[p].OutflowRate += reheated_mass;    
  }

  if(Gal[p].DiscGas[i] < -1e-10)
    printf("DiscGas in update_feedback...%e\n", Gal[p].DiscGas[i]);
  if(Gal[p].DiscGas[i] < 0){
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;}

  if(Gal[p].ColdGas < -1e-10)
    printf("ColdGas in update_feedback...%e\n", Gal[p].ColdGas);
  if(Gal[p].ColdGas < 0){
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;}

  if(Gal[p].HotGas < -1e-10)
    printf("HotGas in update_feedback...%e\n", Gal[p].HotGas);
  if(Gal[p].HotGas < 0){
	Gal[p].HotGas=0.0;
	Gal[p].MetalsHotGas=0.0;}

}


