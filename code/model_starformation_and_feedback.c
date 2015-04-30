#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "core_allvars.h"
#include "core_proto.h"



void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step)
{
  double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, SFE_H2, f_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscDiff, DiscPre, ColdPre;
  int i;

  DiscGasSum = 0.0;
  for(i=0; i<30; i++)
	DiscGasSum += Gal[p].DiscGas[i];
  //printf("Disc gas initial SF....%e\t%e\n", DiscGasSum, Gal[p].ColdGas);
  if(DiscGasSum > 1.001*Gal[p].ColdGas || DiscGasSum < Gal[p].ColdGas/1.001)
  {
	printf("Gas incorrect at start of SF....%e\t%e\n", DiscGasSum, Gal[p].ColdGas);
	//ABORT(1);
  }

  DiscDiff = (DiscGasSum - Gal[p].ColdGas);


  f_H2_const = 1.306e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*0.92);
  SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
  {
    printf("HotGas initial starformation_and_feedback...%e\n", Gal[p].HotGas);
    //ABORT(1);
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
		  f_H2 = 0.75 * 1.0/(1.0/f_H2 + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i]); //Changes f_H2 from being H2/HI to H2/Cold Gas
		else
		  f_H2 = 0.0;
	
		strdot = SfrEfficiency * SFE_H2 * f_H2 * Gal[p].DiscGas[i]/1.3;
	
		if(area!=area || area<0)
		{
			printf("DiscBinEdge, DiscBinEdge, Vvir, area...%e\t%e\t%e\t%e\n", DiscBinEdge[i], DiscBinEdge[i+1], Gal[p].Vvir, area);
			printf("Mvir, ColdGas...%e\t%e\n", Gal[p].Mvir, Gal[p].ColdGas);
			//ABORT(1);
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
		//ABORT(1);
	}
	
	if(stars < 0.0)
	  stars = 0.0;

    if(stars > Gal[p].DiscGas[i])
      stars = Gal[p].DiscGas[i];

    if(SupernovaRecipeOn == 1 && Gal[p].DiscGas[i] > 0.0 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
		area = M_PI * (pow(DiscBinEdge[i+1]/Gal[p].Vvir, 2.0) - pow(DiscBinEdge[i]/Gal[p].Vvir, 2.0));
		Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
        reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area/1.3);

		// can't use more cold gas than is available! so balance SF and feedback 
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
	
	    if(RecycleFraction*stars+reheated_mass > 1.001*Gal[p].DiscGas[i])
	    {
		  printf("Too much cold gas trying to be converted...stars/reheated/gas avail....%e\t%e\t%e\n", stars, reheated_mass, Gal[p].DiscGas[i]);
		  printf("ratio...%e\n", (RecycleFraction*stars+reheated_mass)/Gal[p].DiscGas[i]);
		  //ABORT(1);
		}
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




	//printf("stars, reheated, ejected.........%e\t%e\t%e\n", stars, reheated_mass, ejected_mass);

	stars_sum += stars;

    if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
    {
      printf("HotGas starformation_and_feedback (2)...%e\n", Gal[p].HotGas);
      //ABORT(1);
    }

	if(Gal[p].DiscGasMetals[i]>Gal[p].DiscGas[i])
	{
		printf("More metals than total gas before updating.......%e\t%e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);
		ABORT(1);
	}

	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

    // update for star formation
    metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    update_from_star_formation(p, stars, metallicity, i);

	if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
	  reheated_mass = Gal[p].DiscGas[i];

	if((Gal[p].ColdGas-ColdPre)/(Gal[p].DiscGas[i]-DiscPre) > 1.001 || (Gal[p].ColdGas-ColdPre)/(Gal[p].DiscGas[i]-DiscPre) < 0.999 || (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)<0.0)
	{
		printf("Update SF the cause of problem.....%e\n", (Gal[p].ColdGas-ColdPre)/(Gal[p].DiscGas[i]-DiscPre));
		//printf("ColdDiff, DiscDiff......%e\t%e\n", (Gal[p].ColdGas-ColdPre), (Gal[p].DiscGas[i]-DiscPre));
		//ABORT(1);
	}

    if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
    {
      printf("HotGas starformation_and_feedback (3)...%e\n", Gal[p].HotGas);
      //ABORT(1);
    }

	if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass > 1.0e-8)
	  printf("reheated mass too high in standard SF recipe....%e\n", reheated_mass/Gal[p].DiscGas[i]);

	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

	if(Gal[p].DiscGasMetals[i]>Gal[p].DiscGas[i])
	{
		printf("More metals than total gas between updates.......%e\t%e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);
		ABORT(1);
	}

    // update from SN feedback
	metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    update_from_feedback(p, centralgal, reheated_mass, ejected_mass, metallicity, i);

	if((Gal[p].ColdGas-ColdPre)/(Gal[p].DiscGas[i]-DiscPre) > 1.001 || (Gal[p].ColdGas-ColdPre)/(Gal[p].DiscGas[i]-DiscPre) < 0.999 || (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)<0.0)
	{
		printf("Update feedback the cause of problem.....%e\n", (Gal[p].ColdGas-ColdPre)/(Gal[p].DiscGas[i]-DiscPre));
		//printf("Cold, Disc, reheated.........%e\t%e\t%e\n", Gal[p].ColdGas, Gal[p].DiscGas[i], reheated_mass);
		//printf("Differences.....%e\t%e\n", (Gal[p].ColdGas-ColdPre), (Gal[p].DiscGas[i]-DiscPre));
		//ABORT(1);
	}

    if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
    {
      printf("HotGas starformation_and_feedback (4)...%e\n", Gal[p].HotGas);
	  printf("reheated, ejected...%e\t%e\n", reheated_mass, ejected_mass);
      //ABORT(1);
    }

	if(Gal[p].DiscGasMetals[i]>Gal[p].DiscGas[i])
	{
		printf("More metals than total gas after updating.......%e\t%e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);
		ABORT(1);
	}

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
	
	if(Gal[p].DiscGasMetals[i]>Gal[p].DiscGas[i])
	{
		printf("More metals than total after injecting.......%e\t%e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);
		ABORT(1);
	}
  }

  // update the star formation rate 
  Gal[p].SfrDisk[step] += stars_sum / dt;

  DiscGasSum = 0.0;
  for(i=0; i<30; i++)
	DiscGasSum += Gal[p].DiscGas[i];
	
  //if(abs(DiscDiff) > 1.001*abs(DiscGasSum-Gal[p].ColdGas) || abs(DiscDiff) < abs(DiscGasSum-Gal[p].ColdGas)/1.001 || DiscDiff*(DiscGasSum-Gal[p].ColdGas)<0.0)
  //{
	//printf("SF inducing absolute differences between ColdGas and DiscGas......%e\n", DiscDiff/(DiscGasSum-Gal[p].ColdGas));
	//ABORT(1);
  //}

  if(DiscGasSum > 1.001*Gal[p].ColdGas || DiscGasSum < Gal[p].ColdGas/1.001)
  {
	printf("Gas incorrectly dealt with in SF PRE INSTABILITY....%e\t%e\n", DiscGasSum, Gal[p].ColdGas);
	//ABORT(1);
  }




  // check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, halonr, time, dt, step);

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0)
  {
    printf("HotGas final starformation_and_feedback...%e\n", Gal[p].HotGas);
    //ABORT(1);
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
  double metallicityHot, ColdInit, DiscInit;

  // check first just to be sure 
  if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass > 1.0e-8)
  {
    printf("Something strange here (SF2)....%d\t%e\t%e\n", i, reheated_mass, Gal[p].ColdGas);
    ABORT(19);
    //reheated_mass = Gal[p].ColdGas;
  }

  ColdInit = Gal[p].ColdGas;
  DiscInit = Gal[p].DiscGas[i];

  if(SupernovaRecipeOn == 1)
  {
    Gal[p].ColdGas -= reheated_mass;
    Gal[p].DiscGas[i] -= reheated_mass;
    Gal[centralgal].HotGas += reheated_mass;

    //if((Gal[p].ColdGas-ColdInit) != (Gal[p].DiscGas[i]-DiscInit))
    //{
	  //printf("Update feedback CONFIRMED cause of problem.....%e\n", (Gal[p].ColdGas-ColdInit)/(Gal[p].DiscGas[i]-DiscInit));
	  //printf("Cold, Disc, reheated.........%e\t%e\t%e\n", Gal[p].ColdGas, Gal[p].DiscGas[i], reheated_mass);
	  //printf("Differences.....%e\t%e\n", (Gal[p].ColdGas-ColdInit), (Gal[p].DiscGas[i]-DiscInit));
	  //ABORT(1);
    //}

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


