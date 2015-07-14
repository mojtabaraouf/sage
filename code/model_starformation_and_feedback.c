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
    double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, SFE_H2, f_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscPre, ColdPre;//, cos_theta_gas_stars;
  double NewStars[30], NewStarsMetals[30];
  int i;

  // Checks that the deconstructed disc is being treated properly and not generating NaNs
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  f_H2_const = 1.306e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*0.92);
  SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR;

  // Initialise variables
  strdot = 0.0;
  stars_sum = 0.0;

  Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
  Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas;
    
  //cos_theta_gas_stars = Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2];

  for(i=0; i<30; i++)
  {
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		
	if(Gal[p].Vvir>0) //There can be galaxies passed through they don't appear to be real
	{
		area = M_PI * (pow(DiscBinEdge[i+1]/Gal[p].Vvir, 2.0) - pow(DiscBinEdge[i]/Gal[p].Vvir, 2.0));
        //if(cos_theta_gas_stars >= 0.9397)
            f_H2 = f_H2_const * pow(pow(Gal[p].DiscGas[i]/area, 2.0) + 0.1*Gal[p].DiscGas[i]/area * pow(Gal[p].DiscStars[i]*Gal[p].DiscStars[0], 0.5)/area, 0.92);
        //else
            //f_H2 = f_H2_const * pow(Gal[p].DiscGas[i]/area, 2.0*0.92);
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
    else // I haven't actually dealt with the situation of Supernovae being turned off here.  But do I even want to turn SN off?
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
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	NewStars[i] = (1 - RecycleFraction) * stars;
	NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
	assert(NewStarsMetals[i] <= NewStars[i]);
    update_from_star_formation(p, stars, metallicity, i);

	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);

	if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
	  reheated_mass = Gal[p].DiscGas[i];

	// These checks ensure numerical uncertainties don't blow up	
    assert(abs(Gal[p].ColdGas-ColdPre) <= 1.001*abs(Gal[p].DiscGas[i]-DiscPre) && abs(Gal[p].ColdGas-ColdPre) >= 0.999*abs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
 
	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

    // Update from SN feedback
	metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
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
	  {
		if(Gal[centralgal].HotGas>0.0)
			Gal[centralgal].MetalsHotGas += Yield * stars;
		else
			Gal[centralgal].MetalsEjectedMass += Yield * stars;
	  }
	}
	assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
  }

  // Sum stellar discs together
  combine_stellar_discs(p, NewStars, NewStarsMetals);
    
  for(i=0; i<30; i++){
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}

  // Update the star formation rate 
  Gal[p].SfrDisk[step] += stars_sum / dt;

  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);

  for(i=0; i<30; i++){
	metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}

  // Check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, time, dt, step);

  for(i=0; i<30; i++){
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);}

  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}



void update_from_star_formation(int p, double stars, double metallicity, int i)
{
  // In older SAGE, this updated the gas and stellar components.  It only does the gas component now due to the way in which discs are handled.
  
  //double DiscNewSpin[3], OldDisc[30], OldDiscMetals[30];
	
  // update gas and metals from star formation 
  Gal[p].DiscGas[i] -= (1 - RecycleFraction) * stars;
  Gal[p].DiscGasMetals[i] -= metallicity * (1 - RecycleFraction) * stars;

  if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i])
  	printf("update_from_star_formation report -- gas metals, gas mass = %e, %e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);

 // Here's where I need to consider the offset of the gas and star discs
  //Gal[p].DiscStars[i] += (1 - RecycleFraction) * stars;
  //Gal[p].DiscStarsMetals[i] += metallicity * (1 - RecycleFraction) * stars;
  
  Gal[p].ColdGas -= (1 - RecycleFraction) * stars;
  Gal[p].MetalsColdGas -= metallicity * (1 - RecycleFraction) * stars;
  //Gal[p].StellarMass += (1 - RecycleFraction) * stars;
  //Gal[p].MetalsStellarMass += metallicity * (1 - RecycleFraction) * stars;

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
  metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);

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
	  assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
	}
	else
	{
	  Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[i];
      Gal[centralgal].MetalsHotGas += Gal[p].DiscGasMetals[i];
	  Gal[p].DiscGasMetals[i] = 0.0;
    }

    if(ejected_mass > Gal[centralgal].HotGas)
      ejected_mass = Gal[centralgal].HotGas;

	metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    Gal[centralgal].HotGas -= ejected_mass;
    Gal[centralgal].EjectedMass += ejected_mass;

	if(Gal[centralgal].HotGas>0.0)
	{
      Gal[centralgal].MetalsHotGas -= metallicityHot * ejected_mass;
      Gal[centralgal].MetalsEjectedMass += metallicityHot * ejected_mass;
	}
	else
	{
      Gal[centralgal].MetalsEjectedMass += Gal[centralgal].MetalsHotGas;
      Gal[centralgal].MetalsHotGas = 0.0;
	}
	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);

    Gal[p].OutflowRate += reheated_mass;    
  }

  metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);

  if(Gal[p].DiscGas[i] < 0.0)
  {
    printf("DiscGas in update_feedback...%e\n", Gal[p].DiscGas[i]);
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;
  }

  if(Gal[p].ColdGas < 0.0)
  {
    printf("ColdGas in update_feedback...%d, %e\n", i, Gal[p].ColdGas);
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


void combine_stellar_discs(int p, double NewStars[30], double NewStarsMetals[30])
{
	double sdisc_spin_mag, J_sdisc, J_new, J_retro, J_sum, cos_angle_sdisc_comb, cos_angle_new_comb, DiscStarSum;
	double SDiscNewSpin[3];
	double Disc1[30], Disc1Metals[30], Disc2[30], Disc2Metals[30];
	int i;
	
	for(i=0; i<30; i++){
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
		double metallicity = get_metallicity(NewStars[i], NewStarsMetals[i]);
		assert(NewStarsMetals[i] <= NewStars[i]);}
	
	// Try not to get confused, where "new" here implies the newly formed stars.  In the cooling recipe, "new" meant the combined disc, here instead denoted "comb".
	
	J_new = 0.0;
	for(i=0; i<3; i++)
		J_new += NewStars[i] * pow((pow(DiscBinEdge[i],2.0) + pow(DiscBinEdge[i+1],2.0))/2.0, 0.5);
	
	// Determine projection angles for combining discs
	if(Gal[p].StellarMass > Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass)
	{
		// Ensure the stellar disc spin magnitude is normalised
		sdisc_spin_mag = pow(pow(Gal[p].SpinStars[0], 2.0) + pow(Gal[p].SpinStars[1], 2.0) + pow(Gal[p].SpinStars[2], 2.0), 0.5);
		assert(sdisc_spin_mag==sdisc_spin_mag);
		if(sdisc_spin_mag>0.0)
		{
			for(i=0; i<3; i++)
			{
				Gal[p].SpinStars[i] /= sdisc_spin_mag; 
                assert(Gal[p].SpinStars[i]==Gal[p].SpinStars[i]);
			}
		}
	
		J_sdisc = get_disc_ang_mom(p, 1);
	
		// Obtain new spin vector of stellar disc
        for(i=0; i<3; i++){
			SDiscNewSpin[i] = Gal[p].SpinGas[i]*J_new + Gal[p].SpinStars[i]*J_sdisc;
            assert(SDiscNewSpin[i]==SDiscNewSpin[i]);}
        
		// Normalise the new spin
		sdisc_spin_mag = pow(pow(SDiscNewSpin[0], 2.0) + pow(SDiscNewSpin[1], 2.0) + pow(SDiscNewSpin[2], 2.0), 0.5);
        if(sdisc_spin_mag>0.0)
        {
            for(i=0; i<3; i++)
                SDiscNewSpin[i] /= sdisc_spin_mag;
        }
		
		cos_angle_sdisc_comb = Gal[p].SpinStars[0]*SDiscNewSpin[0] + Gal[p].SpinStars[1]*SDiscNewSpin[1] + Gal[p].SpinStars[2]*SDiscNewSpin[2];
		cos_angle_new_comb = Gal[p].SpinGas[0]*SDiscNewSpin[0] + Gal[p].SpinGas[1]*SDiscNewSpin[1] + Gal[p].SpinGas[2]*SDiscNewSpin[2];
	}
	else
	{
		cos_angle_sdisc_comb = 1.0;
		cos_angle_new_comb = 1.0;
		J_sdisc = 0.0;
        for(i=0; i<3; i++)
            SDiscNewSpin[i] = Gal[p].SpinGas[i];
	}
	
	// Combine the discs
	if(cos_angle_sdisc_comb<1.0)
    {
		project_disc(Gal[p].DiscStars, cos_angle_sdisc_comb, p, Disc1);
		project_disc(Gal[p].DiscStarsMetals, cos_angle_sdisc_comb, p, Disc1Metals);
		project_disc(NewStars, cos_angle_new_comb, p, Disc2);
		project_disc(NewStarsMetals, cos_angle_new_comb, p, Disc2Metals);
		
        Gal[p].StellarMass = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
        Gal[p].MetalsStellarMass = Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass;
		for(i=0; i<30; i++)
		{
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			Gal[p].DiscStars[i] = Disc1[i] + Disc2[i];
			Gal[p].DiscStarsMetals[i] = Disc1Metals[i] + Disc2Metals[i];
			if(Gal[p].DiscStars[i]==0.0 && Gal[p].DiscStarsMetals[i] < 1e-20) Gal[p].DiscStarsMetals[i] = 0.0;
            Gal[p].StellarMass += Gal[p].DiscStars[i];
            Gal[p].MetalsStellarMass += Gal[p].DiscStarsMetals[i];
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
	else if(Gal[p].StellarMass == Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass)
	{
		for(i=0; i<30; i++)
		{
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			Gal[p].DiscStars[i] = NewStars[i];
			Gal[p].DiscStarsMetals[i] = NewStarsMetals[i];
            Gal[p].StellarMass += Gal[p].DiscStars[i];
            Gal[p].MetalsStellarMass += Gal[p].DiscStarsMetals[i];
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
	
	// Readjust disc to deal with any retrograde stars
	if(cos_angle_sdisc_comb<0.0)
		J_retro = J_sdisc*fabs(cos_angle_sdisc_comb);
	else if(cos_angle_new_comb<0.0)
		J_retro = J_new*fabs(cos_angle_new_comb);
	else
		J_retro = 0.0;
	J_sum = J_sdisc*fabs(cos_angle_sdisc_comb) + J_new*fabs(cos_angle_new_comb);
		
	if(J_retro>0.0)
	{
		project_disc(Gal[p].DiscStars, (J_sum - 2.0*J_retro)/J_sum, p, Disc1);
		project_disc(Gal[p].DiscStarsMetals, (J_sum - 2.0*J_retro)/J_sum, p, Disc1Metals);
		for(i=0; i<30; i++)
		{
			Gal[p].DiscStars[i] = Disc1[i];
			Gal[p].DiscStarsMetals[i] = Disc1Metals[i];
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
	
	// Set the new spin direction of the stellar disc
	for(i=0; i<3; i++){
		Gal[p].SpinStars[i] = SDiscNewSpin[i];
		assert(Gal[p].SpinStars[i]==Gal[p].SpinStars[i]);}
    
    DiscStarSum = get_disc_stars(p);
    if(DiscStarSum>0.0) assert(DiscStarSum <= 1.001*(Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass) && DiscStarSum >= (Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass)/1.001);

	for(i=0; i<30; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);}

}


void project_disc(double DiscMass[30], double cos_angle, int p, double *NewDisc)
{
	double high_bound, ratio_last_bin;
	int i, j, j_old, k;
    
	cos_angle = fabs(cos_angle); // This function will not deal with retrograde motion so needs an angle less than pi/2
	
	j_old = 0;

	for(i=0; i<30; i++)
	{
		high_bound = DiscBinEdge[i+1] / cos_angle;
		j = j_old;
		
		while(DiscBinEdge[j]<=high_bound)
		{
			j++;
			if(j==30) break;
		} 
		j -= 1;
		
		NewDisc[i] = 0.0;
		for(k=j_old; k<j; k++) 
		{
			NewDisc[i] += DiscMass[k];
			DiscMass[k] = 0.0;
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
			NewDisc[i] += ratio_last_bin * DiscMass[j];
			DiscMass[j] -= ratio_last_bin * DiscMass[j];
		}
		else
		{
			NewDisc[i] = DiscMass[i];
		}
		assert(NewDisc[i]>=0.0);

		j_old = j;
	}
}






