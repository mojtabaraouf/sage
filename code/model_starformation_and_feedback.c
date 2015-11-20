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
    double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, SFE_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscPre, ColdPre;//, cos_theta_gas_stars;
    double r_inner, r_outer;
    double reff, tdyn, cold_crit, strdotfull, H2sum; // For SFprescription==3

  double NewStars[N_BINS], NewStarsMetals[N_BINS];
  int i;

    double StarsPre = Gal[p].StellarMass;
    assert(Gal[p].StellarMass >= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)/1.01 && Gal[p].StellarMass <= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)*1.01);

    
  // Checks that the deconstructed disc is being treated properly and not generating NaNs
  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  f_H2_const = 1.38e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*H2FractionExponent);
  SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR;

  // Initialise variables
  strdot = 0.0;
  stars_sum = 0.0;

  Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
  Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas; // Do I ever use these or even know what they represent?
    
  //cos_theta_gas_stars = Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2];
    
  update_HI_H2(p);
    
  if(SFprescription==3) // Prescription based on SAGE
  {
      reff = 3.0 * Gal[p].DiskScaleRadius;
      tdyn = reff / Gal[p].Vvir;
      cold_crit = 0.19 * Gal[p].Vvir * reff;
      if(Gal[p].ColdGas > cold_crit && tdyn > 0.0)
          strdotfull = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
      else
          strdotfull = 0.0;
      
      H2sum = 0.0;
      for(i=0; i<N_BINS; i++) H2sum += Gal[p].DiscH2[i];
  }

  for(i=0; i<N_BINS; i++)
  {
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
      
    //r_inner = get_annulus_radius(p, i);
    //r_outer = get_annulus_radius(p, i+1);
      r_inner = Gal[p].DiscRadii[i];
      r_outer = Gal[p].DiscRadii[i+1];
      
      area = M_PI * (r_outer*r_outer - r_inner*r_inner);
		
	if(Gal[p].Vvir>0) // These galaxies (which aren't useful for science) won't have H2 to form stars
	{
        if(SFprescription==1 && Gal[p].DiscH2[i]<0.5*Gal[p].DiscHI[i])
        {
            double bb = pow(Gal[p].DiscStars[i]*Gal[p].DiscStars[0], 0.5)/area; // quadratic b term
            double cc = -pow(0.5/f_H2_const, 1.0/0.92);
            double Sig_gas_half = 0.5*(-bb + pow(bb*bb-4.0*cc, 0.5));
            double SFE_gas = SFE_H2 * 0.75 * 1.0/(1.0/0.5 + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i])/1.3 / Sig_gas_half;
            strdot = SfrEfficiency * SFE_gas * pow(Gal[p].DiscGas[i], 2.0) / area;
        }
        else if(SFprescription==3)
            strdot = strdotfull * Gal[p].DiscH2[i] / H2sum;
        else
            strdot = SfrEfficiency * SFE_H2 * Gal[p].DiscH2[i];
    }
    else
        strdot = 0.0;

	stars = strdot * dt;
	
	if(stars < 0.0)
	  stars = 0.0;

    if(stars > Gal[p].DiscGas[i])
      stars = Gal[p].DiscGas[i];

    if(SupernovaRecipeOn == 1 && Gal[p].DiscGas[i] > 0.0 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
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
		  reheated_mass = Gal[p].DiscGas[i] - stars; // Used to have (1-RecycleFraction)* in front of stars here, but changed philosophy
	    }
	
	    ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) - FeedbackReheatingEpsilon) * stars;
	    if(ejected_mass < 0.0)
	        ejected_mass = 0.0;
	
		assert(stars+reheated_mass < 1.01*Gal[p].DiscGas[i]);
	  }

	  else
	  {
		//reheated_mass = RecycleFraction * stars;
          reheated_mass = 0.0;
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
    assert(abs(Gal[p].ColdGas-ColdPre) <= 1.01*abs(Gal[p].DiscGas[i]-DiscPre) && abs(Gal[p].ColdGas-ColdPre) >= 0.999*abs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
 
	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

    // Update from SN feedback
	metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
    update_from_feedback(p, centralgal, reheated_mass, ejected_mass, metallicity, i);

	assert(abs(Gal[p].ColdGas-ColdPre) <= 1.01*abs(Gal[p].DiscGas[i]-DiscPre) && abs(Gal[p].ColdGas-ColdPre) >= 0.999*abs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);

	// Inject new metals from SN II
	if(SupernovaRecipeOn == 1 && stars>1e-9)
	{
	  if(stars>=1e-8)
	  {
	    Gal[p].DiscGasMetals[i] += Yield * stars*(1.0 - get_metallicity(NewStars[i],NewStarsMetals[i]));
	    Gal[p].MetalsColdGas += Yield * stars*(1.0 - get_metallicity(NewStars[i],NewStarsMetals[i]));
  	  }
//	  else
//	  {
//		if(Gal[centralgal].HotGas>0.0)
//			Gal[centralgal].MetalsHotGas += Yield * stars;
//		else
//			Gal[centralgal].MetalsEjectedMass += Yield * stars;
//	  }
	}
    if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i]) printf("DiscGas, Metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
  }

    double NewStarSum = 0.0;
    for(i=0; i<N_BINS; i++) NewStarSum += NewStars[i];
    //printf("NewStarSum, stars_sum*(1-R) = %e, %e\n", NewStarSum, (1-RecycleFraction)*stars_sum);
    
  // Sum stellar discs together
  if(NewStarSum>0.0)
    combine_stellar_discs(p, NewStars, NewStarsMetals);
    
  for(i=0; i<N_BINS; i++){
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}

  // Update the star formation rate 
  Gal[p].SfrDisk[step] += stars_sum / dt;
  Gal[p].StarsInSitu += (1-RecycleFraction)*stars_sum;
    
    //printf("StarsPre, stars_sum(1-R), StarsPre+formed, StellarMass = %e, %e, %e, %e\n\n", StarsPre, (1-RecycleFraction)*stars_sum, StarsPre+(1-RecycleFraction)*stars_sum, Gal[p].StellarMass);
  if(Gal[p].StellarMass > 1e-8)
  {
      assert(Gal[p].StellarMass >= (StarsPre + (1-RecycleFraction)*stars_sum)/1.01 && Gal[p].StellarMass <= (StarsPre + (1-RecycleFraction)*stars_sum)*1.01);
      assert(Gal[p].StellarMass >= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)/1.01 && Gal[p].StellarMass <= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)*1.01);
  }


  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);

  for(i=0; i<N_BINS; i++){
	metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}

  // Check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, time, dt, step);

  for(i=0; i<N_BINS; i++){
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


void combine_stellar_discs(int p, double NewStars[N_BINS], double NewStarsMetals[N_BINS])
{
	double sdisc_spin_mag, J_sdisc, J_new, J_retro, J_sum, cos_angle_sdisc_comb, cos_angle_new_comb, DiscStarSum;
	double SDiscNewSpin[3];
	double Disc1[N_BINS], Disc1Metals[N_BINS], Disc2[N_BINS], Disc2Metals[N_BINS];
	int i;
	
	for(i=0; i<N_BINS; i++){
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
		//double metallicity = get_metallicity(NewStars[i], NewStarsMetals[i]);
		assert(NewStarsMetals[i] <= NewStars[i]);}
	
	// Try not to get confused, where "new" here implies the newly formed stars.  In the cooling recipe, "new" meant the combined disc, here instead denoted "comb".
	
	J_new = 0.0;
	for(i=0; i<3; i++)
		J_new += NewStars[i] * pow((pow(DiscBinEdge[i],2.0) + pow(DiscBinEdge[i+1],2.0))/2.0, 0.5); // The assumption that DiscBinEdge is now proportional to radius has broken down
	
	// Determine projection angles for combining discs
	if(Gal[p].StellarMass > Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass)
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
		
        sdisc_spin_mag = pow(pow(SDiscNewSpin[0], 2.0) + pow(SDiscNewSpin[1], 2.0) + pow(SDiscNewSpin[2], 2.0), 0.5);
        if(sdisc_spin_mag<0.99 || sdisc_spin_mag>1.01)
        {
            printf("SpinStars somehow became %e\n", sdisc_spin_mag);
            printf("with J_sdisc, J_new = %e, %e\n", J_sdisc, J_new);
        }
        assert(sdisc_spin_mag >= 0.99 && sdisc_spin_mag <= 1.01);
        
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
        //printf("Combining two discs\n");
		project_disc(Gal[p].DiscStars, cos_angle_sdisc_comb, p, Disc1);
		project_disc(Gal[p].DiscStarsMetals, cos_angle_sdisc_comb, p, Disc1Metals);
		project_disc(NewStars, cos_angle_new_comb, p, Disc2);
		project_disc(NewStarsMetals, cos_angle_new_comb, p, Disc2Metals);
		
        Gal[p].StellarMass = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
        Gal[p].MetalsStellarMass = Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass;
		for(i=0; i<N_BINS; i++)
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
	//else if(Gal[p].StellarMass == Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass)
    else
	{
        //printf("Should be 1.0 and is %e\n", cos_angle_sdisc_comb);
		for(i=0; i<N_BINS; i++)
		{
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			Gal[p].DiscStars[i] += NewStars[i];
			Gal[p].DiscStarsMetals[i] += NewStarsMetals[i];
            Gal[p].StellarMass += NewStars[i];
            Gal[p].MetalsStellarMass += NewStarsMetals[i];
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
    //else
    //printf("Nothing actually updated\n");
	
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
        //printf("Dealing with retrograde stars\n");
		project_disc(Gal[p].DiscStars, (J_sum - 2.0*J_retro)/J_sum, p, Disc1);
		project_disc(Gal[p].DiscStarsMetals, (J_sum - 2.0*J_retro)/J_sum, p, Disc1Metals);
		for(i=0; i<N_BINS; i++)
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
    if(DiscStarSum > 1.01*(Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass) || DiscStarSum < (Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass)/1.01)
        printf("Stellar Disc, bulge, total = %e, %e, %e\n", DiscStarSum, Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass, Gal[p].StellarMass);
    
    if(DiscStarSum>0.0) assert(DiscStarSum+Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass <= 1.01*Gal[p].StellarMass && DiscStarSum+Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass >= Gal[p].StellarMass/1.01);

	for(i=0; i<N_BINS; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);}

}


void project_disc(double DiscMass[N_BINS], double cos_angle, int p, double *NewDisc)
{
	double high_bound, ratio_last_bin;
	int i, j, j_old, k;
    
	cos_angle = fabs(cos_angle); // This function will not deal with retrograde motion so needs an angle less than pi/2
	
	j_old = 0;

	for(i=0; i<N_BINS; i++)
	{
		high_bound = DiscBinEdge[i+1] / cos_angle;
		j = j_old;
		
		while(DiscBinEdge[j]<=high_bound)
		{
			j++;
			if(j==N_BINS) break;
		} 
		j -= 1;
		
		NewDisc[i] = 0.0;
		for(k=j_old; k<j; k++) 
		{
			NewDisc[i] += DiscMass[k];
			DiscMass[k] = 0.0;
		}
		if(i!=N_BINS-1)
		{
			if(j!=N_BINS-1){
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


void update_HI_H2(int p)
{
    double area, f_H2, f_H2_HI, Pressure, f_sigma;
    //double M_B_tot, M_B_inf, a_B, M_B, M_B_0;
    int i;
    
    //f_H2_const = 1.38e-3 * H2FractionFactor * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*H2FractionExponent);
    double P_0 = 5.93e-12 / UnitMass_in_g * UnitLength_in_cm * UnitTime_in_s * UnitTime_in_s;
    
    if(Gal[p].Vvir>0.0)
    {
        for(i=0; i<N_BINS; i++)
        {
            //area = M_PI * (pow(get_annulus_radius(p,i+1), 2.0) - pow(get_annulus_radius(p,i), 2.0));
            area = M_PI * (pow(Gal[p].DiscRadii[i+1],2.0) - pow(Gal[p].DiscRadii[i],2.0));
            
            if(SFprescription!=2)
            {
                // Determine bulge stars in each annulus
//                M_B_tot = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
//                if(M_B_tot > 0.0)
//                {
//                    a_B = ((Gal[p].ClassicalBulgeMass * Gal[p].ClassicalBulgeRadius) + (Gal[p].SecularBulgeMass * 0.2*Gal[p].DiskScaleRadius)) / M_B_tot / (1.0 + sqrt(0.5));
//                    M_B_inf = M_B_tot * pow((Gal[p].Rvir+a_B)/Gal[p].Rvir, 2.0);
//                    M_B = M_B_inf * (pow(Gal[p].DiscRadii[i+1]/(Gal[p].DiscRadii[i+1] + a_B), 2.0) - pow(Gal[p].DiscRadii[i]/(Gal[p].DiscRadii[i] + a_B), 2.0));
//                }
//                else
//                    M_B = 0.0;
//                
//                if(i==0)
//                    M_B_0 = M_B;
                
                
                f_sigma =  1.1e6/UnitVelocity_in_cm_per_s / (0.5*Gal[p].Vvir*exp(-(Gal[p].DiscRadii[i]+Gal[p].DiscRadii[i+1])/4.0/Gal[p].DiskScaleRadius)); // Ratio of gas vel dispersion to stars', assuming gas is always 11 km/s
                Pressure = 0.5*M_PI*G * Gal[p].DiscGas[i] * (Gal[p].DiscGas[i] + f_sigma*Gal[p].DiscStars[i]) / pow(area,2.0);
                f_H2_HI = H2FractionFactor * pow(Pressure/P_0, H2FractionExponent);
                
                //f_H2_HI = f_H2_const * pow(pow(Gal[p].DiscGas[i]/area, 2.0) + 0.1*Gal[p].DiscGas[i]/area * pow((Gal[p].DiscStars[i]+M_B)*(Gal[p].DiscStars[0]+M_B_0), 0.5)/area, H2FractionExponent);
                
                
                
            }
            else
            {
                double s, Zp, chi, c_f, Sigma_comp0, Tau_c;
//                c_f = 1.5; // Ideally this should be a free parameter
//                Sigma_comp0 = c_f * (Gal[p].DiscGas[i]*1e10/Hubble_h)/(area*1e6/pow(Hubble_h,2.0));
                Zp = Gal[p].DiscGasMetals[i] / Gal[p].DiscGas[i] / 0.02; // Might also want solar metal fraction to be variable too
//                if(Zp>10.0) printf("Zp too high for recipe = %e\n", Zp);
//                assert(Gal[p].DiscGas[i] >= Gal[p].DiscGasMetals[i]);
//                chi = 2.3 * (1 + 3.1*pow(Zp,0.365)) / 3.0;
//                double psi = chi * (2.5+chi) / (2.5+chi*2.71828);
//                s = Sigma_comp0 * Zp / psi;
//                f_H2_HI = pow(1+pow(s/11.0,3.0)*pow((125.0+s)/(96.0+s),3.0), (1.0/3.0)) - 1.0;
                
                if(Zp>0.01 && Zp<1)
                    c_f = pow(Zp, -0.7);
                else if(Zp>=1)
                    c_f = 1.0;
                else
                    pow(0.01, -0.7);
                
                Sigma_comp0 = c_f * Gal[p].DiscGas[i]/area;
                Tau_c = 0.66 * Zp *(Sigma_comp0 * (CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm));
                chi = 3.1 * (1+ 3.1*pow(Zp,0.365)) / 4.1;
                s = log(1 + 0.6*chi + 0.01*chi*chi) / (0.6*Tau_c);
                f_H2 = 1.0 - 0.75*s/(1+0.25*s); // Not actual H2/cold, but rather H2/(H2+HI)
                f_H2_HI = 1.0 / (1.0/f_H2 - 1.0);
            }
            
            if(f_H2_HI > 0.0)
            {
                assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
                f_H2 = 0.75 * 1.0/(1.0/f_H2_HI + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i]) / 1.3; //Changes f_H2 from being H2/HI to H2/Cold Gas
                Gal[p].DiscH2[i] = f_H2 * Gal[p].DiscGas[i];
                Gal[p].DiscHI[i] = Gal[p].DiscH2[i] / f_H2_HI;
            }
            else
            {
                Gal[p].DiscH2[i] = 0.0;
                Gal[p].DiscHI[i] = 0.0;
            }
        }
    }
}



