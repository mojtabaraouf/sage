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
	double Q_star, Q_gas, V_rot, Q_gas_min;
	double unstable_gas, unstable_stars, metallicity, stars, stars_sum, gas_sink;
    double r_inner, r_outer, r_av, Omega, Kappa, omega_R, c_s;
    double GG;
	double NewStars[N_BINS], NewStarsMetals[N_BINS];
	int i;
    int first, first_gas, first_star;
	
    GG = GRAVITY * UnitMass_in_g * UnitTime_in_s * UnitTime_in_s / pow(UnitLength_in_cm,3.0);
    c_s = 1.7e6 / UnitVelocity_in_cm_per_s; // Speed of sound assumed for cold gas
    
	for(i=0; i<N_BINS; i++){
		metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}
	
	if(Gal[p].Vvir>0.0)
		V_rot = Gal[p].Vvir;
	else
		V_rot = Gal[p].Vmax;
	
	// Deal with gaseous instabilities
	stars_sum = 0.0;
	gas_sink = -Gal[p].BlackHoleMass;
	
    // If first is 1 then that means the first unstable annulus hasn't been found yet
    first = 1;
    first_gas = 1;
    first_star = 1;
    
	for(i=N_BINS-1; i>=0; i--)
	{
        //r_inner = get_annulus_radius(p, i);
        //r_outer = get_annulus_radius(p, i+1);
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        
//        if(r_inner<r0)
//            Omega = V_rot / r0;
        if(i>0)
            Omega = DiscBinEdge[i] / r_inner / r_inner;
        else
            Omega = DiscBinEdge[i+1] / r_outer / r_outer;
        
        // Q_gas assumes c_s=17 km/s.  Would have UnitLength_in_cm once more on both denominator and numerator for human-reading, but redundant
        //Q_gas = 1.7e6 * Omega*UnitVelocity_in_cm_per_s * (r_outer*r_outer - r_inner*r_inner)*UnitLength_in_cm / (GRAVITY * Gal[p].DiscGas[i]*UnitMass_in_g);
        
        Q_gas = c_s * Omega * (r_outer*r_outer - r_inner*r_inner) / GG / Gal[p].DiscGas[i];
        
//        if(Gal[p].StellarMass>0.0)
//            Q_gas_min = QGasMin * (1.0 - (Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass)/Gal[p].StellarMass);
//        else
            Q_gas_min = QGasMin;
		
		if(Q_gas<Q_gas_min)
		{
            if(first==1)
                Gal[p].TotInstabEvents += 1;
            
            if(first_gas==1)
            {
                Gal[p].FirstUnstableGas += i;
                Gal[p].TotInstabEventsGas += 1;
            }
            
            first = 0;
            first_gas = 0;
            Gal[p].TotInstabAnnuliGas +=1;
            
            unstable_gas = Gal[p].DiscGas[i]*(1.0 - Q_gas/Q_gas_min);
			metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
			
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
			
			stars = deal_with_unstable_gas(unstable_gas, p, i, V_rot, metallicity, centralgal, 0, r_inner, r_outer);
			
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
			
			stars_sum += stars;
			NewStars[i] = (1 - RecycleFraction) * stars;
			NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
			assert(NewStarsMetals[i] <= NewStars[i]);
		}
		else
		{
			NewStars[i] = 0.0;
			NewStarsMetals[i] = 0.0;
		}
	}
	
	gas_sink += Gal[p].BlackHoleMass;
	if(gas_sink>0.0 && AGNrecipeOn > 0)
		quasar_mode_wind(p, gas_sink);
	
	// Merge new-star disc with previous stellar disc
	if(stars_sum>0.0)
	{
		for(i=0; i<N_BINS; i++) assert(NewStarsMetals[i] <= NewStars[i]);
		combine_stellar_discs(p, NewStars, NewStarsMetals);
		Gal[p].SfrDisk[step] += stars_sum / dt; // Some of these stars may quickly be transferred to the bulge, so simply updating SfrDisk might be crude
        Gal[p].StarsInstability += (1-RecycleFraction)*stars_sum;
        assert(Gal[p].StellarMass >= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)/1.001 && Gal[p].StellarMass <= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)*1.001);

	}
    
	for(i=0; i<N_BINS; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}
	
	// Deal with stellar instabilities
	for(i=N_BINS-1; i>=0; i--)
	{
        //r_inner = get_annulus_radius(p, i);
        //r_outer = get_annulus_radius(p, i+1);
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        r_av = pow((r_inner*r_inner+r_outer*r_outer)/2.0, 0.5);
        
//        if(r_inner<r0)
//            Omega = sqrt(2.0) * V_rot / r0; // Not the actual value of Omega, but serves the equation below to calculate the right Q.
        if(i>0)
            Kappa = sqrt(2.0*DiscBinEdge[i]/pow(r_inner,3.0) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        else
            Kappa = sqrt(2.0*DiscBinEdge[i+1]/pow(r_outer,3.0) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        
        //Q_star = 0.66 * V_rot*UnitVelocity_in_cm_per_s * (r_outer*r_outer - r_inner*r_inner)*UnitLength_in_cm * exp(-r_av/(2*Gal[p].DiskScaleRadius)) * Omega*UnitVelocity_in_cm_per_s / (GRAVITY * Gal[p].DiscStars[i]*UnitMass_in_g);
        
        omega_R = 0.5*Gal[p].Vvir*exp(-r_av/2.0/Gal[p].DiskScaleRadius);
        Q_star = Kappa * omega_R * 0.935 * (r_outer*r_outer - r_inner*r_inner) / GG / Gal[p].DiscStars[i];
        
		if(Q_star<QStarMin)
		{
            if(first==1)
                Gal[p].TotInstabEvents += 1;
            
            if(first_star==1)
            {
                Gal[p].FirstUnstableStar += i;
                Gal[p].TotInstabEventsStar += 1;
            }
            
            first = 0;
            first_star = 0;
            Gal[p].TotInstabAnnuliStar +=1;
            
			unstable_stars = Gal[p].DiscStars[i] * (1.0 - Q_star/QStarMin);
			metallicity = get_metallicity(Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i]<=Gal[p].DiscStars[i]);
			Gal[p].DiscStars[i] -= unstable_stars;
			Gal[p].DiscStarsMetals[i] = metallicity * Gal[p].DiscStars[i];
            Gal[p].TotSinkStar[i] = unstable_stars;
			if(i!=0)
			{
				Gal[p].DiscStars[i-1] += unstable_stars;
				Gal[p].DiscStarsMetals[i-1] += metallicity * unstable_stars;
				assert(Gal[p].DiscStarsMetals[i-1]<=Gal[p].DiscStars[i-1]);
			}
			else
			{
				Gal[p].SecularBulgeMass += unstable_stars;
				Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars;
			}
		}
	}
	
	for(i=0; i<N_BINS; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);}
}

double deal_with_unstable_gas(double unstable_gas, int p, int i, double V_rot, double metallicity, int centralgal, int direct_to_BH, double r_inner, double r_outer)
{
	double gas_sink, gas_sf;
	double stars, reheated_mass, ejected_mass, Sigma_0gas, fac, area;
    double metallicity_new;
	
    if(unstable_gas > Gal[p].DiscGas[i])
        unstable_gas = Gal[p].DiscGas[i];

    double GasOrig = Gal[p].DiscGas[i];
    double GasMetalsOrig = Gal[p].DiscGasMetals[i];
    
	// Let gas sink -- I may well want to change this formula
    gas_sink = GasSinkRate * unstable_gas;
    if(Gal[p].StellarMass > 0.0)
        gas_sink *= (1.0 - (Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass)/Gal[p].StellarMass); // / (1.0 + pow(280.0 / V_rot, 2.0));
    Gal[p].DiscGas[i] -= gas_sink;
    Gal[p].DiscGasMetals[i] -= metallicity * gas_sink;
    
    Gal[p].TotSinkGas[i] += gas_sink;

    if(direct_to_BH>0 || i==0)
	{
		Gal[p].BlackHoleMass += gas_sink;
		Gal[p].ColdGas -= gas_sink;
		Gal[p].MetalsColdGas -= metallicity * gas_sink;
	}
	else
	{
		Gal[p].DiscGas[i-1] += gas_sink;
		Gal[p].DiscGasMetals[i-1] += metallicity * gas_sink;
		assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
	}

	// Calculate new stars formed in that annulus
	gas_sf = unstable_gas - gas_sink;
	stars = unstable_gas - gas_sink;
	if(Gal[p].DiscGas[i] > 0.0 && stars > 0.0) // Quasar feedback could blow out the unstable gas
	{
		if(SupernovaRecipeOn == 1)
		{
			area = M_PI * (r_outer*r_outer - r_inner*r_inner);
			Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
            reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area/1.3);
						
			// Can't use more cold gas than is available, so balance SF and feedback 
		    if((stars + reheated_mass) > gas_sf && (stars + reheated_mass) > 0.0)
		    {
		    	fac = gas_sf / (stars + reheated_mass);
		    	stars *= fac;
		    	reheated_mass *= fac;
		    }
		
			if(stars<1e-8)
		    {
				if(gas_sf >= 1e-8)
				{
		    		stars = 1e-8;
					reheated_mass = gas_sf - stars; // Previously had (1-RecycleFration)* in front of stars, which would have ensured all the unstable gas was removed in some way, but this would be inconsistent with what's done for the case that stars>1e-8.
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
				ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (V_rot * V_rot) - FeedbackReheatingEpsilon) * stars;  // Maybe change V_rot to standard Vvir?
			    if(ejected_mass < 0.0)
			        ejected_mass = 0.0;
			}
		}
		else
		{
			reheated_mass = 0.0;// not a great treatment right now
			ejected_mass = 0.0;
		}
				
	    update_from_star_formation(p, stars, metallicity, i);
	
		if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
		  reheated_mass = Gal[p].DiscGas[i];
		
		metallicity_new = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	    update_from_feedback(p, centralgal, reheated_mass, ejected_mass, metallicity_new, i);
	
		// Update metals from SN II feedback
//		if(stars <= 1e-8)
//		{
//			if(Gal[centralgal].HotGas > 0.0)
//				Gal[centralgal].MetalsHotGas += Yield * stars;
//			else
//				Gal[centralgal].MetalsEjectedMass += Yield * stars;
//		}
//		else
        if(stars > 1e-8)
        {
			Gal[p].DiscGasMetals[i] += Yield * stars*(1-metallicity);
	    	Gal[p].MetalsColdGas += Yield * stars*(1-metallicity);
		}

        if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i])
        {
            printf("i, DiscGasOrig, DiscGasMetalsOrig = %d, %e, %e\n", i, GasOrig, GasMetalsOrig);
            printf("DiscGas, Metals, %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
            printf("unstable_gas, gas_sf = %e, %e\n", unstable_gas, gas_sf);
            printf("stars formed, reheated_mass = %e, %e\n", stars, reheated_mass);
        }
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	}
	
	return stars;
		
}


void precess_gas(int p, double dt, int halonr)
{
    int i;
    double tdyn, deg_ann, deg, DiscGasSum, DiscStarSum, NewDisc[N_BINS], NewDiscMetals[N_BINS];
    
    double cos_angle_gas_stars = Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2];
        
    DiscGasSum = get_disc_gas(p);
    assert(DiscGasSum <= 1.001*Gal[p].ColdGas || DiscGasSum >= Gal[p].ColdGas/1.001);
    
    DiscStarSum = get_disc_stars(p);
    
    if(cos_angle_gas_stars<1.0 && DiscGasSum>0.0 && DiscStarSum>0.0)
    {
        deg = 0.0;
        for(i=0; i<N_BINS; i++)
        {
            tdyn = pow((pow(DiscBinEdge[i],2.0)+pow(DiscBinEdge[i+1],2.0))/2.0, 0.5) / Gal[p].Vvir / Gal[p].Vvir;
            if(tdyn!=tdyn) printf("tdyn = %e\n", tdyn);
            deg_ann = DegPerTdyn * dt / tdyn; // degrees this annulus wants to precess
            deg += deg_ann * Gal[p].DiscGas[i] / DiscGasSum;
        }
        
        double cos_angle_precess = cos(deg*M_PI/180.0);
        
        if(cos_angle_precess < fabs(cos_angle_gas_stars))
            cos_angle_precess = fabs(cos_angle_gas_stars); // Gas stops precessing once it aligns or counter-aligns with stars
        
        project_disc(Gal[p].DiscGas, cos_angle_precess, p, NewDisc);
        project_disc(Gal[p].DiscGasMetals, cos_angle_precess, p, NewDiscMetals);
        
        for(i=0; i<N_BINS; i++)
        {
            Gal[p].DiscGas[i] = NewDisc[i];
            Gal[p].DiscGasMetals[i] = NewDiscMetals[i];
        }
        
        if(cos_angle_precess == cos_angle_gas_stars && cos_angle_gas_stars >= 0.0)
            for(i=0; i<3; i++) Gal[p].SpinGas[i] = Gal[p].SpinStars[i];
        else if(cos_angle_precess == cos_angle_gas_stars && cos_angle_gas_stars < 0.0)
            for(i=0; i<3; i++) Gal[p].SpinGas[i] = -Gal[p].SpinStars[i];
        else
        {
            double axis[3], axis_mag, NewSpin[3];
            double sin_angle_precess = sin(acos(cos_angle_precess));
            axis[0] = Gal[p].SpinGas[1]*Gal[p].SpinStars[2] - Gal[p].SpinGas[2]*Gal[p].SpinStars[1];
            axis[1] = Gal[p].SpinGas[2]*Gal[p].SpinStars[0] - Gal[p].SpinGas[0]*Gal[p].SpinStars[2];
            axis[2] = Gal[p].SpinGas[0]*Gal[p].SpinStars[1] - Gal[p].SpinGas[1]*Gal[p].SpinStars[0];
            axis_mag = pow(pow(axis[0],2.0)+pow(axis[1],2.0)+pow(axis[2],2.0),0.5);
            for(i=0; i<3; i++) axis[i] /= axis_mag;
            double dot = axis[0]*Gal[p].SpinGas[0] + axis[1]*Gal[p].SpinGas[1] + axis[2]*Gal[p].SpinGas[2];
            NewSpin[0] = axis[0]*dot*(1.0-cos_angle_precess) + Gal[p].SpinGas[0]*cos_angle_precess + (axis[1]*Gal[p].SpinGas[2] - axis[2]*Gal[p].SpinGas[1])*sin_angle_precess;
            NewSpin[1] = axis[1]*dot*(1.0-cos_angle_precess) + Gal[p].SpinGas[1]*cos_angle_precess + (axis[2]*Gal[p].SpinGas[0] - axis[0]*Gal[p].SpinGas[2])*sin_angle_precess;
            NewSpin[2] = axis[2]*dot*(1.0-cos_angle_precess) + Gal[p].SpinGas[2]*cos_angle_precess + (axis[0]*Gal[p].SpinGas[1] - axis[1]*Gal[p].SpinGas[0])*sin_angle_precess;
            for(i=0; i<3; i++)
            {
                if(NewSpin[i]!=NewSpin[i])
                {
                    printf("angle, cos_angle_precess, cos_angle_gas_stars = %e, %e, %e\n", deg, cos_angle_precess, cos_angle_gas_stars);
                    printf("SpinStars = %e, %e, %e\n", Gal[p].SpinStars[0], Gal[p].SpinStars[1], Gal[p].SpinStars[2]);
                    printf("HaloSpin = %e, %e, %e\n", Halo[halonr].Spin[0], Halo[halonr].Spin[1], Halo[halonr].Spin[2]);
                    printf("axis = %e, %e, %e\n", axis[0], axis[1], axis[2]);
                    printf("NewSpin = %e, %e, %e \n", NewSpin[0], NewSpin[1], NewSpin[2]);
                }
                assert(NewSpin[i]==NewSpin[i]);
                Gal[p].SpinGas[i] = NewSpin[i];
                
            }
        }
        
        // check instability here
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
