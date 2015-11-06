#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>

#include "core_allvars.h"
#include "core_proto.h"



void init_galaxy(int p, int halonr)
{
  int j, step;

  if(halonr != Halo[halonr].FirstHaloInFOFgroup)
  {
    printf("Hah?\n");
    ABORT(1);
  }

  Gal[p].Type = 0;

  Gal[p].GalaxyNr = GalaxyCounter;
  GalaxyCounter++;
  
  Gal[p].HaloNr = halonr;
  Gal[p].MostBoundID = Halo[halonr].MostBoundID;
  Gal[p].SnapNum = Halo[halonr].SnapNum - 1;

  Gal[p].mergeType = 0;
  Gal[p].mergeIntoID = -1;
  Gal[p].mergeIntoSnapNum = -1;
  Gal[p].dT = -1.0;

  for(j = 0; j < 3; j++)
  {
    Gal[p].Pos[j] = Halo[halonr].Pos[j];
    Gal[p].Vel[j] = Halo[halonr].Vel[j];
	Gal[p].SpinStars[j] = Halo[halonr].Spin[j] / pow(pow(Halo[halonr].Spin[0], 2.0) + pow(Halo[halonr].Spin[1], 2.0) + pow(Halo[halonr].Spin[2], 2.0), 0.5);
	Gal[p].SpinGas[j] = Halo[halonr].Spin[j] / pow(pow(Halo[halonr].Spin[0], 2.0) + pow(Halo[halonr].Spin[1], 2.0) + pow(Halo[halonr].Spin[2], 2.0), 0.5);
  }

  Gal[p].Len = Halo[halonr].Len;
  Gal[p].Vmax = Halo[halonr].Vmax;
  Gal[p].Vvir = get_virial_velocity(halonr);
  Gal[p].Mvir = get_virial_mass(halonr);
  Gal[p].Rvir = get_virial_radius(halonr);

  Gal[p].deltaMvir = 0.0;

  Gal[p].ColdGas = 0.0;
  Gal[p].StellarMass = 0.0;
  Gal[p].ClassicalBulgeMass = 0.0;
  Gal[p].SecularBulgeMass = 0.0;
  Gal[p].HotGas = 0.0;
  Gal[p].EjectedMass = 0.0;
  Gal[p].BlackHoleMass = 0.0;
  Gal[p].ICS = 0.0;
    
  Gal[p].StarsInSitu = 0.0;
  Gal[p].StarsInstability = 0.0;
  Gal[p].StarsMergeBurst = 0.0;

  Gal[p].MetalsColdGas = 0.0;
  Gal[p].MetalsStellarMass = 0.0;
  Gal[p].ClassicalMetalsBulgeMass = 0.0;
  Gal[p].SecularMetalsBulgeMass = 0.0;
  Gal[p].MetalsHotGas = 0.0;
  Gal[p].MetalsEjectedMass = 0.0;
  Gal[p].MetalsICS = 0.0;
    
  Gal[p].TotInstabEvents = 0;
  Gal[p].TotInstabEventsGas = 0;
  Gal[p].TotInstabEventsStar = 0;
  Gal[p].TotInstabAnnuliGas = 0;
  Gal[p].TotInstabAnnuliStar = 0;
  Gal[p].FirstUnstableGas = 0;
  Gal[p].FirstUnstableStar = 0;
  
  Gal[p].DiscRadii[0] = 0.0;
  for(j=0; j<30; j++)
  {
    Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1] / Gal[p].Vvir;
	Gal[p].DiscGas[j] = 0.0;
	Gal[p].DiscStars[j] = 0.0;
	Gal[p].DiscGasMetals[j] = 0.0;
	Gal[p].DiscStarsMetals[j] = 0.0;
    Gal[p].TotSinkGas[j] = 0.0;
    Gal[p].TotSinkStar[j] = 0.0;
  }

  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrBulge[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
    Gal[p].SfrBulgeColdGas[step] = 0.0;
    Gal[p].SfrBulgeColdGasMetals[step] = 0.0;
  }

  Gal[p].DiskScaleRadius = get_disk_radius(halonr, p);
  Gal[p].ClassicalBulgeRadius = 0.0;
  Gal[p].MergTime = 999.9;
  Gal[p].Cooling = 0.0;
  Gal[p].Heating = 0.0;
  Gal[p].r_heat = 0.0;
  Gal[p].LastMajorMerger = -1.0;
  Gal[p].OutflowRate = 0.0;
  Gal[p].TotalSatelliteBaryons = 0.0;

  Gal[p].infallMvir = -1.0;  //infall properties
  Gal[p].infallVvir = -1.0;
  Gal[p].infallVmax = -1.0;
  
}



double get_disk_radius(int halonr, int p)
{
  // See Mo, Shude & White (1998) eq12, and using a Bullock style lambda.
  double SpinMagnitude, SpinParameter;
  
	SpinMagnitude = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] + 
		Halo[halonr].Spin[1] * Halo[halonr].Spin[1] + Halo[halonr].Spin[2] * Halo[halonr].Spin[2]);
  
  // trim the extreme tail of the spin distribution for more a realistic r_s
  //if(SpinMagnitude > 1.5) SpinMagnitude = 1.5;
  
  SpinParameter = SpinMagnitude / (1.414 * Gal[p].Vvir * Gal[p].Rvir);
    
  return (SpinParameter / 1.414) * Gal[p].Rvir;

}



double get_metallicity(double gas, double metals)
{
  double metallicity;

  if(metals>gas)
	printf("get_metallicity report: metals, gas/stars = %e, %e\n", metals, gas);
  //assert(gas>metals);

  if(gas > 0.0 && metals > 0.0)
  {
    metallicity = metals / gas;
    if(metallicity < 1.0)
      return metallicity;
    else
      return 1.0;
  }
  else
    return 0.0;

}



double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}



double get_virial_mass(int halonr)
{
  if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir >= 0.0)
    return Halo[halonr].Mvir;   /* take spherical overdensity mass estimate */ 
  else
    return Halo[halonr].Len * PartMass;
}



double get_virial_velocity(int halonr)
{
	double Rvir;
	
	Rvir = get_virial_radius(halonr);
	
  if(Rvir > 0.0)
		return sqrt(G * get_virial_mass(halonr) / Rvir);
	else
		return 0.0;
}



double get_virial_radius(int halonr)
{
  // return Halo[halonr].Rvir;  // Used for Bolshoi

  double zplus1, hubble_of_z_sq, rhocrit, fac;
  
  zplus1 = 1 + ZZ[Halo[halonr].SnapNum];
  hubble_of_z_sq =
    Hubble * Hubble *(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
    OmegaLambda);
  
  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
  
  return cbrt(get_virial_mass(halonr) * fac);
}


double get_disc_gas(int p)
{
	double DiscGasSum, DiscMetalsSum;
	int l;
	
	DiscGasSum = DiscMetalsSum = 0.0;
	for(l=0; l<30; l++)
	{
		DiscGasSum += Gal[p].DiscGas[l];
		DiscMetalsSum += Gal[p].DiscGasMetals[l];
	}
    
    if((Gal[p].ColdGas < 1e-15 && Gal[p].ColdGas!=0.0) || (DiscGasSum < 1e-15 && DiscGasSum!=0.0))
    {
        printf("get_disc_gas initial DiscGasSum, ColdGas = %e, %e\n", DiscGasSum, Gal[p].ColdGas);
        Gal[p].ColdGas = 0.0;
        Gal[p].MetalsColdGas = 0.0;
        for(l=0; l<30; l++)
        {
            Gal[p].DiscGas[l] = 0.0;
            Gal[p].DiscGasMetals[l] = 0.0;
        }
        DiscGasSum = 0.0;
    }

	if(DiscGasSum>1.001*Gal[p].ColdGas || DiscGasSum<Gal[p].ColdGas/1.001)
	{
		printf("get_disc_gas report ... DiscSum, ColdGas =  %e, %e\n", DiscGasSum, Gal[p].ColdGas);
		printf("get_disc_gas report ... MetalsSum, ColdMetals =  %e, %e\n", DiscMetalsSum, Gal[p].MetalsColdGas);
		
		if(Gal[p].ColdGas<=1e-15)
		{
			for(l=0; l<30; l++)
				Gal[p].DiscGas[l] = 0.0; // Sometimes a tiny non-zero difference can creep in (probably due to projecting discs).  This just takes care of that.
			DiscGasSum = 0.0;
			Gal[p].ColdGas = 0.0;
		}
        
        if(Gal[p].ColdGas<=1e-15 || Gal[p].MetalsColdGas<=0.0)
        {
            for(l=0; l<30; l++)
                Gal[p].DiscGasMetals[l] = 0.0;
            DiscMetalsSum = 0.0;
            Gal[p].MetalsColdGas = 0.0;
        }
		
		if((DiscGasSum<1.01*Gal[p].ColdGas && DiscGasSum>Gal[p].ColdGas/1.01) || (DiscGasSum==0.0 && Gal[p].ColdGas<1e-10) || DiscGasSum < 1e-10)
		{
			Gal[p].ColdGas = DiscGasSum; // If difference is small, just set the numbers to be the same to prevent small errors from blowing up
			Gal[p].MetalsColdGas = DiscMetalsSum;
		}
		
  	}
	return DiscGasSum;
}

double get_disc_stars(int p)
{
    double DiscStarSum, DiscAndBulge;
    int l;
    
    DiscStarSum = 0.0;
    for(l=0; l<30; l++)
    {
        if(Gal[p].DiscStars[l] < 1e-11)
        {
            Gal[p].DiscStars[l] = 0.0;
            Gal[p].DiscStarsMetals[l] = 0.0;
        }
        DiscStarSum += Gal[p].DiscStars[l];
    }
    
    DiscAndBulge = DiscStarSum + Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
    
    if(DiscAndBulge>1.001*Gal[p].StellarMass || DiscAndBulge<Gal[p].StellarMass/1.001)
    {
        printf("get_disc_stars report %e\t%e\n", DiscAndBulge, Gal[p].StellarMass);
        if(DiscAndBulge<1.01*Gal[p].StellarMass && DiscAndBulge>Gal[p].StellarMass/1.01)
            Gal[p].StellarMass = DiscAndBulge; // If difference is small, just set the numbers to be the same to prevent small errors from blowing up
        if(Gal[p].StellarMass <= Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass)
        {
            for(l=0; l<30; l++)
                Gal[p].DiscStars[l] = 0.0;
            DiscStarSum = 0.0;
        }
    }
    return DiscStarSum;
}

double get_disc_ang_mom(int p, int type)
{
	// type=0 for gas, type=1 for stars
	double J_sum;
	int l;
	
	J_sum = 0.0;
	if(type==0)
	{
		for(l=0; l<30; l++) // Calculation of average J_sum in bin no longer so straight forward
			J_sum += Gal[p].DiscGas[l] * pow((pow(DiscBinEdge[l],2.0) + pow(DiscBinEdge[l+1],2.0))/2.0, 0.5);
	}
	else if(type==1)
	{
		for(l=0; l<30; l++)
			J_sum += Gal[p].DiscStars[l] * pow((pow(DiscBinEdge[l],2.0) + pow(DiscBinEdge[l+1],2.0))/2.0, 0.5);
	}
	
	return J_sum;
}


double get_annulus_radius(int p, int i)
{
    double vel, radius;
    // if i=0, radius=0 --- could add that explicity
    
    //printf("Bessel %e, %e\n", gsl_sf_bessel_K0(1e-3), gsl_sf_bessel_K0(1e1));
    
    if(Gal[p].Vvir > 0.0)
        vel = Gal[p].Vvir;
    else
        vel = Gal[p].Vmax;
    
    if(vel>0.0)
    {
        if(DiscBinEdge[i] >= vel*r0)
            radius = DiscBinEdge[i]/vel;
        else
            radius = sqrt(DiscBinEdge[i] * r0 / vel);
    }
    else
    {
        printf("Annulus radius set as 0 for i=%d", i);
        radius = 0.0;
    }

    return radius;
}


void update_disc_radii(int p)
{
    // Calculate the radius corresponding to an annulus edge for a given galaxy.  Calculation is iterative given a generic rotation curve format.
    int i, j, j_max;
    double left, right, tol, r_try, j_try, dif;
    double M_D, R_D, M_B, GG;
    double rbar, rtilde;
    double r_0, rho_0, a_H;
    double v2disc, v2halo, v2bulge;
    
    //printf("Entering update_disc_radii\n");
    
    tol = 1e-3;
    j_max = 1000;
    
    M_B = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
    M_D = Gal[p].StellarMass + Gal[p].ColdGas - M_B;
    R_D = Gal[p].DiskScaleRadius;
    GG = GRAVITY * UnitMass_in_g * UnitTime_in_s * UnitTime_in_s / pow(UnitLength_in_cm,3.0);
    
    r_0 = pow(10.0, 0.66+0.58*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e11*Hubble_h))) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h;
    rho_0 = pow(10.0, -23.515 - 0.964*pow(M_D*UnitMass_in_g/(SOLAR_MASS*1e11*Hubble_h),0.31)) / UnitDensity_in_cgs / pow(Hubble_h,2.0);
    a_H = pow(10.0, (log(M_B*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h;
    if(a_H > Gal[p].Vvir/40.0) a_H = Gal[p].Vvir/40.0;

    left = 0.0;
    if(Gal[p].Mvir>0.0)
    {
        for(i=1; i<31; i++)
        {
            right = 2.0*DiscBinEdge[i] / Gal[p].Vvir;
            if(right<Gal[p].Rvir) right = Gal[p].Rvir;
            if(right<8.0*left) right = 8.0*left;
            
            for(j=0; j<j_max; j++)
            {
                r_try = (left+right)/2.0;
                
                // Disc contribution to rotation curve
                rtilde = r_try / (2.0*R_D);
                if(rtilde<100) // Need this to prevent underflow errors from Bessel functions (the velocity contribution becomes totally negligible)
                    v2disc = 0.5*GG*M_D/R_D * pow(rtilde, 2.0) * (gsl_sf_bessel_K0(rtilde)*gsl_sf_bessel_I0(rtilde) - gsl_sf_bessel_K1(rtilde)*gsl_sf_bessel_I1(rtilde));
                else
                    v2disc = 0.0;
                
                // Halo contribution to rotation curve
                rbar = r_try / r_0;
                v2halo = 6.4*GG*rho_0*r_0*r_0/rbar * (log(1+rbar) - atan(rbar) + 0.5*log(1+rbar*rbar));
                if(v2halo<0 || v2halo!=v2halo) v2halo = 0.0; // Can happen due to precision in the above calculation
                
                // Bulge contribution to rotation curve
                v2bulge = GG*M_B*r_try / pow(r_try+a_H, 2.0);
                
                j_try = r_try * sqrt(v2disc+v2halo+v2bulge);
                if(j_try!=j_try)
                {
                    printf("v2disc, v2halo, v2bulge, r_try = %e, %e, %e, %e\n", v2disc, v2halo, v2bulge, r_try);
                    ABORT(1);
                }
                
                dif = j_try/DiscBinEdge[i] - 1.0;
                
                // Found correct r
                if(fabs(dif) <= tol)
                {
                    //printf("completed radius calculation in %d iterations\n", j+1);
                    break;
                }
                
                // Reset boundaries for next guess
                if(dif>0)
                    right = r_try;
                else
                    left = r_try;
                
                if(j==j_max-1) printf("Max iterations hit for radius calculation with dif = %e\n", dif);
                
                if(j==j_max-1 || r_try!=r_try)
                {
                    printf("i, rtry, left, right, r[i-1] = %d, %e, %e, %e, %e\n", i, r_try*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h, left*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h, right*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h, Gal[p].DiscRadii[i-1]*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h);
                    printf("j_try, DiscBinEdge = %e, %e\n", j_try*UnitLength_in_cm/(CM_PER_MPC/1e3)*UnitVelocity_in_cm_per_s/1e5/Hubble_h, DiscBinEdge[i]*UnitLength_in_cm/(CM_PER_MPC/1e3)*UnitVelocity_in_cm_per_s/1e5/Hubble_h);
                    printf("v2disc, v2halo, v2bulge = %e, %e, %e\n", v2disc, v2halo, v2bulge);
                    printf("M_D, M_B, M_vir, R_D, R_vir, V_vir = %e, %e, %e, %e, %e, %e\n\n", M_D/Hubble_h, M_B/Hubble_h, Gal[p].Mvir/Hubble_h, Gal[p].DiskScaleRadius/Hubble_h, Gal[p].Rvir/Hubble_h, Gal[p].Vvir);
                }
                
            }
            Gal[p].DiscRadii[i] = r_try;
            left = r_try;
        }
    }
    //printf("Exiting update_disc_radii\n");

}
