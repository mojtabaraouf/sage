#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_hyperg.h>
#include <assert.h>

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
    double M_D, M_int, M_DM, M_B, M_ICS, M_hot;
    double z, a, b, c_DM, c, r_s, r_2, alpha, beta, gamma, X, M_DM_tot, rho_s;
    double a_B, M_B_inf, M_B_tot, a_ICS, M_ICS_inf;
    double ia, ib, ic, iz, hyper_2F1;
    
    clock_t begin, end, h1, h2;
    double time_update, time_hyper;
    
    begin = clock();
    
    double GG = GRAVITY * UnitMass_in_g * UnitTime_in_s * UnitTime_in_s / pow(UnitLength_in_cm,3.0);

    // Try to stably set discs up first
    M_D = Gal[p].StellarMass + Gal[p].ColdGas - Gal[p].SecularBulgeMass - Gal[p].ClassicalBulgeMass;
    if(M_D = 0.0)
    {
        for(i=1; i<31; i++)
            Gal[p].DiscRadii[i] = DiscBinEdge[i] / Gal[p].Vvir;
        return;
    }
    
    
    
    tol = 1e-3;
    j_max = 100;
    
    // Determine the distribution of dark matter in the halo =====
    M_DM_tot = Gal[p].Mvir - Gal[p].HotGas - Gal[p].ColdGas - Gal[p].StellarMass - Gal[p].ICS - Gal[p].BlackHoleMass;
    if(M_DM_tot < 0)
    {
        printf("\nBaryon mass exceeds virial mass\n");
        printf("M_DM_tot, M_vir = %e, %e\n", M_DM_tot, Gal[p].Mvir);
        printf("Hot, Cold, Stars, ICS, BH, SatBary = %e, %e, %e, %e, %e, %e\n\n", Gal[p].HotGas, Gal[p].ColdGas, Gal[p].StellarMass, Gal[p].ICS, Gal[p].BlackHoleMass, Gal[p].TotalSatelliteBaryons);
        //ABORT(1);
    }
    
    if(M_DM_tot < 0.0) M_DM_tot = 0.0;

    
    X = log10(Gal[p].StellarMass/Gal[p].Mvir);
    if(X>-4.1 && X<-1.3) // Di Cintio et al 2014b
    {
        alpha = 2.94 - log10(pow(pow(10.0,X+2.33),-1.08) + pow(pow(10.0,X+2.33),2.29));
        beta = 4.23 + 1.34*X + 0.26*X*X;
        gamma = -0.06 + log10(pow(pow(10.0,X+2.56),-0.68) + pow(10.0,X+2.56));
    }
    else
    {
        alpha = 1.0;
        beta = 3.0;
        gamma = 1.0;
    }
    
    z = ZZ[Gal[p].SnapNum];
    if(z>5.0) z=5.0;
    a = 0.520 + (0.905-0.520)*exp(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
    b = -1.01 + 0.026*z; // Dutton & Maccio 2014
    c_DM = pow(10.0, a+b*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
    c = c_DM * (1.0 + 3e-5*exp(3.4*(X+4.5))); // Di Cintio et al 2014b
    r_2 = Gal[p].Rvir / c; // Di Cintio et al 2014b
    r_s = pow((2.0-gamma)/(beta-2.0), -1.0/alpha) * r_2; // Di Cintio et al 2014b
    
    if(r_s==0)
    {
        printf("\nz, a, b = %e, %e, %e\n", z, a, b);
        printf("M_vir, c_DM, c, r_2, r_s = %e, %e, %e, %e, %e\n", Gal[p].Mvir, c_DM, c, r_2, r_s);
    }
    
    ia = (3.0-gamma)/alpha;
    ib = (beta-gamma)/alpha;
    ic = (alpha-gamma+3.0)/alpha;
    iz = -pow(Gal[p].Rvir/r_s,alpha);
    
    if(iz!=iz || fabs(iz)==INFINITY) printf("\niz, Rvir, r_s, alpha = %e, %e, %e, %e\n", iz, Gal[p].Rvir, r_s, alpha);
    
    h1 = clock();
    if(iz>=-1 && iz<1)
    {
        //printf("trying -0 %e, %e, %e, %e\n", ia, ib, ic, iz);
        hyper_2F1 = gsl_sf_hyperg_2F1(ia, ib, ic, iz);
    }
    else if(iz<-1)
    {
        if(ic-ia<0 || (ic-ia>=0 && ic-ib>=0 && ia*(ic-ib)<ib*(ic-ia)))
        {
            //printf("trying -1 %e, %e, %e, %e\n", ia, ib, ic, iz);
            hyper_2F1 = pow(1.0-iz, -ia) * gsl_sf_hyperg_2F1(ia, ic-ib, ic, iz/(iz-1.0));
        }
        else if(ic-ib<0 || (ic-ia>=0 && ic-ib>=0 && ia*(ic-ib)>=ib*(ic-ia)))
        {
            //printf("trying -2 %e, %e, %e, %e\n", ia, ib, ic, iz);
            hyper_2F1 = pow(1.0-iz, -ib) * gsl_sf_hyperg_2F1(ic-ia, ib, ic, iz/(iz-1.0));
        }
        else
            printf("Potential issue");
    }
    else
        printf("Problem with input for hypergeometric function");
    h2 = clock();
    time_hyper = (double)(h2-h1)/CLOCKS_PER_SEC;
    
    rho_s = M_DM_tot / (4.0*M_PI * pow(Gal[p].Rvir,3.0-gamma) * pow(r_s,gamma) / (3.0-gamma) * hyper_2F1);
    // ===========================================================
    
    //r_0 = pow(10.0, 0.66+0.58*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e11*Hubble_h))) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h;
    //rho_0 = pow(10.0, -23.515 - 0.964*pow(M_D*UnitMass_in_g/(SOLAR_MASS*1e11*Hubble_h),0.31)) / UnitDensity_in_cgs / pow(Hubble_h,2.0);

    
    
    // Determine distribution for bulge and ICS ==================
    M_B_tot = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
    a_B = pow(10.0, (log10(M_B_tot*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h; // Sofue 2015
    if(a_B > Gal[p].Rvir/40.0) a_B = Gal[p].Rvir/40.0; // Arbitrary upper limit.  May want to motivate/change later.
    M_B_inf = M_B_tot * pow((Gal[p].Rvir+a_B)/Gal[p].Rvir, 2.0);
    
    if(a_B>0.0)
        a_ICS = 13.0 * a_B; // Gonzalez et al (2005)
    else if(Gal[p].DiskScaleRadius>0.0)
        a_ICS = 3.0 * Gal[p].DiskScaleRadius; // This is totally arbitrary and essentially a placeholder
    else
        printf("Issue with ICS size");
    M_ICS_inf = Gal[p].ICS * pow((Gal[p].Rvir+a_ICS)/Gal[p].Rvir, 2.0);
    // ===========================================================

    M_D = 0.0;
    left = 0.0;
    if(Gal[p].Mvir>0.0)
    {
        for(i=1; i<31; i++)
        {
            right = 2.0*DiscBinEdge[i] / Gal[p].Vvir;
            if(right<Gal[p].Rvir) right = Gal[p].Rvir;
            if(right<8.0*left) right = 8.0*left;
            M_D += Gal[p].DiscStars[i-1] + Gal[p].DiscGas[i-1];
            
            for(j=0; j<j_max; j++)
            {
                r_try = (left+right)/2.0;
                
                h1 = clock();
                iz = -pow(r_try/r_s,alpha);
                if(iz>=-1 && iz<1)
                {
                    //printf("trying 0 %e, %e, %e, %e\n", ia, ib, ic, iz);
                    hyper_2F1 = gsl_sf_hyperg_2F1(ia, ib, ic, iz);
                }
                else if(iz<-1)
                {
                    if(ic-ia<0 || (ic-ia>=0 && ic-ib>=0 && ia*(ic-ib)<ib*(ic-ia)))
                    {
                        //printf("trying 1 %e, %e, %e, %e\n", ia, ib, ic, iz);
                        hyper_2F1 = pow(1.0-iz, -ia) * gsl_sf_hyperg_2F1(ia, ic-ib, ic, iz/(iz-1.0));
                    }
                    else if(ic-ib<0 || (ic-ia>=0 && ic-ib>=0 && ia*(ic-ib)>=ib*(ic-ia)))
                    {
                        //printf("trying 2 %e, %e, %e, %e\n", ia, ib, ic, iz);
                        hyper_2F1 = pow(1.0-iz, -ib) * gsl_sf_hyperg_2F1(ic-ia, ib, ic, iz/(iz-1.0));
                    }
                    else
                        printf("Potential issue");
                }
                else
                    printf("Problem with input for hypergeometric function");
                h2 = clock();
                time_hyper += (double)(h2-h1)/CLOCKS_PER_SEC;
                
                // Determine mass contributions and hence total mass internal to r_try
                M_DM = 4*M_PI * rho_s * pow(r_try,3.0-gamma) * pow(r_s,gamma) / (3.0-gamma) * hyper_2F1;
                M_B = M_B_inf * pow(r_try/(r_try + a_B), 2.0);
                M_ICS = M_ICS_inf * pow(r_try/(r_try + a_ICS), 2.0);
                M_hot = Gal[p].HotGas * r_try / Gal[p].Rvir;
                M_int = M_DM + M_D + M_B + M_ICS + M_hot + Gal[p].BlackHoleMass;
                
                j_try = sqrt(GG*M_int*r_try);
                dif = j_try/DiscBinEdge[i] - 1.0;
                
                if(r_try<1e-6)
                {
                    printf("\nLow radius being fitted\n");
                    printf("i, j, r_try, j_try, BinEdge = %d, %d, %e, %e, %e\n", i, j, r_try, j_try, DiscBinEdge[i]);
                    printf("M_int, M_DM, M_D, M_B = %e, %e, %e, %e\n", M_int, M_DM, M_D, M_B);
                    printf("M_ICS, M_hot, M_BH = %e, %e, %e\n", M_ICS, M_hot, Gal[p].BlackHoleMass);
                    printf("Vcirc, Vvir = %e, %e\n", DiscBinEdge[i]/r_try, Gal[p].Vvir);
                    printf("Mvir, SM = %e, %e\n", Gal[p].Mvir, Gal[p].StellarMass);
                    
                    //if(fabs(dif) <= tol)
                    //ABORT(1);
                }
                
                if(j_try!=j_try || j_try<=0 || j_try==INFINITY)
                {
                    printf("\nj_try has illogical value\n");
                    printf("j_try, M_int, r_try = %e, %e, %e\n", j_try, M_int, r_try);
                    printf("M_DM, M_D, M_B, M_ICS, M_Hot, M_BH = %e, %e, %e, %e, %e, %e\n", M_DM, M_D, M_B, M_ICS, M_hot, Gal[p].BlackHoleMass);
                    printf("M_B_tot, M_B_inf, M_ICS_tot, M_ICS_inf = %e, %e, %e, %e\n", M_B_tot, M_B_inf, Gal[p].ICS, M_ICS_inf);
                    printf("a_B, a_ICS = %e, %e\n", a_B, a_ICS);
                    printf("R_vir = %e\n", Gal[p].Rvir);
                }
                assert(j_try==j_try && j_try>0.0);
                
                
                // Found correct r
                if(fabs(dif) <= tol)
                {
                    //printf("completed radius calculation in %d iterations\n", j+1);
                    break;
                }
                
                if((right-left)/right <= tol)
                {
                    printf("\nBoundaries for finding r(j) hit tolerance\n");
                    printf("i, rtry, left, right, r[i-1] = %d, %e, %e, %e, %e\n", i, r_try*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h, left*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h, right*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h, Gal[p].DiscRadii[i-1]*UnitLength_in_cm/(CM_PER_MPC/1e3)/Hubble_h);
                    printf("j_try, DiscBinEdge = %e, %e\n", j_try*UnitLength_in_cm/(CM_PER_MPC/1e3)*UnitVelocity_in_cm_per_s/1e5/Hubble_h, DiscBinEdge[i]*UnitLength_in_cm/(CM_PER_MPC/1e3)*UnitVelocity_in_cm_per_s/1e5/Hubble_h);
                    printf("M_D, M_B, M_vir, R_D, R_vir, V_vir = %e, %e, %e, %e, %e, %e\n\n", M_D/Hubble_h, M_B/Hubble_h, Gal[p].Mvir/Hubble_h, Gal[p].DiskScaleRadius/Hubble_h, Gal[p].Rvir/Hubble_h, Gal[p].Vvir);
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
                    //printf("v2disc, v2halo, v2bulge = %e, %e, %e\n", v2disc, v2halo, v2bulge);
                    printf("M_D, M_B, M_vir, R_D, R_vir, V_vir = %e, %e, %e, %e, %e, %e\n\n", M_D/Hubble_h, M_B/Hubble_h, Gal[p].Mvir/Hubble_h, Gal[p].DiskScaleRadius/Hubble_h, Gal[p].Rvir/Hubble_h, Gal[p].Vvir);
                }
                
            }
            Gal[p].DiscRadii[i] = r_try;
            left = r_try;
            //if(DiscBinEdge[i]/r_try > Gal[p].Vvir) printf("Vcirc, Vvir = %e, %e\n", DiscBinEdge[i]/r_try, Gal[p].Vvir);
        }
    }
    //printf("Exiting update_disc_radii\n");
    
    end = clock();
    time_update = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Update time v hyper --- %e, %e, %e\n", time_update, time_hyper, time_hyper/time_update);

}
