#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
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
        Gal[p].SpinClassicalBulge[j] = 0.0;
        Gal[p].SpinSecularBulge[j] = 0.0;
    }
    
    Gal[p].Len = Halo[halonr].Len;
    Gal[p].Vmax = Halo[halonr].Vmax;
    Gal[p].Vvir = get_virial_velocity(halonr, p);
    Gal[p].Mvir = get_virial_mass(halonr, p);
    Gal[p].Rvir = get_virial_radius(halonr, p);
    
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
    for(j=0; j<N_BINS; j++)
    {
        if(Gal[p].Vvir>0.0)
        Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1] / Gal[p].Vvir;
        else
        Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1]; // This is essentially just a place-holder for problem galaxies
        Gal[p].DiscGas[j] = 0.0;
        Gal[p].DiscStars[j] = 0.0;
        Gal[p].DiscGasMetals[j] = 0.0;
        Gal[p].DiscStarsMetals[j] = 0.0;
        Gal[p].DiscSFR[j] = 0.0;
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
    // See Mo, Mao & White (1998) eq12, and using a Bullock style lambda.
    double SpinMagnitude, SpinParameter, radius;
    
    SpinMagnitude = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] +
                         Halo[halonr].Spin[1] * Halo[halonr].Spin[1] + Halo[halonr].Spin[2] * Halo[halonr].Spin[2]);
    
    // trim the extreme tail of the spin distribution for more a realistic r_s
    //if(SpinMagnitude > 1.5) SpinMagnitude = 1.5;
    
    SpinParameter = SpinMagnitude / (1.414 * Gal[p].Vvir * Gal[p].Rvir);
    
    radius = (SpinParameter / 1.414) * Gal[p].Rvir;
    
    if(radius < 0.0 || radius!=radius || radius==INFINITY)
    {
        //printf("Reset DSR from %e to 0\n", radius);
        if(!(SpinMagnitude>0)) printf("SpinMagntiude 0\n");
        radius = 0.0;
    }
    
    return radius;
    
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



double get_virial_mass(int halonr, int p)
{
    if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir >= 0.0)
    return Halo[halonr].Mvir;   /* take spherical overdensity mass estimate */
    else if(Halo[halonr].Len>0)
    return Halo[halonr].Len * PartMass;
    else if(p!=-1)
    return Gal[p].StellarMass + Gal[p].ColdGas + Gal[p].HotGas + Gal[p].BlackHoleMass + Gal[p].ICS;
}



double get_virial_velocity(int halonr, int p)
{
    double Rvir;
    
    Rvir = get_virial_radius(halonr, p);
    
    if(Rvir > 0.0)
    return sqrt(G * get_virial_mass(halonr, p) / Rvir);
    else
    return 0.0;
}



double get_virial_radius(int halonr, int p)
{
    // return Halo[halonr].Rvir;  // Used for Bolshoi
    
    double zplus1, hubble_of_z_sq, rhocrit, fac;
    
    zplus1 = 1 + ZZ[Halo[halonr].SnapNum];
    hubble_of_z_sq =
    Hubble * Hubble *(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
                      OmegaLambda);
    
    rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * G);
    fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
    
    return cbrt(get_virial_mass(halonr, p) * fac);
}


double get_disc_gas(int p)
{
    double DiscGasSum, DiscMetalsSum;
    int l;
    
    DiscGasSum = DiscMetalsSum = 0.0;
    for(l=0; l<N_BINS; l++)
    {
        if(Gal[p].DiscGas[l] < 1e-20) // Negligible gas content is negligible
        {
            Gal[p].DiscGas[l] = 0.0;
            Gal[p].DiscGasMetals[l] = 0.0;
        }
        
        DiscGasSum += Gal[p].DiscGas[l];
        DiscMetalsSum += Gal[p].DiscGasMetals[l];
    }
    
    if((Gal[p].ColdGas < 1e-15 && Gal[p].ColdGas!=0.0) || (DiscGasSum < 1e-15 && DiscGasSum!=0.0) || (Gal[p].ColdGas<=0 && Gal[p].MetalsColdGas>0.0))
    {
        //printf("get_disc_gas initial DiscGasSum, ColdGas = %e, %e\n", DiscGasSum, Gal[p].ColdGas);
        Gal[p].ColdGas = 0.0;
        Gal[p].MetalsColdGas = 0.0;
        for(l=0; l<N_BINS; l++)
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
            for(l=0; l<N_BINS; l++)
            Gal[p].DiscGas[l] = 0.0; // Sometimes a tiny non-zero difference can creep in (probably due to projecting discs).  This just takes care of that.
            DiscGasSum = 0.0;
            Gal[p].ColdGas = 0.0;
        }
        
        if(Gal[p].ColdGas<=1e-15 || Gal[p].MetalsColdGas<=0.0)
        {
            for(l=0; l<N_BINS; l++)
            Gal[p].DiscGasMetals[l] = 0.0;
            DiscMetalsSum = 0.0;
            Gal[p].MetalsColdGas = 0.0;
        }
        
//        if((DiscGasSum<1.01*Gal[p].ColdGas && DiscGasSum>Gal[p].ColdGas/1.01) || (DiscGasSum==0.0 && Gal[p].ColdGas<1e-10) || DiscGasSum < 1e-10)
//        {
            Gal[p].ColdGas = DiscGasSum; // If difference is small, just set the numbers to be the same to prevent small errors from blowing up
            Gal[p].MetalsColdGas = DiscMetalsSum;
//        }
        
    }
    return DiscGasSum;
}

double get_disc_stars(int p)
{
    double DiscStarSum, DiscAndBulge, dumped_mass;
    int l;
    
    DiscStarSum = 0.0;
    dumped_mass = 0.0;
    for(l=0; l<N_BINS; l++)
    {
        if(Gal[p].DiscStars[l] < 1e-11) // This would be less than a single star in an aperture.  Not likely.
        {
            dumped_mass += Gal[p].DiscStars[l];
            Gal[p].DiscStars[l] = 0.0;
            Gal[p].DiscStarsMetals[l] = 0.0;
        }
        DiscStarSum += Gal[p].DiscStars[l];
    }
    
    DiscAndBulge = DiscStarSum + Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
    
    if(DiscAndBulge>1.001*Gal[p].StellarMass || DiscAndBulge<Gal[p].StellarMass/1.001)
    {
        printf("get_disc_stars report: DiscAndBulge, StellarMass, dumped_mass = %e, %e, %e\n", DiscAndBulge, Gal[p].StellarMass, dumped_mass);
//        if(DiscAndBulge<1.01*Gal[p].StellarMass && DiscAndBulge>Gal[p].StellarMass/1.01)
        Gal[p].StellarMass = DiscAndBulge; // If difference is small, just set the numbers to be the same to prevent small errors from blowing up
        if(Gal[p].StellarMass <= Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass)
        {
            for(l=0; l<N_BINS; l++)
                Gal[p].DiscStars[l] = 0.0;
            DiscStarSum = 0.0;
            Gal[p].StellarMass = Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
        }
    }
    
    return DiscStarSum;
}

double get_channel_stars(int p)
{
    double ChannelFrac;
    
    if(Gal[p].StellarMass>0.0)
    {
        ChannelFrac = (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst) / Gal[p].StellarMass;
        if(ChannelFrac>1.01 || ChannelFrac<0.99) printf("1. ChannelFrac, StellarMass = %e, %e\n", ChannelFrac, Gal[p].StellarMass);
        Gal[p].StarsInSitu /= ChannelFrac;
        Gal[p].StarsInstability /= ChannelFrac;
        Gal[p].StarsMergeBurst /= ChannelFrac;
        ChannelFrac = (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst) / Gal[p].StellarMass;
        if(ChannelFrac>1.01 || ChannelFrac<0.99) printf("2. ChannelFrac, StellarMass = %e, %e\n", ChannelFrac, Gal[p].StellarMass);
    }
    else
    {
        Gal[p].StarsInSitu = 0.0;
        Gal[p].StarsInstability = 0.0;
        Gal[p].StarsMergeBurst = 0.0;
    }
    
    return Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst;
}

double get_disc_ang_mom(int p, int type)
{
    // type=0 for gas, type=1 for stars
    double J_sum;
    int l;
    
    J_sum = 0.0;
    if(type==0)
    {
        for(l=0; l<N_BINS; l++) // Calculation of average J_sum in bin no longer so straight forward
            J_sum += Gal[p].DiscGas[l] * (DiscBinEdge[l]+DiscBinEdge[l+1])/2.0;
    }
    else if(type==1)
    {
        for(l=0; l<N_BINS; l++)
            J_sum += Gal[p].DiscStars[l] * (DiscBinEdge[l]+DiscBinEdge[l+1])/2.0;
    }
    
    return J_sum;
}

void check_ejected(int p)
{
    if(!(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass))
    {
        printf("ejected mass, metals = %e, %e\n", Gal[p].EjectedMass, Gal[p].MetalsEjectedMass);
        if(Gal[p].EjectedMass <= 1e-10 || Gal[p].MetalsEjectedMass <= 1e-10)
        {
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
        }
    }
}


void update_disc_radii(int p)
{
    // Calculate the radius corresponding to an annulus edge for a given galaxy.  Calculation is iterative given a generic rotation curve format.
    int i, j, j_max;
    double left, right, tol, r_try, j_try, dif, v_try;
    double M_D, M_int, M_DM, M_CB, M_SB, M_ICS, M_hot;
    double z, a, b, c_DM, c, r_2, X, M_DM_tot, rho_const;
    double a_CB, M_CB_inf, a_SB, M_SB_inf, a_ICS, M_ICS_inf;
    
    // Try to stably set discs up first -- this might no longer be necessary with NFW treatment
    //    M_D = Gal[p].StellarMass + Gal[p].ColdGas - Gal[p].SecularBulgeMass - Gal[p].ClassicalBulgeMass;
    //    if(M_D == 0.0)
    //    {
    //        if(Gal[p].Vvir>0.0)
    //            Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1] / Gal[p].Vvir;
    //        else
    //            Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1]; // This is essentially just a place-holder for problem galaxies
    //        return;
    //    }
    
    tol = 1e-3;
    j_max = 100;
    
    // Determine the distribution of dark matter in the halo =====
    M_DM_tot = Gal[p].Mvir - Gal[p].HotGas - Gal[p].ColdGas - Gal[p].StellarMass - Gal[p].ICS - Gal[p].BlackHoleMass; // Should this include ejected mass?!!!
    
    if(M_DM_tot < 0.0) M_DM_tot = 0.0;
    
    X = log10(Gal[p].StellarMass/Gal[p].Mvir);
    
    z = ZZ[Gal[p].SnapNum];
    if(z>5.0) z=5.0;
    a = 0.520 + (0.905-0.520)*exp(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
    b = -1.01 + 0.026*z; // Dutton & Maccio 2014
    c_DM = pow(10.0, a+b*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
    c = c_DM * (1.0 + 3e-5*exp(3.4*(X+4.5))); // Di Cintio et al 2014b
    r_2 = Gal[p].Rvir / c; // Di Cintio et al 2014b
    rho_const = M_DM_tot / (log((Gal[p].Rvir+r_2)/r_2) - Gal[p].Rvir/(Gal[p].Rvir+r_2));
    // ===========================================================
    
    // Determine distribution for bulge and ICS ==================
    a_SB = 0.2 * Gal[p].DiskScaleRadius / (1.0 + sqrt(0.5)); // Fisher & Drory (2008)
    M_SB_inf = Gal[p].SecularBulgeMass * pow((Gal[p].Rvir+a_SB)/Gal[p].Rvir, 2.0);
    
    a_CB = pow(10.0, (log10(Gal[p].ClassicalBulgeMass*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h; // Sofue 2015
    M_CB_inf = Gal[p].ClassicalBulgeMass * pow((Gal[p].Rvir+a_CB)/Gal[p].Rvir, 2.0);
    
    if(Gal[p].ClassicalBulgeMass>0.0)
    a_ICS = 13.0 * a_CB; // Gonzalez et al (2005)
    else if(Gal[p].DiskScaleRadius>0.0)
    a_ICS = 13.0 * a_SB; // Feeding Fisher & Drory (2008) relation into Gonzalez et al (2005)
    else if(Gal[p].ICS > 0.0)
    printf("Issue with ICS size\n");
    M_ICS_inf = Gal[p].ICS * pow((Gal[p].Rvir+a_ICS)/Gal[p].Rvir, 2.0);
    // ===========================================================
    
    M_D = 0.0;
    left = 0.0;
    
    //    if(Gal[p].Mvir<=0.0)
    //    {
    //        printf("In update_disc_radii\n");
    //        printf("Mvir = %e\n", Gal[p].Mvir);
    //        printf("Cold, Hot, Ejected = %e, %e, %e\n", Gal[p].ColdGas, Gal[p].HotGas, Gal[p].EjectedMass);
    //        printf("Stars total, secular, classical = %e, %e, %e\n", Gal[p].StellarMass, Gal[p].SecularBulgeMass, Gal[p].ClassicalBulgeMass);
    //        printf("Black hole = %e\n", Gal[p].BlackHoleMass);
    //    }
    
    
    if(Gal[p].Mvir>0.0)
    {
        for(i=1; i<N_BINS+1; i++)
        {
            right = 2.0*DiscBinEdge[i] / Gal[p].Vvir;
            if(right<Gal[p].Rvir) right = Gal[p].Rvir;
            if(right<8.0*left) right = 8.0*left;
            M_D += Gal[p].DiscStars[i-1] + Gal[p].DiscGas[i-1];
            
            for(j=0; j<j_max; j++)
            {
                r_try = (left+right)/2.0;
                
                if(Gal[p].Mvir>0.0)
                M_DM = rho_const * (log((r_try+r_2)/r_2) - r_try/(r_try+r_2));
                else
                M_DM = 0.0;
                M_SB = M_SB_inf * pow(r_try/(r_try + a_SB), 2.0);
                M_CB = M_CB_inf * pow(r_try/(r_try + a_CB), 2.0);
                M_ICS = M_ICS_inf * pow(r_try/(r_try + a_ICS), 2.0);
                M_hot = Gal[p].HotGas * r_try / Gal[p].Rvir;
                M_int = M_DM + M_D + M_CB + M_SB + M_ICS + M_hot + Gal[p].BlackHoleMass;
                
                v_try = sqrt(G*M_int/r_try);
                if(v_try<Gal[p].Vmax || Gal[p].Vmax <= 0)
                    j_try = r_try*v_try;
                else
                    j_try = r_try*Gal[p].Vmax;
                dif = j_try/DiscBinEdge[i] - 1.0;
                
                if(j_try!=j_try || j_try<=0 || j_try==INFINITY)
                {
                    printf("\nj_try has illogical value\n");
                    printf("j_try, M_int, r_try = %e, %e, %e\n", j_try, M_int, r_try);
                    printf("M_DM, M_D, M_CB, M_SB, M_ICS, M_Hot, M_BH = %e, %e, %e, %e, %e, %e, %e\n", M_DM, M_D, M_CB, M_SB, M_ICS, M_hot, Gal[p].BlackHoleMass);
                    printf("M_CB_inf, M_SB_inf, M_ICS_tot, M_ICS_inf = %e, %e, %e, %e\n", M_CB_inf, M_SB_inf, Gal[p].ICS, M_ICS_inf);
                    printf("a_CB, a_SB, a_ICS = %e, %e, %e\n", a_CB, a_SB, a_ICS);
                    printf("R_vir = %e\n", Gal[p].Rvir);
                }
                assert(j_try==j_try && j_try>0.0 && j_try!=INFINITY);
                
                
                // Found correct r (within tolerance)
                if(fabs(dif) <= tol || (right-left)/right <= tol)
                break;
                
                // Reset boundaries for next guess
                if(dif>0)
                right = r_try;
                else
                left = r_try;
            }
            assert(r_try==r_try && r_try!=INFINITY);
            Gal[p].DiscRadii[i] = r_try;
            left = r_try;
        }
    }
}
