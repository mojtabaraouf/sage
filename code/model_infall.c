#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double infall_recipe(int centralgal, int ngal, double Zcurr)
{
  int i;
  double tot_stellarMass, tot_BHMass, tot_coldMass, tot_hotMass, tot_hotMetals, tot_ejected, tot_ejectedMetals;
  double tot_ICS, tot_ICSMetals;
  double infallingMass, reionization_modifier, DiscGasSum, Rsat, ExpFac;

  ExpFac = AA[Gal[centralgal].SnapNum]; // Expansion factor
    
  // take care of any potential numerical issues regarding hot and cold gas
  DiscGasSum = get_disc_gas(centralgal);
  assert(DiscGasSum <= 1.001*Gal[centralgal].ColdGas && DiscGasSum >= Gal[centralgal].ColdGas/1.001);
  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // need to add up all the baryonic mass asociated with the full halo 
  tot_stellarMass = tot_coldMass = tot_hotMass = tot_hotMetals = tot_ejected = tot_BHMass = tot_ejectedMetals = tot_ICS = tot_ICSMetals = 0.0;

  for(i = 0; i < ngal; i++)      // Loop over all galaxies in the FoF-halo 
  {
    Rsat = sqrt(sqr(Gal[i].Pos[0]-Gal[centralgal].Pos[0]) + sqr(Gal[i].Pos[1]-Gal[centralgal].Pos[1]) + sqr(Gal[i].Pos[2]-Gal[centralgal].Pos[2])) * ExpFac;
      
    // "Total baryons" is only concerned within Rvir
    if(Rsat <= Gal[centralgal].Rvir)
    {
      tot_stellarMass += Gal[i].StellarMass;
      tot_BHMass += Gal[i].BlackHoleMass;
      tot_coldMass += Gal[i].ColdGas;
      tot_hotMass += Gal[i].HotGas;
      tot_hotMetals += Gal[i].MetalsHotGas;
    }
      
    tot_ejected += Gal[i].EjectedMass;
    tot_ejectedMetals += Gal[i].MetalsEjectedMass;
    tot_ICS += Gal[i].ICS;
    tot_ICSMetals += Gal[i].MetalsICS;

    // satellite ejected gas goes to central ejected reservior
    if(i != centralgal)
      Gal[i].EjectedMass = Gal[i].MetalsEjectedMass = 0.0;

    // satellite ICS goes to central ICS
    if(i != centralgal) 
      Gal[i].ICS = Gal[i].MetalsICS = 0.0; 
  }

  // include reionization if necessary 
  if(ReionizationOn)
    reionization_modifier = do_reionization(centralgal, Zcurr);
  else
    reionization_modifier = 1.0;

  infallingMass = reionization_modifier * BaryonFrac * Gal[centralgal].Mvir - (tot_stellarMass + tot_coldMass + tot_hotMass + tot_ejected + tot_BHMass + tot_ICS);

  // the central galaxy keeps all the ejected mass
  Gal[centralgal].EjectedMass = tot_ejected;
  Gal[centralgal].MetalsEjectedMass = tot_ejectedMetals;

  // take care of any potential numerical issues regarding ejected mass
  if(Gal[centralgal].MetalsEjectedMass > Gal[centralgal].EjectedMass)
    Gal[centralgal].MetalsEjectedMass = Gal[centralgal].EjectedMass;
  if(Gal[centralgal].EjectedMass < 0.0)
    Gal[centralgal].EjectedMass = Gal[centralgal].MetalsEjectedMass = 0.0;
  if(Gal[centralgal].MetalsEjectedMass < 0.0)
    Gal[centralgal].MetalsEjectedMass = 0.0;

  // the central galaxy keeps all the ICS (mostly for numerical convenience)
  Gal[centralgal].ICS = tot_ICS;
  Gal[centralgal].MetalsICS = tot_ICSMetals;

  // take care of any potential numerical issues regarding intracluster stars
  if(Gal[centralgal].MetalsICS > Gal[centralgal].ICS)
    Gal[centralgal].MetalsICS = Gal[centralgal].ICS;
  if(Gal[centralgal].ICS < 0.0)
    Gal[centralgal].ICS = Gal[centralgal].MetalsICS = 0.0;
  if(Gal[centralgal].MetalsICS < 0.0)
    Gal[centralgal].MetalsICS = 0.0;

  return infallingMass;
}



double strip_from_satellite(int halonr, int centralgal, int gal, double max_strippedGas)
{
  double reionization_modifier, strippedGas, strippedGasMetals, metallicity;
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
    
  // Intialise
  strippedGas = 0.0;
    
  if(HotStripOn==1)
  {
    if(ReionizationOn)
      reionization_modifier = do_reionization(gal, ZZ[Halo[halonr].SnapNum]);
    else
      reionization_modifier = 1.0;
  
    strippedGas = -1.0 * (reionization_modifier * BaryonFrac * Gal[gal].Mvir - (Gal[gal].StellarMass + Gal[gal].ColdGas + Gal[gal].HotGas + Gal[gal].EjectedMass + Gal[gal].BlackHoleMass + Gal[gal].ICS) ) / STEPS;
      
      if(strippedGas<0)
      strippedGas = 0.0;
      
      assert(Gal[gal].EjectedMass>=Gal[gal].MetalsEjectedMass);
      assert(Gal[centralgal].EjectedMass>=Gal[centralgal].MetalsEjectedMass);
      
    if(Gal[gal].EjectedMass > strippedGas)
    {
          metallicity = get_metallicity(Gal[gal].EjectedMass, Gal[gal].MetalsEjectedMass);
          Gal[gal].EjectedMass -= strippedGas;
          Gal[gal].MetalsEjectedMass -= strippedGas * metallicity;
          Gal[centralgal].EjectedMass += strippedGas;
          Gal[centralgal].MetalsEjectedMass += strippedGas * metallicity;
          strippedGas = 0.0;
    }
    else
    {
          Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
          Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
          strippedGas -= Gal[gal].EjectedMass;
          Gal[gal].EjectedMass = 0.0;
          Gal[gal].MetalsEjectedMass = 0.0;
    }
      
      assert(Gal[gal].EjectedMass>=Gal[gal].MetalsEjectedMass);
      assert(Gal[centralgal].EjectedMass>=Gal[centralgal].MetalsEjectedMass);
      
    metallicity = get_metallicity(Gal[gal].HotGas, Gal[gal].MetalsHotGas);
	assert(Gal[gal].MetalsHotGas <= Gal[gal].HotGas);
    strippedGasMetals = strippedGas * metallicity;
  
    if(strippedGas > Gal[gal].HotGas) strippedGas = Gal[gal].HotGas;
    if(strippedGasMetals > Gal[gal].MetalsHotGas) strippedGasMetals = Gal[gal].MetalsHotGas;

    if(strippedGas>0)
    {
        Gal[gal].HotGas -= strippedGas;
        Gal[gal].MetalsHotGas -= strippedGasMetals;

        Gal[centralgal].HotGas += strippedGas;
        Gal[centralgal].MetalsHotGas += strippedGas * metallicity;
    }
  }
  else if(HotStripOn==2)
  {
      Gal[centralgal].HotGas += Gal[gal].HotGas;
      Gal[centralgal].MetalsHotGas += Gal[gal].MetalsHotGas;
      Gal[gal].HotGas = 0.0;
      Gal[gal].MetalsHotGas = 0.0;
      Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
      Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
      Gal[gal].EjectedMass = 0.0;
      Gal[gal].MetalsEjectedMass = 0.0;
  }
  else if(HotStripOn==3)
  {
      double r_gal2, v_gal2, rho_IGM, Pram, Pgrav, left, right, r_try, dif;
      int i, ii;
      
      r_gal2 = (sqr(Gal[gal].Pos[0]-Gal[centralgal].Pos[0]) + sqr(Gal[gal].Pos[1]-Gal[centralgal].Pos[1]) + sqr(Gal[gal].Pos[2]-Gal[centralgal].Pos[2])) * sqr(AA[Gal[centralgal].SnapNum]);
      v_gal2 = (sqr(Gal[gal].Vel[0]-Gal[centralgal].Vel[0]) + sqr(Gal[gal].Vel[1]-Gal[centralgal].Vel[1]) + sqr(Gal[gal].Vel[2]-Gal[centralgal].Vel[2]));
      rho_IGM = Gal[centralgal].HotGas/ (4 * M_PI * Gal[centralgal].Rvir * r_gal2);
      Pram = rho_IGM*v_gal2;
      
      Gal[gal].Mvir = get_virial_mass(halonr); // Do this need to be updated here?  When else is it updated?
      Gal[gal].Rvir = get_virial_radius(halonr);
      Pgrav = G * Gal[gal].Mvir * Gal[gal].HotGas / 8.0 / sqr(sqr(Gal[gal].Rvir)); // First calculate restoring force at the virial radius of the subhalo
      if(Pram >= Pgrav)
      {
          double M_int, M_DM, M_CB, M_SB, M_hot;
          double z, a, b, c_DM, c, r_2, X, M_DM_tot, rho_const;
          double a_CB, M_CB_inf, a_SB, M_SB_inf;
          
          // Gets rid of any ejected gas immediately if ram pressure is strong enough
          Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
          Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
          Gal[gal].EjectedMass = 0.0;
          Gal[gal].MetalsEjectedMass = 0.0;
          
          M_DM_tot = Gal[gal].Mvir - Gal[gal].HotGas - Gal[gal].ColdGas - Gal[gal].StellarMass - Gal[gal].BlackHoleMass;
          if(M_DM_tot < 0.0) M_DM_tot = 0.0;
          
          X = log10(Gal[gal].StellarMass/Gal[gal].Mvir);
          z = ZZ[Gal[gal].SnapNum];
          if(z>5.0) z=5.0;
          a = 0.520 + (0.905-0.520)*exp(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
          b = -0.101 + 0.026*z; // Dutton & Maccio 2014
          c_DM = pow(10.0, a+b*log10(Gal[gal].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
          c = c_DM * (1.0 + 3e-5*exp(3.4*(X+4.5))); // Di Cintio et al 2014b
          r_2 = Gal[gal].Rvir / c; // Di Cintio et al 2014b
          rho_const = M_DM_tot / (log((Gal[gal].Rvir+r_2)/r_2) - Gal[gal].Rvir/(Gal[gal].Rvir+r_2));
          
          a_SB = 0.2 * Gal[gal].DiskScaleRadius / (1.0 + sqrt(0.5)); // Fisher & Drory (2008)
          M_SB_inf = Gal[gal].SecularBulgeMass * sqr((Gal[gal].Rvir+a_SB)/Gal[gal].Rvir);
          
          a_CB = pow(10.0, (log10(Gal[gal].ClassicalBulgeMass*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h; // Sofue 2015
          M_CB_inf = Gal[gal].ClassicalBulgeMass * sqr((Gal[gal].Rvir+a_CB)/Gal[gal].Rvir);
          
          left = 0.0;
          right = Gal[gal].Rvir;
          r_try = Gal[gal].Rvir; // initialise to suppress warning
          for(ii=0; ii<100; ii++)
          {
              r_try = (left+right)/2.0;
              M_DM = rho_const * (log((r_try+r_2)/r_2) - r_try/(r_try+r_2));
              M_SB = M_SB_inf * sqr(r_try/(r_try + a_SB));
              M_CB = M_CB_inf * sqr(r_try/(r_try + a_CB));
              M_hot = Gal[gal].HotGas * r_try / Gal[gal].Rvir;
              M_int = M_DM + M_CB + M_SB + M_hot + Gal[gal].BlackHoleMass;
              
              // Add mass from the disc
              for(i=0; i<N_BINS; i++)
              {
                  if(Gal[gal].DiscRadii[i+1] <= r_try)
                      M_int += (Gal[gal].DiscGas[i] + Gal[gal].DiscStars[i]);
                  else
                  {
                      M_int += ((Gal[gal].DiscGas[i] + Gal[gal].DiscStars[i]) * sqr((r_try - Gal[gal].DiscRadii[i])/(Gal[gal].DiscRadii[i+1]-Gal[gal].DiscRadii[i])));
                      break;
                  }
              }
              Pgrav = G * M_int * Gal[gal].HotGas / (8.0 * cube(r_try) * Gal[gal].Rvir);
              dif = fabs(Pram-Pgrav)/Pram;
              if(dif <= 1e-3)
                  break;
              else if(Pgrav>Pram)
                  left = 1.0*r_try;
              else
                  right = 1.0*r_try;
              
              if(ii==99) printf("ii is 99 in hot RPS\n");
          }
          
          // Actually strip the gas
          strippedGas = (1.0 - r_try/Gal[gal].Rvir) * Gal[gal].HotGas / STEPS;
          if(strippedGas>max_strippedGas) strippedGas = max_strippedGas;
          metallicity = get_metallicity(Gal[gal].HotGas, Gal[gal].MetalsHotGas);
          strippedGasMetals = metallicity * strippedGas;
          
          if(strippedGas > Gal[gal].HotGas) strippedGas = Gal[gal].HotGas;
          if(strippedGasMetals > Gal[gal].MetalsHotGas) strippedGasMetals = Gal[gal].MetalsHotGas;
          
          if(strippedGas>0.0)
          {
              Gal[gal].HotGas -= strippedGas;
              Gal[gal].MetalsHotGas -= strippedGasMetals;
              
              Gal[centralgal].HotGas += strippedGas;
              Gal[centralgal].MetalsHotGas += strippedGas * metallicity;
          }
          else
              printf("Stripped gas is = %e\n", strippedGas);
      }
  }
    
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
  return strippedGas;
}


void ram_pressure_stripping(int centralgal, int gal)
{
    double r_gal2, v_gal2, rho_IGM, Sigma_gas, area;
    double ExpFac = AA[Gal[centralgal].SnapNum];
    double angle = acos(Gal[gal].SpinStars[0]*Gal[gal].SpinGas[0] + Gal[gal].SpinStars[1]*Gal[gal].SpinGas[1] + Gal[gal].SpinStars[2]*Gal[gal].SpinGas[2])*180.0/M_PI;
    double Sigma_disc;
    double Pram, Pgrav, Mstrip, MstripZ;
    int i, j;
    
    r_gal2 = (sqr(Gal[gal].Pos[0]-Gal[centralgal].Pos[0]) + sqr(Gal[gal].Pos[1]-Gal[centralgal].Pos[1]) + sqr(Gal[gal].Pos[2]-Gal[centralgal].Pos[2])) * sqr(ExpFac);
    v_gal2 = (sqr(Gal[gal].Vel[0]-Gal[centralgal].Vel[0]) + sqr(Gal[gal].Vel[1]-Gal[centralgal].Vel[1]) + sqr(Gal[gal].Vel[2]-Gal[centralgal].Vel[2]));
    rho_IGM = Gal[centralgal].HotGas/ (4 * M_PI * Gal[centralgal].Rvir * r_gal2);
    
    for(i=0; i<N_BINS; i++)
    {
        area = M_PI * (sqr(Gal[gal].DiscRadii[i+1]) - sqr(Gal[gal].DiscRadii[i]));
        Sigma_gas = Gal[gal].DiscGas[i] / area;
        
        if(angle<=ThetaThresh)
            Sigma_disc = Sigma_gas + Gal[gal].DiscStars[i] / area;
        else
            Sigma_disc = Sigma_gas;
        
        Pram = rho_IGM*v_gal2;
        Pgrav = 2*M_PI*G*Sigma_disc*Sigma_gas;
        Mstrip = Gal[gal].DiscGas[i]*Pram/Pgrav / STEPS;
        MstripZ = Gal[gal].DiscGasMetals[i]*Pram/Pgrav / STEPS;
        
        if(Pram >= Pgrav && i==0 && Sigma_gas>0.0 && RamPressureOn==1) // If the central gas is stripped, assume all gas will be stripped.
        {
            if(HeatedToCentral)
            {
                Gal[centralgal].HotGas += Gal[gal].ColdGas;
                Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas;
            }
            else
            {
                Gal[gal].HotGas += Gal[gal].ColdGas;
                Gal[gal].MetalsHotGas += Gal[gal].MetalsColdGas;
            }
            Gal[gal].ColdGas = 0.0;
            Gal[gal].MetalsColdGas = 0.0;
            for(j=0; j<N_BINS; j++)
            {
                Gal[gal].DiscGas[j] = 0.0;
                Gal[gal].DiscGasMetals[j] = 0.0;
            }
            break;
        }
        else if(((Pram >= Pgrav && RamPressureOn==1) || ((Mstrip>=Gal[gal].DiscGas[i] || MstripZ>=Gal[gal].DiscGasMetals[i]) && RamPressureOn==2)) && Sigma_gas>0.0)
        {
            if(HeatedToCentral)
            {
                Gal[centralgal].HotGas += Gal[gal].DiscGas[i];
                Gal[centralgal].MetalsHotGas += Gal[gal].DiscGasMetals[i];
            }
            else
            {
                Gal[gal].HotGas += Gal[gal].DiscGas[i];
                Gal[gal].MetalsHotGas += Gal[gal].DiscGasMetals[i];
            }
            Gal[gal].ColdGas -= Gal[gal].DiscGas[i];
            Gal[gal].MetalsColdGas -= Gal[gal].DiscGasMetals[i];
            Gal[gal].DiscGas[i] = 0.0;
            Gal[gal].DiscGasMetals[i] = 0.0;
        }
        else if(Pram >= Pgrav && RamPressureOn==2 && Sigma_gas>0.0)
        {
            if(HeatedToCentral)
            {
                Gal[centralgal].HotGas += Mstrip;
                Gal[centralgal].MetalsHotGas += MstripZ;
            }
            else
            {
                Gal[gal].HotGas += Mstrip;
                Gal[gal].MetalsHotGas += MstripZ;
            }
            Gal[gal].ColdGas -= Mstrip;
            Gal[gal].MetalsColdGas -= MstripZ;
            Gal[gal].DiscGas[i] -= Mstrip;
            Gal[gal].DiscGasMetals[i] -= MstripZ;
            assert(Gal[centralgal].MetalsHotGas<=Gal[centralgal].HotGas);
        }
    }
}



double do_reionization(int gal, double Zcurr)
{
  double alpha, a, f_of_a, a_on_a0, a_on_ar, Mfiltering, Mjeans, Mchar, mass_to_use, modifier;
  double Tvir, Vchar, omegaZ, xZ, deltacritZ, HubbleZ;

  // we employ the reionization recipie described in Gnedin (2000), however use the fitting 
  // formulas given by Kravtsov et al (2004) Appendix B 

  // here are two parameters that Kravtsov et al keep fixed, alpha gives the best fit to the Gnedin data 
  alpha = 6.0;
  Tvir = 1e4;

  // calculate the filtering mass 
  a = 1.0 / (1.0 + Zcurr);
  a_on_a0 = a / a0;
  a_on_ar = a / ar;

  if(a <= a0)
    f_of_a = 3.0 * a / ((2.0 + alpha) * (5.0 + 2.0 * alpha)) * pow(a_on_a0, alpha);
  else if((a > a0) && (a < ar))
    f_of_a =
    (3.0 / a) * a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
    a * a / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5));
  else
    f_of_a =
    (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
    (ar * ar / 10.0) * (5.0 - 4.0 * pow(a_on_ar, -0.5)) - (a0 * a0 / 10.0) * (5.0 -
    4.0 *
    pow(a_on_a0,
    -0.5)) +
    a * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * pow(a_on_ar, -0.5)));

  // this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21 
  Mjeans = 25.0 * pow(Omega, -0.5) * 2.21;
  Mfiltering = Mjeans * pow(f_of_a, 1.5);

  // calculate the characteristic mass coresponding to a halo temperature of 10^4K 
  Vchar = sqrt(Tvir / 36.0);
  omegaZ = Omega * (cube(1.0 + Zcurr) / (Omega * cube(1.0 + Zcurr) + OmegaLambda));
  xZ = omegaZ - 1.0;
  deltacritZ = 18.0 * M_PI * M_PI + 82.0 * xZ - 39.0 * xZ * xZ;
  HubbleZ = Hubble * sqrt(Omega * cube(1.0 + Zcurr) + OmegaLambda);

  Mchar = Vchar * Vchar * Vchar / (G * HubbleZ * sqrt(0.5 * deltacritZ));

  // we use the maximum of Mfiltering and Mchar 
  mass_to_use = dmax(Mfiltering, Mchar);
  modifier = 1.0 / cube(1.0 + 0.26 * (mass_to_use / Gal[gal].Mvir));

  return modifier;

}



void add_infall_to_hot(int centralgal, double infallingGas)
{
  double metallicity;

  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // if the halo has lost mass, subtract baryons from the ejected mass first, then the hot gas
  if(infallingGas < 0.0 && Gal[centralgal].EjectedMass > 0.0)
  {
    check_ejected(centralgal);
    metallicity = get_metallicity(Gal[centralgal].EjectedMass, Gal[centralgal].MetalsEjectedMass);
	assert(Gal[centralgal].MetalsEjectedMass <= Gal[centralgal].EjectedMass);
    Gal[centralgal].MetalsEjectedMass += infallingGas*metallicity;
    if(Gal[centralgal].MetalsEjectedMass < 0.0) Gal[centralgal].MetalsEjectedMass = 0.0;

    Gal[centralgal].EjectedMass += infallingGas;
    if(Gal[centralgal].EjectedMass < 0.0)
    {
      infallingGas = Gal[centralgal].EjectedMass;
      Gal[centralgal].EjectedMass = Gal[centralgal].MetalsEjectedMass = 0.0;
    }
    else
      infallingGas = 0.0;
  }

  // if the halo has lost mass, subtract hot metals mass next, then the hot gas
  if(infallingGas < 0.0 && Gal[centralgal].MetalsHotGas > 0.0)
  {
    metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    Gal[centralgal].MetalsHotGas += infallingGas*metallicity;
	Gal[centralgal].HotGas += infallingGas;
	if(Gal[centralgal].HotGas < 0.0) Gal[centralgal].HotGas = Gal[centralgal].MetalsHotGas = 0.0;
    if(Gal[centralgal].MetalsHotGas < 0.0) Gal[centralgal].MetalsHotGas = 0.0;
	metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
  }

  // limit the infalling gas so that the hot halo alone doesn't exceed the baryon fraction
  if(infallingGas > 0.0 && (Gal[centralgal].HotGas + infallingGas) / Gal[centralgal].Mvir > BaryonFrac)
    infallingGas = BaryonFrac * Gal[centralgal].Mvir - Gal[centralgal].HotGas;

  metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // add (subtract) the ambient (enriched) infalling gas to the central galaxy hot component 
  if(infallingGas > 0.0)
    Gal[centralgal].HotGas += infallingGas;

  metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}

