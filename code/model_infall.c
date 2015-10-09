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
  double infallingMass, reionization_modifier, DiscGasSum;
  double newSatBaryons, tot_satBaryons;

  DiscGasSum = get_disc_gas(centralgal);
  assert(DiscGasSum <= 1.001*Gal[centralgal].ColdGas && DiscGasSum >= Gal[centralgal].ColdGas/1.001);
  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // need to add up all the baryonic mass asociated with the full halo 
  tot_stellarMass = tot_coldMass = tot_hotMass = tot_hotMetals = tot_ejected = tot_BHMass = tot_ejectedMetals = tot_ICS = tot_ICSMetals = tot_satBaryons = 0.0;

  for(i = 0; i < ngal; i++)      // Loop over all galaxies in the FoF-halo 
  {
    tot_stellarMass += Gal[i].StellarMass;
    tot_BHMass += Gal[i].BlackHoleMass;
    tot_coldMass += Gal[i].ColdGas;
    tot_hotMass += Gal[i].HotGas;
    tot_hotMetals += Gal[i].MetalsHotGas;
    tot_ejected += Gal[i].EjectedMass;
    tot_ejectedMetals += Gal[i].MetalsEjectedMass;
    tot_ICS += Gal[i].ICS;
    tot_ICSMetals += Gal[i].MetalsICS;

		// record the current baryons in satellites only
    if(i != centralgal)
			tot_satBaryons += Gal[i].StellarMass + Gal[i].BlackHoleMass + Gal[i].ColdGas + Gal[i].HotGas;

    // satellite ejected gas goes to central ejected reservior
    if(i != centralgal)
      Gal[i].EjectedMass = Gal[i].MetalsEjectedMass = 0.0;

    // satellite ICS goes to central ICS
    if(i != centralgal) 
      Gal[i].ICS = Gal[i].MetalsICS = 0.0; 
  }

	// the existing baryons that have fallen in with substructure since the last timestep
	newSatBaryons = tot_satBaryons - Gal[centralgal].TotalSatelliteBaryons;

  // include reionization if necessary 
  if(ReionizationOn)
    reionization_modifier = do_reionization(centralgal, Zcurr);
  else
    reionization_modifier = 1.0;

  infallingMass =
    // reionization_modifier * BaryonFrac * Gal[centralgal].Mvir - (tot_stellarMass + tot_coldMass + tot_hotMass + tot_ejected + tot_BHMass + tot_ICS);
    reionization_modifier * BaryonFrac * Gal[centralgal].deltaMvir - newSatBaryons;

  // the central galaxy keeps all the ejected mass
  Gal[centralgal].EjectedMass = tot_ejected;
  Gal[centralgal].MetalsEjectedMass = tot_ejectedMetals;

  if(Gal[centralgal].MetalsEjectedMass > Gal[centralgal].EjectedMass)
    Gal[centralgal].MetalsEjectedMass = Gal[centralgal].EjectedMass;
  if(Gal[centralgal].EjectedMass < 0.0)
    Gal[centralgal].EjectedMass = Gal[centralgal].MetalsEjectedMass = 0.0;
  if(Gal[centralgal].MetalsEjectedMass < 0.0)
    Gal[centralgal].MetalsEjectedMass = 0.0;

  // the central galaxy keeps all the ICS (mostly for numerical convenience)
  Gal[centralgal].ICS = tot_ICS;
  Gal[centralgal].MetalsICS = tot_ICSMetals;

  if(Gal[centralgal].MetalsICS > Gal[centralgal].ICS)
    Gal[centralgal].MetalsICS = Gal[centralgal].ICS;
  if(Gal[centralgal].ICS < 0.0)
    Gal[centralgal].ICS = Gal[centralgal].MetalsICS = 0.0;
  if(Gal[centralgal].MetalsICS < 0.0)
    Gal[centralgal].MetalsICS = 0.0;

  DiscGasSum = get_disc_gas(centralgal);
  assert(DiscGasSum <= 1.001*Gal[centralgal].ColdGas && DiscGasSum >= Gal[centralgal].ColdGas/1.001);
  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
  return infallingMass;
}



void strip_from_satellite(int halonr, int centralgal, int gal)
{
  double reionization_modifier, strippedGas, strippedGasMetals, metallicity;
  double r_gal2, v_gal2, rho_IGM, Sigma_gas, area;//, Sigma_star
  int i, j;
  
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  if(ReionizationOn)
    reionization_modifier = do_reionization(gal, ZZ[Halo[halonr].SnapNum]);
  else
    reionization_modifier = 1.0;
  
  strippedGas = -1.0 *
    // (reionization_modifier * BaryonFrac * Gal[gal].Mvir - (Gal[gal].StellarMass + Gal[gal].ColdGas + Gal[gal].HotGas + Gal[gal].EjectedMass + Gal[gal].BlackHoleMass + Gal[gal].ICS) ) / STEPS;
    ( reionization_modifier * BaryonFrac * Gal[gal].deltaMvir ) / STEPS;

  if(strippedGas > 0.0)
  {
    metallicity = get_metallicity(Gal[gal].HotGas, Gal[gal].MetalsHotGas);
	assert(Gal[gal].MetalsHotGas <= Gal[gal].HotGas);
    strippedGasMetals = strippedGas * metallicity;
  
    if(strippedGas > Gal[gal].HotGas) strippedGas = Gal[gal].HotGas;
    if(strippedGasMetals > Gal[gal].MetalsHotGas) strippedGasMetals = Gal[gal].MetalsHotGas;

    Gal[gal].HotGas -= strippedGas;
    Gal[gal].MetalsHotGas -= strippedGasMetals;

    Gal[centralgal].HotGas += strippedGas;
    Gal[centralgal].MetalsHotGas += strippedGas * metallicity;
  }
  
	// Ram pressure stripping of cold gas
    if(RamPressureOn)
    {
        r_gal2 = (pow(Gal[gal].Pos[0]-Gal[centralgal].Pos[0], 2.0) + pow(Gal[gal].Pos[1]-Gal[centralgal].Pos[1], 2.0) + pow(Gal[gal].Pos[2]-Gal[centralgal].Pos[2], 2.0)) * pow(UnitLength_in_cm, 2.0);
        v_gal2 = (pow(Gal[gal].Vel[0]-Gal[centralgal].Vel[0], 2.0) + pow(Gal[gal].Vel[1]-Gal[centralgal].Vel[1], 2.0) + pow(Gal[gal].Vel[2]-Gal[centralgal].Vel[2], 2.0)) * pow(UnitVelocity_in_cm_per_s, 2.0);
        rho_IGM = Gal[centralgal].HotGas*UnitMass_in_g / (4 * M_PI * Gal[centralgal].Rvir*UnitLength_in_cm * r_gal2);
        
        area = M_PI * (pow(get_annulus_radius(gal,1), 2.0) - pow(get_annulus_radius(gal,0), 2.0)) * pow(UnitLength_in_cm, 2.0);
        Sigma_gas = Gal[gal].DiscGas[0]*UnitMass_in_g / area;
        //Sigma_star = Gal[gal].DiscStars[0]*UnitMass_in_g / area;
        
        if(rho_IGM*v_gal2 >= 2*M_PI*GRAVITY*Sigma_gas*Sigma_gas && Gal[gal].DiscGas[0] > 0.0) // I don't think this assumption the Sigma gas will always be less in outer annuli is right, especially after mergers
        {
            //printf("LHS, RHS of RPS = %e, %e at z = %e\n", rho_IGM*v_gal2, 2*M_PI*GRAVITY*Sigma_gas*Sigma_gas, ZZ[Halo[halonr].SnapNum]);
            Gal[centralgal].HotGas += Gal[gal].ColdGas;
            Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas;
            Gal[gal].ColdGas = 0.0;
            Gal[gal].MetalsColdGas = 0.0;
            for(i=0; i<30; i++)
            {
                Gal[gal].DiscGas[i] = 0.0;
                Gal[gal].DiscGasMetals[i] = 0.0;
            }
        }
        else if(Gal[gal].ColdGas > 0.0)
        {
            for(i=1; i<30; i++)
            {
                area = M_PI * (pow(get_annulus_radius(gal,i+1), 2.0) - pow(get_annulus_radius(gal,i), 2.0)) * pow(UnitLength_in_cm, 2.0);
                Sigma_gas = Gal[gal].DiscGas[i]*UnitMass_in_g / area;
            
                if(rho_IGM*v_gal2 >= 2*M_PI*GRAVITY*Sigma_gas*Sigma_gas) // Currently no accounting for gravity of stellar disc
                {
                    //printf("LHS, RHS of RPS = %e, %e\n", rho_IGM*v_gal2, 2*M_PI*GRAVITY*Sigma_gas*Sigma_gas);
                    for(j=i; j<30; j++)
                    {
                        Gal[centralgal].HotGas += Gal[gal].DiscGas[i];
                        Gal[centralgal].MetalsHotGas += Gal[gal].DiscGasMetals[i];
                        Gal[gal].ColdGas -= Gal[gal].DiscGas[i];
                        Gal[gal].MetalsColdGas -= Gal[gal].DiscGasMetals[i];
                        Gal[gal].DiscGas[i] = 0.0;
                        Gal[gal].DiscGasMetals[i] = 0.0;
                    }
                    break;
                }
            }
        }
    }
	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
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
  omegaZ = Omega * (pow(1.0 + Zcurr, 3.0) / (Omega * pow(1.0 + Zcurr, 3.0) + OmegaLambda));
  xZ = omegaZ - 1.0;
  deltacritZ = 18.0 * M_PI * M_PI + 82.0 * xZ - 39.0 * xZ * xZ;
  HubbleZ = Hubble * sqrt(Omega * pow(1.0 + Zcurr, 3.0) + OmegaLambda);

  Mchar = Vchar * Vchar * Vchar / (G * HubbleZ * sqrt(0.5 * deltacritZ));

  // we use the maximum of Mfiltering and Mchar 
  mass_to_use = dmax(Mfiltering, Mchar);
  modifier = 1.0 / pow(1.0 + 0.26 * (mass_to_use / Gal[gal].Mvir), 3.0);

  return modifier;

}



void add_infall_to_hot(int centralgal, double infallingGas)
{
  float metallicity;

  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // if the halo has lost mass, subtract baryons from the ejected mass first, then the hot gas
  if(infallingGas < 0.0 && Gal[centralgal].EjectedMass > 0.0)
  {  
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
  if(infallingGas<0.0)
  	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
  else
  	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}

