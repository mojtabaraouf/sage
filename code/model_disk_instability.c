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
	double Q_star, Q_gas, V_rot, Q_gas_min, Q_star_min, Q_tot, W, Q_stable;
	double unstable_gas, unstable_stars, metallicity, stars, stars_sum, gas_sink;
    double r_inner, r_outer, r_av, Omega, Kappa, sigma_R, c_s;
	double NewStars[N_BINS], NewStarsMetals[N_BINS], SNgas[N_BINS], spinmag, angle, DiscGasSum, DiscStarSum, StarChannelSum;
    double old_spin[3], SNgas_copy[N_BINS], SNgas_proj[N_BINS], cos_angle;
	int i, s;
    int first, first_gas, first_star;
	
    double star_init = Gal[p].StellarMass;
    
    DiscStarSum = get_disc_stars(p);
    StarChannelSum = get_channel_stars(p);
    assert(Gal[p].StellarMass >= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)/1.01 && Gal[p].StellarMass <= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)*1.01);
    
    DiscGasSum = get_disc_gas(p);
    assert(Gal[p].ColdGas >= DiscGasSum/1.01 && Gal[p].ColdGas <= DiscGasSum*1.01);
    
    c_s = 1.1e6 / UnitVelocity_in_cm_per_s; // Speed of sound assumed for cold gas, now set to be the same as vel disp of gas at 11 km/s
    
    angle = acos(Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2])*180.0/M_PI;
    
	for(i=0; i<N_BINS; i++)
    {
		metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
        NewStars[i] = 0.0;
        NewStarsMetals[i] = 0.0;
        SNgas[i] = 0.0;
    }
	
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
        r_av = pow((r_inner*r_inner+r_outer*r_outer)/2.0, 0.5);
        
        if(Gal[p].DiscGas[i]==0.0)
            continue;
        
        if(i>0)
            Kappa = sqrt(2.0*DiscBinEdge[i]/pow(r_inner,3.0) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        else
            Kappa = sqrt(2.0*DiscBinEdge[i+1]/pow(r_outer,3.0) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        
        sigma_R = 0.5*Gal[p].Vvir*exp(-r_av/2.0/Gal[p].DiskScaleRadius);

        Q_gas = c_s * Kappa * (r_outer*r_outer - r_inner*r_inner) / G / Gal[p].DiscGas[i];
        
        if(Gal[p].DiscStars[i]>0.0 && angle<=10.0)
        {
            Q_star = Kappa * sigma_R * 0.935 * (r_outer*r_outer - r_inner*r_inner) / G / Gal[p].DiscStars[i];
            
            W = 2.0*sigma_R*c_s / (sigma_R*sigma_R + c_s*c_s);
            if(Q_gas >= Q_star)
                Q_tot = 1.0 / (W/Q_gas + 1.0/Q_star);
            else
                Q_tot = 1.0 / (1.0/Q_gas + W/Q_star);
            
            if(Q_tot>=QTotMin)
                continue;
            
            Q_stable = QTotMin + W; // This would be the quantity to make both Q_s and Q_g if they're both lower than this.
            if(Q_gas<Q_stable && Q_star<Q_stable)
                Q_gas_min = Q_stable;
            else if(Q_gas>=Q_stable && Q_star<Q_stable)
                continue; // The stars' responsibility to sort out the instability
            else if(Q_gas<Q_stable && Q_star>=Q_stable)
                Q_gas_min = 1.0 / (1.0/QTotMin - W/Q_star);
        }
        else
            Q_gas_min = QTotMin;
        
        assert(Q_gas_min >= QTotMin);
        
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
            
            if(unstable_gas>1e-10)
            {
                metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
                
                if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
                assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
                assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
                
                double before = Gal[p].DiscGas[i];
                assert(r_inner!=INFINITY && r_inner==r_inner);
                stars = deal_with_unstable_gas(unstable_gas, p, i, V_rot, metallicity, centralgal, 0, r_inner, r_outer);
                if(stars>=MIN_STARS_FOR_SN)
                    SNgas[i] = RecycleFraction * stars;
                else
                    SNgas[i] = 0.0;

                if(before-(Gal[p].DiscGas[i]-SNgas[i])<0.99*unstable_gas || before-(Gal[p].DiscGas[i]-SNgas[i])>1.01*unstable_gas)
                {
                    printf("before, after, diff, unstable_gas = %e, %e, %e, %e\n", before, Gal[p].DiscGas[i]-SNgas[i], before-(Gal[p].DiscGas[i]-SNgas[i]), unstable_gas);
                    ABORT(0);
                }
                
                if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
                assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
                assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
                
                stars_sum += stars;
                Gal[p].DiscSFR[i] += stars / dt;
                if(stars>=MIN_STARS_FOR_SN)
                {
                    NewStars[i] = (1 - RecycleFraction) * stars;
                    NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
                }
                else
                {
                    NewStars[i] = stars;
                    NewStarsMetals[i] = metallicity * stars;
                }
                assert(NewStarsMetals[i] <= NewStars[i]);
            }
            else
            {
                NewStars[i] = 0.0;
                NewStarsMetals[i] = 0.0;
            }
		}
        
        Q_gas = c_s * Kappa * (r_outer*r_outer - r_inner*r_inner) / G / (Gal[p].DiscGas[i]-SNgas[i]);
        assert(Q_gas>0);
        if(Q_gas < 0.99*Q_gas_min) printf("Q_gas final, min = %e, %e\n", Q_gas, Q_gas_min);
        assert(Q_gas >= 0.99*Q_gas_min);
        assert(Q_gas >= 0.99*QTotMin);
        
//        Q_star = Kappa * sigma_R * 0.935 * (r_outer*r_outer - r_inner*r_inner) / G / Gal[p].DiscStars[i];
//        printf("\nAt end of Q_gas check\n");
//        printf("i, Q_gas, Q_star = %d, %e, %e\n", i, Q_gas, Q_star);
//        printf("SNgas, DiscGas, DiscStars = %e, %e, %e\n", SNgas[i], Gal[p].DiscGas[i], Gal[p].DiscStars[i]);
//        printf("Kappa, r_inner, r_outer, c_s = %e, %e, %e, %e\n", Kappa, r_inner, r_outer, c_s);
	}
	
	gas_sink += Gal[p].BlackHoleMass; // Because this was set as -BHMass at the start, this is actually the increase in BH mass from the instab.
	if(gas_sink>0.0 && AGNrecipeOn > 0)
		quasar_mode_wind(p, gas_sink);
	
    for(i=0; i<N_BINS; i++) SNgas_copy[i] = SNgas[i];
    
	// Merge new-star disc with previous stellar disc
	if(stars_sum>0.0)
	{
        assert(Gal[p].StellarMass==star_init);
        
        double NewStarsSum = 0.0;
		for(i=0; i<N_BINS; i++)
        {
            assert(NewStarsMetals[i] <= NewStars[i]);
            NewStarsSum += NewStars[i];
        }
        //double StarPre = get_disc_stars(p);
        for(i=0; i<3; i++) old_spin[i] = Gal[p].SpinStars[i];
		combine_stellar_discs(p, NewStars, NewStarsMetals);
        cos_angle = Gal[p].SpinStars[0]*old_spin[0] + Gal[p].SpinStars[1]*old_spin[1] + Gal[p].SpinStars[2]*old_spin[2];
        project_disc(SNgas_copy, cos_angle, p, SNgas_proj);
        
        //double StarPost = get_disc_stars(p);
        
        //if(!(StarPost-StarPre <= 1.01*NewStarsSum && StarPost-StarPre >= NewStarsSum/1.01))
            //printf("StarPre, StarPost, StarPost-StarPre, NewStarsSum, star_init = %e, %e, %e, %e, %e\n", StarPre, StarPost, StarPost-StarPre, NewStarsSum, star_init);
        //else
            //printf("working\n");
        //assert(StarPost-StarPre <= 1.01*NewStarsSum && StarPost-StarPre >= NewStarsSum/1.01);
        
        
		Gal[p].SfrDisk[step] += stars_sum / dt; // Some of these stars may quickly be transferred to the bulge, so simply updating SfrDisk might be crude
        Gal[p].StarsInstability += NewStarsSum;
        
        if(!(NewStarsSum<=1.001*stars_sum)) printf("NewStarsSum, stars_sum = %e, %e\n", NewStarsSum, stars_sum);
        assert(NewStarsSum<=1.001*stars_sum);
        StarChannelSum = get_channel_stars(p);
        
        if(!(Gal[p].StellarMass >= StarChannelSum/1.01 && Gal[p].StellarMass <= StarChannelSum*1.01))
            //printf("stars, insitu, instab, mergeburst = %e, %e, %e, %e\n", Gal[p].StellarMass, Gal[p].StarsInSitu, Gal[p].StarsInstability, Gal[p].StarsMergeBurst);
            printf("stars actual, sum channels = %e, %e\n", Gal[p].StellarMass, Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst);
        
        assert(Gal[p].StellarMass >= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)/1.01 && Gal[p].StellarMass <= (Gal[p].StarsInSitu+Gal[p].StarsInstability+Gal[p].StarsMergeBurst)*1.01);

	}
    else
    {
        for(i=0; i<N_BINS; i++) SNgas_proj[i] = 0.0;
    }
    
	for(i=0; i<N_BINS; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);}
	
	// Deal with stellar instabilities
	for(i=N_BINS-1; i>=0; i--)
	{
        if(Gal[p].DiscStars[i]==0.0)
            continue;

        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        r_av = pow((r_inner*r_inner+r_outer*r_outer)/2.0, 0.5);
        
        if(i>0)
            Kappa = sqrt(2.0*DiscBinEdge[i]/pow(r_inner,3.0) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        else
            Kappa = sqrt(2.0*DiscBinEdge[i+1]/pow(r_outer,3.0) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        
        sigma_R = 0.5*Gal[p].Vvir*exp(-r_av/2.0/Gal[p].DiskScaleRadius);
        
        if(Gal[p].DiscGas[i]-SNgas[i]>0.0 && angle<=10.0)
        {
            Q_star = Kappa * sigma_R * 0.935 * (r_outer*r_outer - r_inner*r_inner) / G / (Gal[p].DiscStars[i]+SNgas_proj[i]);
            Q_gas = c_s * Kappa * (r_outer*r_outer - r_inner*r_inner) / G / (Gal[p].DiscGas[i]-SNgas[i]);
            
            W = 2.0*sigma_R*c_s / (sigma_R*sigma_R + c_s*c_s);
            if(Q_gas >= Q_star)
            {
                Q_tot = 1.0 / (W/Q_gas + 1.0/Q_star);
                Q_star_min = 1.0 / (1.0/QTotMin - W/Q_gas);
            }
            else
            {
                Q_tot = 1.0 / (1.0/Q_gas + W/Q_star);
                Q_star_min = W / (1.0/QTotMin - 1.0/Q_gas);
            }
            
            if(Q_tot>=QTotMin || Q_star>=Q_star_min)
                continue;
            
//            if(Q_gas < 0.9*QTotMin)
//            {
//                printf("\nWarning, Q_gas too low after dealing with gas instabilities\n");
//                printf("i, Q_gas, Q_star, Q_tot = %d, %e, %e, %e\n", i, Q_gas, Q_star, Q_tot);
//                printf("Annulus gas, gas-SN, stars+SN = %e, %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGas[i]-SNgas[i], Gal[p].DiscStars[i]+SNgas_proj[i]);
//                printf("Kappa, r_inner, r_outer, c_s = %e, %e, %e, %e\n", Kappa, r_inner, r_outer, c_s);
//                printf("Angle = %e\n", angle);
//                ABORT(0);
//            }
        }
        else
        {
            Q_star = Kappa * sigma_R * 0.935 * (r_outer*r_outer - r_inner*r_inner) / G / Gal[p].DiscStars[i];
            Q_star_min = QTotMin;
        }
        

        
		if(Q_star<Q_star_min)
		{
            double j_lose, j_gain, m_up, m_down;
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
            
			unstable_stars = (Gal[p].DiscStars[i]+SNgas[i]) * (1.0 - Q_star/Q_star_min);
            if(unstable_stars > Gal[p].DiscStars[i]) unstable_stars = Gal[p].DiscStars[i];
            
            if(unstable_stars>1e-10)
            {
                metallicity = get_metallicity(Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
                assert(Gal[p].DiscStarsMetals[i]<=Gal[p].DiscStars[i]);
                Gal[p].DiscStars[i] -= unstable_stars;
                Gal[p].DiscStarsMetals[i] = metallicity * Gal[p].DiscStars[i];
                Gal[p].TotSinkStar[i] += unstable_stars;
                
                
                if(i==N_BINS-1)
                {
                    Gal[p].DiscStars[i-1] += unstable_stars;
                    Gal[p].DiscStarsMetals[i-1] += metallicity * unstable_stars;
                    assert(Gal[p].DiscStarsMetals[i-1] <= Gal[p].DiscStars[i-1]);
                }
                else if(r_inner > 0.2*Gal[p].DiskScaleRadius || DiskInstabilityOn<2) // Conserve angular momentum while moving stars to restore stability
                {
                    j_gain = (DiscBinEdge[i+2]-DiscBinEdge[i])/2.0;
                    if(i!=0)
                    {
                        j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i-1])/2.0;
                        m_up = j_lose / (j_gain + j_lose) * unstable_stars;
                        m_down = m_up * j_gain / j_lose;
                        assert((m_up+m_down)<=1.01*unstable_stars && (m_up+m_down)>=0.99*unstable_stars);
                        
                        Gal[p].DiscStars[i-1] += m_down;
                        Gal[p].DiscStarsMetals[i-1] += metallicity * m_down;
                        assert(Gal[p].DiscStarsMetals[i-1]<=Gal[p].DiscStars[i-1]);
                    }
                    else // In principle, this shouldn't happen, but I guess there could be a really tiny DiskScaleRadius on the odd occasion
                    {
                        j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i])/2.0;
                        m_up = j_lose / (j_gain + j_lose) * unstable_stars;
                        m_down = m_up * j_gain / j_lose;
                        assert((m_up+m_down)<=1.01*unstable_stars && (m_up+m_down)>=0.99*unstable_stars);
                        
                        for(s=0; s<3; s++)
                        {
                            assert(Gal[p].SpinSecularBulge[s]==Gal[p].SpinSecularBulge[s] && Gal[p].SpinSecularBulge[s]!=INFINITY);
                            Gal[p].SpinSecularBulge[s] = Gal[p].SpinSecularBulge[s]*Gal[p].SecularBulgeMass/(Gal[p].SecularBulgeMass+m_down);
                            assert(Gal[p].SpinSecularBulge[s]==Gal[p].SpinSecularBulge[s] && Gal[p].SpinSecularBulge[s]!=INFINITY);
                        }
                        Gal[p].SecularBulgeMass += m_down;
                        Gal[p].SecularMetalsBulgeMass += metallicity * m_down;
                        
                    }
                    
                    Gal[p].DiscStars[i+1] += m_up;
                    Gal[p].DiscStarsMetals[i+1] += metallicity * m_up;
                    assert(Gal[p].DiscStarsMetals[i+1]<=Gal[p].DiscStars[i+1]);
                }
                else // Transfer unstable stars directly into the pseudobulge.  The annuli are already within it!
                {
                    j_lose = (DiscBinEdge[i+1]+DiscBinEdge[i])/2.0;
                    for(s=0; s<3; s++)
                    {
                        assert(Gal[p].SpinSecularBulge[s]==Gal[p].SpinSecularBulge[s] && Gal[p].SpinSecularBulge[s]!=INFINITY);
                        Gal[p].SpinSecularBulge[s] = (Gal[p].SpinSecularBulge[s]*Gal[p].SecularBulgeMass + Gal[p].SpinStars[s]*unstable_stars*j_lose) / (Gal[p].SecularBulgeMass + unstable_stars);
                        if(!(Gal[p].SpinSecularBulge[s]==Gal[p].SpinSecularBulge[s] && Gal[p].SpinSecularBulge[s]!=INFINITY))
                        printf("SecBulgeMass, unstable_stars = %e, %e\n", Gal[p].SecularBulgeMass, unstable_stars);
                        assert(Gal[p].SpinSecularBulge[s]==Gal[p].SpinSecularBulge[s] && Gal[p].SpinSecularBulge[s]!=INFINITY);
                    }
                    Gal[p].SecularBulgeMass += unstable_stars;
                    Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars;

                }

            }
		}
	}
	
	for(i=0; i<N_BINS; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);}
    
    update_disc_radii(p);
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
    
    double ejected_sum = 0.0;
    
    double j_lose, j_gain, m_up, m_down;
    
	// Let gas sink -- I may well want to change this formula
    gas_sink = GasSinkRate * unstable_gas;
    
    if(unstable_gas - gas_sink < MIN_STARFORMATION) // Not enough unstable gas to form stars
        gas_sink = unstable_gas;
    
//    if(Gal[p].StellarMass > 0.0)
//        gas_sink *= (1.0 - (Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass)/Gal[p].StellarMass); // / (1.0 + pow(280.0 / V_rot, 2.0));
    Gal[p].DiscGas[i] -= gas_sink;
    Gal[p].DiscGasMetals[i] -= metallicity * gas_sink;
    
    Gal[p].TotSinkGas[i] += gas_sink;
    
    if(i==N_BINS-1)
    {
        Gal[p].DiscGas[i-1] += gas_sink;
        Gal[p].DiscGasMetals[i-1] += metallicity * gas_sink;
        assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
    }
    else // Conserve angular momentum while moving gas to restore stability
    {
        j_gain = (DiscBinEdge[i+2]-DiscBinEdge[i])/2.0;
        if(i==0)
        {
            j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i])/2.0;
            m_up = j_lose / (j_gain + j_lose) * gas_sink;
            m_down = m_up * j_gain / j_lose;
            Gal[p].BlackHoleMass += m_down;
            Gal[p].ColdGas -= m_down;
            Gal[p].MetalsColdGas -= metallicity * m_down;
            assert(Gal[p].MetalsColdGas<=Gal[p].ColdGas);
        }
        else
        {
            j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i-1])/2.0;
            m_up = j_lose / (j_gain + j_lose) * gas_sink;
            m_down = m_up * j_gain / j_lose;
            Gal[p].DiscGas[i-1] += m_down;
            Gal[p].DiscGasMetals[i-1] += metallicity * m_down;
            assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
        }
        
        Gal[p].DiscGas[i+1] += m_up;
        Gal[p].DiscGasMetals[i+1] += metallicity * m_up;
        assert(Gal[p].DiscGasMetals[i+1] <= Gal[p].DiscGas[i+1]);
    }

    
    
//    if(direct_to_BH>0 || i==0)
//	{
//		Gal[p].BlackHoleMass += gas_sink;
//        assert(Gal[p].BlackHoleMass>=0);
//		Gal[p].ColdGas -= gas_sink;
//		Gal[p].MetalsColdGas -= metallicity * gas_sink;
//	}
//	else
//	{
//		Gal[p].DiscGas[i-1] += gas_sink;
//		Gal[p].DiscGasMetals[i-1] += metallicity * gas_sink;
//		assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
//	}

	// Calculate new stars formed in that annulus
	gas_sf = unstable_gas - gas_sink;
	stars = unstable_gas - gas_sink;
	if(Gal[p].DiscGas[i] > 0.0 && stars > 0.0) // Quasar feedback could blow out the unstable gas
	{
		if(SupernovaRecipeOn == 1)
		{
			area = M_PI * (r_outer*r_outer - r_inner*r_inner);
			Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
            reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area/1.3);
            
            if(!(reheated_mass==reheated_mass && reheated_mass!=INFINITY))
            {
                printf("area, reheated = %e, %e\n", area, reheated_mass);
                printf("r_inner, r_outer = %e, %e\n", r_inner, r_outer);
            }
            assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);

			// Can't use more cold gas than is available, so balance SF and feedback
		    if((stars + reheated_mass) > gas_sf && (stars + reheated_mass) > 0.0)
		    {
		    	fac = gas_sf / (stars + reheated_mass);
		    	stars *= fac;
		    	reheated_mass *= fac;
                assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);

		    }
		
			if(stars<MIN_STARS_FOR_SN)
		    {
				if(gas_sf >= MIN_STARS_FOR_SN)
				{
		    		stars = MIN_STARS_FOR_SN;
					reheated_mass = gas_sf - stars; // Previously had (1-RecycleFration)* in front of stars, which would have ensured all the unstable gas was removed in some way, but this would be inconsistent with what's done for the case that stars>MIN_STARS_FOR_SN.
                    assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);

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
                ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (V_rot * V_rot) - FeedbackReheatingEpsilon) * stars;
			    if(ejected_mass < 0.0)
			        ejected_mass = 0.0;
			}
		}
		else
		{
			reheated_mass = 0.0;// not a great treatment right now
			ejected_mass = 0.0;
		}
				
        ejected_sum += ejected_mass;
        
	    update_from_star_formation(p, stars, metallicity, i);
	
		if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
		  reheated_mass = Gal[p].DiscGas[i];
		
		metallicity_new = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
        
        if(!(reheated_mass==reheated_mass && reheated_mass!=INFINITY))
        {
            printf("stars, reheated, ejected = %e, %e, %e\n", stars, reheated_mass, ejected_mass);
            printf("gas_sf, gas_sink = %e, %e\n", gas_sf, gas_sink);
            printf("annulus gas = %e\n", Gal[p].DiscGas[i]);
        }
        
        assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
	    update_from_feedback(p, centralgal, reheated_mass, metallicity_new, i);

        if(SupernovaRecipeOn == 1 && stars>=MIN_STARS_FOR_SN)
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
    
    update_from_ejection(centralgal, ejected_sum);
	
	return stars;
		
}


void precess_gas(int p, double dt, int halonr)
{
    int i;
    double tdyn, deg_ann, deg, DiscGasSum, DiscStarSum, NewDisc[N_BINS], NewDiscMetals[N_BINS], cos_angle_gas_stars;
    double StarSpin[3];
    
    // Axis of symmetry assumed to be the bulge in a bulge-dominated system, else it's the disc
    for(i=0; i<3; i++)
    {
        if(Gal[p].ClassicalBulgeMass>0.5*Gal[p].StellarMass && fabs(Gal[p].SpinClassicalBulge[0]+Gal[p].SpinClassicalBulge[1]+Gal[p].SpinClassicalBulge[2]) > 0.0){
            StarSpin[i] = Gal[p].SpinClassicalBulge[i] / pow(pow(Gal[p].SpinClassicalBulge[0],2.0)+pow(Gal[p].SpinClassicalBulge[1],2.0)+pow(Gal[p].SpinClassicalBulge[2],2.0),0.5);
            assert(fabs(StarSpin[0]+StarSpin[1]+StarSpin[2]) > 0.0);}
        else{
            StarSpin[i] = Gal[p].SpinStars[i];
            assert(fabs(StarSpin[0]+StarSpin[1]+StarSpin[2]) > 0.0);}
    }
    
    assert(fabs(StarSpin[0]+StarSpin[1]+StarSpin[2]) > 0.0);
    
    cos_angle_gas_stars = StarSpin[0]*Gal[p].SpinGas[0] + StarSpin[1]*Gal[p].SpinGas[1] + StarSpin[2]*Gal[p].SpinGas[2];
        
    DiscGasSum = get_disc_gas(p);
    assert(DiscGasSum <= 1.001*Gal[p].ColdGas || DiscGasSum >= Gal[p].ColdGas/1.001);
    
    //printf("disc stars from instability\n");
    DiscStarSum = get_disc_stars(p);
    
    if(cos_angle_gas_stars==0)
        printf("Spin of gas and stars orthogonal -- no precession\n");
    
    if(fabs(cos_angle_gas_stars)<1.0 && DiscGasSum>0.0 && DiscStarSum>0.0 & cos_angle_gas_stars!=0.0)
    {
        deg = 0.0;
        for(i=0; i<N_BINS; i++)
        {
            tdyn = pow(Gal[p].DiscRadii[i+1],2.0) / DiscBinEdge[i+1];
            if(tdyn!=tdyn) printf("tdyn = %e\n", tdyn);
            deg_ann = DegPerTdyn * dt / tdyn; // degrees this annulus wants to precess
            deg += deg_ann * Gal[p].DiscGas[i] / DiscGasSum;
        }
        
        double cos_angle_precess = cos(deg*M_PI/180.0);
        
        if(cos_angle_precess < fabs(cos_angle_gas_stars))
            cos_angle_precess = fabs(cos_angle_gas_stars); // Gas stops precessing once it aligns or counter-aligns with stars
        
        assert(cos_angle_precess > 0.0);
        
        project_disc(Gal[p].DiscGas, cos_angle_precess, p, NewDisc);
        project_disc(Gal[p].DiscGasMetals, cos_angle_precess, p, NewDiscMetals);
        
        for(i=0; i<N_BINS; i++)
        {
            Gal[p].DiscGas[i] = NewDisc[i];
            Gal[p].DiscGasMetals[i] = NewDiscMetals[i];
        }
        
        if(cos_angle_precess == fabs(cos_angle_gas_stars) && cos_angle_gas_stars >= 0.0)
            for(i=0; i<3; i++) Gal[p].SpinGas[i] = StarSpin[i];
        else if(cos_angle_precess == fabs(cos_angle_gas_stars) && cos_angle_gas_stars < 0.0)
            for(i=0; i<3; i++) Gal[p].SpinGas[i] = -StarSpin[i];
        else
        {
            double axis[3], axis_mag, NewSpin[3];
            double sin_angle_precess = sin(acos(cos_angle_precess));
            axis[0] = Gal[p].SpinGas[1]*StarSpin[2] - Gal[p].SpinGas[2]*StarSpin[1];
            axis[1] = Gal[p].SpinGas[2]*StarSpin[0] - Gal[p].SpinGas[0]*StarSpin[2];
            axis[2] = Gal[p].SpinGas[0]*StarSpin[1] - Gal[p].SpinGas[1]*StarSpin[0];
            if(cos_angle_gas_stars < 0.0)
                for(i=0; i<3; i++) axis[i] *= -1.0;
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
                    printf("SpinStars = %e, %e, %e\n", StarSpin[0], StarSpin[1], StarSpin[2]);
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