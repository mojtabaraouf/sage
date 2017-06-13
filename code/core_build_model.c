#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <mpi.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void construct_galaxies(int halonr, int tree)
{
  static int halosdone = 0;
  int prog, fofhalo, ngal;

  HaloAux[halonr].DoneFlag = 1;
  halosdone++;

  prog = Halo[halonr].FirstProgenitor;
  while(prog >= 0)
  {
    if(HaloAux[prog].DoneFlag == 0)
      construct_galaxies(prog, tree);
    prog = Halo[prog].NextProgenitor;
  }

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 0)
  {
    HaloAux[fofhalo].HaloFlag = 1;
    while(fofhalo >= 0)
    {
      prog = Halo[fofhalo].FirstProgenitor;
      while(prog >= 0)
      {
        if(HaloAux[prog].DoneFlag == 0)
          construct_galaxies(prog, tree);
        prog = Halo[prog].NextProgenitor;
      }

      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }
  }

  // At this point, the galaxies for all progenitors of this halo have been
  // properly constructed. Also, the galaxies of the progenitors of all other 
  // halos in the same FOF-group have been constructed as well. We can hence go
  // ahead and construct all galaxies for the subhalos in this FOF halo, and
  // evolve them in time. 

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 1)
  {
    ngal = 0;
    HaloAux[fofhalo].HaloFlag = 2;

    while(fofhalo >= 0)
    {
      ngal = join_galaxies_of_progenitors(fofhalo, ngal);
      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }

    evolve_galaxies(Halo[halonr].FirstHaloInFOFgroup, ngal, tree);
  }

}



int join_galaxies_of_progenitors(int halonr, int ngalstart)
{
  int ngal, prog, mother_halo=-1, i, j, first_occupied, lenmax, lenoccmax, centralgal;
  double previousMvir, previousVvir, previousVmax, SpinMag;
  int step;

  lenmax = 0;
  lenoccmax = 0;
  first_occupied = Halo[halonr].FirstProgenitor;
  prog = Halo[halonr].FirstProgenitor;

  if(prog >=0)
    if(HaloAux[prog].NGalaxies > 0)
    lenoccmax = -1;

  // Find most massive progenitor that contains an actual galaxy
  // Maybe FirstProgenitor never was FirstHaloInFOFGroup and thus got no galaxy

  while(prog >= 0)
  {
    if(Halo[prog].Len > lenmax)
    {
      lenmax = Halo[prog].Len;
      mother_halo = prog;
    }
    if(lenoccmax != -1 && Halo[prog].Len > lenoccmax && HaloAux[prog].NGalaxies > 0)
    {
      lenoccmax = Halo[prog].Len;
      first_occupied = prog;
    }
    prog = Halo[prog].NextProgenitor;
  }

  ngal = ngalstart;
  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0)
  {
    for(i = 0; i < HaloAux[prog].NGalaxies; i++)
    {
      if(ngal >= FoF_MaxGals)
      {
        printf("Opps. We reached the maximum number FoF_MaxGals=%d of galaxies in FoF... exiting.\n", FoF_MaxGals);
        ABORT(1);
      }

      // This is the cruical line in which the properties of the progenitor galaxies 
      // are copied over (as a whole) to the (temporary) galaxies Gal[xxx] in the current snapshot 
      // After updating their properties and evolving them 
      // they are copied to the end of the list of permanent galaxies HaloGal[xxx] 

      Gal[ngal] = HaloGal[HaloAux[prog].FirstGalaxy + i];
      Gal[ngal].HaloNr = halonr;      
      Gal[ngal].dT = -1.0;

      // this deals with the central galaxies of (sub)halos 
      if(Gal[ngal].Type == 0 || Gal[ngal].Type == 1)
      {

        // remember properties from the last snapshot
        previousMvir = Gal[ngal].Mvir;
        previousVvir = Gal[ngal].Vvir;
        previousVmax = Gal[ngal].Vmax;

        if(prog == first_occupied)
        {
          // update properties of this galaxy with physical properties of halo 
          Gal[ngal].MostBoundID = Halo[halonr].MostBoundID;

          for(j = 0; j < 3; j++)
          {
            Gal[ngal].Pos[j] = Halo[halonr].Pos[j];
            Gal[ngal].Vel[j] = Halo[halonr].Vel[j];
          }
  
          if(Halo[halonr].Len > Gal[ngal].LenMax) Gal[ngal].LenMax = Halo[halonr].Len;
          Gal[ngal].Len = Halo[halonr].Len;
          Gal[ngal].Vmax = Halo[halonr].Vmax;
            assert(Gal[ngal].LenMax>=Gal[ngal].Len);

          Gal[ngal].deltaMvir = get_virial_mass(halonr, ngal) - Gal[ngal].Mvir;

          //if(get_virial_mass(halonr, ngal) > Gal[ngal].Mvir)
          //{
            Gal[ngal].Rvir = get_virial_radius(halonr, ngal);  //use the maximum Rvir in model
              //use the maximum Vvir in model
          //}

          Gal[ngal].Mvir = get_virial_mass(halonr, ngal);

          Gal[ngal].Cooling = 0.0;
          Gal[ngal].Heating = 0.0;
          Gal[ngal].OutflowRate = 0.0;

          for(step = 0; step < STEPS; step++)
          {
            Gal[ngal].SfrDisk[step] = Gal[ngal].SfrBulge[step] = 0.0;
            Gal[ngal].SfrDiskColdGas[step] = Gal[ngal].SfrDiskColdGasMetals[step] = 0.0;
            Gal[ngal].SfrBulgeColdGas[step] = Gal[ngal].SfrBulgeColdGasMetals[step] = 0.0;
          }

          if(halonr == Halo[halonr].FirstHaloInFOFgroup)
          {
            // a central galaxy
            Gal[ngal].mergeType = 0;
            Gal[ngal].mergeIntoID = -1;
            Gal[ngal].MergTime = 999.9;            

              Gal[ngal].Vvir = get_virial_velocity(halonr, ngal); // Keeping Vvir fixed for subhaloes as a new thing
              
            Gal[ngal].DiskScaleRadius = get_disk_radius(halonr, ngal); // Only update the scale radius for centrals.  Satellites' spin will be too variable and untrustworthy for new cooling.
              
            SpinMag = pow(pow(Halo[halonr].Spin[0], 2.0) + pow(Halo[halonr].Spin[1], 2.0) + pow(Halo[halonr].Spin[2], 2.0), 0.5);
            if(SpinMag>0)
            {
                for(j = 0; j < 3; j++) // Also only update the spin direction of the hot gas for centrals
                    Gal[ngal].SpinHot[j] = Halo[halonr].Spin[j] / SpinMag;
            }

            Gal[ngal].Type = 0;
          }
          else
          {
            // a satellite with subhalo
              Gal[ngal].mergeType = 0;
              Gal[ngal].mergeIntoID = -1;
              
            if(Gal[ngal].Type == 0)  // remember the infall properties before becoming a subhalo
            {
              Gal[ngal].infallMvir = previousMvir;
              Gal[ngal].infallVvir = previousVvir;
              Gal[ngal].infallVmax = previousVmax;
            }

            if(Gal[ngal].Type == 0 || Gal[ngal].MergTime > 999.0)
			{
              // here the galaxy has gone from type 1 to type 2 or otherwise doesn't have a merging time.
              Gal[ngal].MergTime = estimate_merging_time(halonr, Halo[halonr].FirstHaloInFOFgroup, ngal);
			}
            Gal[ngal].Type = 1;
          }
        }
        else
        {
          // an orhpan satellite galaxy
          Gal[ngal].deltaMvir = -1.0*Gal[ngal].Mvir;
          Gal[ngal].Mvir = 0.0;

          if(Gal[ngal].MergTime > 999.0 || Gal[ngal].Type == 0)
          {
            // Here the galaxy has gone from type 0 to type 2. Merge it!
            Gal[ngal].MergTime = 0.0;
          
            Gal[ngal].infallMvir = previousMvir;
            Gal[ngal].infallVvir = previousVvir;
            Gal[ngal].infallVmax = previousVmax;
          }

          Gal[ngal].Type = 2;
        }
      }

      // Note: Galaxies that are already type 2 do not need special treatment at this point 
      if(Gal[ngal].Type < 0 || Gal[ngal].Type > 2)
      {
        printf("what's that????\n");
        ABORT(88);
      }

      ngal++;

    }

    prog = Halo[prog].NextProgenitor;
  }

  if(ngal == 0)
  {
    // We have no progenitors with galaxies. This means we create a new galaxy. 
    init_galaxy(ngal, halonr);
    ngal++;
  }

  // Per Halo there can be only one Type 0 or 1 galaxy, all others are Type 2 
  // In fact this galaxy is very likely to be the first galaxy in the halo if first_occupied==FirstProgenitor 
  // and the Type0/1 galaxy in FirstProgenitor was also the first one 
  // This cannot be guaranteed though for the pathological first_occupied != FirstProgenitor case 

  for(i = ngalstart, centralgal = -1; i < ngal; i++)
  {
    if(Gal[i].Type == 0 || Gal[i].Type == 1)
    {
      if(centralgal != -1)
      {
        printf("can't be\n");
        ABORT(8);
      }
      centralgal = i;
    }
  }

  for(i = ngalstart; i < ngal; i++)
    Gal[i].CentralGal = centralgal;

  return ngal;

}



void evolve_galaxies(int halonr, int ngal, int tree)	// note: halonr is here the FOF-background subhalo (i.e. main halo) 
{
  int p, i, step, centralgal, merger_centralgal, currenthalo, offset;
  double infallingGas, coolingGas, deltaT, time, galaxyBaryons, currentMvir, DiscGasSum;
    double max_hotstrip;

  centralgal = Gal[0].CentralGal;
  if(Gal[centralgal].Type != 0 || Gal[centralgal].HaloNr != halonr)
  {
    printf("Something wrong here ..... \n");
    ABORT(54);
  }

  // basically compute the diff between the hot gas obtained at the end of the 
  // previous snapshot and the one obtained at the beginning of the new snapshot
  // using the conservation of baryons 

  infallingGas = infall_recipe(centralgal, ngal, ZZ[Halo[halonr].SnapNum]);
 
  // Reset SFRs for each galaxy
  for(p = 0; p < ngal; p++)
  {
      for(i=0; i<N_BINS; i++)
          Gal[p].DiscSFR[i] = 0.0;
  }
    
  // we integrate things forward by using a number of intervals equal to STEPS 
  for(step = 0; step < STEPS; step++)
  {

    // Loop over all galaxies in the halo 
    for(p = 0; p < ngal; p++)
    {
	  DiscGasSum = get_disc_gas(p);
	  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
	  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
	  assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);

        if(step==0 || HotStripOn<2) max_hotstrip = Gal[p].HotGas / STEPS;
        
      // don't treat galaxies that have already merged 
      if(Gal[p].mergeType > 0)
        continue;

      deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
      time = Age[Gal[p].SnapNum] - (step + 0.5) * (deltaT / STEPS);
      
      if(Gal[p].dT < 0.0)
        Gal[p].dT = deltaT;
        
	  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
        assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);

      // for central galaxy only
      if(p == centralgal)
      {
        add_infall_to_hot(centralgal, infallingGas / STEPS);
        assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
          assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
      }
      else if(HotStripOn>0 && Gal[p].Type == 1 && Gal[p].HotGas > 0.0 && max_hotstrip>0.0)
      {
            max_hotstrip = strip_from_satellite(halonr, centralgal, p, max_hotstrip);
            assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
          assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
      }
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);

      assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
      if(ReIncorporationFactor > 0.0)
        reincorporate_gas(p, deltaT / STEPS);
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);
        
      // Ram pressure stripping of cold gas from satellites
      if(RamPressureOn>0 && Gal[p].Type == 1 && Gal[p].ColdGas>0.0 && (Gal[p].ColdGas+Gal[p].StellarMass)>Gal[p].HotGas)
          ram_pressure_stripping(centralgal, p);
	
	  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
        assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);
        
      // determine cooling gas given halo properties
      // Preventing cooling when there's no angular momentum, as the code isn't built to handle that.  This only cropped up for Vishnu haloes with 2 particles, which clearly weren't interesting/physical.  Haloes always have some spin otherwise.
        if(!(Gal[p].SpinHot[0]==0 & Gal[p].SpinHot[1]==0 && Gal[p].SpinHot[2]==0))
        {
          coolingGas = cooling_recipe(p, deltaT / STEPS);
          Gal[p].AccretedGasMass += coolingGas;
          cool_gas_onto_galaxy(p, centralgal, coolingGas, deltaT / STEPS, step);
        }
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);

      DiscGasSum = get_disc_gas(p);
      assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
	  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
	  assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);
	  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
        assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);

        // Update radii of the annuli
        if(Gal[p].Mvir > 0 && Gal[p].Rvir > 0) update_disc_radii(p);
	
	  // stars form and then explode!
      //if(p==centralgal)
        starformation_and_feedback(p, centralgal, time, deltaT / STEPS, halonr, step);

      DiscGasSum = get_disc_gas(p);
	  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
	  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
	  assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);
	  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
        assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
        assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);
        
      // precess gas disc
      if(GasPrecessionOn && Gal[p].StellarMass>0.0 && Gal[p].ColdGas>0.0)
        precess_gas(p, deltaT / STEPS, halonr);
      assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);
	
    }

    // check for satellite disruption and merger events 
    for(p = 0; p < ngal; p++)
    {

      if((Gal[p].Type == 1 || Gal[p].Type == 2) && Gal[p].mergeType == 0)  // satellite galaxy!
      {
        if(Gal[p].MergTime > 999.0)
        {
          printf("satellite doesn't have a merging time! %f\n", Gal[p].MergTime);
          //ABORT(77);
        }

        deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
        Gal[p].MergTime -= deltaT / STEPS;
        
        // only consider mergers or disruption for halo-to-baryonic mass ratios below the threshold
        // or for satellites with no baryonic mass (they don't grow and will otherwise hang around forever)
        currentMvir = Gal[p].Mvir - Gal[p].deltaMvir * (1.0 - ((double)step + 1.0) / (double)STEPS);
        galaxyBaryons = Gal[p].StellarMass + Gal[p].ColdGas;
        if((galaxyBaryons == 0.0) || (galaxyBaryons > 0.0 && (currentMvir / galaxyBaryons <= ThresholdSatDisruption)))        
        {
          if(Gal[p].Type==1) 
            merger_centralgal = centralgal;
          else
            merger_centralgal = Gal[p].CentralGal;

          if(Gal[merger_centralgal].mergeType > 0) 
            merger_centralgal = Gal[merger_centralgal].CentralGal;

          Gal[p].mergeIntoID = NumGals + merger_centralgal;  // position in output 

          if(Gal[p].MergTime > 0.0)  // disruption has occured!
            disrupt_satellite_to_ICS(merger_centralgal, p);
          else
          {
            time = Age[Gal[p].SnapNum] - (step + 0.5) * (deltaT / STEPS);
            deal_with_galaxy_merger(p, merger_centralgal, centralgal, time, deltaT / STEPS, halonr, step);
          }
 
        }
     
      }
	  DiscGasSum = get_disc_gas(p);
	  assert(DiscGasSum <= 1.001*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.001);
	  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
	  assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);
    }

  } // end move forward in interval STEPS 


  // extra miscellaneous stuff 
  deltaT = Age[Gal[0].SnapNum] - Age[Halo[halonr].SnapNum];
  Gal[centralgal].TotalSatelliteBaryons = 0.0;
  for(p = 0; p < ngal; p++)
  {
    Gal[p].Cooling /= deltaT;
    Gal[p].Heating /= deltaT;
    Gal[p].OutflowRate /= deltaT;
      
    if(Gal[p].Mvir > 0 && Gal[p].Rvir > 0)
    {
        if(Gal[p].Type==0) update_disc_radii(p);
        update_HI_H2(p);
    }

    if(p != centralgal)
		Gal[centralgal].TotalSatelliteBaryons += (Gal[p].StellarMass + Gal[p].BlackHoleMass + Gal[p].ColdGas + Gal[p].HotGas);
  }


  // attach final galaxy list to halo 
  offset = 0;
  for(p = 0, currenthalo = -1; p < ngal; p++)
  {
    
	assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);
    if(Gal[p].HaloNr != currenthalo)
    {
      currenthalo = Gal[p].HaloNr;
      HaloAux[currenthalo].FirstGalaxy = NumGals;
      HaloAux[currenthalo].NGalaxies = 0;
    }
      
    // Merged galaxies won't be output. So go back through its history and find it
    // in the previous timestep. Then copy the current merger info there.
    offset = 0;
    i = p-1;
    while(i >= 0)
    {
     if(Gal[i].mergeType > 0) 
       if(Gal[p].mergeIntoID > Gal[i].mergeIntoID)
         offset++;  // these galaxies won't be kept so offset mergeIntoID below
     i--;
    }
    
    i = -1;
    if(Gal[p].mergeType > 0)
    {
      i = HaloAux[currenthalo].FirstGalaxy - 1;
      while(i >= 0)
      {
        if(HaloGal[i].GalaxyNr == Gal[p].GalaxyNr)
          break;
        else
          i--;
      }
      
      if(i < 0)
      {
        printf("Ran over the end of HaloGal looking for progenitor!\n");
        ABORT(12);
      }
      
      HaloGal[i].mergeType = Gal[p].mergeType;
      HaloGal[i].mergeIntoID = Gal[p].mergeIntoID - offset;
      HaloGal[i].mergeIntoSnapNum = Halo[currenthalo].SnapNum;
    }
    
    if(Gal[p].mergeType == 0)
    {
      if(NumGals >= MaxGals)
      {
        printf("maximum number of galaxies reached... exiting.\n");
        ABORT(1);
      }

      Gal[p].SnapNum = Halo[currenthalo].SnapNum;
      HaloGal[NumGals++] = Gal[p];
      HaloAux[currenthalo].NGalaxies++;
    }
	assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);
  }

  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
}






