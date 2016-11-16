#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>

#include "core_allvars.h"
#include "core_proto.h"


#define TREE_MUL_FAC        (1000000000LL)
#define FILENR_MUL_FAC      (1000000000000000LL)

// keep a static file handle to remove the need to do constant seeking.
FILE* save_fp[ABSOLUTEMAXSNAPS] = { NULL };
int save_fd[ABSOLUTEMAXSNAPS] = {-1};


void save_galaxies(int filenr, int tree)
{
#ifndef MINIMIZE_IO
    char buf[1000];
#endif
    int i, n;
    struct GALAXY_OUTPUT galaxy_output;
    int OutputGalCount[MAXSNAPS], *OutputGalOrder;
    
    OutputGalOrder = (int*)malloc( NumGals*sizeof(int) );
    assert( OutputGalOrder );
    
    // reset the output galaxy count and order
    for(i = 0; i < MAXSNAPS; i++)
    OutputGalCount[i] = 0;
    for(i = 0; i < NumGals; i++)
    OutputGalOrder[i] = -1;
    
    // first update mergeIntoID to point to the correct galaxy in the output
    for(n = 0; n < NOUT; n++)
    {
        for(i = 0; i < NumGals; i++)
        {
            if(HaloGal[i].SnapNum == ListOutputSnaps[n])
            {
                OutputGalOrder[i] = OutputGalCount[n];
                OutputGalCount[n]++;
            }
        }
    }
    
    for(i = 0; i < NumGals; i++)
    if(HaloGal[i].mergeIntoID > -1)
    HaloGal[i].mergeIntoID = OutputGalOrder[HaloGal[i].mergeIntoID];
    
    // now prepare and write galaxies
    for(n = 0; n < NOUT; n++)
    {
#ifndef MINIMIZE_IO
        // only open the file if it is not already open.
        if( !save_fp[n] )
        {
            sprintf(buf, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);
            
            if(!(save_fp[n] = fopen(buf, "r+")))
            {
                printf("can't open file `%s'\n", buf);
                ABORT(0);
            }
            
            save_fd[n] = fileno(save_fp[n]);
            
            // write out placeholders for the header data.
            size_t size = (Ntrees + 2)*sizeof(int);
            int* tmp_buf = (int*)malloc( size );
            memset( tmp_buf, 0, size );
            pwrite( tmp_buf, sizeof(int), Ntrees + 2, save_fp[n] );
            free( tmp_buf );
        }
#endif
        
        for(i = 0; i < NumGals; i++)
        {
            if(HaloGal[i].SnapNum == ListOutputSnaps[n])
            {
                prepare_galaxy_for_output(filenr, tree, &HaloGal[i], &galaxy_output);
//                myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, save_fp[n]);
                ssize_t pout = mypwrite(save_fd[n], &galaxy_output, sizeof(struct GALAXY_OUTPUT)*1, fd_offsets[n]);
                fd_offsets[n] += sizeof(struct GALAXY_OUTPUT)*1;
                
                TotGalaxies[n]++;
                TreeNgals[n][tree]++;
            }
        }
        
    }
    
    // don't forget to free the workspace.
    free( OutputGalOrder );
    
}



void prepare_galaxy_for_output(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o)
{
  int j, step;
    
  o->SnapNum = g->SnapNum;
  o->Type = g->Type;
    
    assert( g->GalaxyNr < TREE_MUL_FAC ); // breaking tree size assumption
    assert(tree < FILENR_MUL_FAC/TREE_MUL_FAC);
    o->GalaxyIndex = g->GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
    assert( (o->GalaxyIndex - g->GalaxyNr - TREE_MUL_FAC*tree)/FILENR_MUL_FAC == filenr );
    assert( (o->GalaxyIndex - g->GalaxyNr -FILENR_MUL_FAC*filenr) / TREE_MUL_FAC == tree );
    assert( o->GalaxyIndex - TREE_MUL_FAC*tree - FILENR_MUL_FAC*filenr == g->GalaxyNr );
    
    o->CentralGalaxyIndex = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
    
    o->HaloIndex = g->HaloNr;
    o->TreeIndex = tree;
    o->SimulationHaloIndex = Halo[g->HaloNr].SubhaloIndex;


  o->mergeType = g->mergeType;
  o->mergeIntoID = g->mergeIntoID;
  o->mergeIntoSnapNum = g->mergeIntoSnapNum;
  o->dT = g->dT * UnitTime_in_s / SEC_PER_MEGAYEAR;

  for(j = 0; j < 3; j++)
  {
    o->Pos[j] = g->Pos[j];
    o->Vel[j] = g->Vel[j];
    o->Spin[j] = Halo[g->HaloNr].Spin[j];
    o->SpinStars[j] = g->SpinStars[j];
    o->SpinGas[j] = g->SpinGas[j];
//      o->SpinSecularBulge[j] = g->SpinSecularBulge[j];
      o->SpinClassicalBulge[j] = g->SpinClassicalBulge[j];
  }

  o->Len = g->Len;
    o->LenMax = g->LenMax;
  o->Mvir = g->Mvir;
  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup, -1);
  o->Rvir = get_virial_radius(g->HaloNr, -1);  //output the actual Rvir, not the maximum Rvir
  o->Vvir = get_virial_velocity(g->HaloNr, -1);  //output the actual Vvir, not the maximum Vvir
  o->Vmax = g->Vmax;
  o->VelDisp = Halo[g->HaloNr].VelDisp;
    
    for(j=0; j<N_BINS+1; j++)
    o->DiscRadii[j] = g->DiscRadii[j];
    
  o->ColdGas = g->ColdGas;
  o->StellarMass = g->StellarMass;
  o->ClassicalBulgeMass = g->ClassicalBulgeMass;
  o->SecularBulgeMass = g->SecularBulgeMass;
  o->HotGas = g->HotGas;
  o->EjectedMass = g->EjectedMass;
  o->BlackHoleMass = g->BlackHoleMass;
  o->ICS = g->ICS;

  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsStellarMass = g->MetalsStellarMass;
  o->ClassicalMetalsBulgeMass = g->ClassicalMetalsBulgeMass;
  o->SecularMetalsBulgeMass = g->SecularMetalsBulgeMass;
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;
  o->MetalsICS = g->MetalsICS;
    
//  o->StarsInSitu = g->StarsInSitu;
//  o->StarsInstability = g->StarsInstability;
//  o->StarsMergeBurst = g->StarsMergeBurst;
    
//    o->AccretedGasMass = g->AccretedGasMass;
//    o->EjectedSNGasMass = g->EjectedSNGasMass;
//    o->EjectedQuasarGasMass = g->EjectedQuasarGasMass;

//    o->TotInstabEvents = g->TotInstabEvents;
//    o->TotInstabEventsGas = g->TotInstabEventsGas;
//    o->TotInstabEventsStar = g->TotInstabEventsStar;
//    o->TotInstabAnnuliGas = g->TotInstabAnnuliGas;
//    o->TotInstabAnnuliStar = g->TotInstabAnnuliStar;
//    
//    if(g->TotInstabEventsGas)
//        o->FirstUnstableAvGas = (1.0*g->FirstUnstableGas) / (1.0*g->TotInstabEventsGas);
//    else
//        o->FirstUnstableAvGas = 0.0;
//    
//    if(g->TotInstabEventsStar)
//        o->FirstUnstableAvStar = (1.0*g->FirstUnstableStar) / (1.0*g->TotInstabEventsStar);
//    else
//        o->FirstUnstableAvStar = 0.0;
  
  for(j=0; j<N_BINS; j++)
  {
	o->DiscGas[j] = g->DiscGas[j];
	o->DiscStars[j] = g->DiscStars[j];
	o->DiscGasMetals[j] = g->DiscGasMetals[j];
	o->DiscStarsMetals[j] = g->DiscStarsMetals[j];
      o->DiscHI[j] = g->DiscHI[j];
      o->DiscH2[j] = g->DiscH2[j];
      o->DiscSFR[j] = g->DiscSFR[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
//      o->TotSinkGas[j] = g->TotSinkGas[j];
//      o->TotSinkStar[j] = g->TotSinkStar[j];
  }

  o->SfrDisk = 0.0;
  o->SfrBulge = 0.0;
  o->SfrDiskZ = 0.0;
  o->SfrBulgeZ = 0.0;
  
  // NOTE: in Msun/yr 
  for(step = 0; step < STEPS; step++)
  {
    o->SfrDisk += g->SfrDisk[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    o->SfrBulge += g->SfrBulge[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    
    if(g->SfrDiskColdGas[step] > 0.0)
      o->SfrDiskZ += g->SfrDiskColdGasMetals[step] / g->SfrDiskColdGas[step] / STEPS;

    if(g->SfrBulgeColdGas[step] > 0.0)
      o->SfrBulgeZ += g->SfrBulgeColdGasMetals[step] / g->SfrBulgeColdGas[step] / STEPS;
  }

  o->DiskScaleRadius = g->DiskScaleRadius;
  
  if(g->ClassicalBulgeMass + g->SecularBulgeMass > 0.0)
    o->BulgeEffectiveRadius = 
      ((g->ClassicalBulgeMass * g->ClassicalBulgeRadius) + (g->SecularBulgeMass * 0.2*g->DiskScaleRadius)) /
         (g->ClassicalBulgeMass + g->SecularBulgeMass);
  else 
    o->BulgeEffectiveRadius = 0.0;
    

  if (g->Cooling > 0.0)
    o->Cooling = log10(g->Cooling * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Cooling = 0.0;
  if (g->Heating > 0.0)
    o->Heating = log10(g->Heating * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Heating = 0.0;

  o->LastMajorMerger = g->LastMajorMerger * UnitTime_in_Megayears;
  o->LastMinorMerger = g->LastMinorMerger * UnitTime_in_Megayears;
  o->OutflowRate = g->OutflowRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;

  //infall properties
  if(g->Type != 0)
  {
    o->infallMvir = g->infallMvir;
    o->infallVvir = g->infallVvir;
    o->infallVmax = g->infallVmax;
  }
  else
  {
    o->infallMvir = 0.0;
    o->infallVvir = 0.0;
    o->infallVmax = 0.0;
  }

}



void finalize_galaxy_file(int filenr)
{
    int n;
    
    for(n = 0; n < NOUT; n++)
    {
#ifndef MINIMIZE_IO
        // file must already be open.
        assert( save_fp[n] );
        fsync(save_fp[n]);
        
        // seek to the beginning.
        fseek( save_fp[n], 0, SEEK_SET );
#endif
        myfwrite(&Ntrees, sizeof(int), 1, save_fp[n]);
        myfwrite(&TotGalaxies[n], sizeof(int), 1, save_fp[n]);
        myfwrite(TreeNgals[n], sizeof(int), Ntrees, save_fp[n]);
        
#ifndef MINIMIZE_IO
        // close the file and clear handle after everything has been written
        fclose( save_fp[n] );
        save_fp[n] = NULL;
#else
        write_galaxy_data_snap(n, filenr);
#endif
    }
    
}

#undef TREE_MUL_FAC
#undef FILENR_MUL_FAC

