//****************************************************************************
// CDICE
//     C/C++ VERSION OF DICE
//             Economic - Climate Projections
//             Based on GAMS version of DICE2013R
//				(DICE2013Rv2_102213_vanilla_v24b.gms)
//				
//            William D. Nordhaus, Yale University
//				The Climate Casino (2013), Yale University Press
//              http://www.econ.yale.edu/~nordhaus/homepage/index.html
//			 for previous versions, see also
//				DICE 2007:
//				William D. Nordhaus
//				A Question of Balance (2008), Yale University Press
//				DICE1999:
//				William D. Norhaus and Joseph Boyer
//				Warming the World (2000), The MIT Press
//				DICE 1994:
//				William D. Nordhaus,
//				Managing the Global Commons (1994), The MIT Press
//            
//           
//          Translated to C/C++ by Gregory Garner
//            V 1.0  2013 Penn State University
//
//*****************************************************************************/
/* Copyright (C) 2013 Gregory Garner, Klaus Keller, Patrick Reed, Martha Butler, and others

  CDICE is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at your
  option) any later version.
 
  CDICE is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.
 
  You should have received a copy of the GNU Lesser General Public License
  along with the CDICE code.  If not, see <http://www.gnu.org/licenses/>.
 */


//#ifdef CDICEMOD


#include "CDICE.h"

DICE dice;
DICE *dicePtr = &dice;




//=================================================================================
//
// calc_CDICE(DICE *dice)
//
// Run the full model.  
//
// Return()
// - 'dice' structure variables are updated
//=================================================================================
void calc_CDICE(DICE *dice)
{
	// Run the full model
	processModel(dice);
	
	return;
}


//=================================================================================
//
// optim_CDICE(double* vars, double* objs, double* consts)
//
// Function over which BORG optimizes.
//
// Return()
// - Values of objectives and constraint returned in 'objs' and 'consts'.
//=================================================================================
void optim_CDICE(double* vars, double* objs, double* consts)
{
  // Initialize for loop iteration variable
  int i;
  int pind;
  
  // Initialize constraint values
  int dvaroutofbounds = 0;
  
  // Assign the decision variable vector to the control rate and savings rate
  for(i=0; i < dicePtr->ctrl.nvars; i++) {
    if (i < 59) {
      if (dicePtr->ctrl.policy_inertia <= 0.0) {
      	dicePtr->dvars.miu[i+1] = vars[i];
      }
      else {
	      dicePtr->dvars.miu[i+1] = dicePtr->dvars.miu[i] + (vars[i]*dicePtr->ctrl.policy_inertia);
	    }
      if(dicePtr->dvars.miu[i+1] < 0.0 || dicePtr->dvars.miu[i+1] > dicePtr->lims.miu_up[i+1]) {
      	dvaroutofbounds = 1;
      }
    }
    else {
      dicePtr->dvars.s[i-59] = vars[i];
      if(vars[i] < 0.0 || vars[i] > 1.0) {
      	dvaroutofbounds = 1;
      }
    }
  }
  
  // Initialize some variables to hold information for later
  // calculations ------------------------------------------------
  
  // Sum of utility for expected utility calculation (Obj 1)
  double sum_util = 0.0;
  
  // Max, sum, and expectation of temperature (Obj 2, Constraint)
  double maxTemp = -9999.99;
  
  // Max, sum, and expectation of carbon emissions (Constraint)
  double max_cca = -1.0;

  // Sum of damages and abatement (Obj 3,4)
  double sum_dam = 0.0;
  double sum_abate = 0.0;
  double *dam_vec = new double[dicePtr->ctrl.csSamples];
  double *abate_vec = new double[dicePtr->ctrl.csSamples];
  
  // Number of times temperature goes above 2C (Obj 2)
  int nTempUnder2C = 0;
  
  // Number of times the fossil fuel limit is exceeded (Constraint)
  int nCCAOverFosslim = 0;
  
  // Loop over the states of climate sensitivity (if needed)
  for (i=0; i < dicePtr->ctrl.csSamples; i++) {
    
    // Reset the maxtemp and maxcca variables
    maxTemp = -9999.99;
    max_cca = -1.0;
    
    // Are we using DOEclim?
    if(dicePtr->ctrl.use_doeclim == 1) {
      
      // Assign the cs sample and ocean diff
      dicePtr->clim.t2co = dicePtr->clim.cs_sample_vec[i];
      dicePtr->clim.kappa = dicePtr->clim.ocean_diff_vec[i];
      dicePtr->clim.alpha = dicePtr->clim.alpha_vec[i];
      doeclim_DICE_init(dicePtr);
      
    }
    else {
      
      // Assign this cs sample
      dicePtr->clim.t2xco2 = dicePtr->clim.cs_sample_vec[i];
      
    }
    
    // Run the model with this sample
    calc_CDICE(dicePtr);
    
    // Store the discounted utility in the utility vector
    sum_util += dicePtr->econ.utility;
    
    // Find the maximum temperature increase and cumulative 
    // carbon emisions reached and store it
    int iter = 0;
    if(dicePtr->ctrl.use_doeclim == 1) {
      iter = dicePtr->config.tstep * dicePtr->config.nPeriods;
      for (pind=0; pind < iter; pind++) {
				if (dicePtr->clim.temp[pind+dicePtr->clim.hind_ns] > maxTemp) {
	  			maxTemp = dicePtr->clim.temp[pind+dicePtr->clim.hind_ns];
				}
      }
    }
    else {
      iter = dicePtr->config.nPeriods;
      for (pind=0; pind < iter; pind++) {
				if (dicePtr->clim.tatm[pind] > maxTemp) {
					maxTemp = dicePtr->clim.tatm[pind];
				}
      }
    }
    for(pind=0; pind < dicePtr->config.nPeriods; pind++) {
      if (dicePtr->econ.cca[pind] > max_cca) {
				max_cca = dicePtr->econ.cca[pind];
      }
    }
	   
    // Determine if this max temp is over 2C
    if (maxTemp < 2.0) {
      nTempUnder2C++;
    }
	   
    // Determine if these emissions exceed the fossil fuel limit
    if (max_cca > dicePtr->carb.fosslim) {
      nCCAOverFosslim++;
    }
	   
    // Populate the damage and abatement cost vector
    dam_vec[i] = dicePtr->econ.pv_damages[dicePtr->config.nPeriods-1];
    abate_vec[i] = dicePtr->econ.pv_abatecost[dicePtr->config.nPeriods-1];
  }
	
			
  // Calculate and assign the objective values
  // Maximize Expectation of discounted utility
  objs[0] = -1.0 * sum_util / (double) dicePtr->ctrl.csSamples;

  // Probability of T < 2C
  objs[1] = -1.0 * (double) nTempUnder2C / (double) dicePtr->ctrl.csSamples;


  // Minimize the net present value of damages and abatement
  for(i=0; i<dicePtr->ctrl.csSamples; i++) {
    sum_dam += dam_vec[i];
    sum_abate += abate_vec[i];
  }
  objs[2] = sum_dam / (double) dicePtr->ctrl.csSamples;
  objs[3] = sum_abate / (double) dicePtr->ctrl.csSamples;

	
	// Assign Constraints ------------------------------------------------
	
	// Decision variable was out of bounds
	consts[0] = (double) dvaroutofbounds;
	
  // Evaluate the fossil fuel constraint
  consts[1] = (double) nCCAOverFosslim;

  // Clean up -----------------------------------------------------------
  delete dam_vec;
  delete abate_vec;


  return;
}


//=================================================================================
//
// main()
//
// Main function.
//
// Return(0) if successful
//=================================================================================
int main(int argc, char* argv[])
{
	// DEBUG - Remove the test_out.txt file if it exists
//	int ret = remove("test_out.txt");

  // Initialize a CDICE model run
  allocateConfig(dicePtr);
  allocateCtrl(dicePtr);
  allocateDvars(dicePtr);
  ctrl_CDICE(dicePtr, "cdice_control.txt");
  allocateCarb(dicePtr);
  allocateClim(dicePtr);
  allocateEcon(dicePtr);
  allocateLims(dicePtr);
  initTS(dicePtr);

  // If the user wants to use the DOEclim model,
  // load the forcings and calibration for the hindcast.
  if(dicePtr->ctrl.use_doeclim == 1) {
    doeclim_load_hind_forc(dicePtr);
    doeclim_load_calibration(dicePtr);
  }
	
  // Initialize variables for this run
  int i, j, n;
  int rank;
  double segsize, cdf_val;
  char runtime[256];
  char borgresults[256];
  char ncfile[256];
	
  // Sample from the CS distribution to create the states of the world (sow).
  // These sow will be used to introduce uncertainty into DICE and 
  // allow us to estimate the objective of P(T>=2K). 
  if (dicePtr->ctrl.csSamples == 1)
    {
      if(dicePtr->ctrl.use_doeclim == 1) {
				dicePtr->clim.cs_sample_vec[0] =  dicePtr->clim.t2co;
      }
      else {
				dicePtr->clim.cs_sample_vec[0] =  dicePtr->clim.t2xco2;
      }
    }
	
  else
    {
      segsize = 1.0 / (double) dicePtr->ctrl.csSamples;
      for (i=0; i < dicePtr->ctrl.csSamples; i++) {
	  		cdf_val = ((double) i + 0.5) * segsize;
	  		dicePtr->clim.cs_sample_vec[i] =  quantile(dicePtr->config.dist, cdf_val);
			}
    }
  
  // Apply the adjustment to ocean diffusivity and the aerosol forcing factor
  // using the calibration data
  if(dicePtr->ctrl.use_doeclim == 1) {
	  apply_calibration(dicePtr);
	}
  
  // Create a BORG problem
  BORG_Problem problem = BORG_Problem_create(dicePtr->ctrl.nvars, dicePtr->ctrl.nobjs, dicePtr->ctrl.nconsts, optim_CDICE);
  
  // Set the bounds on the decision variables
  for (j=0; j < dicePtr->ctrl.nvars; j++) {
    if (j < 59) {
    	if(dicePtr->ctrl.policy_inertia <= 0.0) {
	      BORG_Problem_set_bounds(problem, j, 0.0, dicePtr->lims.miu_up[j+1]);
	    }
	    else {
	    	BORG_Problem_set_bounds(problem, j, -1.0, 1.0);
	    }
    }
    else {
      BORG_Problem_set_bounds(problem, j, 0.0, 1.0);
    }
  }
	
  // Set the epsilon value for each objective
  BORG_Problem_set_epsilon(problem, 0, 0.1);		// Discounted Utility
  BORG_Problem_set_epsilon(problem, 1, 0.002);		// P(T<2K)
  BORG_Problem_set_epsilon(problem, 2, 0.05);		// PV Damages
  BORG_Problem_set_epsilon(problem, 3, 0.05);		// PV Abatement
	
  // This sets up an experiment with multiple runs and seeds. As of now,
  // it's simply set to run once.
  for(n=0; n<1; n++) {
    
    // Set the seed for BORG
    BORG_Random_seed(13*n*(rank+1));
    
    // Run the optimization
		BORG_Archive result = BORG_Algorithm_run(problem, dicePtr->ctrl.nevals);
      
    // Open a file stream for the BORG results
    FILE* borg_results;
    sprintf(borgresults, "borg_results_%d.txt", n);
    borg_results = fopen(borgresults, "w");
    if(borg_results == NULL) {
			perror ("Failed to open borg results file. Exiting...");
			exit(4);
    }
      
    // Print the results in the archive to the file stream
    // and, if set in the control file, dump the model
    // output to a netCDF file
    BORG_Archive_print(result, borg_results);
    fclose(borg_results);      
    if(dicePtr->ctrl.writeModelDump == 1) {
			sprintf(ncfile, "dice_model_output_%d.nc", n);
			WriteNetCDF(dicePtr, result, ncfile);
    }
    BORG_Archive_destroy(result);
      
  }
  
  // Free the BORG objects
  BORG_Problem_destroy(problem);
  

  return 0;
}


//#endif //CDICEMOD
