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

#include "CDICE.h"

extern DICE dice;
extern DICE* dicePtr;

//=========================================================================
//
// invCDFSample(lndist_type &dist, rgen_type &gen, double lb, double ub)
//
// Returns a value of the inverted CDF of dist using a random probability
// produced by gen between the bounds lb and ub.
//
// Return(x) The value of the inverted CDF.
//=========================================================================

double invCDFSample(lndist_type &dist, rgen_type &gen, double lb, double ub)
{
  boost::uniform_real<> runif(lb, ub);
  double x = quantile(dist, runif(gen));
  return(x);
}


//========================================================================
//
// ncError(int val)
//
// Prints the netcdf error to stdout and exits with an error code
//
//========================================================================
void ncError(int val)
{
  printf("Error: %s\n", nc_strerror(val));
  exit(2);
}
	

//=============================================================================
// calcSCC(DICE *dice1, BORG_Solution solution)
//
// Calculate the social cost of carbon in *dice1 using the provided
// BORG solution
//
//=============================================================================
void calcSCC(DICE *dice1, BORG_Solution solution)
{	
  // Make a temporary dice object
  DICE dice2;
  DICE *dice2ptr = &dice2;

  // Initialize a CDICE model run
  allocateConfig(dice2ptr);
  allocateCtrl(dice2ptr);
  allocateDvars(dice2ptr);
  ctrl_CDICE(dice2ptr, "cdice_control.txt");
  allocateCarb(dice2ptr);
  allocateClim(dice2ptr);
  allocateEcon(dice2ptr);
  allocateLims(dice2ptr);
  initTS(dice2ptr);
  
  if(dice1->ctrl.use_doeclim == 1) {
    doeclim_load_hind_forc(dice2ptr);
  }
  
  // Make temporary variables to hold the changes in utility
  double dutil_e;
  double dutil_c;
  int i;
  int t;
  int tsub;
  
  // Assign this control policy to the temp dice object
  for(i=0; i<dice2ptr->ctrl.nvars; i++) {
    if(i < 59) {
    	if(dice1->ctrl.policy_inertia > 0.0) {
    		dice1->dvars.miu[i+1] = dice1->dvars.miu[i] + (BORG_Solution_get_variable(solution, i) * dice1->ctrl.policy_inertia);
	      dice2ptr->dvars.miu[i+1] = dice2ptr->dvars.miu[i] + (BORG_Solution_get_variable(solution, i) * dice2ptr->ctrl.policy_inertia);
	    }
	    else {
      dice1->dvars.miu[i+1] = BORG_Solution_get_variable(solution, i);
      dice2ptr->dvars.miu[i+1] = BORG_Solution_get_variable(solution, i);
      }
    }
    else {
      dice1->dvars.s[i-59] = BORG_Solution_get_variable(solution, i);
      dice2ptr->dvars.s[i-59] = BORG_Solution_get_variable(solution, i);
    }
  }
  
  // Run the model over the dice object
  calc_CDICE(dice1);
  
  // Pass the updated configuration parameters to
  // the temp dice object
  if(dice1->ctrl.use_doeclim == 1) {
    dice2ptr->clim.t2co = dice1->clim.t2co;
    dice2ptr->clim.kappa = dice1->clim.kappa;
    dice2ptr->clim.alpha = dice1->clim.alpha;
    doeclim_DICE_init(dice2ptr);
  }
  else {
    dice2ptr->clim.t2xco2 = dice1->clim.t2xco2;
  }
  
  // Loop over the available time steps
  for(t=1; t<dice1->config.nPeriods; t++) {
    
    // Reset the temp dice object
    calc_CDICE(dice2ptr);
    
    // Increment the consumption at this time step
    dice2ptr->econ.c[t-1] += 1.0;
    
    // Do the post processing of this run
    postProcess(dice2ptr);
    
    // Store the utility from this run
    dutil_c = dice2ptr->econ.utility - dice1->econ.utility;
    
    // Reset the consumption
    dice2ptr->econ.c[t-1] = dice1->econ.c[t-1];        
    
    // Increment the emissions at this time step
    dice2ptr->econ.e[t-1] += 1.0;
    
    // Run the model on the temporary dice object
    for(tsub=t; tsub<dice1->config.nPeriods; tsub++) {
      calcCarb(dice2ptr, tsub);
      if(dice2ptr->ctrl.use_doeclim) {
	doeclim_dice_ts(dice2ptr, tsub);
      }
      else {
	calcClim(dice2ptr, tsub);
      }
      calcEcon(dice2ptr, tsub);
    }
    
    // Do the post processing of this run
    postProcess(dice2ptr);
    
    // Store the utility from this run
    dutil_e = dice2ptr->econ.utility - dice1->econ.utility;
    
    // Calculate the SCC from these differences in utility
    dice1->econ.scc[t-1] = -1000.0 * (dutil_e / dutil_c);
  }
  
  // Delete the temporary DICE object
  free_CDICE(dice2ptr);
  
  return;
}


//=========================================================================
//
// WriteNetCDF(DICE *dice, BORG_Archive result, string ncfile)
//
// Writes the full model output for a set of solutions and (if applicable)
// climate sentitivities.
//
//=========================================================================
void WriteNetCDF(DICE *dice, BORG_Archive result, const char *ncfile)
{
  // Loop control variables
  int iInd, jInd, tInd;
  
  // Was doeclim used?
  int use_doeclim = dice->ctrl.use_doeclim;  
  
  // Define the dimensions for the vectors
  int nPeriods = dice->config.nPeriods;
  int nSamples = dice->ctrl.csSamples;
  int nDoeclimts = dice->clim.ns;
  int nSolutions = BORG_Archive_get_size(result);
  
  int ndims = 3;

  // Dimension ids
  int nPeriodsid, nSamplesid, nSolutionsid, nDoeclimtsid;
  int fulldimids[ndims];
  int extfulldimids[ndims];
  int per_sol_dimids[2];
  int sol_sam_dimids[2];
  int sam_dimids[1];
  int per_dimids[1];
  
  // netcdf file id
  int ncid;
  
  // Variable ids
  int tatmid, matid, miuid, sid, utilityid, csid, oceandiffid, alphaid;
  int yearid, forcid, toceanid, muid, mlid, ccaid, kid, cid;
  int cpcid, iid, riid, yid, ygrossid, ynetid, damagesid, damfracid;
  int abatecostid, perioduid, sccid, lid, gaid, alid, gsigid;
  int sigmaid, cost1id, etreeid, eindid, rrid, forcothid, partfractid;
  int pv_abatecostid, pv_damagesid, totalcostid, pv_totalcostid;
  int pbacktimeid, eid, mcabateid, cemutotperid;
  
  int tempid, doeclim_forcid, heat_mixedid;
  
  // Initialize the variable vectors
  float tatm[nPeriods][nSolutions][nSamples]; 
  float mat[nPeriods][nSolutions][nSamples];
  float forc[nPeriods][nSolutions][nSamples];  
  float tocean[nPeriods][nSolutions][nSamples]; 
  float mu[nPeriods][nSolutions][nSamples];
  float ml[nPeriods][nSolutions][nSamples];
  float cca[nPeriods][nSolutions][nSamples];
  float k[nPeriods][nSolutions][nSamples];
  float c[nPeriods][nSolutions][nSamples];
	float cpc[nPeriods][nSolutions][nSamples];
  float i[nPeriods][nSolutions][nSamples];
  float ri[nPeriods][nSolutions][nSamples];
  float y[nPeriods][nSolutions][nSamples];
  float ygross[nPeriods][nSolutions][nSamples];
	float ynet[nPeriods][nSolutions][nSamples];
  float damages[nPeriods][nSolutions][nSamples];
  float damfrac[nPeriods][nSolutions][nSamples];
  float abatecost[nPeriods][nSolutions][nSamples];
  float periodu[nPeriods][nSolutions][nSamples];
  float scc[nPeriods][nSolutions][nSamples];
  float eind[nPeriods][nSolutions][nSamples];
  float pv_abatecost[nPeriods][nSolutions][nSamples];
  float pv_damages[nPeriods][nSolutions][nSamples];
  float totalcost[nPeriods][nSolutions][nSamples];
  float pv_totalcost[nPeriods][nSolutions][nSamples];
  float pbacktime[nPeriods][nSolutions][nSamples];
  float e[nPeriods][nSolutions][nSamples];
  float mcabate[nPeriods][nSolutions][nSamples]; 
  float cemutotper[nPeriods][nSolutions][nSamples];
  
  float temp[nDoeclimts][nSolutions][nSamples];
	float doeclim_forc[nDoeclimts][nSolutions][nSamples];
  float heat_mixed[nDoeclimts][nSolutions][nSamples];

  float miu[nPeriods][nSolutions];
	float utility[nSolutions][nSamples];
     
  float cs[nSamples];
  float oceandiff[nSamples];
  float alpha[nSamples];

  short year[nPeriods];
  float l[nPeriods];
  float ga[nPeriods];
  float al[nPeriods];
  float gsig[nPeriods];
  float sigma[nPeriods];
  float cost1[nPeriods];
  float etree[nPeriods];
  float rr[nPeriods];
  float forcoth[nPeriods];
  float partfract[nPeriods];
  
  // If the savings rate is set from the control file,
  // we will only need to provide a single dimension
  // for the savings rate output
  float s1[nPeriods];				// Savings rate if set in control file
  float s2[nPeriods][nSolutions];	// Savings rate if used as decision variable
  
  // Setup the netcdf file
  int retval;
  if((retval = nc_create(ncfile, NC_64BIT_OFFSET, &ncid))) { ncError(retval); }
  
  // Define the dimensions in the netcdf file
  if((retval = nc_def_dim(ncid, "nPeriods", nPeriods, &nPeriodsid))) { ncError(retval); }
  if((retval = nc_def_dim(ncid, "nSolutions", nSolutions, &nSolutionsid))) { ncError(retval); }
  if((retval = nc_def_dim(ncid, "nSamples", nSamples, &nSamplesid))) { ncError(retval); }
  
  if(use_doeclim == 1) {
    if((retval = nc_def_dim(ncid, "nDoeclimts", nDoeclimts, &nDoeclimtsid))) { ncError(retval); }
  }

  // Gather the dimension IDs into a vector
  fulldimids[0] = nPeriodsid;
  fulldimids[1] = nSolutionsid;
  fulldimids[2] = nSamplesid;
  
  extfulldimids[0] = nDoeclimtsid;
  extfulldimids[1] = nSolutionsid;
  extfulldimids[2] = nSamplesid;
  
  per_sol_dimids[0] = nPeriodsid;
  per_sol_dimids[1] = nSolutionsid;
  
  sol_sam_dimids[0] = nSolutionsid;
  sol_sam_dimids[1] = nSamplesid;
  
  sam_dimids[0] = nSamplesid;
  
  per_dimids[0] = nPeriodsid;
  
  // Create the netcdf variables
  if((retval = nc_def_var(ncid, "mat", NC_FLOAT, 3, fulldimids, &matid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "forc", NC_FLOAT, 3, fulldimids, &forcid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "mu", NC_FLOAT, 3, fulldimids, &muid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "ml", NC_FLOAT, 3, fulldimids, &mlid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "cca", NC_FLOAT, 3, fulldimids, &ccaid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "k", NC_FLOAT, 3, fulldimids, &kid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "c", NC_FLOAT, 3, fulldimids, &cid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "cpc", NC_FLOAT, 3, fulldimids, &cpcid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "i", NC_FLOAT, 3, fulldimids, &iid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "ri", NC_FLOAT, 3, fulldimids, &riid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "y", NC_FLOAT, 3, fulldimids, &yid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "ygross", NC_FLOAT, 3, fulldimids, &ygrossid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "ynet", NC_FLOAT, 3, fulldimids, &ynetid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "damages", NC_FLOAT, 3, fulldimids, &damagesid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "damfrac", NC_FLOAT, 3, fulldimids, &damfracid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "abatecost", NC_FLOAT, 3, fulldimids, &abatecostid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "periodu", NC_FLOAT, 3, fulldimids, &perioduid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "eind", NC_FLOAT, 3, fulldimids, &eindid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "pv_abatecost", NC_FLOAT, 3, fulldimids, &pv_abatecostid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "pv_damages", NC_FLOAT, 3, fulldimids, &pv_damagesid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "totalcost", NC_FLOAT, 3, fulldimids, &totalcostid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "pv_totalcost", NC_FLOAT, 3, fulldimids, &pv_totalcostid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "e", NC_FLOAT, 3, fulldimids, &eid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "mcabate", NC_FLOAT, 3, fulldimids, &mcabateid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "cemutotper", NC_FLOAT, 3, fulldimids, &cemutotperid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "scc", NC_FLOAT, 3, fulldimids, &sccid))) { ncError(retval); }

  if(use_doeclim == 1) {
    if((retval = nc_def_var(ncid, "temp", NC_FLOAT, 3, extfulldimids, &tempid))) { ncError(retval); }
    if((retval = nc_def_var(ncid, "doeclim_forc", NC_FLOAT, 3, extfulldimids, &doeclim_forcid))) { ncError(retval); }
    if((retval = nc_def_var(ncid, "heat_mixed", NC_FLOAT, 3, extfulldimids, &heat_mixedid))) { ncError(retval); }
  }
  else {
    if((retval = nc_def_var(ncid, "tatm", NC_FLOAT, 3, fulldimids, &tatmid))) { ncError(retval); }
    if((retval = nc_def_var(ncid, "tocean", NC_FLOAT, 3, fulldimids, &toceanid))) { ncError(retval); }
  }
  
  
  if((retval = nc_def_var(ncid, "miu", NC_FLOAT, 2, per_sol_dimids, &miuid))) { ncError(retval); }
  
  if((retval = nc_def_var(ncid, "utility", NC_FLOAT, 2, sol_sam_dimids, &utilityid))) { ncError(retval); }
  
  if((retval = nc_def_var(ncid, "cs", NC_FLOAT, 1, sam_dimids, &csid))) { ncError(retval); }  
  if(use_doeclim == 1) {
  	if((retval = nc_def_var(ncid, "oceandiff", NC_FLOAT, 1, sam_dimids, &oceandiffid))) { ncError(retval); }
  	if((retval = nc_def_var(ncid, "alpha", NC_FLOAT, 1, sam_dimids, &alphaid))) { ncError(retval); }
  }
  
  if((retval = nc_def_var(ncid, "Year", NC_SHORT, 1, per_dimids, &yearid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "l", NC_FLOAT, 1, per_dimids, &lid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "ga", NC_FLOAT, 1, per_dimids, &gaid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "al", NC_FLOAT, 1, per_dimids, &alid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "gsig", NC_FLOAT, 1, per_dimids, &gsigid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "sigma", NC_FLOAT, 1, per_dimids, &sigmaid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "cost1", NC_FLOAT, 1, per_dimids, &cost1id))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "etree", NC_FLOAT, 1, per_dimids, &etreeid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "rr", NC_FLOAT, 1, per_dimids, &rrid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "forcoth", NC_FLOAT, 1, per_dimids, &forcothid))) { ncError(retval); }
  if((retval = nc_def_var(ncid, "partfract", NC_FLOAT, 1, per_dimids, &partfractid))) { ncError(retval); }
  
  // If the savings rate is set in the control file,
  // the savings rate will only need one dimension
  if (dice->ctrl.setSavings < 0) {
    if((retval = nc_def_var(ncid, "s", NC_FLOAT, 2, per_sol_dimids, &sid))) { ncError(retval); }
  }
  else {
    if((retval = nc_def_var(ncid, "s", NC_FLOAT, 1, per_dimids, &sid))) { ncError(retval); }
  }
  
  // End "Metadata" mode
  if((retval = nc_enddef(ncid))) { ncError(retval); }
  
  // Populate the variable vectors with the appropriate data
  // First, loop over the BORG solutions
  for(iInd=0; iInd<nSolutions; iInd++) {

    // Get this BORG solution
    BORG_Solution this_solution = BORG_Archive_get(result, iInd);
    
    // Apply the decision vector to the dice object
    for(tInd=0; tInd<dice->ctrl.nvars; tInd++) {
      if(tInd < 59) {
      	if(dice->ctrl.policy_inertia <= 0.0) {
					dice->dvars.miu[tInd+1] = BORG_Solution_get_variable(this_solution, tInd);
				}
				else {
					dice->dvars.miu[tInd+1] = dice->dvars.miu[tInd] + (BORG_Solution_get_variable(this_solution, tInd) * dice->ctrl.policy_inertia);
				}
      }
      else {
	dice->dvars.s[tInd-59] = BORG_Solution_get_variable(this_solution, tInd);
      }
    }
    
    // Store the decision variables from this solution
    for(tInd=0; tInd<nPeriods; tInd++) {
      miu[tInd][iInd] = (float) dice->dvars.miu[tInd];
      if (dice->ctrl.setSavings < 0) {
	s2[tInd][iInd] = (float) dice->dvars.s[tInd];
      }
    }			
    
    // Loop over the climate sensitivity samples
    for(jInd=0; jInd<nSamples; jInd++) {
      
      // If this is the first BORG solution, put
      // the climate sensitivity in its vector
      if(iInd==0) {
				cs[jInd] = (float) dice->clim.cs_sample_vec[jInd];
				if(use_doeclim == 1) {
					oceandiff[jInd] = (float) dice->clim.ocean_diff_vec[jInd];
					alpha[jInd] = (float) dice->clim.alpha_vec[jInd];
				}
      }
      
      // Apply this climate sensitivity to the dice object
      if(use_doeclim == 1) {
	dice->clim.t2co = dice->clim.cs_sample_vec[jInd];
	dice->clim.kappa = dice->clim.ocean_diff_vec[jInd];
	dice->clim.alpha = dice->clim.alpha_vec[jInd];
	doeclim_DICE_init(dice);
      }
      else {
	dice->clim.t2xco2 = dice->clim.cs_sample_vec[jInd];
      }
            
      // Run CDICE with the SCC calculation
      calcSCC(dice, this_solution);
      //calc_CDICE(dice);
      
      // Capture the full time series of the model variables
      // in the appropriate vector
      utility[iInd][jInd] = (float) dice->econ.utility;
      for(tInd=0; tInd<nPeriods; tInd++) {
	
	if(use_doeclim == 0) {
	  tatm[tInd][iInd][jInd] = (float) dice->clim.tatm[tInd];
	  tocean[tInd][iInd][jInd] = (float) dice->clim.tocean[tInd];
	}
	mat[tInd][iInd][jInd] = (float) dice->carb.mat[tInd];
	forc[tInd][iInd][jInd] = (float) dice->carb.forc[tInd];
	mu[tInd][iInd][jInd] = (float) dice->carb.mu[tInd];
	ml[tInd][iInd][jInd] = (float) dice->carb.ml[tInd];
	cca[tInd][iInd][jInd] = (float) dice->econ.cca[tInd];
	k[tInd][iInd][jInd] = (float) dice->econ.k[tInd];
	c[tInd][iInd][jInd] = (float) dice->econ.c[tInd];
	cpc[tInd][iInd][jInd] = (float) dice->econ.cpc[tInd];
	i[tInd][iInd][jInd] = (float) dice->econ.i[tInd];
	ri[tInd][iInd][jInd] = (float) dice->econ.ri[tInd];
	y[tInd][iInd][jInd] = (float) dice->econ.y[tInd];
	ygross[tInd][iInd][jInd] = (float) dice->econ.ygross[tInd];
	ynet[tInd][iInd][jInd] = (float) dice->econ.ynet[tInd];
	damages[tInd][iInd][jInd] = (float) dice->econ.damages[tInd];
	damfrac[tInd][iInd][jInd] = (float) dice->econ.damfrac[tInd];
	abatecost[tInd][iInd][jInd] = (float) dice->econ.abatecost[tInd];
	periodu[tInd][iInd][jInd] = (float) dice->econ.periodu[tInd];
	eind[tInd][iInd][jInd] = (float) dice->econ.eind[tInd];
	pv_abatecost[tInd][iInd][jInd] = (float) dice->econ.pv_abatecost[tInd];
	pv_damages[tInd][iInd][jInd] = (float) dice->econ.pv_damages[tInd];
	totalcost[tInd][iInd][jInd] = (float) dice->econ.totalcost[tInd];
	pv_totalcost[tInd][iInd][jInd] = (float) dice->econ.pv_totalcost[tInd];
	e[tInd][iInd][jInd] = (float) dice->econ.e[tInd];
	mcabate[tInd][iInd][jInd] = (float) dice->econ.mcabate[tInd];
	cemutotper[tInd][iInd][jInd] = (float) dice->econ.cemutotper[tInd];
	scc[tInd][iInd][jInd] = (float) dice->econ.scc[tInd];
	
	
	if((iInd == 0) && (jInd == 0)) {
	  year[tInd] = (short) dice->config.dateSeries[tInd];
	  l[tInd] = (float) dice->econ.l[tInd];
	  ga[tInd] = (float) dice->econ.ga[tInd];
	  al[tInd] = (float) dice->econ.al[tInd];
	  gsig[tInd] = (float) dice->econ.gsig[tInd];
	  sigma[tInd] = (float) dice->econ.sigma[tInd];
	  cost1[tInd] = (float) dice->econ.cost1[tInd];
	  etree[tInd] = (float) dice->econ.etree[tInd];
	  rr[tInd] = (float) dice->econ.rr[tInd];
	  forcoth[tInd] = (float) dice->carb.forcoth[tInd];
	  partfract[tInd] = (float) dice->econ.partfract[tInd];
	  if (dice->ctrl.setSavings >= 0) {
	    s1[tInd] = (float) dice->dvars.s[tInd];
	  }
	}
	
	
      } // end loop over nPeriods
      
      // If DOEclim was used, store the variables with appropriate
      // dimensions.
      if(use_doeclim == 1) {
	for(tInd=0; tInd<nDoeclimts; tInd++) {
	  temp[tInd][iInd][jInd] = (float) dice->clim.temp[tInd];
	  doeclim_forc[tInd][iInd][jInd] = (float) dice->clim.forc[tInd];
	  heat_mixed[tInd][iInd][jInd] = (float) dice->clim.heat_mixed[tInd];
	}
      }
      
    } // end loop over nSamples
    
  } // end loop over nSolutions

  // Write these variable vectors to the netcdf file
  if(use_doeclim == 0) {
    if((retval = nc_put_var(ncid, tatmid, &tatm[0][0][0]))) { ncError(retval); }
    if((retval = nc_put_var(ncid, toceanid, &tocean[0][0][0]))) { ncError(retval); }
  }
  if((retval = nc_put_var(ncid, matid, &mat[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, forcid, &forc[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, muid, &mu[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, mlid, &ml[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, ccaid, &cca[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, kid, &k[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, cid, &c[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, cpcid, &cpc[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, iid, &i[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, riid, &ri[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, yid, &y[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, ygrossid, &ygross[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, ynetid, &ynet[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, damagesid, &damages[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, damfracid, &damfrac[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, abatecostid, &abatecost[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, perioduid, &periodu[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, eindid, &eind[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, pv_abatecostid, &pv_abatecost[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, pv_damagesid, &pv_damages[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, totalcostid, &totalcost[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, pv_totalcostid, &pv_totalcost[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, eid, &e[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, mcabateid, &mcabate[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, cemutotperid, &cemutotper[0][0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, sccid, &scc[0][0][0]))) { ncError(retval); }

  if(use_doeclim == 1) {
    if((retval = nc_put_var(ncid, tempid, &temp[0][0][0]))) { ncError(retval); }
    if((retval = nc_put_var(ncid, doeclim_forcid, &doeclim_forc[0][0][0]))) { ncError(retval); }
    if((retval = nc_put_var(ncid, heat_mixedid, &heat_mixed[0][0][0]))) { ncError(retval); }
    if((retval = nc_put_var(ncid, oceandiffid, &oceandiff[0]))) { ncError(retval); }
	  if((retval = nc_put_var(ncid, alphaid, &alpha[0]))) { ncError(retval); }
  }
  
  
  if((retval = nc_put_var(ncid, miuid, &miu[0][0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, utilityid, &utility[0][0]))) { ncError(retval); }
  
  if((retval = nc_put_var(ncid, csid, &cs[0]))) { ncError(retval); }
  
  if((retval = nc_put_var(ncid, yearid, &year[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, lid, &l[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, gaid, &ga[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, alid, &al[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, gsigid, &gsig[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, sigmaid, &sigma[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, cost1id, &cost1[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, etreeid, &etree[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, rrid, &rr[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, forcothid, &forcoth[0]))) { ncError(retval); }
  if((retval = nc_put_var(ncid, partfractid, &partfract[0]))) { ncError(retval); }
  
  if (dice->ctrl.setSavings < 0) {
    if((retval = nc_put_var(ncid, sid, &s2[0][0]))) { ncError(retval); }
  }
  else {
    if((retval = nc_put_var(ncid, sid, &s1[0]))) { ncError(retval); }
  }
  
  
  // Close the netcdf file
  if((retval = nc_close(ncid))) { ncError(retval); }
  
  return;
}
