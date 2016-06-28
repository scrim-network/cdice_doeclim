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


//===========================================================================================
// allocateConfig(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Configuration structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateConfig(DICE *dice)
{
    // Allocate the configuration structure----------------------------------------------------

  dice->config.nPeriods = 60;	// Number of time steps
  dice->config.tstep = 5.0;	// Number of years per time step
  dice->config.startYear = 2010;	// Model start year	
  dice->config.dateSeries = new double[dice->config.nPeriods];	// Human-readable year
	
  dice->config.distMean = 1.098001424;	// Mean of lognormal distribution of CS (Olson et al.)
  dice->config.distSD = 0.265206276;	// SD of lognormal distribution of CS (Olson et al.)
  lndist_type tempdist(dice->config.distMean, dice->config.distSD);
  dice->config.dist = tempdist;			// Distribution of CS
  rgen_type tempgen(time(NULL));
  dice->config.gen = tempgen;			// Random number generator for CS sampling

  return;
}


//===========================================================================================
// allocateCtrl(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Control structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateCtrl(DICE *dice)
{
  // Allocate the control variable structure-----------------------------------------------------
  // DEFAULTS
  dice->ctrl.writeModelDump = 0;			// Write model output to a file
  dice->ctrl.nobjs = 1;				// Number of objectives to calculate
  dice->ctrl.nvars = 2*dice->config.nPeriods-11;	// Number of input variables to read from user
  dice->ctrl.nconsts = 1;				// Number of constraints for optimization
  dice->ctrl.nevals = 10000;			// Number of function evaluations for optimization
  dice->ctrl.csSamples = 1;			// Number of samples to pull from CS distribution
  dice->ctrl.maxtime = 0.1;			// Number of hours to cap the optimization
  dice->ctrl.runtimeFreq = 0;			// Runtime dynamics output frequency (NFE)
  dice->ctrl.setSavings = -1.0;			// Set the savings rate (<0 for decision var)
  dice->ctrl.bau = 0;				// Run BAU scenario (default is no)
  dice->ctrl.use_doeclim = 0;			// Use the DOEClim model (1 = DOEClim, 0 = DICE)
  dice->ctrl.policy_inertia = -1.0;		// Policy Inertia ( <= 0.0 for no policy inertia)

  return;
}


//===========================================================================================
// allocateCarb(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Carbon structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateCarb(DICE *dice)
{
  
  dice->carb.mat0 = 830.4;		// Initial concentration in atmosphere 2010 [GtC]
  dice->carb.mu0 = 1527.0;		// Initial concentration in upper strata [GtC]
  dice->carb.ml0 = 10010.0;		// Initial concentration in lower strata [GtC]
  dice->carb.mateq = 588.0;		// Equilibrium concentration in atmosphere [GtC]
  dice->carb.mueq = 1350.0;		// Equilibrium concentration in upper strata [GtC]
  dice->carb.mleq = 10000.0;		// Equilibrium concentration in lower strata [GtC]
  
  dice->carb.b12 = 0.088;
  dice->carb.b23 = 0.00250;
  dice->carb.b11 = 1 - dice->carb.b12;
  dice->carb.b21 = dice->carb.b12 * dice->carb.mateq / dice->carb.mueq;
  dice->carb.b22 = 1 - dice->carb.b21 - dice->carb.b23;
  dice->carb.b32 = dice->carb.b23 * dice->carb.mueq / dice->carb.mleq;
  dice->carb.b33 = 1 - dice->carb.b32;

  dice->carb.fex0 = 0.25;	// 2010 forcings of non-CO2 greenhouse gases (GHG) [Wm-2]
  dice->carb.fex1 = 0.70;	// 2100 forcings of non-CO2 GHG [Wm-2]
  dice->carb.fco22x = 3.8;	// Forcings of equilibrium CO2 doubling [Wm-2]
  dice->carb.fosslim = 6000.0;	// Maximum cummulative extraction fossil fuels [GtC]

  dice->carb.forcoth = new double[dice->config.nPeriods];	// Exogenous forcing for other GHG
  dice->carb.forc = new double[dice->config.nPeriods];		// Increase in radiative forcing [Wm-2 from 1900]
  dice->carb.mat = new double[dice->config.nPeriods];		// Carbon concentration increase in atmosphere [GtC from 1750]
  dice->carb.mu = new double[dice->config.nPeriods];		// Carbon concentration increase in shallow oceans [GtC from 1750]
  dice->carb.ml = new double[dice->config.nPeriods];		// Carbon concentration increase in lower oceans [GtC from 1750]


  return;
}

//===========================================================================================
// allocateClim(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Climate structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateClim(DICE *dice)
{
  
  if(dice->ctrl.use_doeclim==1) {
    // DOEclim model variables
    dice->clim.t2co = 3.0;
    dice->clim.kappa = 0.55;
    dice->clim.alpha = 0.44;
    dice->clim.hind_ns = 110; // Number of hindcast timesteps (1900 - 2010)
    dice->clim.i1900 = 0;  // Index of year 1900
    dice->clim.dt = 1;
    dice->clim.ns = dice->clim.hind_ns+(dice->config.tstep*(dice->config.nPeriods-1)); // 1-yr res.
    dice->clim.forc = new double[dice->clim.ns];
    dice->clim.IB = new double[4];
    dice->clim.A = new double[4];
    dice->clim.Ker = new double[dice->clim.ns];
    dice->clim.n_calpoints = 0;
    
    dice->clim.temp = new double[dice->clim.ns];
    dice->clim.temp_landair = new double[dice->clim.ns];
    dice->clim.temp_sst = new double[dice->clim.ns];
    dice->clim.heat_mixed = new double[dice->clim.ns];
    dice->clim.heat_interior = new double[dice->clim.ns];
    dice->clim.heatflux_mixed = new double[dice->clim.ns];
    dice->clim.heatflux_interior = new double[dice->clim.ns];
    
   	dice->clim.ocean_diff_vec = new double [dice->ctrl.csSamples];
		dice->clim.alpha_vec = new double [dice->ctrl.csSamples];
    dice->clim.cs_cal = new double [dice->ctrl.csSamples];
    dice->clim.kv_cal = new double [dice->ctrl.csSamples];
    dice->clim.alpha_cal = new double [dice->ctrl.csSamples];
    
		doeclim_DICE_init(dice);

  }

  else {
    // DICE Climate model variables
    dice->clim.t2xco2 = 2.9;	// Equilibrium temperature impact [dC per doubling CO2]
    dice->clim.tocean0 = 0.0068;	// Initial lower stratum temperature change [dC from 1900]
    dice->clim.tatm0 = 0.80;	// Initial atmospheric temperature change [dC from 1900]
    dice->clim.c10 = 0.098;	// Initial climate equation coefficient for upper level
    dice->clim.c1beta = 0.01243;	// Regression slope coefficient (SoA~Equil TSC)
    dice->clim.c1 = 0.098;	// Climate equation coefficient for upper level
    dice->clim.c3 = 0.088;	// Transfer coefficient upper to lower stratum
    dice->clim.c4 = 0.025;	// Transfer coefficient for lower level
    dice->clim.lam = dice->carb.fco22x / dice->clim.t2xco2;	// Climate model parameter
    dice->clim.tatm = new double[dice->config.nPeriods];		// Increase temperature of atmosphere [dC from 1900]
    dice->clim.tocean = new double[dice->config.nPeriods];	// Increase temperature of lower oceans [dC from 1900]
    
  }

  dice->clim.cs_sample_vec = new double [dice->ctrl.csSamples];

  return;
}


//===========================================================================================
// allocateEcon(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Economics structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateEcon(DICE *dice)
{
  // Allocate the economic model structure---------------------------------------------

  dice->econ.elasmu = 1.45;		// Elasticity of marginal utility of consumption
  dice->econ.prstp = 0.015;		// Initial rate of social time preference (per year)

  dice->econ.gama = 0.300;		// Capital elasticity in production function
  dice->econ.pop0 = 6838.0;		// Initial world population [Millions]
  dice->econ.popadj = 0.134;		// Growth rate to calibrate to 2050 population projection
  dice->econ.popasym = 10500.0;	// Asymptotic world population [Millions]
  dice->econ.dk = 0.100;			// Depreciation rate on capital (per year)
  dice->econ.q0 = 63.69;			// Initial world gross output [Trillions 2005 US $]
  dice->econ.k0 = 135.0;			// Initial capital value [Trillions 2005 US $]
  dice->econ.a0 = 3.80;			// Initial level of total factor productivity (TFP)
  dice->econ.ga0 = 0.079;		// Initial growth rate for TFP (per 5 years)
  dice->econ.dela = 0.006;		// Decline rate of TFP (per 5 years)

  dice->econ.a10 = 0.0;			// Initial damage intercept
  dice->econ.a20 = 0.00267;		// Initial damage quadratic term
  dice->econ.a1 = 0.0;			// Damage intercept
  dice->econ.a2 = 0.00267;		// Damage quadratic term
  dice->econ.a3 = 2.00;			// Damage exponent

  dice->econ.expcost2 = 2.8;		// Exponent of control cost function
  dice->econ.pback = 344.0;		// Cost of backstop [2005$ per tCO2 2010]
  dice->econ.gback = 0.025;		// Initial cost decline backstop [cost per period]
  dice->econ.tnopol = 45;		// Period before which no emissions controls base
  dice->econ.cprice0 = 1.0;		// Initial base carbon price [2005$ per tCO2]
  dice->econ.gcprice = 0.02;		// Growth rate of base carbon price (per year)
		
  dice->econ.periodfullpart = 21;	// Period at which have full participation
  dice->econ.partfract2010 = 1.0;	// Fraction of emissions under control in 2010
  dice->econ.partfractfull = 1.0;	// Fraction of emissions under control at full time

  dice->econ.gsigma1 = -0.01;		// Initial growth of sigma (per year)
  dice->econ.dsig = -0.001;		// Decline rate of decarbonization (per period)
//dice->econ.sig0 = dice->carb.e0 / (dice->econ.q0 * (1 - dice->econ.miu0));
						// Carbon Intensity 2010 [kgCO2 per output 2005 USD 2010]
  dice->econ.eland0 = 3.3;		// Carbon emissions from land 2010 [GtCO2 per year]
  dice->econ.deland = 0.2;		// Decline rate of land emissions (per period)
  dice->econ.e0 = 33.61;		// Industrial emissions 2010 [GtCO2 per year]

  dice->econ.scale1 = 0.016408662;	// Multiplicitive scaling coefficient
  dice->econ.scale2 = -3855.106895;     // Additive scaling coefficient

  dice->econ.utility = 0.0;		// Welfare function (Sum of discounted utility of per capita consumption)

  dice->econ.l = new double[dice->config.nPeriods];		// Level of population and labor (millions)
  dice->econ.al = new double[dice->config.nPeriods];		// Level of total factor productivity
  dice->econ.rr  = new double [dice->config.nPeriods];		// Average utility social discount rate
  dice->econ.ga = new double[dice->config.nPeriods];		// Growth rate of productivity from (sic)
  dice->econ.gl = new double[dice->config.nPeriods];		// Growth rate of labor
  dice->econ.gcost = new double[dice->config.nPeriods];		// Growth of cost factor
  dice->econ.cost1 = new double[dice->config.nPeriods];		// Adjusted cost for backup
  dice->econ.partfract = new double[dice->config.nPeriods];	// Fraction of emissions in control regime
  //dice->econ.gfacpop = new double[dice->config.nPeriods];	// Growth factor population (Appears in GAMS version, but never used
  dice->econ.pbacktime = new double[dice->config.nPeriods];	// Backstop price
  dice->econ.scc = new double[dice->config.nPeriods];		// Social cost of carbon
  dice->econ.cpricebase = new double[dice->config.nPeriods];	// Carbon price in base case
  dice->econ.photel = new double[dice->config.nPeriods];	// Carbon price under no damages (Hotelling rent condition)
  dice->econ.gsig = new double[dice->config.nPeriods];		// Change in sigma (cumulative improvement of energy efficiency)
  dice->econ.sigma = new double[dice->config.nPeriods];		// CO2-equivalent-emissions output ratio 
  dice->econ.etree = new double[dice->config.nPeriods];		// Emissions from deforestation

  dice->econ.c = new double[dice->config.nPeriods];		// Consumption [Trillions 2005 US$ per year]
  dice->econ.k = new double[dice->config.nPeriods];		// Capital stock [Trillions 2005 US$]
  dice->econ.cpc = new double[dice->config.nPeriods];		// Per capita consumption [Thousands 2005 US$ per year]
  dice->econ.i = new double[dice->config.nPeriods];		// Investment [trillions 2005 US$ per year]
  dice->econ.ri = new double[dice->config.nPeriods];		// Real interest rate (per annum)
  dice->econ.y = new double[dice->config.nPeriods];		// Gross world product net of abatement and damages [Trillions 2005 US$ per year]
  dice->econ.ygross = new double[dice->config.nPeriods];	// Gross world product GROSS of abatement and damages [Trillions 2005 US$ per year]
  dice->econ.ynet = new double[dice->config.nPeriods];		// Output net of damages equation [Trillions of 2005 US$ per year]
  dice->econ.damages = new double[dice->config.nPeriods];	// Damages [Trillions 2005 US$ per year]
  dice->econ.damfrac = new double[dice->config.nPeriods];	// Damages as fraction of gross output
  dice->econ.abatecost = new double[dice->config.nPeriods];	// Cost of emissions reductions [Trillions 2005 US$ per year]
  dice->econ.mcabate = new double[dice->config.nPeriods];	// Marginal cost of abatement [2005 US$ per ton CO2]
  dice->econ.periodu = new double[dice->config.nPeriods];	// One period utility function
  dice->econ.cprice = new double[dice->config.nPeriods];	// Carbon price [2005 US$ per ton CO2]
  dice->econ.cemutotper = new double[dice->config.nPeriods];	// Period utility
  dice->econ.ri_disc = new double[dice->config.nPeriods];	// Real interest rate (discounted) for present value calculations
  dice->econ.pv_damages = new double[dice->config.nPeriods];	// Present value of damages
  dice->econ.pv_abatecost = new double[dice->config.nPeriods];	// Present value of abatement costs
  dice->econ.totalcost = new double[dice->config.nPeriods];	// Total costs (abatement + damages)
  dice->econ.pv_totalcost = new double[dice->config.nPeriods];	// Present value of total costs (abatement + damages)
  dice->econ.e = new double[dice->config.nPeriods];		// Total CO2 emissions [GtCO2 per year]
  dice->econ.eind = new double[dice->config.nPeriods];		// Industrial emissions [GtCO2 per year]
  dice->econ.cca = new double[dice->config.nPeriods];		// Cumulative industrial carbon emissions [GtC]


  return;
}

//===========================================================================================
// allocateDvars(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Decision Variable structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateDvars(DICE *dice)
{

  dice->dvars.limmiu = 1.2;		// Upper limit on control rate after 2150
  dice->dvars.miu0 = 0.039;		// Initial emissions control rate for base case 2010

  dice->dvars.miu = new double[dice->config.nPeriods];	// Emission control rate GHGs
  dice->dvars.s = new double[dice->config.nPeriods];	// Gross savings rate as function of gross world production
  

  return;
}


//===========================================================================================
// allocateLims(DICE *dice) 
//
// Allocates space and initializes the bounds in the Lims structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateLims(DICE *dice)
{
  dice->lims.optlrsav = (dice->econ.dk + 0.004) / (dice->econ.dk + 0.004 * dice->econ.elasmu + dice->econ.prstp) * dice->econ.gama;
				// Optimal long-run savings rate used for transversality
  dice->lims.miu_up = new double[dice->config.nPeriods];	// Upper limit of control policy for each time period
  dice->lims.k_lo = 1.0;
  dice->lims.mat_lo = 10.0;
  dice->lims.mu_lo = 100.0;
  dice->lims.ml_lo = 1000.0;
  dice->lims.c_lo = 2.0;
  dice->lims.tocean_up = 20.0;
  dice->lims.tocean_lo = -1.0;
  dice->lims.tatm_lo = 0.0;
  dice->lims.tatm_up = 40.0;
  dice->lims.cpc_lo = 0.01;
  dice->lims.y_lo = 0.0;
  dice->lims.ygross_lo = 0.0;
  dice->lims.i_lo = 0.0;
  dice->lims.cca_up = dice->carb.fosslim;
}


//===========================================================================================
// initTS(DICE *dice) 
//
// Initializes the time series variables in the DICE structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void initTS(DICE *dice)
{

  for (int i=0; i<dice->config.nPeriods; i++) {
    
    dice->config.dateSeries[i] = dice->config.startYear + (dice->config.tstep * i);
    dice->econ.l[i] = 0.0;
    dice->econ.al[i] = 0.0;		
    dice->econ.sigma[i] = 0.0;	 
    dice->econ.rr[i] = 0.0;		
    dice->econ.ga[i] = 0.0;		
    dice->carb.forcoth[i] = 0.0;
    dice->econ.gl[i] = 0.0;		
    dice->econ.gcost[i] = 0.0;
    dice->econ.gsig[i] = 0.0;		
    dice->econ.etree[i] = 0.0;
    dice->econ.cost1[i] = 0.0;
    dice->econ.partfract[i] = 0.0;
    //dice->econ.gfacpop[i] = 0.0;	// Appears in GAMS version, but is never used
    dice->econ.pbacktime[i] = 0.0;
    dice->econ.scc[i] = 0.0;
    dice->econ.cpricebase[i] = 0.0;
    dice->econ.photel[i] = 0.0;

    dice->dvars.miu[i] = dice->dvars.miu0;
    dice->carb.forc[i] = 0.0;
    dice->carb.mat[i] = 0.0;
    dice->carb.mu[i] = 0.0;
    dice->carb.ml[i] = 0.0;
    dice->econ.e[i] = 0.0;
    dice->econ.eind[i] = 0.0;
    dice->econ.c[i] = 0.0;
    dice->econ.k[i] = 0.0;
    dice->econ.cpc[i] = 0.0;
    dice->econ.i[i] = 0.0;
    dice->econ.ri[i] = 0.0;
    dice->econ.y[i] = 0.0;
    dice->econ.ygross[i] = 0.0;
    dice->econ.ynet[i] = 0.0;
    dice->econ.damages[i] = 0.0;
    dice->econ.damfrac[i] = 0.0;
    dice->econ.abatecost[i] = 0.0;
    dice->econ.mcabate[i] = 0.0;
    dice->econ.cca[i] = 0.0;
    dice->econ.periodu[i] = 0.0;
    dice->econ.cprice[i] = 0.0;
    dice->econ.cemutotper[i] = 0.0;

    dice->econ.ri_disc[i] = 0.0;
    dice->econ.pv_damages[i] = 0.0;
    dice->econ.pv_abatecost[i] = 0.0;
    dice->econ.totalcost[i] = 0.0;
    dice->econ.pv_totalcost[i] = 0.0;

    // If not using DOEclim, then initialize DICE clim time series
    if(dice->ctrl.use_doeclim == 0) {
      dice->clim.tatm[i] = 0.0;
      dice->clim.tocean[i] = 0.0;
    }
		
    // Set the last 10 periods of the savings rate to optlrsav
    if(i < dice->config.nPeriods-10) {
      dice->dvars.s[i] = 0.0;
    }
    else {
      dice->dvars.s[i] = dice->lims.optlrsav;
    }
    
    // Override the savings rate if prescribed from the control file
    if(dice->ctrl.setSavings > 0.0) {
    	dice->dvars.s[i] = dice->ctrl.setSavings;
    }
    
  }
	
	// If using DOEclim, then initialize the appropriate time series
	if(dice->ctrl.use_doeclim == 1) {
		for(int i=0; i<dice->clim.ns; i++) {
  	  dice->clim.temp[i] = 0.0;
    	dice->clim.temp_landair[i] = 0.0;
	    dice->clim.temp_sst[i] = 0.0;
  	  dice->clim.heat_mixed[i] = 0.0;
    	dice->clim.heat_interior[i] = 0.0;
	    dice->clim.heatflux_mixed[i] = 0.0;
 	   dice->clim.heatflux_interior[i] = 0.0;
  	}
  }
	
  // Set the first control policy timestep to the initial value
  dice->dvars.miu[0] = dice->dvars.miu0;
	

  // Calculate values for the exogenous variables -------------------------
  // Information from the exogenous variables is used
  // to help define some of the bounds initialized later
  // in this function.
  calcExog(dice);

	
  // Time series
  for(int i=0; i<dice->config.nPeriods; i++) {

    // Control rate limit
    if (i < 29) {
      dice->lims.miu_up[i] = 1.0;
    }
    else	{
      dice->lims.miu_up[i] = dice->dvars.limmiu * dice->econ.partfract[i];
    }
    
  }


  // Initialize the vectors that contain values for the stochastic runs
  for (int i=0; i < dice->ctrl.csSamples; i++) {
    dice->clim.cs_sample_vec[i] = 0.0;
    if(dice->ctrl.use_doeclim == 1) {
    	dice->clim.ocean_diff_vec[i] = 0.0;
    	dice->clim.alpha_vec[i] = 0.0;
    	dice->clim.cs_cal[i] = 0.0;
    	dice->clim.kv_cal[i] = 0.0;
    	dice->clim.alpha_cal[i] = 0.0;
    }

  }
  
  return;

}


//=================================================================================
//
// void ctrl_CDICE(DICE *dice, string controlFileName)
//
// Reads in and applies values from dice control file.
//
// Return: Update DICE object
//=================================================================================
void ctrl_CDICE(DICE *dice, string controlFileName)
{
	
  // Temporary variable to hold info from control file
  char *junk;
  junk = new char[1000];
  
  //cout << "init_CDICE: reading control file\n";   //debug
  //cin >> dummychar;
  
  // Open the control file for reading
  ifstream control_stream;
  control_stream.open(controlFileName.c_str());
  if (!control_stream.is_open()) {
    cout << "Failed to open control file. Exiting ...";
    exit(1);
  }

  // Read a line, parsing for the specific line content, discard rest of line
  long unsigned int tempint;
  double tempdouble;
  while (!control_stream.eof()) {
    control_stream >> junk;
    
    // Should we write the model output to a file?
    if (!strcmp(junk,"<model_dump>")) 	{
      control_stream >> tempint; 
      control_stream.ignore(2000,'\n');
      if (tempint == 0 || tempint == 1) {
	dice->ctrl.writeModelDump = (int) tempint;
      }
      else	{
	cout << "model_dump must be 0 or 1. Exiting ...";
	exit(1);
      }
    }
    
    // What's the start year (should be 2010)
    else if (!strcmp(junk,"<start_year>")) {
      control_stream >> tempint;
      control_stream.ignore(2000,'\n');
      if (tempint > 1985 && tempint < 2305) {
	dice->config.startYear = (int) tempint;
      }
      else	{
	cout << "start year must be reasonable. Exiting ...";
	exit(1);
      }
    }
    
    // How many objectives are we calculating? (Must be at least one)
    else if (!strcmp(junk,"<nobjs>")) {
      control_stream >> tempint;
      control_stream.ignore(2000,'\n');
      if (tempint > 0) {
	dice->ctrl.nobjs = (int) tempint;
      }
      else	{
	cout << "Number of objectives must be greater than 0. Exiting ...";
	exit(1);
      }
    }
    
    // How many constraints are there? (One by default...fossil fuel limit)
    else if (!strcmp(junk,"<nconsts>")) {
      control_stream >> tempint;
      control_stream.ignore(2000,'\n');
      if (tempint > 0) {
	dice->ctrl.nconsts = (int) tempint;
      }
      else	{
	cout << "Number of constraints must be greater than 0. Exiting ...";
	exit(1);
      }
    }
    
    // How many function evaluations are needed to optimize? (Should be more than 0, 10,000+ is reasonable)
    else if (!strcmp(junk,"<nevals>")) {
      control_stream >> tempint;
      control_stream.ignore(2000,'\n');
      if (tempint >= 0) {
	dice->ctrl.nevals = tempint;
      }
      else	{
	cout << "Number of function evaluations must be >= 0. Exiting ...";
	exit(1);
      }
    }
    
    // How many samples should we pull from the CS distribution?
    // Note: A value of 1 will use the default settings in cdice
    else if (!strcmp(junk,"<csSamples>")) {
      control_stream >> tempint;
      control_stream.ignore(2000,'\n');
      if (tempint > 0) {
	dice->ctrl.csSamples = tempint;
      }
      else	{
	cout << "Number of csSamples must be greater than 0. Exiting ...";
	exit(1);
      }
    }
    
    // How much time (in hours) is needed to optimize? (Should be more than 0)
    else if (!strcmp(junk,"<maxtime>")) {
      control_stream >> tempdouble;
      control_stream.ignore(2000,'\n');
      if (tempdouble >= 0.0) {
	dice->ctrl.maxtime = tempdouble;
      }
      else	{
	cout << "Number of hours for optimization must be >= 0. Exiting ...";
	exit(1);
      }
    }
    
    // How often should the runtime dynamics be output? (NFE per output)
    else if (!strcmp(junk,"<runtimeFreq>")) {
      control_stream >> tempint;
      control_stream.ignore(2000,'\n');
      if (tempint >= 0) {
	dice->ctrl.runtimeFreq = tempint;
      }
      else	{
	cout << "NFE for runtime dynamics must be 0 or positive. Exiting ...";
	exit(1);
      }
    }
    
    // Set the savings rate if necessary (set < 0 to use savings rate as decision variable)
    else if (!strcmp(junk,"<setSavings>")) {
      control_stream >> tempdouble;
      control_stream.ignore(2000,'\n');
      if (tempdouble <= 1.0) {
	dice->ctrl.setSavings = tempdouble;
      }
      else	{
	cout << "Savings rate cannot be greater than 1.0. Exiting ...";
	exit(1);
      }
    }
    
    // Should we do a BAU optimization run?
    else if (!strcmp(junk,"<bau>")) 	{
      control_stream >> tempint; 
      control_stream.ignore(2000,'\n');
      if (tempint == 0 || tempint == 1) {
	dice->ctrl.bau = (int) tempint;
      }
      else	{
	cout << "BAU must be 0 or 1. Exiting ...";
	exit(1);
      }
    }
    

    // Should we use the DOECLIM climate model?
    else if (!strcmp(junk,"<use_doeclim>")) 	{
      control_stream >> tempint; 
      control_stream.ignore(2000,'\n');
      if (tempint == 0 || tempint == 1) {
	dice->ctrl.use_doeclim = (int) tempint;
      }
      else	{
	cout << "use_doeclim must be 0 or 1. Exiting ...";
	exit(1);
      }
    }
    
    // Set the policy inertia (set < 0 lets BORG pass control policy w/o inertia)
    else if (!strcmp(junk,"<policyInertia>")) {
      control_stream >> tempdouble;
      control_stream.ignore(2000,'\n');
      if (tempdouble <= 1.0) {
	dice->ctrl.policy_inertia = tempdouble;
      }
      else	{
	cout << "Policy Inertia cannot be greater than 1.0. Exiting ...";
	exit(1);
      }
    }
    
    
  }
  
  control_stream.close();

  // If the savings rate is set from the control file, then
  // re-initialize the savings rate and reset the number
  // of decision variables
  if (dice->ctrl.setSavings >= 0.0) {
    dice->ctrl.nvars = dice->config.nPeriods - 1;
    for (int i=0; i < dice->config.nPeriods; i++) {
      dice->dvars.s[i] = dice->ctrl.setSavings;
    }
  }
  


  delete junk;
  

  return;
}


//=================================================================================
//
// void free_CDICE(DICE *dice)
//
// Frees up memory by destroying a DICE object
//
//=================================================================================
void free_CDICE(DICE *dice)
{
	// Delete the allocated arrays
	delete dice->config.dateSeries;	
	delete dice->carb.forcoth;	
  delete dice->carb.forc;		
  delete dice->carb.mat;		
  delete dice->carb.mu;		
  delete dice->carb.ml;		
  delete dice->clim.cs_sample_vec;
  delete dice->econ.l;		
  delete dice->econ.al;		
  delete dice->econ.rr ;		
  delete dice->econ.ga;		
  delete dice->econ.gl;		
  delete dice->econ.gcost;		
  delete dice->econ.cost1;		
  delete dice->econ.partfract;	
  delete dice->econ.pbacktime;	
  delete dice->econ.scc;		
  delete dice->econ.cpricebase;	
  delete dice->econ.photel;	
  delete dice->econ.gsig;		
  delete dice->econ.sigma;		
  delete dice->econ.etree;		
  delete dice->econ.c;		
  delete dice->econ.k;		
  delete dice->econ.cpc;		
  delete dice->econ.i;		
  delete dice->econ.ri;		
  delete dice->econ.y;		
  delete dice->econ.ygross;	
  delete dice->econ.ynet;		
  delete dice->econ.damages;	
  delete dice->econ.damfrac;	
  delete dice->econ.abatecost;	
  delete dice->econ.mcabate;	
  delete dice->econ.periodu;	
  delete dice->econ.cprice;	
  delete dice->econ.cemutotper;	
  delete dice->econ.ri_disc;	
  delete dice->econ.pv_damages;	
  delete dice->econ.pv_abatecost;	
  delete dice->econ.totalcost;	
  delete dice->econ.pv_totalcost;	
  delete dice->econ.e;		
  delete dice->econ.eind;		
  delete dice->econ.cca;		
  delete dice->dvars.miu;	
  delete dice->dvars.s;	
  delete dice->lims.miu_up;	

	if(dice->ctrl.use_doeclim == 1) {	
		delete dice->clim.forc;
    delete dice->clim.IB;
    delete dice->clim.A;
    delete dice->clim.Ker;
    delete dice->clim.temp;
    delete dice->clim.temp_landair;
    delete dice->clim.temp_sst;
    delete dice->clim.heat_mixed;
    delete dice->clim.heat_interior;
    delete dice->clim.heatflux_mixed;
    delete dice->clim.heatflux_interior;
   	delete dice->clim.ocean_diff_vec;
		delete dice->clim.alpha_vec;
    delete dice->clim.cs_cal;
    delete dice->clim.kv_cal;
    delete dice->clim.alpha_cal;
  } else {
    delete dice->clim.tatm;		
    delete dice->clim.tocean;	
  }
  
/*  // Delete the structures in the DICE object
  delete dice->config;
  delete dice->dvars;
  delete dice->ctrl;
  delete dice->clim;
  delete dice->econ;
  delete dice->carb;
  delete dice->lims;
*/  
  // Delete the parent structure of the DICE object
//  delete dice;
  
}

//#endif  //CDICEMOD
