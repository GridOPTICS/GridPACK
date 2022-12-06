/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wtdta1_cblock.cpp
 * @author Shuangshuang Jin 
 * @Created on:    Dec 05, 2022
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <string>
#include <cstdio>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_mechanical_model.hpp"
#include "wtdta1_cblock.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wtdta1Model::Wtdta1Model(void)
{
    Pmech = 0.0;
    Pgen = 0.0;
    wt = 1.0;
    //dx0 = 0;
    //dx1 = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Wtdta1Model::~Wtdta1Model(void)
{
}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * Wtdta1Model
 */
void gridpack::dynamic_simulation::Wtdta1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
    if (!data->getValue(WIND_DT_H, &Ht, idx)) Ht = 0.0; // Ht
    if (!data->getValue(WIND_DT_DAMP,&damp,idx)) damp = 0.0;
    if (!data->getValue(WIND_DT_HFRAC,&hfrac,idx)) hfrac = 0.0;
    if (!data->getValue(WIND_DT_FREQ1,&freq1,idx)) freq1 = 0.0;
    if (!data->getValue(WIND_DT_DSHAFT,&dshaft,idx)) dshaft = 0.0;

    // Using 1-Mass model since 2-Mass model parameters are not provided
    itblock1.setparams(2*Ht);
    itblock2.setparams(1); 
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wtdta1Model::init(double mag, double ang, double ts)
{
    double pi = 4.0*atan(1.0);
    w0 = 2*pi*60.0;//1.0;
    //x0 = wt - w0; // x0 = dOut - WTGT.w0; dOut = wt;
    
    double u1, u2; 
    u1 = itblock2.init_given_y(wt - w0);
    u2 = itblock1.init_given_y(u1); 
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtdta1Model::predictor(double t_inc, bool flag)
{
    /*// Step-1: update predictor state variables using corrector
    if (!flag) {
        x0 = x1;
    }
    // Step-2: calculate predictor dx0 (dx/dt)
    dx0 = (Pmech - Pgen)/((x0 + w0)*2*Ht); // In = Pmech - Pgen;
    // Step-3: integrate
    x1 = x0 + dx0 * t_inc;
    // Step-4: update outputs
    wt = x1 + w0;*/

    double u1, y1;
    u1 = (Pmech - Pgen - dshaft) / wt;
    y1 = itblock1.getoutput(u1, t_inc, PREDICTOR, true);
    double u2, y2;
    u2 = y1;
    y2 = itblock2.getoutput(u2, t_inc, PREDICTOR, true);
    wt = y2 + w0;
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtdta1Model::corrector(double t_inc, bool flag)
{
    /*// Step-1: calculate corrector dx'/dt
    dx1 = (Pmech - Pgen)/((x0 + w0)*2*Ht);
    
    // Step-2: integrate
    x1 = x0 + (dx0 + dx1) / 2.0 * t_inc;

    // Step-3: Update output
    wt = x1 + w0;*/

    double u1, y1;
    u1 = (Pmech - Pgen - dshaft) / wt;
    y1 = itblock1.getoutput(u1, t_inc, CORRECTOR, true);
    double u2, y2;
    u2 = y1;
    y2 = itblock2.getoutput(u2, t_inc, CORRECTOR, true);
    wt = y2 + w0;
}

/**
 * Set the Pgen parameter inside the Drive-Train model
 * @param Pgen value of the electrical power
 */
void  gridpack::dynamic_simulation::Wtdta1Model::setPgen(double Pgen1)
{
    Pgen = Pgen1;
}

/**
 * Set the Pmech parameter inside the Drive-Train model
 * @param Pmech value of the mechanical power
 */
void  gridpack::dynamic_simulation::Wtdta1Model::setPmech(double Pmech1)
{
    Pmech = Pmech1;
}

/**
 * Get the value of the rotational speed
 * @return value of rotational speed
 */
double  gridpack::dynamic_simulation::Wtdta1Model::getRotationalSpeed()
{
  return wt;
}
