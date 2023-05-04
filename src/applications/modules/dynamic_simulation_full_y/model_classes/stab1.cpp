/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   stab1.cpp
 * @author Shuangshuang Jin
 * @Created on: Apr 03, 2023
 * @Last modified:   Apr 20, 2023
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_pss_model.hpp"
#include "stab1.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Stab1Model::Stab1Model(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Stab1Model::~Stab1Model(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * Stab1Model
 */
void gridpack::dynamic_simulation::Stab1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(STAB1_J, &j, idx)) j = 1.0; 
  if (!data->getValue(STAB1_J1, &j1, idx)) j1 = 1.0; 
  if (!data->getValue(STAB1_J2, &j2, idx)) j2 = 1.0; 
  if (!data->getValue(STAB1_J3, &j3, idx)) j3 = 1.0; 
  if (!data->getValue(STAB1_J4, &j4, idx)) j4 = 1.0; 
  if (!data->getValue(STAB1_J5, &j5, idx)) j5 = 1.0; 
  if (!data->getValue(STAB1_J6, &j6, idx)) j6 = 1.0; 

    T = j1;
    K = T * j; // J: K/T
    T3 = j3;
    T1 = T3 * j2; // J+2: T1/T3
    T4 = j5;
    T2 = T4 * j4; // J+4: T2/T4
    Hlim = j6;
    
  // Set up blocks
    // Feedback block
    double a[2], b[2];
    a[0] = T; a[1] = 1.0;
    b[0] = K;  b[1] = 0.0;
    Feedback_blk.setcoeffs(a,b);
    
    // Lead lag blocks
    Leadlag_blk1.setparams(T1, T3);
    Leadlag_blk2.setparams(T2, T4);
    
    // Gain limit blocks
    Hlim_blk.setparams(1.0, -Hlim, Hlim); 

    //printf("%f %f %f %f %f %f %f\n", j, j1, j2, j3, j4, j5, j6);
  
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Stab1Model::init(double mag, double ang, double ts)
{
    double u1, u2, u3, u4;
    u1 = vothsg;
    u2 = Leadlag_blk2.init_given_y(u1);
    u3 = Leadlag_blk1.init_given_y(u2);
    u4 = Feedback_blk.init_given_y(u3);
    speed = u4;
     
    /* Alternatively
    u4 = Feedback_blk.init_given_u(speed);
    u3 = Leadlag_blk1.init_given_u(u4);
    u2 = Leadlag_blk2.init_given_u(u3);
    vothsg = u2;*/
}

void gridpack::dynamic_simulation::Stab1Model::computeModel(double t_inc, IntegrationStage int_flag)
{
    double u1, y1, u2, y2, u3, y3, u4, y4;
    u1 = speed;
    y1 = Feedback_blk.getoutput(u1, t_inc, int_flag, true);
    u2 = y1;
    y2 = Leadlag_blk1.getoutput(u2, t_inc, int_flag, true);
    u3 = y2;
    y3 = Leadlag_blk2.getoutput(u3, t_inc, int_flag, true);
    u4 = y3;
    y4 = Hlim_blk.getoutput(u4);
    vothsg = y4;
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Stab1Model::predictor(double t_inc, bool flag)
{
    computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Stab1Model::corrector(double t_inc, bool flag)
{
    computeModel(t_inc,CORRECTOR);
}

double gridpack::dynamic_simulation::Stab1Model::getVothsg( )
{
		return vothsg;
}

void gridpack::dynamic_simulation::Stab1Model::setSpeed(double input_speed)
{
   speed = input_speed;
}


