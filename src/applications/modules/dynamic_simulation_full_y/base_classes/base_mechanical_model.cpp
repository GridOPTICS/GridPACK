/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_mechanical_model.cpp
 * @author Shuangshuang Jin
 * @Created on:   Nov 15, 2022
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BaseMechanicalModel::BaseMechanicalModel(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BaseMechanicalModel::~BaseMechanicalModel(void)
{
}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * BaseMechanicalModel
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
}

/**
 * Initialize mechanical model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::init(
    double mag, double ang, double ts)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::corrector(
    double t_inc, bool flag)
{
}

