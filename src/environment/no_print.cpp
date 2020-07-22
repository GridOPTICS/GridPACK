/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include <stdio.h>
#include "gridpack/environment/no_print.hpp"

gridpack::NoPrint
         *gridpack::NoPrint::p_instance = NULL;

/**
 * Retrieve instance of the NoPrint object
 */
gridpack::NoPrint *gridpack::NoPrint::instance()
{
  if (p_instance == NULL) {
    p_instance = new NoPrint();
  }
  return p_instance;
}

/**
 * set status of print object
 * @param flag status of print object
 */
void gridpack::NoPrint::setStatus(bool flag)
{
  p_status = flag;
}

/**
 * return status of NoPrint object
 * @return true if no print desired
 */
bool gridpack::NoPrint::status()
{
  return p_status;
}

/**
 * Constructor
 */
gridpack::NoPrint::NoPrint()
{
  p_status = false;
}

/**
 * Destructor
 */
gridpack::NoPrint::~NoPrint()
{
}
