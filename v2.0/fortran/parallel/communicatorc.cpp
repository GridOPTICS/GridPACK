/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   communicatorc.cpp
 * @author William A. Perkins
 * @date   2014-05-14 10:25:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May 12, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <gridpack/parallel/communicator.hpp>

extern "C" void
communicator_initialize(gridpack::parallel::Communicator **c)
{
  *c = new gridpack::parallel::Communicator;
}

extern "C" void
communicator_finalize(gridpack::parallel::Communicator **c)
{
  delete (*c);
  *c = NULL;
}

extern "C" int
communicator_size(const gridpack::parallel::Communicator *c)
{
  return c->size();
}

extern "C" int
communicator_rank(const gridpack::parallel::Communicator *c)
{
  return c->rank();
}

extern "C" int
communicator_world_rank(const gridpack::parallel::Communicator *c)
{
  return c->worldRank();
}

extern "C" void
communicator_barrier(const gridpack::parallel::Communicator *c)
{
  c->barrier();
}

extern "C" void
communicator_sync(const gridpack::parallel::Communicator *c)
{
  c->sync();
}

