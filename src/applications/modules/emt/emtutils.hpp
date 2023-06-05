/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emt.hpp
 *
 * @brief  Header file for the dynamic simulation application utilities
 *    -- EmtProfiler -- Code Profiler
 *    -- EmtParams   -- Simulation Parameters
 *    -- EmtConfig   -- Configuration Parameters
 *
 *
 */
// -------------------------------------------------------------

#ifndef _emtutils_h_
#define _emtutils_h_

#include <gridpack/include/gridpack.hpp>

/**
   Class for code profiling
*/
class EmtProfiler
{ 
  public:
  
  EmtProfiler() 
  {
    p_timer = gridpack::utility::CoarseTimer::instance();
    // Set up timers
    p_ttotal = p_timer->createCategory("Dynamic Simulation: Total Application");
    p_tsetup = p_timer->createCategory("Dynamic Simulation: Set up");
    p_tsolve = p_timer->createCategory("Dynamic Simulation: Solve");
    p_tpostprocess = p_timer->createCategory("Dynamic Simulation: Post Processing");

    p_timer->start(p_ttotal);
  }

  void startdatareadtimer(void) 
  {
    p_treaddata = p_timer->createCategory("Dynamic Simulation: Data read");
    p_timer->start(p_treaddata);
  }

  void stopdatareadtimer(void) 
  {
    p_timer->stop(p_treaddata);
  }

  void startsetuptimer(void) 
  {
    p_tsetup = p_timer->createCategory("Dynamic Simulation: Set Up");
    p_timer->start(p_tsetup);
  }

  void stopsetuptimer(void) 
  {
    p_timer->stop(p_tsetup);
  }

  void startsolvetimer(void) 
  {
    p_tsolve = p_timer->createCategory("Dynamic Simulation: Solve");
    p_timer->start(p_tsolve);
  }

  void stopsolvetimer(void) 
  {
    p_timer->stop(p_tsolve);
  }


  ~EmtProfiler()
  {
    p_timer->stop(p_ttotal);
    p_timer->dump();
  }

  private:
  
  gridpack::utility::CoarseTimer *p_timer;

  // Some basic code profiling categories
  int p_ttotal;    // Total Application Time
  int p_treaddata; // Read data
  int p_tsetup;   // Set Up time (after reading data)
  int p_tsolve;   // Solve time
  int p_tpostprocess;  // Post processing time
};

/**
   Class for holding simulation parameters
*/
class EmtParams
{
  public:
  /** 
      Basic constructor
  */
  EmtParams(void) {};

  /** 
      Basic constructor with simulation time and initial time-step
  */
  EmtParams(double dt0, double final_time)
  {
    p_dt0 = dt0;
    p_final_time = final_time;
  }
  private:
  double p_dt0; // Initial time-step
  double p_final_time; // Final time
};

#endif
