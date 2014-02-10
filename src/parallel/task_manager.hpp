// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   task_manager.hpp
 * @author Bruce Palmer
 * @date   February 10, 2014
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _task_manager_hpp_
#define _task_manager_hpp_

#include "gridpack/parallel/communicator"
#include <ga.h>

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class TaskManager
// -------------------------------------------------------------
class TaskManager {
public:

  // Constructor on world communicator
  TaskManager(void)
  {
    int grp = GA_Pgroup_get_world();

    p_GAcounter = GA_Create_handle();
    int one = 1;
    GA_Set_data(p_GAcounter,one,&one,C_INT);
    GA_Set_pgroup(p_GAcounter,grp);
    if (!GA_Allocate(p_GAcounter)) {
      // TODO: some kind of error
    }
    GA_zero(p_GAcounter);
    p_ntasks = 0;
  }

  // Constructor that uses user-specified communicator
  // @param comm user-specified communicator
  TaskManager(Communicator &comm)
  {
    int grp = comm.getGroup();

    p_GAcounter = GA_Create_handle();
    int one = 1;
    GA_Set_data(p_GAcounter,one,&one,C_INT);
    GA_Set_pgroup(p_GAcounter,grp);
    if (!GA_Allocate(p_GAcounter)) {
      // TODO: some kind of error
    }
    GA_zero(p_GAcounter);
    p_ntasks = 0;
  }

  // Destructor
  ~TaskManager(void)
  {
    GA_Destroy(p_GAcounter);
  }

  // Specify total number of tasks and set task manager to zero
  // @param ntasks total number of tasks
  void set(int ntasks)
  {
    GA_Zero(p_GAcounter);
    p_ntasks = ntasks;
  }
  
  // Get the next task from the task manager. If the manager finds a task it
  // returns true and next is set to the index of the task, otherwise it returns
  // false and next is set to -1
  // @param next index of next task
  // @return false if no other tasks are found
  bool nextTask(int *next) {
    int zero = 0;
    long one = 1;
    *next = static_cast<int>(GA_Read_inc(p_GAcounter,&zero,one));
    if (*next < p_ntasks) {
      return true;
    } else {
      *next = -1;
      return false;
    }
  }

  // Get the next task for the whole communicator. The same value of next is
  // returned for all processors in the communicator comm. If the manager finds
  // a task it returns true and next is set to the index of the task, otherwise
  // it returns false and next is set to -1

  bool nextTask(Communicator &comm, int *next) {
    int zero = 0;
    long one = 1;
    int me = comm.rank();
    if (me == 0) {
      *next = static_cast<int>(GA_Read_inc(p_GAcounter,&zero,one));
    } else {
      *next = 0;
    }
    GA_Pgroup_igop(comm.getGroup(),next,one,"+");
    if (*next < p_ntasks) {
      return true;
    } else {
      *next = -1;
      return false;
    }
  }

protected:
  
  int p_GAcounter;
  int p_ntasks;
};


} // namespace gridpack
} // namespace parallel

#endif

