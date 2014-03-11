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

#include "gridpack/parallel/communicator.hpp"
#include <ga.h>

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class TaskManager
// -------------------------------------------------------------
class TaskManager {
public:

  /**
   * Constructor on world communicator
   */
  TaskManager(void)
  {
    p_grp = GA_Pgroup_get_world();

    p_GAcounter = GA_Create_handle();
    int one = 1;
    GA_Set_data(p_GAcounter,one,&one,C_INT);
    GA_Set_pgroup(p_GAcounter,p_grp);
    if (!GA_Allocate(p_GAcounter)) {
      // TODO: some kind of error
    }
    p_task_list = NULL;
    GA_Zero(p_GAcounter);
    p_ntasks = 0;
  }

  /**
   * Constructor that uses user-specified communicator
   * @param comm user-specified communicator
   */
  TaskManager(Communicator &comm)
  {
    p_grp = comm.getGroup();

    p_GAcounter = GA_Create_handle();
    int one = 1;
    GA_Set_data(p_GAcounter,one,&one,C_INT);
    GA_Set_pgroup(p_GAcounter,p_grp);
    if (!GA_Allocate(p_GAcounter)) {
      // TODO: some kind of error
    }
    p_task_list = NULL;
    GA_Zero(p_GAcounter);
    p_ntasks = 0;
  }

  /**
   * Destructor
   */
  ~TaskManager(void)
  {
    delete [] p_task_list;
    GA_Destroy(p_GAcounter);
  }

  /**
   * Specify total number of tasks and set task manager to zero
   * @param ntasks total number of tasks
   */
  void set(int ntasks)
  {
    GA_Zero(p_GAcounter);
    p_ntasks = ntasks;
    p_task_list = new int[ntasks];
    for (int i=0; i<ntasks; i++) p_task_list[i] = 0;
  }
  
  /**
   * Get the next task from the task manager. If the manager finds a task it
   * returns true and next is set to the index of the task, otherwise it returns
   * false and next is set to -1
   * @param next index of next task
   * @return false if no other tasks are found
   */
  bool nextTask(int *next) {
    int zero = 0;
    long one = 1;
    *next = static_cast<int>(NGA_Read_inc(p_GAcounter,&zero,one));
    if (*next < p_ntasks) {
      p_task_list[*next] = GA_Pgroup_nodeid(GA_Pgroup_get_world());
      return true;
    } else {
      *next = -1;
      GA_Pgroup_sync(p_grp);
      return false;
    }
  }

  /**
   * Get the next task for the whole communicator. The same value of next is
   * returned for all processors in the communicator comm. If the manager finds
   * a task it returns true and next is set to the index of the task, otherwise
   * it returns false and next is set to -1
   * @param comm communicator for next task
   * @param next index of next task
   * @return false if no other tasks are found
   */

  bool nextTask(Communicator &comm, int *next) {
    int zero = 0;
    long one = 1;
    int me = comm.rank();
    if (me == 0) {
      *next = static_cast<int>(NGA_Read_inc(p_GAcounter,&zero,one));
    } else {
      *next = 0;
    }
    GA_Pgroup_igop(comm.getGroup(),next,one,"+");
    if (*next < p_ntasks) {
      if (me == 0) p_task_list[*next] = GA_Pgroup_nodeid(p_grp);
      return true;
    } else {
      *next = -1;
      GA_Pgroup_sync(p_grp);
      return false;
    }
  }

  /**
   * Print out statistics on how task are distributed on processors
   */
  void printStats() {
    GA_Igop(p_task_list,p_ntasks,"+");
    int nprocs = GA_Pgroup_nnodes(p_grp);
    int me = GA_Pgroup_nodeid(p_grp);
    int procs[nprocs];
    int i;
    for (i=0; i<nprocs; i++) procs[i] = 0;
    for (i=0; i<p_ntasks; i++) procs[p_task_list[i]]++;
    // print out number of tasks evaluated on each processor
    if (me == 0) {
      printf("\nNumber of tasks per processors\n");
      for (i=0; i<nprocs; i++) {
        printf("  Number of tasks on process %6d: %6d\n",i,procs[i]);
      }
    }
  }

protected:
  
  int p_GAcounter;
  int p_ntasks;
  int p_grp;
  int *p_task_list;
};


} // namespace gridpack
} // namespace parallel

#endif

