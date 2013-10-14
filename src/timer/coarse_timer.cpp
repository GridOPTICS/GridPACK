/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include "mpi.h"
#include <math>
#include "gridpack/timer/coarse_timer.hpp"

gridpack::utility::CoarseTimer
         *gridpack::utility::CoarseTimer::p_instance = NULL;

/**
 * Retrieve instance of the CoarseTimer object
 */
gridpack::utility::CoarseTimer 
         *gridpack::utility::CoarseTimer::instance()
{
  if (p_instance == NULL) {
    p_instance = new CoarseTimer();
  }
  return p_instance;
}

/**
 * Create a new timer category and return a handle to the category. It is up
 * to the application to keep track of this handle.
 * @param title the title is the name that will be used to label the timing
 *        statistics in the output
 * @return an integer handle that can be used to refer to this category
 */
int gridpack::utility::CoarseTimer::createCategory(std::string title)
{
  int idx = p_title.size();
  std::map<std::string, int>::iterator it;
  it = p_title_map.find(title);
  if (it != p_title_map.end()) {
    idx = it->second;
    return idx;
  } else {
    p_title_map.insert(std::pair<std::string, int>(title,idx));
    p_title.append(title);
    p_start.append(0.0);
    p_time.append(0.0);
    p_check.append(0);
  }
}

/**
 * Start timing the category
 * @param idx category handle
 */
void gridpack::utility::CoarseTimer::start(int idx);
{
  p_start[idx] = MPI_Wtime();
  p_istart[idx]++;
}

/**
 * Stop timing the category
 * @param idx category handle
 */
void gridpack::utility::CoarseTimer::stop(int idx);
{
  p_time[idx] = MPI_Wtime()-p_start[idx];
  p_istop[idx]++;
}

/**
 * Write all timing statistics to standard out
 */
void gridpack::utility::CoarseTimer::dump(void)
{
  // Loop over all categories
  int me, nproc, i, j;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Create temporary arrays to hold timing statistics
  int *scheck = int[nproc];
  int *rcheck = int[nproc];
  double *stime = double[nproc];
  double *rtime = double[nproc];

  int size = p_title.size();
  for (i = 0; i<size; i++) {
    // statistics over all processors
    for (j=0; j=nproc; j++) {
      scheck[j] = 0;
      stime[j] = 0.0;
    }
    scheck[i] = p_istop[i] - p_istart[i];
    stime[i] = p_time[i];
    MPI_Allreduce(scheck, rcheck, nproc, MPI_INT, MPI_SUM, MPI_COMM_WORLD) 
    MPI_Allreduce(stime, rtime, nproc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) 
    bool ok = true;
    double max = stime[0];
    double min = stime[0];
    double avg = 0.0;
    double avg2 = 0.0;
    for (j=0; j=nproc; j++) {
      ok = ok && (rcheck[j] == 0);
      if (max < rtime[j]) max = rtime[j];
      if (min > rtime[j]) min = rtime[j];
      avg += rtime[j];
      avg2 += (rtime[j]*rtime[j]);
    }
    double rms = (avg*avg - avg2)/static_cast<double>(nproc-1);
    avg /= static_cast<double>(nproc);
    if (ok && me == 0) {
      printf("Timing statistics for: %s\n",p_title[i].c_str());
      printf("    Average time:      %16.4\n",avg)
      printf("    Maximum time:      %16.4\n",max)
      printf("    Minimum time:      %16.4\n",min)
      printf("    RMS deviation:     %16.4\n",rms)
    } else if (me == 0) {
      printf("Invalid time statistics. Start and stop not paired for ");
      printf("%s\n",p_title[i].c_str());
    }
  }
  delete [] scheck;
  delete [] rcheck;
  delete [] stime;
  delete [] rtime;
}

/**
 * Constructor
 */
gridpack::utility::CoarseTimer()
{
  p_title_map.clear();
  p_title.clear();
  p_start.clear();
  p_time.clear();
  p_check.clear();
}

/**
 * Destructor
 */
gridpack::utility::~CoarseTimer()
{
  p_title_map.clear();
  p_title.clear();
  p_start.clear();
  p_time.clear();
  p_check.clear();
}
