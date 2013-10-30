/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include "mpi.h"
#include <math.h>
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
int gridpack::utility::CoarseTimer::createCategory(const std::string title)
{
  int idx = p_title.size();
  std::map<std::string, int>::iterator it;
  it = p_title_map.find(title);
  if (it != p_title_map.end()) {
    idx = it->second;
  } else {
    p_title_map.insert(std::pair<std::string, int>(title,idx));
    p_title.push_back(title);
    p_start.push_back(0.0);
    p_time.push_back(0.0);
    p_istart.push_back(0);
    p_istop.push_back(0);
  }
  return idx;
}

/**
 * Start timing the category
 * @param idx category handle
 */
void gridpack::utility::CoarseTimer::start(const int idx)
{
  p_start[idx] = MPI_Wtime();
  p_istart[idx]++;
}

/**
 * Stop timing the category
 * @param idx category handle
 */
void gridpack::utility::CoarseTimer::stop(const int idx)
{
  p_time[idx] += MPI_Wtime()-p_start[idx];
  p_istop[idx]++;
}

/**
 * Write all timing statistics to standard out
 */
void gridpack::utility::CoarseTimer::dump(void) const
{
  // Loop over all categories
  int me, nproc, i, j;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Create temporary arrays to hold timing statistics
  int *scheck = new int[nproc];
  int *rcheck = new int[nproc];
  double *stime = new double[nproc];
  double *rtime = new double[nproc];

  int size = p_title.size();
  for (i = 0; i<size; i++) {
    // statistics over all processors
    for (j=0; j<nproc; j++) {
      scheck[j] = 0;
      stime[j] = 0.0;
    }
    scheck[me] = p_istop[i] - p_istart[i];
    stime[me] = p_time[i];
    MPI_Allreduce(scheck, rcheck, nproc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(stime, rtime, nproc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    bool ok = true;
    double max = stime[0];
    double min = stime[0];
    double avg = 0.0;
    double avg2 = 0.0;
    for (j=0; j<nproc; j++) {
      ok = ok && (rcheck[j] == 0);
      if (max < rtime[j]) max = rtime[j];
      if (min > rtime[j]) min = rtime[j];
      avg += rtime[j];
      avg2 += (rtime[j]*rtime[j]);
    }
    avg /= static_cast<double>(nproc);
    double rms = avg2-static_cast<double>(nproc)*avg*avg;
    if (nproc > 1) {
      rms = rms/static_cast<double>(nproc-1);
      if (rms > 0.0) {
        rms = sqrt(rms);
      }
    } else {
      rms = -1.0;
    }
    if (ok && me == 0) {
      printf("Timing statistics for: %s\n",p_title[i].c_str());
      printf("    Average time:      %16.4f\n",avg);
      printf("    Maximum time:      %16.4f\n",max);
      printf("    Minimum time:      %16.4f\n",min);
      if (rms > 0.0) {
        printf("    RMS deviation:     %16.4f\n",rms);
      }
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
gridpack::utility::CoarseTimer::CoarseTimer()
{
  p_title_map.clear();
  p_title.clear();
  p_start.clear();
  p_time.clear();
  p_istart.clear();
  p_istop.clear();
}

/**
 * Destructor
 */
gridpack::utility::CoarseTimer::~CoarseTimer()
{
  p_title_map.clear();
  p_title.clear();
  p_start.clear();
  p_time.clear();
  p_istart.clear();
  p_istop.clear();
}
