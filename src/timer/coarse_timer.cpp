/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include "gridpack/parallel/communicator.hpp"
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
  if (!p_profile) return;
  p_start[idx] = MPI_Wtime();
  p_istart[idx]++;
}

/**
 * Stop timing the category
 * @param idx category handle
 */
void gridpack::utility::CoarseTimer::stop(const int idx)
{
  if (!p_profile) return;
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
  gridpack::parallel::Communicator comm;
  MPI_Comm world = static_cast<MPI_Comm>(comm);
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nproc);

  // Create temporary arrays to hold timing statistics
  int *scheck = new int[nproc];
  int *rcheck = new int[nproc];
  double *stime = new double[nproc];
  double *rtime = new double[nproc];

  int size = p_title.size();
  int size_min,size_max;
  i = size;
  MPI_Allreduce(&i,&size_min,1,MPI_INT,MPI_MIN, world);
  i = size;
  MPI_Allreduce(&i,&size_max,1,MPI_INT,MPI_MAX, world);
  if (size_max != size_min) {
    if (me == 0) {
      printf ("Different numbers of timing catagories on\n");
      printf ("different processors min: %d max: %d\n",size_min,size_max);
      return;
    }
  }
  for (i = 0; i<size; i++) {
    // statistics over all processors
    for (j=0; j<nproc; j++) {
      scheck[j] = 0;
      stime[j] = 0.0;
    }
    scheck[me] = p_istop[i] - p_istart[i];
    int sncheck = 0;
    int rncheck = 0;
    if (p_istop[i] > 0 || p_start[i] > 0) sncheck = 1;
    stime[me] = p_time[i];
    MPI_Allreduce(scheck, rcheck, nproc, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&sncheck, &rncheck, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(stime, rtime, nproc, MPI_DOUBLE, MPI_SUM, world);
    bool ok = true;
    double max = rtime[0];
    double min = rtime[0];
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
    if (ok && me == 0 && rncheck > 0) {
      printf("Timing statistics for: %s\n",p_title[i].c_str());
      printf("    Average time:      %16.4f\n",avg);
      printf("    Maximum time:      %16.4f\n",max);
      printf("    Minimum time:      %16.4f\n",min);
      if (rms > 0.0) {
        printf("    RMS deviation:     %16.4f\n",rms);
      }
    } else if (me == 0 && rncheck > 0) {
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
 * Write out profile for all processors for requested handle
 * @param idx category handle
 */
void gridpack::utility::CoarseTimer::dumpProfile(const int idx) const
{
  // Loop over all categories
  int me, nproc, j;
  gridpack::parallel::Communicator comm;
  MPI_Comm world = static_cast<MPI_Comm>(comm);
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nproc);

  // Create temporary arrays to hold timing statistics
  int *scheck = new int[nproc];
  int *rcheck = new int[nproc];
  double *stime = new double[nproc];
  double *rtime = new double[nproc];

  int size = p_title.size();
  // statistics over all processors
  for (j=0; j<nproc; j++) {
    scheck[j] = 0;
    stime[j] = 0.0;
  }
  scheck[me] = p_istop[idx] - p_istart[idx];
  int sncheck = 0;
  int rncheck = 0;
  if (p_istop[idx] > 0 || p_start[idx] > 0) sncheck = 1;
  stime[me] = p_time[idx];
  MPI_Allreduce(scheck, rcheck, nproc, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&sncheck, &rncheck, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(stime, rtime, nproc, MPI_DOUBLE, MPI_SUM, world);
  bool ok = true;
  if (ok && me == 0 && rncheck > 0) {
    printf("Timing profile for: %s\n",p_title[idx].c_str());
    for (j=0; j<nproc; j++) {
      printf("    Time on process[%6d]: %16.4f\n",j,rtime[j]);
    }
  } else if (me == 0 && rncheck > 0) {
    printf("Invalid time statistics. Start and stop not paired for ");
    printf("%s\n",p_title[idx].c_str());
  }
  delete [] scheck;
  delete [] rcheck;
  delete [] stime;
  delete [] rtime;
}

/**
 * Write out profile for all processors for requested title
 * @param idx category title
 */
void gridpack::utility::CoarseTimer::dumpProfile(const std::string title) const
{
  std::map<std::string, int>::const_iterator it;
  gridpack::parallel::Communicator comm;
  MPI_Comm world = static_cast<MPI_Comm>(comm);
  it = p_title_map.find(title);
  if (it != p_title_map.end()) {
    int idx = it->second;
    dumpProfile(idx);
  } else {
    int me;
    MPI_Comm_rank(world, &me);
    if (me == 0) {
      printf("No statistics for:%s\n",title.c_str()); 
    }
  }
}

/**
 * Turn timing on and off. If timing is off, no data is collected.
 * @param flag turn timer on (true) or off (false)
 */
void gridpack::utility::CoarseTimer::configTimer(bool flag)
{
  p_profile = flag;
}

/**
 * Return current time. Can be used to solve timing problems that can't be
 * handled using the regular timing capabilities
 * @return current time in seconds according to internal clock
 */
double gridpack::utility::CoarseTimer::currentTime()
{
  return MPI_Wtime();
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
  p_profile = true;
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
