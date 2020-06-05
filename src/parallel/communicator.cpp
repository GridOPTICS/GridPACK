// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   communicator.cpp
 * @author William A. Perkins
 * @date   2019-12-05 08:34:28 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <functional>

#if USE_PROGRESS_RANKS
#include <ga-mpi.h>
#endif
#include <ga++.h>
#include "gridpack/utilities/uncopyable.hpp"
#include "communicator.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class CommunicatorPrivate
// -------------------------------------------------------------
class CommunicatorPrivate {
public:

  /// How many things were created?
  static int created;
  /// How many things were destroyed?
  static int destroyed;

  /// Default constructor.
  CommunicatorPrivate(void)
    : p_handle(GA_Pgroup_get_world()), p_created(false)
  {
  }

  /// Construct on the processes in the specified communicators
  CommunicatorPrivate(const boost::mpi::communicator& comm)
    : p_handle(), p_created(true)
  {
    int me(comm.rank()), nprocs(comm.size());
    int gaWrld = GA_Pgroup_get_world();
    int defGrp = GA_Pgroup_get_default();
    GA_Pgroup_set_default(gaWrld);
    std::vector<int> gaSrc(nprocs, 0), gaDest(nprocs, 0);
    gaSrc[me] = GA_Nodeid();
    boost::mpi::all_reduce(comm, &gaSrc[0], nprocs, &gaDest[0], std::plus<int>());
    p_handle = GA_Pgroup_create(&gaDest[0],nprocs);
    // GA_Pgroup_sync(p_handle);
    GA_Pgroup_set_default(defGrp);
    created++;
  }

  /// Destructor
  ~CommunicatorPrivate(void)
  {
    if (p_created) {
      GA_Pgroup_destroy(p_handle);
      destroyed--;
    }
  }

  /// Get GA process group handle
  const int& handle(void) const
  {
    return p_handle;
  }

private:

  /// The GA process group handle
  int p_handle;
  
  /// Was the GA process group created for this instance
  const bool p_created;
};

int CommunicatorPrivate::created(0);
int CommunicatorPrivate::destroyed(0);

// -------------------------------------------------------------
//  class Communicator
// -------------------------------------------------------------

// -------------------------------------------------------------
// Communicator::report
//
// Just used to make sure Communicator creation and destruction does
// not leak memory.
// -------------------------------------------------------------
void
Communicator::report(void)
{
  std::cout << "CommunicatorPrivate instances created:   " << CommunicatorPrivate::created
            << std::endl
            << "CommunicatorPrivate instances destroyed: " << CommunicatorPrivate::destroyed
            << std::endl;
}

// -------------------------------------------------------------
// Communicator:: constructors / destructor
// -------------------------------------------------------------
Communicator::Communicator(void)
#if USE_PROGRESS_RANKS
  : p_comm(GA_MPI_Comm(),boost::mpi::comm_duplicate),
    p_private(new CommunicatorPrivate())
#else
  : p_comm(), p_private(new CommunicatorPrivate())
#endif
{
  
}

Communicator::Communicator(const Communicator& old)
  : p_comm(old.p_comm), p_private(old.p_private)
{
}

Communicator::Communicator(const boost::mpi::communicator& comm)
  : p_comm(comm), p_private(new CommunicatorPrivate(p_comm))
{
}


Communicator::Communicator(const MPI_Comm& comm)
  : p_comm(comm, boost::mpi::comm_duplicate),
    p_private(new CommunicatorPrivate(p_comm))
{
}


Communicator::~Communicator(void)
{
}

// -------------------------------------------------------------
// Communicator::swap
// -------------------------------------------------------------
void
Communicator::swap(Communicator& other)
{
  std::swap(p_comm, other.p_comm);
  std::swap(p_private, other.p_private);
}

// -------------------------------------------------------------
// Communicator::operator=
// -------------------------------------------------------------
Communicator & 
Communicator::operator= (const Communicator & rhs) 
{
  Communicator temp(rhs);
  this->swap(temp);
  return (*this);
}

// -------------------------------------------------------------
// Communicator::worldRank
// -------------------------------------------------------------
int 
Communicator::worldRank(void) const
{
  return GA::SERVICES.nodeid();
}

// -------------------------------------------------------------
// Communicator::getGroup
// -------------------------------------------------------------
int
Communicator::getGroup(void) const
{
  return p_private->handle();
}

// -------------------------------------------------------------
// Communicator::divide
// -------------------------------------------------------------
Communicator 
Communicator::divide(int nsize) const
{
  int nprocs(size());
  int me(rank());
  // find out how many communicators need to be created
  int ngrp = nprocs/nsize;
  if (ngrp*nsize < nprocs) ngrp++;
  // evaluate how many communicators are of size nsize
  int nlarge = nprocs;
  if (ngrp*nsize > nprocs) {
    int nsmall = ngrp*nsize - nprocs;
    nlarge = nprocs - nsmall;
  }
  nlarge = nlarge/nsize;
  if (nprocs - nlarge*nsize  > nsize) nlarge++;
  // figure out what color this process is and use split functionality to
  // create new communicators
  int i;
  int color = 0;
  for (i=0; i<=me; i++) {
    if (i<nlarge*nsize) {
      if (i%nsize == 0) color++;
    } else {
      if ((i-nlarge*nsize)%(nsize-1) == 0) color++;
    }
  }
  return this->split(color);
}

// -------------------------------------------------------------
// Communicator::sync
// -------------------------------------------------------------
void 
Communicator::sync(void) const
{
  GA_Pgroup_sync(p_private->handle());
}    

/**
 * Sum vector over all processors in the communicator
 * @param x vector of values to be summed
 * @param nvals number of values in vector
 */
void Communicator::sum(float *x, int nvals) const
{
  int i;
  float *src = new float[nvals];
  float *dest = new float[nvals];
  for (i=0; i<nvals; i++) {
    src[i] = x[i];
  }
  boost::mpi::all_reduce(p_comm, &src[0], nvals, &dest[0], std::plus<float>());
  for (i=0; i<nvals; i++) {
    x[i] = dest[i];
  }
  delete [] src;
  delete [] dest;
}

void Communicator::sum(double *x, int nvals) const
{
  int i;
  double *src = new double[nvals];
  double *dest = new double[nvals];
  for (i=0; i<nvals; i++) {
    src[i] = x[i];
  }
  boost::mpi::all_reduce(p_comm, &src[0], nvals, &dest[0], std::plus<double>());
  for (i=0; i<nvals; i++) {
    x[i] = dest[i];
  }
  delete [] src;
  delete [] dest;
}

void Communicator::sum(int *x, int nvals) const
{
  int i;
  int *src = new int[nvals];
  int *dest = new int[nvals];
  for (i=0; i<nvals; i++) {
    src[i] = x[i];
  }
  boost::mpi::all_reduce(p_comm, &src[0], nvals, &dest[0], std::plus<int>());
  for (i=0; i<nvals; i++) {
    x[i] = dest[i];
  }
  delete [] src;
  delete [] dest;
}

void Communicator::sum(long *x, int nvals) const
{
  int i;
  long *src = new long[nvals];
  long *dest = new long[nvals];
  for (i=0; i<nvals; i++) {
    src[i] = x[i];
  }
  boost::mpi::all_reduce(p_comm, &src[0], nvals, &dest[0], std::plus<long>());
  for (i=0; i<nvals; i++) {
    x[i] = dest[i];
  }
  delete [] src;
  delete [] dest;
}

void Communicator::sum(gridpack::ComplexType *x, int nvals) const
{
  int i;
  gridpack::ComplexType *src = new gridpack::ComplexType[nvals];
  gridpack::ComplexType *dest = new gridpack::ComplexType[nvals];
  for (i=0; i<nvals; i++) {
    src[i] = x[i];
  }
  boost::mpi::all_reduce(p_comm, &src[0], nvals, &dest[0],
      std::plus<gridpack::ComplexType>());
  for (i=0; i<nvals; i++) {
    x[i] = dest[i];
  }
  delete [] src;
  delete [] dest;
}

/**
 * Find maximum of vector components over all processors
 * in the communicator
 * @param x vector of values to be evaluated
 * @param nvals number of values in vector
 */
void Communicator::max(float *x, int nvals) const
{
  char cmax[4];
  strcpy(cmax,"max");
  GA_Pgroup_fgop(p_private->handle(),x,nvals,cmax);
}

void Communicator::max(double *x, int nvals) const
{
  char cmax[4];
  strcpy(cmax,"max");
  GA_Pgroup_dgop(p_private->handle(),x,nvals,cmax);
}

void Communicator::max(int *x, int nvals) const
{
  char cmax[4];
  strcpy(cmax,"max");
  GA_Pgroup_igop(p_private->handle(),x,nvals,cmax);
}

void Communicator::max(long *x, int nvals) const
{
  char cmax[4];
  strcpy(cmax,"max");
  GA_Pgroup_lgop(p_private->handle(),x,nvals,cmax);
}

/**
 * Find minimum of vector components over all processors
 * in the communicator
 * @param x vector of values to be evaluated
 * @param nvals number of values in vector
 */
void Communicator::min(float *x, int nvals) const
{
  char cmin[4];
  strcpy(cmin,"min");
  GA_Pgroup_fgop(p_private->handle(),x,nvals,cmin);
}

void Communicator::min(double *x, int nvals) const
{
  char cmin[4];
  strcpy(cmin,"min");
  GA_Pgroup_dgop(p_private->handle(),x,nvals,cmin);
}

void Communicator::min(int *x, int nvals) const
{
  char cmin[4];
  strcpy(cmin,"min");
  GA_Pgroup_igop(p_private->handle(),x,nvals,cmin);
}

void Communicator::min(long *x, int nvals) const
{
  char cmin[4];
  strcpy(cmin,"min");
  GA_Pgroup_lgop(p_private->handle(),x,nvals,cmin);
}

bool Communicator::any(const bool& lflag) const
{
  bool result(false);
  boost::mpi::all_reduce(p_comm, lflag, result, std::logical_or<bool>());
  return result;
}

bool Communicator::all(const bool& lflag) const
{
  bool result(false);
  boost::mpi::all_reduce(p_comm, lflag, result, std::logical_and<bool>());
  return result;
}

} // namespace parallel
} // namespace gridpack
