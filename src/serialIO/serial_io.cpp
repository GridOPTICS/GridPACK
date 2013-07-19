// -------------------------------------------------------------
/**
 * @file   serial_io.cpp
 * @author Bruce Palmer
 * @date   2013-07-18 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/serialIO/serial_io.hpp"

#include <iostream>

/**
 * Simple constructor
 * @param max_str_len: the maximum string length written out by any bus
 * @param network: the network for which output is desired
 */
gridpack::serial_io::SerialBusIO::SerialBusIO(int max_str_len,
                     boost::shared_ptr<Network> network){
  p_GA_type = NGA_Register_type(max_str_len);
  p_network = network;
  p_size = max_str_len;
  
  // Find total number of buses in network and create GAs for moving strings
  // around
  int nbus = p_network->totalBuses();
  int one = 1;
  p_stringGA = GA_Create_handle();
  GA_Set_data(p_stringGA,one,&nbus,p_GA_type);
  GA_Allocate(p_stringGA);
  p_maskGA = GA_Create_handle();
  GA_Set_data(p_maskGA,one,&nbus,C_INT);
  GA_Allocate(p_maskGA);
}

/**
 * Simple Destructor
 */
gridpack::serial_io::SerialBusIO::~SerialBusIO(void)
{
  NGA_Deregister_type(p_GA_type);
  GA_Destroy(p_stringGA);
  GA_Destroy(p_maskGA);
}

/**
 * Write output from buses to standard out
 * @param signal: an optional character string used to control contents of
 *                output
 */
void gridpack::serial_io::SerialBusIO::write(char *signal)
{
  int nBus = p_network->totalBuses();
  char string[p_size];
  int nwrites = 0;
  int i;
  int one = 1;

  // Count up total strings being written from this processor
  for (i=0; i<nBus; i++) {
    if (p_network->getActiveBus(i) &&
        p_network->getBus(i)->serialWrite(string,signal)) nwrites++;
  }

  // Set up buffers to scatter strings to global buffer
  int **index;
  index = new int*[nwrites];
  int ones[nwrites];
  char strbuf[nwrites*p_size];
  char *ptr = strbuf;
  nwrites = 0;
  for (i=0; i<nBus; i++) {
    if (p_network->getActiveBus(i) &&
        p_network->getBus(i)->serialWrite(ptr,signal)) {
      index[nwrites] = new int;
      *(index[nwrites]) = p_network->getGlobalBusIndex(i);
      ones[nwrites] = 1;
      nwrites++;
      ptr += p_size;
    }
  }

  // Scatter data to global buffer and set mask array
  GA_Zero(p_maskGA);
  if (nwrites > 1) {
    NGA_Scatter(p_stringGA,strbuf,index,nwrites);
    NGA_Scatter(p_maskGA,ones,index,nwrites);
  }
  GA_Sync();
  for (i=0; i<nwrites; i++) {
    delete index;
  }
  delete [] index;

  // String data is now stored on global array. Process 0 now retrieves data
  // from each successive processor and writes it to standard out  
  if (GA_Nodeid() == 0) {
    int nprocs = GA_Nnodes();
    int lo, hi;
    for (i=0; i<nprocs; i++) {
      NGA_Distribution(p_maskGA, i, &lo, &hi);
      int ld = hi - lo + 1;
      // Figure out how many strings are coming from process i
      int imask[ld];
      NGA_Get(p_maskGA,&lo,&hi,imask,&one);
      int j;
      nwrites = 0;
      for (j=0; j<ld; j++) {
        if (imask[i] == 1) {
          nwrites++;
        }
      }
      // Create buffers to retrieve strings from process i
      if (nwrites > 0) {
        nwrites = 0;
	char iobuf[p_size*nwrites];
	index = new int*[nwrites];
	for (j=0; j<ld; j++) {
	  if (imask[i] == 1) {
            index[nwrites] = new int;
            *(index[nwrites]) = j + lo;
            nwrites++;
	  }
	}
        NGA_Gather(p_stringGA,iobuf,index,nwrites);
        ptr = iobuf;
        nwrites++;
	for (j=0; j<ld; j++) {
	  if (imask[i] == 1) {
            std::cout << ptr;
            ptr += p_size;
            nwrites++;
            delete index[j];
	  }
	}
        delete [] index;
      }
    }
  }
  GA_Sync();
}

/**
 * Simple constructor
 * @param max_str_len: the maximum string length written out by any branch
 * @param network: the network for which output is desired
 */
gridpack::serial_io::SerialBranchIO::SerialBranchIO(int max_str_len,
                     boost::shared_ptr<Network> network)
{
  p_GA_type = NGA_Register_type(max_str_len);
  p_network = network;
  p_size = max_str_len;

  // Find total number of branches in network and create GAs for moving strings
  // around
  int nbranch = p_network->totalBranches();
  int one = 1;
  p_stringGA = GA_Create_handle();
  GA_Set_data(p_stringGA,one,&nbranch,p_GA_type);
  GA_Allocate(p_stringGA);
  p_maskGA = GA_Create_handle();
  GA_Set_data(p_maskGA,one,&nbranch,C_INT);
  GA_Allocate(p_maskGA);
}

/**
 * Simple Destructor
 */
gridpack::serial_io::SerialBranchIO::~SerialBranchIO(void)
{
  NGA_Deregister_type(p_GA_type);
  GA_Destroy(p_stringGA);
  GA_Destroy(p_maskGA);
}

/**
 * Write output from branches to standard out
 * @param signal: an optional character string used to control contents of
 *                output
 */
void gridpack::serial_io::SerialBranchIO::write(char *signal)
{
  int nBranch = p_network->totalBranches();
  char string[p_size];
  int nwrites = 0;
  int i;
  int one = 1;

  // Count up total strings being written from this processor
  for (i=0; i<nBranch; i++) {
    if (p_network->getActiveBranch(i) &&
        p_network->getBranch(i)->serialWrite(string,signal)) nwrites++;
  }

  // Set up buffers to scatter strings to global buffer
  int **index;
  index = new int*[nwrites];
  int ones[nwrites];
  char strbuf[nwrites*p_size];
  char *ptr = strbuf;
  nwrites = 0;
  for (i=0; i<nBranch; i++) {
    if (p_network->getActiveBranch(i) &&
        p_network->getBranch(i)->serialWrite(ptr,signal)) {
      index[nwrites] = new int;
      *(index[nwrites]) = p_network->getGlobalBranchIndex(i);
      ones[nwrites] = 1;
      nwrites++;
      ptr += p_size;
    }
  }

  // Scatter data to global buffer and set mask array
  GA_Zero(p_maskGA);
  if (nwrites > 1) {
    NGA_Scatter(p_stringGA,strbuf,index,nwrites);
    NGA_Scatter(p_maskGA,ones,index,nwrites);
  }
  GA_Sync();
  for (i=0; i<nwrites; i++) {
    delete index;
  }
  delete [] index;

  // String data is now stored on global array. Process 0 now retrieves data
  // from each successive processor and writes it to standard out  
  if (GA_Nodeid() == 0) {
    int nprocs = GA_Nnodes();
    int lo, hi;
    for (i=0; i<nprocs; i++) {
      NGA_Distribution(p_maskGA, i, &lo, &hi);
      int ld = hi - lo + 1;
      // Figure out how many strings are coming from process i
      int imask[ld];
      NGA_Get(p_maskGA,&lo,&hi,imask,&one);
      int j;
      nwrites = 0;
      for (j=0; j<ld; j++) {
        if (imask[i] == 1) {
          nwrites++;
        }
      }
      // Create buffers to retrieve strings from process i
      if (nwrites > 0) {
        nwrites = 0;
	char iobuf[p_size*nwrites];
	index = new int*[nwrites];
	for (j=0; j<ld; j++) {
	  if (imask[i] == 1) {
            index[nwrites] = new int;
            *(index[nwrites]) = j + lo;
            nwrites++;
	  }
	}
        NGA_Gather(p_stringGA,iobuf,index,nwrites);
        ptr = iobuf;
        nwrites++;
	for (j=0; j<ld; j++) {
	  if (imask[i] == 1) {
            std::cout << ptr;
            ptr += p_size;
            nwrites++;
            delete index[j];
	  }
	}
        delete [] index;
      }
    }
  }
  GA_Sync();
}
