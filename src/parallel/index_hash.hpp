// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   index_hash.hpp
 * @author Bruce Palmer
 * @date   2014-06-24 09:25:03 d3g293
 * 
 * @brief  
 * This is a utility that is designed to provide a relatively efficient way of
 * mapping between different sets of indexes in a distributed way. Note that all
 * these operations are collective
 * 
 */

// -------------------------------------------------------------

#ifndef _hash_map_hpp_
#define _hash_map_hpp_

#include <boost/unordered_map.hpp>
#include <vector>
#include <parallel/communicator.hpp>

namespace gridpack {
namespace hash_map {

// -------------------------------------------------------------
//  class GlobalIndexHashMap
// -------------------------------------------------------------
class GlobalIndexHashMap {
public:

  // Default constructor
  GlobalIndexHashMap(const parallel::Communicator &comm);

  // Default destructor
  ~GlobalIndexHashMap(void);

  // add key-value pairs to hash map where key is singe integer
  // @param pairs list of key-value pairs where both keys and values are
  //              integers
  void addPairs(std::vector<std::pair<int,int> > &pairs);

  // add key-value pairs to hash map where key is another index pair of integers
  // @param pairs list of key-value pairs where key is a pair of integers and
  //              value is a single integer
  void addPairs(std::vector<std::pair<std::pair<int,int>,int> > &pairs);

  // get values corresponding to a list of keys from the hash map where key is a
  // single integer
  // @param keys list of integer keys
  // @param values returned list of values corresponding to the list of keys
  void getValues(std::vector<int> &keys, std::vector<int> &values);

  // get values corresponding to a list of keys from the hash map where key is a
  // pair of integers
  // @param keys list of integer-pair keys
  // @param values returned list of values corresponding to the list of keys
  void getValues(std::vector<std::pair<int,int> > &keys, std::vector<int> &values);

private:

  // hash function for indices. Maps the value of key into the interval
  // [0,p_nprocs-1] where p_nprocs is the total number of processors
  // @param key input value of key
  // @return number between 0 and p_nprocs-1
  int hashValue(int key);

  // hash function for index pairs. Maps the value of key into the interval
  // [0,p_nprocs-1] where p_nprocs is the total number of processors
  // @param key input value of key
  // @return number between 0 and nprocs-1
  int pairHashValue(std::pair<int,int> key);

  int p_nprocs;
  int p_me;

  MPI_Comm p_comm;

  boost::unordered_map<int, int> p_umap;
  
  boost::unordered_map<std::pair<int,int>, int> p_pmap;
};


} // namespace gridpack
} // namespace hash_map

#endif

