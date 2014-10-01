/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
/**
 * @file   simple_adjacency.hpp
 * @author William A. Perkins
 * @date   2013-07-08 12:08:32 d3g096
 * 
 * @brief  An inline function for creating a simple AdjacencyList instance
 * 
 * 
 */



#ifndef _simple_adjacency_hpp_
#define _simple_adjacency_hpp_

#include "graph_partitioner.hpp"

/// Make (allocates) AdjacencyList that represents a simple linear network
gridpack::network::AdjacencyList *
simple_adjacency_list(const gridpack::parallel::Communicator& comm,
                      const int& global_nodes);

/// Make (allocate) a partitioner for a simple linear network
gridpack::network::GraphPartitioner *
simple_graph_partitioner(const gridpack::parallel::Communicator& comm,
                         const int& global_nodes);


#endif
