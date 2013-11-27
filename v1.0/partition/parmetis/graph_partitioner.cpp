/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   graph_partitioner.cpp
 * @author William A. Perkins
 * @date   2013-11-08 11:51:42 d3g096
 * 
 * @brief  
 * 
 * 
 */


#include "graph_partitioner.hpp"
#include "parmetis/parmetis_graph_partitioner_impl.hpp"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class GraphPartitioner
// -------------------------------------------------------------

// -------------------------------------------------------------
// GraphPartitioner:: constructors / destructor
// -------------------------------------------------------------
GraphPartitioner::GraphPartitioner(const parallel::Communicator& comm)
  : parallel::WrappedDistributed(),
    utility::Uncopyable(),
    p_impl(new ParMETISGraphPartitionerImpl(comm))
{
  p_setDistributed(p_impl.get());
}

GraphPartitioner::GraphPartitioner(const parallel::Communicator& comm,
                                   const int& local_nodes, const int& local_edges)
  : parallel::WrappedDistributed(),
    utility::Uncopyable(),
    p_impl(new ParMETISGraphPartitionerImpl(comm, local_nodes, local_edges))
{
  p_setDistributed(p_impl.get());
}

} // namespace network
} // namespace gridpack
