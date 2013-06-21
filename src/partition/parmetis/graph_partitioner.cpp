/**
 * @file   graph_partitioner.cpp
 * @author William A. Perkins
 * @date   2013-06-21 11:31:50 d3g096
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
  : utility::Uncopyable(),
    p_impl(new ParMETISGraphPartitionerImpl(comm))
{
  // empty
}

GraphPartitioner::GraphPartitioner(const parallel::Communicator& comm,
                                   const int& local_nodes, const int& local_edges)
  : utility::Uncopyable(),
    p_impl(new ParMETISGraphPartitionerImpl(comm, local_nodes, local_edges))
{
  // empty
}

} // namespace network
} // namespace gridpack
