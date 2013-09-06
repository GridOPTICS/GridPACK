/**
 * @file   graph_partitioner.cpp
 * @author William A. Perkins
 * @date   2013-09-06 13:36:16 d3g096
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
  p_set_distributed(p_impl.get());
}

GraphPartitioner::GraphPartitioner(const parallel::Communicator& comm,
                                   const int& local_nodes, const int& local_edges)
  : parallel::WrappedDistributed(),
    utility::Uncopyable(),
    p_impl(new ParMETISGraphPartitionerImpl(comm, local_nodes, local_edges))
{
  p_set_distributed(p_impl.get());
}

} // namespace network
} // namespace gridpack
