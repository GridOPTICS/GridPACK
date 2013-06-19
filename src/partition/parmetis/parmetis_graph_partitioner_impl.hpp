// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   parmetis_graph_partitioner_impl.hpp
 * @author William A. Perkins
 * @date   2013-06-19 11:28:17 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 19, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _parmetis_graph_partitioner_impl_hpp_
#define _parmetis_graph_partitioner_impl_hpp_

#include "graph_partitioner_implementation.hpp"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class ParMETISGraphPartitionerImpl
// -------------------------------------------------------------
class ParMETISGraphPartitionerImpl 
  : public GraphPartitionerImplementation {
public:

  /// Default constructor.
  ParMETISGraphPartitionerImpl(const parallel::Communicator& comm);

    /// Construct w/ known local sizes (guesses to size containers, maybe)
  ParMETISGraphPartitionerImpl(const parallel::Communicator& comm,
                               const int& local_nodes, const int& local_edges);

  /// Destructor
  ~ParMETISGraphPartitionerImpl(void);

protected:

  /// Partition the graph (specialized)
  void p_partition(void);

};



} // namespace network
} // namespace gridpack

#endif
