/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   parmetis_graph_partitioner_impl.cpp
 * @author William A. Perkins
 * @date   2013-07-09 11:27:50 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include <algorithm>
#include <utility>
#include <parmetis.h>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include "parmetis/parmetis_graph_partitioner_impl.hpp"
#include "parmetis/parmetis_graph_wrapper.hpp"


namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class ParMETISGraphPartitionerImpl
// -------------------------------------------------------------

// -------------------------------------------------------------
// ParMETISGraphPartitionerImpl:: constructors / destructor
// -------------------------------------------------------------
ParMETISGraphPartitionerImpl::ParMETISGraphPartitionerImpl(const parallel::Communicator& comm)
  : GraphPartitionerImplementation(comm)
{
  
}

ParMETISGraphPartitionerImpl::ParMETISGraphPartitionerImpl(const parallel::Communicator& comm,
                               const int& local_nodes, const int& local_edges)
  : GraphPartitionerImplementation(comm, local_nodes, local_edges)
{
  
}

ParMETISGraphPartitionerImpl::~ParMETISGraphPartitionerImpl(void)
{
}

// -------------------------------------------------------------
// ParMETISGraphPartitionerImpl::p_partition
// -------------------------------------------------------------
/**
 * This routine is really all about taking the information in
 * p_adjacency_list and translating it into what ParMETIS understands
 * ("Distributed CSR Format").  This format has three parts: (1)
 * vertex/node distribution on processors (vtxdist), (2) an adjacency
 * vector that contains adjacency information for all local
 * vertexes/nodes (adjncy), and (3) an index into the adjacency vector
 * that indicates the adjacency info for individual nodes (xadj).
 * These are named as in the ParMETIS documentation.
 * 
 * (I'm starting to remember why I don't like ParMETIS. Serious
 * problems with ParMETIS:
 * 
 *    It assumes that the vertex/node numbering is continuous
 *    starting with 0 for the vertex/nodes on process 0.
 *
 *    The must be at least one node/vertex per processor (maybe).
 * 
 */

void 
ParMETISGraphPartitionerImpl::p_partition(void)
{
  int me(this->processor_rank());
  std::vector<idx_t> vtxdist;
  std::vector<idx_t> xadj;
  std::vector<idx_t> adjncy;

  ParMETISGraphWrapper wrap(p_adjacency_list);

  wrap.get_csr_local(vtxdist, xadj, adjncy);

  int nnodes(vtxdist[me+1] - vtxdist[me]);

#if 0
  for (int p = 0; p < this->processor_size(); ++p) {
    if (this->processor_rank() == p) {
      std::cout << "Processor " << p << ":     nodes: ";
      for (Index n = 0; n < nnodes; ++n) {
        std::cout << p_adjacency_list.node_index(n) << ",";
      }
      std::cout << std::endl;
      std::cout << "Processor " << p << ":   vtxdist: ";
      std::copy(vtxdist.begin(), vtxdist.end(),
                std::ostream_iterator<idx_t>(std::cout, ","));
      std::cout << std::endl;
      std::cout << "Processor " << p << ":      xadj: ";
      std::copy(xadj.begin(), xadj.end(),
                std::ostream_iterator<idx_t>(std::cout, ","));
      std::cout << std::endl;
      std::cout << "Processor " << p << ":    adjncy: ";
      std::copy(adjncy.begin(), adjncy.end(),
                std::ostream_iterator<idx_t>(std::cout, ","));
      std::cout << std::endl;
    }
    this->communicator().barrier();
  }
#endif

  // Call the partitioner (try to use variable names that match the documentation)

  int status;

  idx_t ncon(1);
  idx_t wgtflag(3), numflag(0);
  idx_t nparts(this->processor_size());
  std::vector<idx_t> vwgt(nnodes, 1);
  std::vector<idx_t> adjwgt(adjncy.size(), 2);
  std::vector<real_t> tpwgts(nparts*ncon, 1.0/static_cast<real_t>(nparts));
  real_t ubvec(1.05);
  std::vector<idx_t> options(3);
  options[0] = 1;
  options[1] = 127;
  options[2] = 14;
  MPI_Comm comm(this->communicator());

  idx_t edgecut;
  std::vector<idx_t> part(nnodes);
  status = ParMETIS_V3_PartKway(&vtxdist[0], 
                                &xadj[0], 
                                &adjncy[0],
                                &vwgt[0],
                                &adjwgt[0],
                                &wgtflag,
                                &numflag,
                                &ncon,
                                &nparts,
                                &tpwgts[0],
                                &ubvec,
                                &options[0],
                                &edgecut, &part[0],
                                &comm);

  // "part" contains the destination processors; transfer this to the
  // local array

  wrap.set_partition(vtxdist, part);
  wrap.get_partition(p_node_destinations);

}


} // namespace network
} // namespace gridpack
