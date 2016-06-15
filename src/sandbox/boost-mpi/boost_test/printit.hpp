// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: printit.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 16, 2013 by William A. Perkins
// Last Change: 2013-12-16 09:13:18 d3g096
// -------------------------------------------------------------


#ifndef _printit_hpp_
#define _printit_hpp_

#include <boost/mpi.hpp>

// -------------------------------------------------------------
// printit
// -------------------------------------------------------------
template <typename T> 
void
printit(const boost::mpi::communicator& comm, 
        const std::vector<T>& things,
        const std::string& caption)
{
  if (!caption.empty() && comm.rank() == 0) {
    std::cout << caption << std::endl;
    std::cout.flush();
  }
  for (int p = 0; p < comm.size(); ++p) {
    if (comm.rank() == p) {
      std::cout << p << ": ";
      std::copy(things.begin(), things.end(),
                std::ostream_iterator<T>(std::cout, ","));
      std::cout << std::endl;
      std::cout.flush();
    }
    comm.barrier();
  }

  comm.barrier();
  size_t global_size;
  boost::mpi::reduce(comm, things.size(), global_size, std::plus<size_t>(), 0);
  if (comm.rank() == 0) {
    if (!caption.empty()) {
      std::cout << caption;
    }
    std::cout << "Number of things: " << global_size << std::endl;
  }
  comm.barrier();
  std::cout << comm.rank() << ": leaving printit()" << std::endl;
}                                      


#endif
