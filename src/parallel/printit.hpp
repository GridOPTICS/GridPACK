// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   printit.hpp
 * @author William A. Perkins
 * @date   2013-08-01 09:38:00 d3g096
 * 
 * @brief  Declaration of prinit<>()
 * 
 * 
 */

// -------------------------------------------------------------

#include <iostream>
#include <iterator>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>

#ifndef _printit_hpp_
#define _printit_hpp_

// -------------------------------------------------------------
// printit
// -------------------------------------------------------------
template <typename T> 
void
printit(const boost::mpi::communicator& comm, 
        const std::vector<T> things,
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

  size_t global_size;
  boost::mpi::reduce(comm, things.size(), global_size, std::plus<size_t>(), 0);
  if (comm.rank() == 0) {
    if (!caption.empty()) {
      std::cout << caption;
    }
    std::cout << "Number of things: " << global_size << std::endl;
  }
}                                      



#endif
