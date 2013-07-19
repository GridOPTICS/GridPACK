// -------------------------------------------------------------
/**
 * @file   serial_io.hpp
 * @author Bruce Palmer
 * @date   2013-07-18 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _serial_io_h_
#define _serial_io_h_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace serial_io {

// -------------------------------------------------------------
// A set of classes to support output of information from buses
// and branches to standard output. Each bus or branch is
// responsible for creating a string that can be written to
// standard out. These modules then organize these sequentially
// and write them from process 0
// -------------------------------------------------------------

typedef gridpack::network::BaseNetwork<gridpack::component::BaseBusComponent,
        gridpack::component::BaseBranchComponent> Network;

class SerialBusIO {
  public:

    /**
     * Simple constructor
     * @param max_str_len: the maximum string length written out by any bus
     * @param network: the network for which output is desired
     */
    SerialBusIO(int max_str_len, boost::shared_ptr<Network>);

    /**
     * Simple Destructor
     */
    ~SerialBusIO(void);

    /**
     * Write output from buses to standard out
     * @param signal: an optional character string used to control contents of
     *                output
     */
    void write(char *signal);

  private:
    int p_GA_type;
    boost::shared_ptr<Network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
};

class SerialBranchIO {
  public:

    /**
     * Simple constructor
     * @param max_str_len: the maximum string length written out by any branch
     * @param network: the network for which output is desired
     */
    SerialBranchIO(int max_str_len, boost::shared_ptr<Network> network);

    /**
     * Simple Destructor
     */
    ~SerialBranchIO(void);

    /**
     * Write output from buses to standard out
     * @param signal: an optional character string used to control contents of
     *                output
     */
    void write(char *signal);

  private:
    int p_GA_type;
    boost::shared_ptr<Network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
};
} // serial_io
} // gridpack
#endif  // _serial_io_h
