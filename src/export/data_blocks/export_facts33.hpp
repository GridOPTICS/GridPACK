/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_facts33.hpp
 *
 * This class exports FACTS data in PSS/E v33 format
 *
 *  Created on: June 6, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTFACTS33_HPP_
#define EXPORTFACTS33_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "gridpack/parallel/communicator.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/export/base_export.hpp"

namespace gridpack {
namespace expnet {

template <class _network>
class ExportFACTS33
{
  public:

    /**
     * Constructor
     */
    explicit ExportFACTS33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportFACTS33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeFACTSBlock(std::ofstream &fout)
    {
      int me = p_comm.rank();
      // BaseExport<_network> exprt(p_comm);
      if (me == 0) {
        fout << "0 / END OWNER DATA, BEGIN FACTS DATA" 
          << std::endl;
      }
      // exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTFACTS33_HPP_ */
