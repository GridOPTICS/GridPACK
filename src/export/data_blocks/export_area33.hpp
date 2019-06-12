/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_area33.hpp
 *
 * This class exports area data in PSS/E v33 format
 *
 *  Created on: April 4, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTAREA33_HPP_
#define EXPORTAREA33_HPP_

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
class ExportArea33
{
  public:

    /**
     * Constructor
     */
    explicit ExportArea33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportArea33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeAreaBlock(std::ofstream &fout)
    {
      BaseExport<_network> exprt(p_comm);
      int me = p_comm.rank();
      // Form vector of text strings based from data collection objects
      int nbus = p_network->numBuses();
      gridpack::component::DataCollection *data;
      int i, j;
      char buf[MAX_STRING_SIZE];
      std::vector<text_line> text_data;
      for (i=0; i<nbus; i++) {
        if (p_network->getActiveBus(i)) {
          data = p_network->getBusData(i).get();
          double rval;
          int ival;
          std::string sval;
          char *ptr = buf;
          if (data->getValue(AREAINTG_NUMBER,&ival)) {
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            ival = p_network->getOriginalBusIndex(i);
            data->getValue(AREAINTG_ISW,&ival);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            rval = 0.0;
            data->getValue(AREAINTG_PDES,&rval);
            sprintf(ptr,"%f,",ival);
            ptr += strlen(ptr);
            rval = 10.0;
            data->getValue(AREAINTG_PTOL,&rval);
            sprintf(ptr,"%f,",ival);
            ptr += strlen(ptr);
            sval = "            ";
            data->getValue(AREAINTG_NAME,&sval);
            sprintf(ptr,"%s",sval.c_str());
            text_line text;
            strcpy(text.text,buf);
            text.global_idx = p_network->getGlobalBusIndex(i);
            text.device_idx = 0;
            text_data.push_back(text);
          }
        }
      }
      if (me == 0) {
        fout << "0 / END TRANSFORMER DATA, BEGIN AREA DATA" << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTAREA33_HPP_ */
