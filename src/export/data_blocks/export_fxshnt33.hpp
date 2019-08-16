/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_fxshnt33.hpp
 *
 * This class exports fixed shunt data in PSS/E v33 format
 *
 *  Created on: April 4, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTFXSHNT33_HPP_
#define EXPORTFXSHNT33_HPP_

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
class ExportFxShnt33
{
  public:

    /**
     * Constructor
     */
    explicit ExportFxShnt33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportFxShnt33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeFxShntBlock(std::ofstream &fout)
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
          char *ptr;
          int nshnt = 0;
          data->getValue(SHUNT_NUMBER,&nshnt);
          for (j=0; j<nshnt; j++) {
            ptr = buf;
            ival = p_network->getOriginalBusIndex(i);
            data->getValue(SHUNT_BUSNUMBER,&ival,j);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            data->getValue(SHUNT_ID,&sval,j);
            sprintf(ptr," \'%s\',",sval.c_str());
            ptr += strlen(ptr);
            ival = 1;
            data->getValue(SHUNT_STATUS,&ival,j); 
            sprintf(ptr," %d,",ival);
            ptr += strlen(ptr);
            rval = 0.0;
            data->getValue(BUS_SHUNT_GL,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 0.0;
            data->getValue(BUS_SHUNT_BL,&rval,j);
            sprintf(ptr," %f\n",rval);
            text_line text;
            strcpy(text.text,buf);
            text.global_idx = p_network->getGlobalBusIndex(i);
            text.device_idx = j;
            text_data.push_back(text);
          }
        }
      }
      if (me == 0) {
        fout << "0 / END LOAD DATA, BEGIN FIXED SHUNT DATA" << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTFXSHNT33_HPP_ */
