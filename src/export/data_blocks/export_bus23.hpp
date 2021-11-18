/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_bus23.hpp
 *
 * This class exports bus data in PSS/E v23 format
 *
 *  Created on: July 9, 2021
 *      Author: Bruce Palmer
 */

#ifndef EXPORTBUS23_HPP_
#define EXPORTBUS23_HPP_

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
class ExportBus23
{
  public:

    /**
     * Constructor
     */
    explicit ExportBus23(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportBus23(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeBusBlock(std::ofstream &fout)
    {
      BaseExport<_network> exprt(p_comm);
      int me = p_comm.rank();
      // Form vector of text strings based from data collection objects
      int nbus = p_network->numBuses();
      gridpack::component::DataCollection *data;
      int i;
      char buf[MAX_STRING_SIZE];
      std::vector<text_line> text_data;
      // Find value of SBASE and CASE_ID
      double sbase = 0.0;
      double case_id = 0;
      if (nbus > 0) {
        data = p_network->getBusData(0).get();
        sbase = 100.0;
        data->getValue(CASE_SBASE,&sbase);
        data->getValue(CASE_ID,&case_id);
      }

      for (i=0; i<nbus; i++) {
        if (p_network->getActiveBus(i)) {
          data = p_network->getBusData(i).get();
          double rval;
          int ival;
          std::string sval;
          char *ptr = buf;
          data->getValue(BUS_NUMBER,&ival);
          sprintf(ptr,"%d,",ival);
          ptr += strlen(ptr);
          ival = 1;
          data->getValue(BUS_TYPE,&ival); 
          sprintf(ptr," %d,",ival);
          ptr += strlen(ptr);
          rval = 0.0;
          data->getValue(LOAD_PL,&rval,0);
          sprintf(ptr," %f,",rval);
          ptr += strlen(ptr);
          rval = 0.0;
          data->getValue(LOAD_QL,&rval,0);
          sprintf(ptr," %f,",rval);
          ptr += strlen(ptr);
          rval = 0.0;
          data->getValue(BUS_SHUNT_GL,&rval,0);
          sprintf(ptr," %f,",rval);
          ptr += strlen(ptr);
          rval = 0.0;
          data->getValue(BUS_SHUNT_BL,&rval,0);
          sprintf(ptr," %f,",rval);
          ptr += strlen(ptr);
          ival = 1;
          data->getValue(BUS_AREA,&ival); 
          sprintf(ptr," %d,",ival);
          ptr += strlen(ptr);
          rval = 1.0;
          if (!data->getValue("BUS_PF_VMAG",&rval)) {
            data->getValue(BUS_VOLTAGE_MAG,&rval);
          }
          sprintf(ptr," %16.12f,",rval);
          ptr += strlen(ptr);
          rval = 0.0;
          if (!data->getValue("BUS_PF_VANG",&rval)) {
            data->getValue(BUS_VOLTAGE_ANG,&rval);
          }
          sprintf(ptr," %16.12f,",rval);
          ptr += strlen(ptr);
          data->getValue(BUS_NAME,&sval);
          if (sval.find_first_of('\'') != std::string::npos) {
            sprintf(ptr," %s,",sval.c_str());
          } else {
            sprintf(ptr," \'%s\',",sval.c_str());
          }
          ptr += strlen(ptr);
          rval = 0.0;
          data->getValue(BUS_BASEKV,&rval);
          sprintf(ptr," %f,",rval);
          ptr += strlen(ptr);
          ival = 1;
          data->getValue(BUS_ZONE,&ival); 
          sprintf(ptr," %d\n",ival);
          text_line text;
          strcpy(text.text,buf);
          text.global_idx = p_network->getGlobalBusIndex(i);
          text.device_idx = 0;
          text_data.push_back(text);
        }
      }
      if (me == 0) {
        fout << case_id<<",  "<<sbase<<std::endl;
        fout << std::endl;
        fout << "/ BEGIN BUS DATA" << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTBUS23_HPP_ */
