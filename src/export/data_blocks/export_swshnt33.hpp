/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_facts33.hpp
 *
 * This class exports switched shunt data in PSS/E v33 format
 *
 *  Created on: June 6, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTSWSHNT33_HPP_
#define EXPORTSWSHNT33_HPP_

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
class ExportSwShnt33
{
  public:

    /**
     * Constructor
     */
    explicit ExportSwShnt33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportSwShnt33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeSwShntBlock(std::ofstream &fout)
    {
      BaseExport<_network> exprt(p_comm);
      int me = p_comm.rank();
      // Form vector of text strings based on data collection
      int nbus = p_network->numBuses();
      gridpack::component::DataCollection *data;
      int i;
      char buf[MAX_STRING_SIZE];
      std::vector<text_line> text_data;
      for (i=0; i<nbus; i++) {
        if (p_network->getActiveBus(i)) {
          data = p_network->getBusData(i).get();
          double rval;
          int ival;
          std::string sval;
          char *ptr = buf;
          if (data->getValue(SWSHUNT_BUSNUMBER,&ival)) {
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            ival = 1;
            data->getValue(SHUNT_MODSW,&ival);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(SHUNT_ADJM,&ival);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            ival = 1;
            data->getValue(SHUNT_SWCH_STAT,&ival);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(SHUNT_VSWHI,&rval);
            sprintf(ptr,"%f,",rval);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(SHUNT_VSWLO,&rval);
            sprintf(ptr,"%f,",rval);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(SHUNT_SWREM,&ival);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            rval = 100.0;
            data->getValue(SHUNT_RMPCT,&rval);
            sprintf(ptr,"%f,",rval);
            ptr += strlen(ptr);
            sval = "\'            \'"; 
            data->getValue(SHUNT_RMIDNT,&sval);
            sprintf(ptr,"%s,",sval.c_str());
            ptr += strlen(ptr);
            rval = 0.0;
            data->getValue(SHUNT_BINIT,&rval);
            sprintf(ptr,"%f",rval);
            ptr += strlen(ptr);
            int icnt = 0;
            while (icnt<9) {
              char ni[32], bi[32];
              sprintf(ni,"SHUNT_N%d",icnt+1);
              sprintf(bi,"SHUNT_B%d",icnt+1);
              if (data->getValue(ni,&ival) && data->getValue(bi,&rval)) {
                sprintf(ptr,",%d,%f",ival,rval);
                ptr += strlen(ptr);
              } else {
                break;
              }
              icnt++;
            }
            sprintf(ptr,"\n");
            text_line text;
            strcpy(text.text,buf);
            text.global_idx = p_network->getGlobalBusIndex(i);
            text.device_idx = 0;
            text_data.push_back(text);
          }
        }
      }
      if (me == 0) {
        fout << "0 / END FACTS DATA, BEGIN SWITCHED SHUNT DATA" 
          << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
      if (me == 0) {
        fout << "0 / END SWITCHED SHUNT DATA" 
          << std::endl;
      }
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTSWSHNT33_HPP_ */
