/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_gen33.hpp
 *
 * This class exports generator data in PSS/E v33 format
 *
 *  Created on: April 15, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTGEN33_HPP_
#define EXPORTGEN33_HPP_

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
class ExportGen33
{
  public:

    /**
     * Constructor
     */
    explicit ExportGen33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportGen33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeGenBlock(std::ofstream &fout)
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
          gridpack::ComplexType zval;
          char *ptr;
          int ngen = 0;
          data->getValue(GENERATOR_NUMBER,&ngen);
          for (j=0; j<ngen; j++) {
            ptr = buf;
            ival = p_network->getOriginalBusIndex(i);
            data->getValue(GENERATOR_BUSNUMBER,&ival,j);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            data->getValue(GENERATOR_ID,&sval,j);
            sprintf(ptr," \'%s\',",sval.c_str());
            ptr += strlen(ptr);
            rval = 0.0;
            if (data->getValue(GENERATOR_PG,&rval,j)) {
              data->getValue("GENERATOR_PF_PG",&rval,j);
            }
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 0.0;
            if (data->getValue(GENERATOR_QG,&rval,j)) {
              data->getValue("GENERATOR_PF_QG",&rval,j);
            }
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 9999.0;
            data->getValue(GENERATOR_QMAX,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = -9999.0;
            data->getValue(GENERATOR_QMIN,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_VS,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(GENERATOR_IREG,&ival,j);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            rval = 100.0;
            if (!data->getValue(GENERATOR_MBASE,&rval,j)) {
              data->getValue(CASE_SBASE,&rval);
            }
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            zval = ComplexType(0.0,1.0);
            data->getValue(GENERATOR_ZSOURCE,&zval,j);
            sprintf(ptr," %f, %f,",real(zval),imag(zval),j);
            ptr += strlen(ptr);
            rval = 0.0;
            data->getValue(GENERATOR_RT,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 0.0;
            data->getValue(GENERATOR_XT,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_GTAP,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = 1;
            data->getValue(GENERATOR_STAT,&ival,j);
            sprintf(ptr,"%d,",ival);
            ptr += strlen(ptr);
            rval = 100.0;
            data->getValue(GENERATOR_RMPCT,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = 9999.0;
            data->getValue(GENERATOR_PMAX,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            rval = -9999.0;
            data->getValue(GENERATOR_PMIN,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = p_network->getOriginalBusIndex(i);
            data->getValue(GENERATOR_OWNER1,&ival,j);
            sprintf(ptr," %d,",ival);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_OFRAC1,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(GENERATOR_OWNER2,&ival,j);
            sprintf(ptr," %d,",ival);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_OFRAC2,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(GENERATOR_OWNER3,&ival,j);
            sprintf(ptr," %d,",ival);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_OFRAC3,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(GENERATOR_OWNER4,&ival,j);
            sprintf(ptr," %d,",ival);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_OFRAC4,&rval,j);
            sprintf(ptr," %f,",rval);
            ptr += strlen(ptr);
            ival = 0;
            data->getValue(GENERATOR_WMOD,&ival,j);
            sprintf(ptr," %d,",ival);
            ptr += strlen(ptr);
            rval = 1.0;
            data->getValue(GENERATOR_WPF,&rval,j);
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
        fout << "0 / END FIXED SHUNT DATA, BEGIN GENERATOR DATA" << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTGEN33_HPP_ */
