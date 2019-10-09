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
      gridpack::component::DataCollection *data = p_network->getNetworkData().get();
      if (me == 0) {
        fout << "0 / END TRANSFORMER DATA, BEGIN AREA DATA" << std::endl;
        int narea;
        data->getValue(AREA_TOTAL,&narea);
        int i, j;
        char buf[MAX_STRING_SIZE];
        for (i=0; i<narea; i++) {
          double rval;
          int ival;
          std::string sval;
          char *ptr;
          ptr = buf;
          data->getValue(AREAINTG_NUMBER,&ival,i);
          sprintf(ptr,"%d,",ival);
          ptr += strlen(ptr);
          data->getValue(AREAINTG_ISW,&ival,i);
          sprintf(ptr,"%d,",ival);
          ptr += strlen(ptr);
          rval = 0.0;
          data->getValue(AREAINTG_PDES,&rval,i);
          sprintf(ptr,"%f,",rval);
          ptr += strlen(ptr);
          rval = 10.0;
          data->getValue(AREAINTG_PTOL,&rval,i);
          sprintf(ptr,"%f,",rval);
          ptr += strlen(ptr);
          sval = "\'            \'";
          data->getValue(AREAINTG_NAME,&sval,i);
          sprintf(ptr,"%s",sval.c_str());
          fout << buf << std::endl;
        }
      }
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTAREA33_HPP_ */
