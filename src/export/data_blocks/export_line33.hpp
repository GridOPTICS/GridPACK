/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_line33.hpp
 *
 * This class exports line data in PSS/E v33 format
 *
 *  Created on: April 15, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTLINE33_HPP_
#define EXPORTLINE33_HPP_

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
class ExportLine33
{
  public:

    /**
     * Constructor
     */
    explicit ExportLine33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportLine33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeLineBlock(std::ofstream &fout)
    {
      BaseExport<_network> exprt(p_comm);
      int me = p_comm.rank();
      // Form vector of text strings based from data collection objects
      int nbranch = p_network->numBranches();
      gridpack::component::DataCollection *data;
      int i, j;
      char buf[MAX_STRING_SIZE];
      std::vector<text_line> text_data;
      for (i=0; i<nbranch; i++) {
        if (p_network->getActiveBranch(i)) {
          data = p_network->getBranchData(i).get();
          double rval;
          int ival;
          int idx1, idx2;
          std::string sval;
          DoubleComplex zval;
          char *ptr;
          int nline = 0;
          double tap_ratio;
          data->getValue(BRANCH_NUM_ELEMENTS,&nline);
          for (j=0; j<nline; j++) {
            ptr = buf;
            data->getValue(BRANCH_TAP,&tap_ratio,j);
            bool is_xform = data->getValue(TRANSFORMER_WINDV1,&rval,j);
            is_xform = is_xform && data->getValue(TRANSFORMER_WINDV2,&rval,j);
            if (!is_xform && (tap_ratio == 0.0 || tap_ratio == 1.0)) {
              data->getValue(BRANCH_FROMBUS,&idx1);
              data->getValue(BRANCH_TOBUS,&idx2);
              ival = 1;
              data->getValue(BRANCH_SWITCHED,&ival,j);
              if (ival == 1) {
                sprintf(ptr,"%d, %d,",idx1,idx2);
              } else {
                sprintf(ptr,"%d, %d,",idx2,idx1);
              }
              ptr += strlen(ptr);
              sval = " 1",
              data->getValue(BRANCH_CKT,&sval,j);
              sprintf(ptr," \'%s\',",sval.c_str());
              ptr += strlen(ptr);
              // No default for this value
              data->getValue(BRANCH_R,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              // No default for this value
              data->getValue(BRANCH_X,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_B,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              // No default for this value
              data->getValue(BRANCH_RATING_A,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_RATING_B,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_RATING_C,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_SHUNT_ADMTTNC_G1,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_SHUNT_ADMTTNC_B1,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_SHUNT_ADMTTNC_G2,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_SHUNT_ADMTTNC_B2,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              ival = 1;
              data->getValue(BRANCH_STATUS,&ival,j);
              sprintf(ptr,"%d,",ival);
              ptr += strlen(ptr);
              ival = 1;
              data->getValue(BRANCH_METER,&ival,j);
              sprintf(ptr,"%d,",ival);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(BRANCH_LENGTH,&rval,j);
              sprintf(ptr," %f",rval);
              ptr += strlen(ptr);
              if (data->getValue(BRANCH_O1,&ival,j) &&
                  data->getValue(BRANCH_F1,&rval,j)) {
                sprintf(ptr,", %d, %f",ival,rval);
                ptr += strlen(ptr);
              }
              if (data->getValue(BRANCH_O2,&ival,j) &&
                  data->getValue(BRANCH_F2,&rval,j)) {
                sprintf(ptr,", %d, %f",ival,rval);
                ptr += strlen(ptr);
              }
              if (data->getValue(BRANCH_O3,&ival,j) &&
                  data->getValue(BRANCH_F3,&rval,j)) {
                sprintf(ptr,", %d, %f",ival,rval);
                ptr += strlen(ptr);
              }
              if (data->getValue(BRANCH_O4,&ival,j) &&
                  data->getValue(BRANCH_F4,&rval,j)) {
                sprintf(ptr,", %d, %f",ival,rval);
                ptr += strlen(ptr);
              }
              sprintf(ptr,"\n");
              text_line text;
              strcpy(text.text,buf);
              text.global_idx = p_network->getGlobalBranchIndex(i);
              text.device_idx = j;
              text_data.push_back(text);
            }
          }
        }
      }
      if (me == 0) {
        fout << "0 / END GENERATOR DATA, BEGIN LINE DATA" << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTLINE33_HPP_ */
