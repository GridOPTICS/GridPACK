/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * export_xform33.hpp
 *
 * This class exports transformer data in PSS/E v33 format
 *
 *  Created on: April 17, 2019
 *      Author: Bruce Palmer
 */

#ifndef EXPORTXFORM33_HPP_
#define EXPORTXFORM33_HPP_

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
class ExportXform33
{
  public:

    /**
     * Constructor
     */
    explicit ExportXform33(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~ExportXform33(){}

    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeXformBlock(std::ofstream &fout)
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
            if (is_xform || (tap_ratio != 0.0 && tap_ratio != 1.0)) {
              // Note that three winding buses have already been converted to
              // regular lines and buses. May want to revisit this issue in
              // the future
              data->getValue(BRANCH_FROMBUS,&idx1);
              data->getValue(BRANCH_TOBUS,&idx2);
              ival = 1;
              data->getValue(BRANCH_SWITCHED,&ival,j);
              // First line
              if (ival == 1) {
                sprintf(ptr,"%d, %d, 0,",idx1,idx2);
              } else {
                sprintf(ptr,"%d, %d, 0,",idx2,idx1);
              }
              ptr += strlen(ptr);
              sval = " 1";
              data->getValue(BRANCH_CKT,&sval,j);
              sprintf(ptr," \'%s\',",sval.c_str());
              ptr += strlen(ptr);
              // Assume default value for CW
              sprintf(ptr," 1,");
              ptr += strlen(ptr);
              // Assume default value for CZ
              sprintf(ptr," 1,");
              ptr += strlen(ptr);
              // Assume default value for CM
              sprintf(ptr," 1,");
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_MAG1,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_MAG2,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              ival = 2;
              data->getValue(TRANSFORMER_NMETR,&ival,j);
              sprintf(ptr," %d,",ival);
              ptr += strlen(ptr);
              ival = 1;
              sprintf(ptr," \'            \',");
              ptr += strlen(ptr);
              data->getValue(BRANCH_STATUS,&ival,j);
              sprintf(ptr," %d,",ival);
              ptr += strlen(ptr);
              ival = 1;
              data->getValue(BRANCH_O1,&ival,j);
              rval = 0.0;
              data->getValue(BRANCH_F1,&rval,j);
              sprintf(ptr," %d, %f,",ival,rval);
              ptr += strlen(ptr);
              ival = 0;
              data->getValue(BRANCH_O2,&ival,j);
              rval = 0.0;
              data->getValue(BRANCH_F2,&rval,j);
              sprintf(ptr," %d, %f,",ival,rval);
              ptr += strlen(ptr);
              ival = 0;
              data->getValue(BRANCH_O3,&ival,j);
              rval = 0.0;
              data->getValue(BRANCH_F3,&rval,j);
              sprintf(ptr," %d, %f,",ival,rval);
              ptr += strlen(ptr);
              ival = 0;
              data->getValue(BRANCH_O4,&ival,j);
              rval = 0.0;
              data->getValue(BRANCH_F4,&rval,j);
              sprintf(ptr," %d, %f,",ival,rval);
              ptr += strlen(ptr);
              // VECGRP just defaults to 12 blanks
              sprintf(ptr,"\'            \'\n");
              ptr += strlen(ptr);

              // Second line
              rval = 0.0;
              if (!data->getValue(TRANSFORMER_R1_2,&rval,j)) {
                data->getValue(BRANCH_R,&rval,j);
              }
              sprintf(ptr,"%f,",rval);
              ptr += strlen(ptr);
              if (!data->getValue(TRANSFORMER_X1_2,&rval,j)) {
                data->getValue(BRANCH_X,&rval,j);
              }
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              if (!data->getValue(TRANSFORMER_SBASE1_2,&rval,j)) {
                data->getValue(CASE_SBASE,&rval);
              }
              sprintf(ptr," %f,\n",rval);
              ptr += strlen(ptr);

              // Third line
              rval = 1.0;
              if (!data->getValue(TRANSFORMER_WINDV1,&rval,j)) {
                data->getValue(BRANCH_TAP,&rval,j);
              } 
              sprintf(ptr,"%f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_NOMV1,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_ANG1,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
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
              ival = 0;
              data->getValue(TRANSFORMER_CODE1,&ival,j);
              sprintf(ptr," %d,",ival);
              ptr += strlen(ptr);
              ival = 0;
              data->getValue(TRANSFORMER_CONT1,&ival,j);
              sprintf(ptr," %d,",ival);
              ptr += strlen(ptr);
              rval = 1.1;
              data->getValue(TRANSFORMER_RMA,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.9;
              data->getValue(TRANSFORMER_RMI,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 1.1;
              data->getValue(TRANSFORMER_VMA,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.9;
              data->getValue(TRANSFORMER_VMI,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              ival = 33;
              data->getValue(TRANSFORMER_NTP,&ival,j);
              sprintf(ptr," %d,",ival);
              ptr += strlen(ptr);
              ival = 0;
              data->getValue(TRANSFORMER_TAB,&ival,j);
              sprintf(ptr," %d,",ival);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_CR,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_CX,&rval,j);
              sprintf(ptr," %f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_CNXA,&rval,j);
              sprintf(ptr," %f\n",rval);
              ptr += strlen(ptr);

              // Fourth line
              rval = 1.0;
              data->getValue(TRANSFORMER_WINDV2,&rval,j);
              sprintf(ptr,"%f,",rval);
              ptr += strlen(ptr);
              rval = 0.0;
              data->getValue(TRANSFORMER_NOMV1,&rval,j);
              sprintf(ptr," %f",rval);
              ptr += strlen(ptr);
              sprintf(ptr,"\n");
              // Finish up
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
        fout << "0 / END LINE DATA, BEGIN TRANSFORMER DATA" << std::endl;
      }
      exprt.writeDataBlock(fout, text_data);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* EXPORTXFORM33_HPP_ */
