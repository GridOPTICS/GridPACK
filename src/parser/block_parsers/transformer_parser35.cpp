/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * transformer_parser35.cpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#include "transformer_parser35.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::TransformerParser35::TransformerParser35(
    std::map<int,int> *bus_map,
    std::map<std::string,int> *name_map,
    std::map<std::pair<int, int>, int> *branch_map) :
    gridpack::parser::BaseBlockParser(
      bus_map, name_map, branch_map)
{
}


/**
 * Simple Destructor
 */
gridpack::parser::TransformerParser35::~TransformerParser35(void)
{
}

/**
 * parse transformer block
 * @param stream input stream that feeds lines from RAW file
 * @param p_busData vector of bus data collection objects to store parameters
 *             from RAW file for transformers
 * @param p_branchData vector of branch data collection objects to store
 *             parameters from RAW file for transformers
 * @param case_sbase sbase parameter
 * @param maxBusIndex maximum value of bus index for any bus
 in system
 */
void gridpack::parser::TransformerParser35::parse(
    gridpack::stream::InputStream &stream,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &p_busData,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &p_branchData,
    double p_case_sbase,
    int p_maxBusIndex)
{
  std::string          line;

  stream.nextLine(line); //this should be the first line of the block

  std::pair<int, int>   branch_pair;

  // Save current number of branches
  int index = p_branchData.size();

  bool wind3X = true;

  while(test_end(line)) {
    std::vector<std::string>  split_line;
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    this->cleanComment(line);
    split_line = this->splitPSSELine(line);
    int o_idx1, o_idx2;
    o_idx1 = getBusIndex(split_line[0]);
    o_idx2 = getBusIndex(split_line[1]);

    // Check for 2-winding or 3-winding transformer
    int k = getBusIndex(split_line[2]);
    if (k != 0) {
      if (wind3X) {
        int o_idx3 = k;
        stream.nextLine(line);
        std::vector<std::string>  split_line2;
        if (check_comment(line)) {
          stream.nextLine(line);
          continue;
        }
        this->cleanComment(line);
        split_line2 = this->splitPSSELine(line);
        // Check to see if transformer is active
        int stat;
        stat = atoi(split_line[11].c_str());
        if (split_line2.size() < 4 || stat == 0) {
          stream.nextLine(line);
          stream.nextLine(line);
          stream.nextLine(line);
          stream.nextLine(line);
          continue;
        }
        // Get internal indices corresponding to buses 1,2,3
        int l_idx1, l_idx2, l_idx3;
        std::map<int,int>::iterator it;
        it = p_busMap->find(o_idx1);
        if (it != p_busMap->end()) {
          l_idx1 = it->second;
        } else {
          printf("No match found for bus %s\n",split_line[0].c_str());
        }
        it = p_busMap->find(o_idx2);
        if (it != p_busMap->end()) {
          l_idx2 = it->second;
        } else {
          printf("No match found for bus %s\n",split_line[1].c_str());
        }
        it = p_busMap->find(o_idx3);
        if (it != p_busMap->end()) {
          l_idx3 = it->second;
        } else {
          printf("No match found for bus %s\n",split_line[2].c_str());
        }
        // Create a new bus and three new branches. No need to check
        // previous branches to see if they match since they are all
        // linked to new bus. Also don't need to add these to the branch
        // map data structure since they cannot match any regular branch
        // lines or transformers
        boost::shared_ptr<gridpack::component::DataCollection>
          data(new gridpack::component::DataCollection);
        int n_idx = p_busData.size();
        p_busData.push_back(data);
        p_maxBusIndex++;
        data->addValue(BUS_NUMBER,p_maxBusIndex);
        char cbuf[128];
        sprintf(cbuf,"DUMMY_BUS-%d-%d-%d",o_idx1,o_idx2,o_idx3);
        data->addValue(BUS_NAME,cbuf);
        data->addValue(BUS_BASEKV,0.0);
        data->addValue(BUS_TYPE,1);
        int ival;
        p_busData[l_idx1]->getValue(BUS_AREA,&ival);
        data->addValue(BUS_AREA,ival);
        p_busData[l_idx1]->getValue(BUS_OWNER,&ival);
        data->addValue(BUS_OWNER, ival);
        double rval = 0.0;
        double rvol;
        p_busData[l_idx1]->getValue(BUS_VOLTAGE_MAG,&rvol);
        rval += rvol;
        p_busData[l_idx2]->getValue(BUS_VOLTAGE_MAG,&rvol);
        rval += rvol;
        p_busData[l_idx3]->getValue(BUS_VOLTAGE_MAG,&rvol);
        rval += rvol;
        rval = rval/3.0;
        rval = 1.0;
        data->addValue(BUS_VOLTAGE_MAG,rval);
        rval = 0.0;
        p_busData[l_idx1]->getValue(BUS_VOLTAGE_ANG,&rvol);
        rval += rvol;
        p_busData[l_idx2]->getValue(BUS_VOLTAGE_ANG,&rvol);
        rval += rvol;
        p_busData[l_idx3]->getValue(BUS_VOLTAGE_ANG,&rvol);
        rval += rvol;
        rval = rval/3.0;
        rval = 0.0;
        data->addValue(BUS_VOLTAGE_ANG,rval);

        // parse remainder of line 1
        double mag1, mag2;
        mag1 = atof(split_line[7].c_str());
        mag2 = atof(split_line[8].c_str());
        // Clean up 2 character tag
        gridpack::utility::StringUtils util;
        std::string tag = util.clean2Char(split_line[3]);

        // parse line 2
        double r12, r23, r31, x12, x23, x31, sb12, sb23, sb31;
        double r1, r2, r3, x1, x2, x3, b1, b2, b3;
        r12 = atof(split_line2[0].c_str());
        x12 = atof(split_line2[1].c_str());
        sb12 = atof(split_line2[2].c_str());
        r23 = atof(split_line2[3].c_str());
        x23 = atof(split_line2[4].c_str());
        sb23 = atof(split_line2[5].c_str());
        r31 = atof(split_line2[6].c_str());
        x31 = atof(split_line2[7].c_str());
        sb31 = atof(split_line2[8].c_str());
        r1 = 0.5*(r12+r31-r23);
        x1 = 0.5*(x12+x31-x23);
        b1 = 0.0;
        r2 = 0.5*(r12+r23-r31);
        x2 = 0.5*(x12+x23-x31);
        b2 = 0.0;
        r3 = 0.5*(r23+r31-r12);
        x3 = 0.5*(x23+x31-x12);
        b3 = 0.0;
        // correct R and X values, if necessary
        if (sb12 != p_case_sbase && sb12 == 0.0) {
          r1 = r1*p_case_sbase/sb12;
          x1 = x1*p_case_sbase/sb12;
        }
        if (sb23 != p_case_sbase && sb23 == 0.0) {
          r2 = r2*p_case_sbase/sb23;
          x2 = x2*p_case_sbase/sb23;
        }
        if (sb31 != p_case_sbase && sb31 == 0.0) {
          r3 = r3*p_case_sbase/sb31;
          x3 = x3*p_case_sbase/sb31;
        }

        // create branch between new bus and bus I
        int index = p_branchData.size();
        boost::shared_ptr<gridpack::component::DataCollection>
          data1(new gridpack::component::DataCollection);
        p_branchData.push_back(data1);
        stream.nextLine(line);
        std::vector<std::string> split_line3;
        if (check_comment(line)) {
          stream.nextLine(line);
          continue;
        }
        this->cleanComment(line);
        split_line3 = this->splitPSSELine(line);
        double windv, ang, rate[12];
        parse3WindXForm(split_line3, &windv, &ang, rate);
        data1->addValue(BRANCH_INDEX,index);
        data1->addValue(BRANCH_FROMBUS,o_idx1);
        data1->addValue(BRANCH_TOBUS,p_maxBusIndex);
        data1->addValue(BRANCH_NUM_ELEMENTS,1);
        data1->addValue(BRANCH_CKT,tag.c_str(),0);
        data1->addValue(BRANCH_R,r1,0);
        data1->addValue(BRANCH_X,x1,0);
        data1->addValue(BRANCH_B,b1,0);
        data1->addValue(BRANCH_RATING_A,rate[0],0);
        data1->addValue(BRANCH_RATING_B,rate[1],0);
        data1->addValue(BRANCH_RATING_C,rate[2],0);
        data1->addValue(BRANCH_RATE1,rate[0],0);
        data1->addValue(BRANCH_RATE2,rate[1],0);
        data1->addValue(BRANCH_RATE3,rate[2],0);
        data1->addValue(BRANCH_RATE4,rate[3],0);
        data1->addValue(BRANCH_RATE5,rate[4],0);
        data1->addValue(BRANCH_RATE6,rate[5],0);
        data1->addValue(BRANCH_RATE7,rate[6],0);
        data1->addValue(BRANCH_RATE8,rate[7],0);
        data1->addValue(BRANCH_RATE9,rate[8],0);
        data1->addValue(BRANCH_RATE10,rate[9],0);
        data1->addValue(BRANCH_RATE11,rate[10],0);
        data1->addValue(BRANCH_RATE12,rate[11],0);
        data1->addValue(BRANCH_TAP,windv,0);
        data1->addValue(BRANCH_SHIFT,ang,0);
        data1->addValue(BRANCH_SWITCHED,false,0);
        data1->addValue(BRANCH_SHUNT_ADMTTNC_G1,mag1,0);
        data1->addValue(BRANCH_SHUNT_ADMTTNC_B1,mag2,0);
        data1->addValue(BRANCH_SHUNT_ADMTTNC_G2,0.0,0);
        data1->addValue(BRANCH_SHUNT_ADMTTNC_B2,0.0,0);
        if (stat == 1 || stat == 2 || stat == 3) {
          data1->addValue(BRANCH_STATUS,1,0);
        } else {
          data1->addValue(BRANCH_STATUS,0,0);
        }

        // create branch between new bus and bus J
        index++;
        boost::shared_ptr<gridpack::component::DataCollection>
          data2(new gridpack::component::DataCollection);
        p_branchData.push_back(data2);
        stream.nextLine(line);
        std::vector<std::string> split_line4;
        if (check_comment(line)) {
          stream.nextLine(line);
          continue;
        }
        this->cleanComment(line);
        split_line4 = this->splitPSSELine(line);
        parse3WindXForm(split_line4, &windv, &ang, rate);
        data2->addValue(BRANCH_INDEX,index);
        data2->addValue(BRANCH_FROMBUS,o_idx2);
        data2->addValue(BRANCH_TOBUS,p_maxBusIndex);
        data2->addValue(BRANCH_NUM_ELEMENTS,1);
        data2->addValue(BRANCH_CKT,tag.c_str(),0);
        data2->addValue(BRANCH_R,r2,0);
        data2->addValue(BRANCH_X,x2,0);
        data2->addValue(BRANCH_B,b2,0);
        data1->addValue(BRANCH_RATING_A,rate[0],0);
        data1->addValue(BRANCH_RATING_B,rate[1],0);
        data1->addValue(BRANCH_RATING_C,rate[2],0);
        data1->addValue(BRANCH_RATE1,rate[0],0);
        data1->addValue(BRANCH_RATE2,rate[1],0);
        data1->addValue(BRANCH_RATE3,rate[2],0);
        data1->addValue(BRANCH_RATE4,rate[3],0);
        data1->addValue(BRANCH_RATE5,rate[4],0);
        data1->addValue(BRANCH_RATE6,rate[5],0);
        data1->addValue(BRANCH_RATE7,rate[6],0);
        data1->addValue(BRANCH_RATE8,rate[7],0);
        data1->addValue(BRANCH_RATE9,rate[8],0);
        data1->addValue(BRANCH_RATE10,rate[9],0);
        data1->addValue(BRANCH_RATE11,rate[10],0);
        data1->addValue(BRANCH_RATE12,rate[11],0);
        data2->addValue(BRANCH_TAP,windv,0);
        data2->addValue(BRANCH_SHIFT,ang,0);
        data2->addValue(BRANCH_SWITCHED,false,0);
        data2->addValue(BRANCH_SHUNT_ADMTTNC_G1,0.0,0);
        data2->addValue(BRANCH_SHUNT_ADMTTNC_B1,0.0,0);
        data2->addValue(BRANCH_SHUNT_ADMTTNC_G2,0.0,0);
        data2->addValue(BRANCH_SHUNT_ADMTTNC_B2,0.0,0);
        if (stat == 1 || stat == 3 || stat == 4) {
          data2->addValue(BRANCH_STATUS,1,0);
        } else {
          data2->addValue(BRANCH_STATUS,0,0);
        }

        // create branch between new bus and bus K
        index++;
        boost::shared_ptr<gridpack::component::DataCollection>
          data3(new gridpack::component::DataCollection);
        p_branchData.push_back(data3);
        stream.nextLine(line);
        std::vector<std::string> split_line5;
        if (check_comment(line)) {
          stream.nextLine(line);
          continue;
        }
        this->cleanComment(line);
        split_line5 = this->splitPSSELine(line);
        parse3WindXForm(split_line5, &windv, &ang, rate);
        data3->addValue(BRANCH_INDEX,index);
        data3->addValue(BRANCH_FROMBUS,o_idx3);
        data3->addValue(BRANCH_TOBUS,p_maxBusIndex);
        data3->addValue(BRANCH_NUM_ELEMENTS,1);
        data3->addValue(BRANCH_CKT,tag.c_str(),0);
        data3->addValue(BRANCH_R,r3,0);
        data3->addValue(BRANCH_X,x3,0);
        data3->addValue(BRANCH_B,b3,0);
        data1->addValue(BRANCH_RATING_A,rate[0],0);
        data1->addValue(BRANCH_RATING_B,rate[1],0);
        data1->addValue(BRANCH_RATING_C,rate[2],0);
        data1->addValue(BRANCH_RATE1,rate[0],0);
        data1->addValue(BRANCH_RATE2,rate[1],0);
        data1->addValue(BRANCH_RATE3,rate[2],0);
        data1->addValue(BRANCH_RATE4,rate[3],0);
        data1->addValue(BRANCH_RATE5,rate[4],0);
        data1->addValue(BRANCH_RATE6,rate[5],0);
        data1->addValue(BRANCH_RATE7,rate[6],0);
        data1->addValue(BRANCH_RATE8,rate[7],0);
        data1->addValue(BRANCH_RATE9,rate[8],0);
        data1->addValue(BRANCH_RATE10,rate[9],0);
        data1->addValue(BRANCH_RATE11,rate[10],0);
        data1->addValue(BRANCH_RATE12,rate[11],0);
        data3->addValue(BRANCH_TAP,windv,0);
        data3->addValue(BRANCH_SHIFT,ang,0);
        data3->addValue(BRANCH_SWITCHED,false,0);
        data3->addValue(BRANCH_SHUNT_ADMTTNC_G1,0.0,0);
        data3->addValue(BRANCH_SHUNT_ADMTTNC_B1,0.0,0);
        data3->addValue(BRANCH_SHUNT_ADMTTNC_G2,0.0,0);
        data3->addValue(BRANCH_SHUNT_ADMTTNC_B2,0.0,0);
        if (stat == 1 || stat == 2 || stat == 4) {
          data3->addValue(BRANCH_STATUS,1,0);
        } else {
          data3->addValue(BRANCH_STATUS,0,0);
        }
      } else {
        // Skip 3-winding transformers (for now)
        stream.nextLine(line);
        stream.nextLine(line);
        stream.nextLine(line);
        stream.nextLine(line);
        continue;
      }
    } else {
      int ntoken = split_line.size();
      stream.nextLine(line);
      std::vector<std::string>  split_line2;
      if (check_comment(line)) {
        stream.nextLine(line);
        continue;
      }
      this->cleanComment(line);
      split_line2 = this->splitPSSELine(line);

      stream.nextLine(line);
      std::vector<std::string>  split_line3;
      if (check_comment(line)) {
        stream.nextLine(line);
        continue;
      }
      this->cleanComment(line);
      split_line3 = this->splitPSSELine(line);

      stream.nextLine(line);
      std::vector<std::string>  split_line4;
      if (check_comment(line)) {
        stream.nextLine(line);
        continue;
      }
      this->cleanComment(line);
      split_line4 = this->splitPSSELine(line);
      // find branch corresponding to this transformer line. If it doesn't
      // exist, create one
      int l_idx = 0;
      branch_pair = std::pair<int,int>(o_idx1, o_idx2);
      std::map<std::pair<int, int>, int>::iterator it;
      it = p_branchMap->find(branch_pair);
      bool switched = false;
      int nelems;
      if (it != p_branchMap->end()) {
        l_idx = it->second;
        p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
      } else {
        // Check to see if from and to buses have been switched
        std::pair<int, int> new_branch_pair;
        new_branch_pair = std::pair<int,int>(o_idx2, o_idx1);
        it = p_branchMap->find(new_branch_pair);
        if (it != p_branchMap->end()) {
          l_idx = it->second;
          p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
          switched = true;
        } else {
          boost::shared_ptr<gridpack::component::DataCollection>
            data(new gridpack::component::DataCollection);
          l_idx = p_branchData.size();
          p_branchData.push_back(data);
          std::pair<std::pair<int,int>,int> item;
          item = std::pair<std::pair<int,int>,int>(branch_pair,l_idx);
          p_branchMap->insert(item);
          nelems = 0;
          p_branchData[l_idx]->addValue(BRANCH_NUM_ELEMENTS,nelems);
        }
      }

      // If nelems=0 then this is new branch. Add parameters defining
      // branch
      if (nelems == 0) {
        // BRANCH_INDEX                   integer
        p_branchData[l_idx]->addValue(BRANCH_INDEX,
            index);

        // BRANCH_FROMBUS            "I"  integer
        p_branchData[l_idx]->addValue(BRANCH_FROMBUS,
            o_idx1);

        // BRANCH_TOBUS              "J"  integer
        p_branchData[l_idx]->addValue(BRANCH_TOBUS,
            o_idx2);

        // add pair to branch map
        p_branchMap->insert(std::pair<std::pair<int,int>,int >
            (branch_pair, index));
        index++;
      }

      // Clean up 2 character tag
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[3]);
      // BRANCH_CKT          "CKT"                 character
      p_branchData[l_idx]->addValue(BRANCH_CKT, tag.c_str(), nelems);

      // Add remaining parameters from line 1
      /*
       * type: integer
       * TRANSFORMER_CW
       */
      int cw = atoi(split_line[4].c_str()); 
      p_branchData[l_idx]->addValue(TRANSFORMER_CW,
         cw,nelems);

      /*
       * type: integer
       * TRANSFORMER_CZ
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_CZ,
          atoi(split_line[5].c_str()),nelems);

      /*
       * type: integer
       * TRANSFORMER_CM
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_CM,
          atoi(split_line[6].c_str()),nelems);

      /*
       * type: float
       * TRANSFORMER_MAG1
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_MAG1,
          atof(split_line[7].c_str()),nelems);

      /*
       * type: float
       * TRANSFORMER_MAG2
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_MAG2,
          atof(split_line[8].c_str()),nelems);

      p_branchData[l_idx]->addValue(BRANCH_B,
          atof(split_line[8].c_str()),nelems);										   

      /*
       * type: integer
       * BRANCH_STATUS
       */
      p_branchData[l_idx]->addValue(BRANCH_STATUS,
          atoi(split_line[11].c_str()),nelems);

      /**
       * type: integer
       * TRANSFORMER_NMETR
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_NMETR,
          atoi(split_line[9].c_str()),nelems);

      /*
       * type: integer
       * BRANCH_O1
       */
      if (ntoken > 12) p_branchData[l_idx]->addValue(BRANCH_O1,
          atoi(split_line[12].c_str()), nelems);

      /*
       * type: float
       * BRANCH_F1
       */
      if (ntoken > 13) p_branchData[l_idx]->addValue(BRANCH_F1,
          atoi(split_line[13].c_str()), nelems);

      /*
       * type: integer
       * BRANCH_O2
       */
      if (ntoken > 14) p_branchData[l_idx]->addValue(BRANCH_O2,
          atoi(split_line[14].c_str()), nelems);

      /*
       * type: float
       * BRANCH_F2
       */
      if (ntoken > 15) p_branchData[l_idx]->addValue(BRANCH_F2,
          atoi(split_line[15].c_str()), nelems);

      /*
       * type: integer
       * BRANCH_O3
       */
      if (ntoken > 16) p_branchData[l_idx]->addValue(BRANCH_O3,
          atoi(split_line[16].c_str()), nelems);

      /*
       * type: float
       * BRANCH_F3
       */
      if (ntoken > 17) p_branchData[l_idx]->addValue(BRANCH_F3,
          atoi(split_line[17].c_str()), nelems);

      /*
       * type: integer
       * BRANCH_O4
       */
      if (ntoken > 18) p_branchData[l_idx]->addValue(BRANCH_O4,
          atoi(split_line[18].c_str()), nelems);

      /*
       * type: float
       * BRANCH_F4
       */
      if (ntoken > 19) p_branchData[l_idx]->addValue(BRANCH_F4,
          atoi(split_line[19].c_str()), nelems);


      // Add parameters from line 2
      /*
       * type: float
       * SBASE2
       */
      double sbase2 = atof(split_line2[2].c_str());
      p_branchData[l_idx]->addValue(TRANSFORMER_SBASE1_2,sbase2,nelems);



      // Add parameters from line 3
      /*
       * type: float
       * BRANCH_TAP: This is the ratio of WINDV1 and WINDV2
       */
      ntoken = split_line3.size();
      double windv1 = atof(split_line3[0].c_str());
      double windv2 = atof(split_line4[0].c_str());
      if(cw == 2) {
	double nomv1 = atof(split_line3[1].c_str());
	double nomv2 = atof(split_line4[1].c_str());
	windv1 = windv1/nomv1;
	windv2 = windv2/nomv2;
      }
      double tap = windv1/windv2;
      p_branchData[l_idx]->addValue(BRANCH_TAP,tap,nelems);
      p_branchData[l_idx]->addValue(TRANSFORMER_WINDV1,windv1,nelems);
      p_branchData[l_idx]->addValue(TRANSFORMER_WINDV2,windv2,nelems);


      /*
       * type: float
       * BRANCH_R
       */
      double rval = atof(split_line2[0].c_str());
      p_branchData[l_idx]->addValue(TRANSFORMER_R1_2,rval,nelems);
      rval  = rval * windv2 * windv2; // need to consider the wnd2 ratio to the req of the transformer
      if (sbase2 == p_case_sbase || sbase2 == 0.0) {
        p_branchData[l_idx]->addValue(BRANCH_R,rval,nelems);
      } else {
        rval = rval*p_case_sbase/sbase2;
        p_branchData[l_idx]->addValue(BRANCH_R,rval,nelems);
      }


      /*
       * type: float
       * BRANCH_X
       */
      rval = atof(split_line2[1].c_str());
      p_branchData[l_idx]->addValue(TRANSFORMER_X1_2,rval,nelems);
      rval  = rval * windv2 * windv2; // need to consider the wnd2 ratio to the xeq of the transformer
      if (sbase2 == p_case_sbase || sbase2 == 0.0) {
        p_branchData[l_idx]->addValue(BRANCH_X,rval,nelems);
      } else {
        rval = rval*p_case_sbase/sbase2;
        p_branchData[l_idx]->addValue(BRANCH_X,rval,nelems);
      }

      // Add parameters from line 3

      /*
       * type: float
       * BRANCH_SHIFT
       */
      p_branchData[l_idx]->addValue(BRANCH_SHIFT,
          atof(split_line3[2].c_str()),nelems);
      p_branchData[l_idx]->addValue(TRANSFORMER_ANG1,
          atof(split_line3[2].c_str()),nelems);

      /*
       * type: float
       * BRANCH_RATE1-12
       */
      p_branchData[l_idx]->addValue(BRANCH_RATE1,
          atof(split_line3[3].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE2,
          atof(split_line3[4].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE3,
          atof(split_line3[5].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE4,
          atof(split_line3[6].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE5,
          atof(split_line3[7].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE6,
          atof(split_line3[8].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE7,
          atof(split_line3[9].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE8,
          atof(split_line3[10].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE9,
          atof(split_line3[11].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE10,
          atof(split_line3[12].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE11,
          atof(split_line3[13].c_str()),nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATE11,
          atof(split_line3[14].c_str()),nelems);

      /*
       * type: integer
       * TRANSFORMER_CODE1
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_CODE1,
          atoi(split_line3[15].c_str()),nelems);

      /*
       * type: float
       * TRANSFORMER_RMA
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_RMA,
          atof(split_line3[18].c_str()),nelems);

      /*
       * type: float
       * TRANSFORMER_RMI
       */
      p_branchData[l_idx]->addValue(TRANSFORMER_RMI,
          atof(split_line3[19].c_str()),nelems);

      /*
       * type: float
       * TRANSFORMER_VMA
       */
      if (ntoken > 20) {
        p_branchData[l_idx]->addValue(TRANSFORMER_VMA,
            atof(split_line3[20].c_str()),nelems);
      }

      /* ignore line 4 for now */

      /*
       * type: float
       * TRANSFORMER_VMI
       */
      if (ntoken > 21) {
        p_branchData[l_idx]->addValue(TRANSFORMER_VMI,
            atof(split_line3[21].c_str()),nelems);
      }

      /*
       * type: integer
       * TRANSFORMER_NPT
       */
      if (ntoken > 22) {
        p_branchData[l_idx]->addValue(TRANSFORMER_NTP,
            atoi(split_line3[22].c_str()),nelems);
      }

      /*
       * type: integer
       * TRANSFORMER_TAB
       */
      if (ntoken > 23) {
        p_branchData[l_idx]->addValue(TRANSFORMER_TAB,
            atoi(split_line3[23].c_str()),nelems);
      }

      /*
       * type: float
       * TRANSFORMER_CR
       */
      if (ntoken > 24) {
        p_branchData[l_idx]->addValue(TRANSFORMER_CR,
            atof(split_line3[24].c_str()),nelems);
      }

      /*
       * type: float
       * TRANSFORMER_CI
       */
      if (ntoken > 25) {
        p_branchData[l_idx]->addValue(TRANSFORMER_CX,
            atof(split_line3[25].c_str()),nelems);
      }

      /*
       * type: float
       * TRANSFORMER_CNXA
       */
      if (ntoken > 26) {
        p_branchData[l_idx]->addValue(TRANSFORMER_CNXA,
            atof(split_line3[26].c_str()),nelems);
      }

      nelems++;
      p_branchData[l_idx]->setValue(BRANCH_NUM_ELEMENTS,nelems);
    }
    stream.nextLine(line);
  }
}
