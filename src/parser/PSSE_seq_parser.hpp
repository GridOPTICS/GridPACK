/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: December 30, 2014
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef PSSE_SEQ_PARSER_HPP_
#define PSSE_SEQ_PARSER_HPP_

#define OLD_MAP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif


#include "gridpack/component/base_component.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/stream/input_stream.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "gridpack/factory/base_factory.hpp"

#define TERM_CHAR '0'

namespace gridpack {
namespace parser {

template <class _network>
class PSSE_seq_parser : public BaseParser<_network>
{
  public:

    /**
     * Constructor
     */
    explicit PSSE_seq_parser(boost::shared_ptr<_network> network)
      : p_network(network)
    {
      this->setNetwork(network);
      p_timer = gridpack::utility::CoarseTimer::instance();
    }

    /**
     * Destructor
     */
    virtual ~PSSE_seq_parser(){}

    /**
     * parse a sequence file
     * @param fileName name of file
     */
    void parse(const std::string &fileName)
    {
      std::string ext = getExtension(fileName);
      if (ext == "seq") {
        getSeqData(fileName);
        sortSeqData();
      } else {
        printf("Unknown file extension %s\n",ext.c_str());
      }
    }

  protected:

    /* ************************************************************************
     **************************************************************************
     ***** PROTECTED SCOPE
     **************************************************************************
     *********************************************************************** */

    /**
     * Assign network to internal network pointer variable
     */
    void setNetwork(boost::shared_ptr<_network> network)
    {
      p_network = network;
      BaseParser<_network>::setNetwork(network);
    }

    /**
     * This routine opens up a .seq file. Assume that a .raw file has already been parsed
     */
    void getSeqData(const std::string & fileName)
    {
      int t_ds = p_timer->createCategory("PSSE_seq_parser:getSeqData");
      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      if (me == 0) {
        p_input_stream.openFile(fileName);
        if (!p_input_stream.isOpen()) {
          p_timer->stop(t_ds);
          return;
        }
        /* read first line and extract change code */
        std::string line;
        p_input_stream.nextLine(line);
        char sep = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&sep);
        p_change_code = atoi(split_line[0].c_str());
        /* extract data for generators */
        findPosSeqImpData();
        findNegSeqImpData();
        findZeroSeqImpData();
        /* extract data for buses */
        findNegSeqShuntData();
        findZeroSeqShuntData();
        /* extract data for branches */
        findZeroSeqNonTransformerData();
        findZeroSeqMutualImpData();
        findZeroSeqTransformerData();
        /* close input */
        p_input_stream.close();
      }
      /* distribute data to correct processors */
      p_timer->stop(t_ds);
    }

    // Data structure to hold generator params
    struct gen_params{
      int bus_id;     // ID of bus that owns device
      char gen_id[3]; // Generator ID
      bool use_pos;
      bool use_neg;
      bool use_zero;
      double zrpos;
      double zxpos;
      double zrneg;
      double zxneg;
      double rzero;
      double xzero;
    };

    // Data structure to bus parameters
    struct bus_params{
      int bus_id; // ID of bus
      double gneg;
      double bneg;
      bool use_neg;
      bool use_zero;
      double gzero;
      double bzero;
    };

    // Data structure to hold branche parameters
    struct branch_params{
      int from_bus; // ID of from bus
      int to_bus;   // ID of to bus
      char branch_id[3]; // Branch identifier
      bool is_xform;
      bool three_winding;
      double rlinz;
      double xlinz;
      double bchz;
      double gi;
      double bi;
      double gj;
      double bj;
      int k; // ID of bus to which another winding of transformer is connected
      int cc; // winding connection code
      double rg;
      double xg;
      double r1;
      double x1;
      double rg2;
      double xg2;
      double r2;
      double x2;
      double r3;
      double x3;
    };


    // Extract extension from file name and convert it to lower case
    std::string getExtension(const std::string file)
    {
      std::string ret;
      std::string line = file;
      int ntok1 = line.find('.',0);
      if (ntok1 == std::string::npos) return ret;
      int nsav = ntok1;
      while(ntok1 != std::string::npos) {
        ntok1 = line.find('.',ntok1+1);
        if (ntok1 != std::string::npos) nsav=ntok1;
      }
      ntok1 = nsav;
      ntok1++;
      int ntok2 = line.find(' ',ntok1);
      if (ntok2 == std::string::npos)
        ntok2 = line.size();
      // get extension
      ret = line.substr(ntok1,ntok2-ntok1);
      // convert all characters to lower case 
      int size = ret.size();
      int i;
      for (i=0; i<size; i++) {
        if (isalpha(ret[i])) {
          ret[i] = tolower(ret[i]);
        }
      }
      return ret;
    }

    /**
     * Extract positive sequence impedence data and put it in gen_params structs
     */
    void findPosSeqImpData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        /* Check if entry already exists for this generator */
        int bus_id = atoi(split_line[0].c_str());
        std::string tag = p_util.clean2Char(split_line[1]);
        std::pair<int,std::string> gen = std::pair<int,std::string>(bus_id,tag);
        std::map<std::pair<int,std::string>,int>::iterator it = p_genIDs.find(gen);
        if (it != p_genIDs.end()) {
          int idx = it->second;
          p_gen_data[idx].use_pos = true;
          p_gen_data[idx].zrpos = atof(split_line[2].c_str());
          p_gen_data[idx].zxpos = atof(split_line[3].c_str());
        } else {
          gen_params data;
          data.bus_id = bus_id;
          data.use_pos = true;
          data.use_neg = false;
          data.use_zero = false;
          strcpy(data.gen_id,tag.c_str());
          data.gen_id[2] = '\0';
          data.zrpos = atof(split_line[2].c_str());
          data.zxpos = atof(split_line[3].c_str());
          int len = p_gen_data.size();
          p_gen_data.push_back(data);
          p_genIDs.insert(std::pair<std::pair<int,std::string>,int>(gen,len));
        }
      }
    }

    /**
     * Extract negative sequence impedence data and put it in gen_params structs
     */
    void findNegSeqImpData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        /* Check if entry already exists for this generator */
        int bus_id = atoi(split_line[0].c_str());
        std::string tag = p_util.clean2Char(split_line[1]);
        std::pair<int,std::string> gen = std::pair<int,std::string>(bus_id,tag);
        std::map<std::pair<int,std::string>,int>::iterator it = p_genIDs.find(gen);
        if (it != p_genIDs.end()) {
          int idx = it->second;
          p_gen_data[idx].use_neg = true;
          p_gen_data[idx].zrneg = atof(split_line[2].c_str());
          p_gen_data[idx].zxneg = atof(split_line[3].c_str());
        } else {
          gen_params data;
          data.use_pos = false;
          data.use_neg = true;
          data.use_zero = false;
          data.bus_id = bus_id;
          strcpy(data.gen_id,tag.c_str());
          data.gen_id[2] = '\0';
          data.zrneg = atof(split_line[2].c_str());
          data.zxneg = atof(split_line[3].c_str());
          int len = p_gen_data.size();
          p_gen_data.push_back(data);
          p_genIDs.insert(std::pair<std::pair<int,std::string>,int>(gen,len));
        }
      }
    }

    /**
     * Extract zero sequence impedence data and put it in gen_params structs
     */
    void findZeroSeqImpData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        /* Check if entry already exists for this generator */
        int bus_id = atoi(split_line[0].c_str());
        std::string tag = p_util.clean2Char(split_line[1]);
        std::pair<int,std::string> gen = std::pair<int,std::string>(bus_id,tag);
        std::map<std::pair<int,std::string>,int>::iterator it = p_genIDs.find(gen);
        if (it != p_genIDs.end()) {
          int idx = it->second;
          p_gen_data[idx].use_zero = true;
          p_gen_data[idx].rzero = atof(split_line[2].c_str());
          p_gen_data[idx].xzero = atof(split_line[3].c_str());
        } else {
          gen_params data;
          data.use_pos = false;
          data.use_neg = false;
          data.use_zero = true;
          data.bus_id = bus_id;
          strcpy(data.gen_id,tag.c_str());
          data.gen_id[2] = '\0';
          data.rzero = atof(split_line[2].c_str());
          data.xzero = atof(split_line[3].c_str());
          int len = p_gen_data.size();
          p_gen_data.push_back(data);
          p_genIDs.insert(std::pair<std::pair<int,std::string>,int>(gen,len));
        }
      }
    }

    /**
     * Extract negative sequence shunt data and put it in bus_params structs
     */
    void findNegSeqShuntData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        /* Check if entry already exists for this bus */
        int bus_id = atoi(split_line[0].c_str());
        std::map<int,int>::iterator it = p_busIDs.find(bus_id);
        if (it != p_busIDs.end()) {
          int idx = it->second;
          p_bus_data[idx].use_neg = true;
          p_bus_data[idx].gneg = atof(split_line[1].c_str());
          p_bus_data[idx].bneg = atof(split_line[2].c_str());
        } else {
          bus_params data;
          data.use_neg = true;
          data.use_zero = false;
          data.bus_id = bus_id;
          data.gneg = atof(split_line[1].c_str());
          data.bneg = atof(split_line[2].c_str());
          int len = p_bus_data.size();
          p_bus_data.push_back(data);
          p_busIDs.insert(std::pair<int,int>(bus_id,len));
        }
      }
    }

    /**
     * Extract zero sequence shunt data and put it in bus_params structs
     */
    void findZeroSeqShuntData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        /* Check if entry already exists for this bus */
        int bus_id = atoi(split_line[0].c_str());
        std::map<int,int>::iterator it = p_busIDs.find(bus_id);
        if (it != p_busIDs.end()) {
          int idx = it->second;
          p_bus_data[idx].use_zero;
          p_bus_data[idx].gzero = atof(split_line[1].c_str());
          p_bus_data[idx].bzero = atof(split_line[2].c_str());
        } else {
          bus_params data;
          data.use_neg = true;
          data.use_zero = false;
          data.bus_id = bus_id;
          data.gzero = atof(split_line[1].c_str());
          data.bzero = atof(split_line[2].c_str());
          int len = p_bus_data.size();
          p_bus_data.push_back(data);
          p_busIDs.insert(std::pair<int,int>(bus_id,len));
        }
      }
    }

    /**
     * Extract zero sequence non-transformer branch data and put it in branch_params structs
     */
    void findZeroSeqNonTransformerData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        /* Check if entry already exists for this branch */
        int from_bus = atoi(split_line[0].c_str());
        int to_bus = atoi(split_line[1].c_str());
        std::string tag = p_util.clean2Char(split_line[2]);
        std::pair<int,int> branch_id = std::pair<int,int>(from_bus,to_bus);
        std::pair<std::pair<int,int>,std::string> branch = std::pair<std::pair<int,int>,std::string>(branch_id,tag);
        std::map<std::pair<std::pair<int,int>,std::string>,int>::iterator it = p_branchIDs.find(branch);
        if (it != p_branchIDs.end()) {
          int idx = it->second;
          p_branch_data[idx].is_xform = false;
          p_branch_data[idx].from_bus = from_bus;
          p_branch_data[idx].to_bus = to_bus;
          strcpy(p_branch_data[idx].branch_id,tag.c_str());
          p_branch_data[idx].branch_id[2] = '\0';
          p_branch_data[idx].rlinz = atof(split_line[3].c_str());
          p_branch_data[idx].xlinz = atof(split_line[4].c_str());
          p_branch_data[idx].bchz = atof(split_line[5].c_str());
          p_branch_data[idx].gi = atof(split_line[6].c_str());
          p_branch_data[idx].bi = atof(split_line[7].c_str());
          p_branch_data[idx].gj = atof(split_line[8].c_str());
          p_branch_data[idx].bj = atof(split_line[9].c_str());
        } else {
          branch_params data;
          data.is_xform = false;
          data.from_bus = from_bus;
          data.to_bus = to_bus;
          strcpy(data.branch_id,tag.c_str());
          data.branch_id[2] = '\0';
          data.rlinz = atof(split_line[3].c_str());
          data.xlinz = atof(split_line[4].c_str());
          data.bchz = atof(split_line[5].c_str());
          data.gi = atof(split_line[6].c_str());
          data.bi = atof(split_line[7].c_str());
          data.gj = atof(split_line[8].c_str());
          data.bj = atof(split_line[9].c_str());
          int len = p_branch_data.size();
          p_branch_data.push_back(data);
          p_branchIDs.insert(std::pair<std::pair<std::pair<int,int>,std::string>,int>(branch,len));
        }
      }
    }

    /**
     * Extract zero sequence mutual impedence data. Currently ignore this data
     */
    void findZeroSeqMutualImpData()
    {
      std::string line;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        // Ignore this data for now.
      }
    }

    /**
     * Extract zero sequence transformer branch data and put it in branch_params structs
     */
    void findZeroSeqTransformerData()
    {
      std::string line;
      bool three_winding;
      while(p_input_stream.nextLine(line) && test_end(line)) {
        char seq = ',';
        std::vector<std::string> split_line = p_util.charTokenizer(line,&seq);
        if (split_line.size() == 11) {
          three_winding = false;
        } else {
          three_winding = true;
        }
        /* Check if entry already exists for this branch */
        int from_bus = atoi(split_line[0].c_str());
        int to_bus = atoi(split_line[1].c_str());
        std::string tag = p_util.clean2Char(split_line[3]);
        std::pair<int,int> branch_id = std::pair<int,int>(from_bus,to_bus);
        std::pair<std::pair<int,int>,std::string> branch = std::pair<std::pair<int,int>,std::string>(branch_id,tag);
        std::map<std::pair<std::pair<int,int>,std::string>,int>::iterator it = p_branchIDs.find(branch);
        if (it != p_branchIDs.end()) {
          int idx = it->second;
          p_branch_data[idx].is_xform = true;
          p_branch_data[idx].three_winding = three_winding;
          p_branch_data[idx].from_bus = from_bus;
          p_branch_data[idx].to_bus = to_bus;
          strcpy(p_branch_data[idx].branch_id,tag.c_str());
          p_branch_data[idx].branch_id[2] = '\0';
          p_branch_data[idx].k = atoi(split_line[2].c_str());
          p_branch_data[idx].cc = atoi(split_line[4].c_str());
          p_branch_data[idx].rg = atof(split_line[5].c_str());
          p_branch_data[idx].xg = atof(split_line[6].c_str());
          p_branch_data[idx].r1 = atof(split_line[7].c_str());
          p_branch_data[idx].x1 = atof(split_line[8].c_str());
          if (!three_winding) {
            p_branch_data[idx].rg2 = atof(split_line[9].c_str());
            p_branch_data[idx].xg2 = atof(split_line[10].c_str());
          } else {
            p_branch_data[idx].r2 = atof(split_line[11].c_str());
            p_branch_data[idx].x2 = atof(split_line[12].c_str());
            p_branch_data[idx].r3 = atof(split_line[13].c_str());
            p_branch_data[idx].x3 = atof(split_line[14].c_str());
          }
        } else {
          branch_params data;
          data.is_xform = true;
          data.three_winding = three_winding;
          data.from_bus = from_bus;
          data.to_bus = to_bus;
          strcpy(data.branch_id,tag.c_str());
          data.branch_id[2] = '\0';
          data.k = atoi(split_line[2].c_str());
          data.cc = atoi(split_line[4].c_str());
          data.rg = atof(split_line[5].c_str());
          data.xg = atof(split_line[6].c_str());
          data.r1 = atof(split_line[7].c_str());
          data.x1 = atof(split_line[8].c_str());
          if (!three_winding) {
            data.rg2 = atof(split_line[9].c_str());
            data.xg2 = atof(split_line[10].c_str());
          } else {
            data.r2 = atof(split_line[11].c_str());
            data.x2 = atof(split_line[12].c_str());
            data.r3 = atof(split_line[13].c_str());
            data.x3 = atof(split_line[14].c_str());
          }
          int len = p_branch_data.size();
          p_branch_data.push_back(data);
          p_branchIDs.insert(std::pair<std::pair<std::pair<int,int>,std::string>,int>(branch,len));
        }
      }
    }

    void sortSeqData()
    {
      // Store generator parameters in data collection
      int nsize = p_gen_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(p_gen_data[i].bus_id);
      }
      // Distribute data to correct processors
      {
        gridpack::hash_distr::HashDistribution<_network,gen_params,gen_params>
          distr(p_network);
        distr.distributeBusValues(buses,p_gen_data);
      }
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());
        // Find out how many generators are already on bus
        int ngen = 0;
        data->getValue(GENERATOR_NUMBER, &ngen);
        // Identify index of generator to which this data applies
        int g_id = -1;
        if (ngen > 0) {
          // Find 2 character tag for generator ID
          std::string tag = p_gen_data[i].gen_id;
          int j;
          for (j=0; j<ngen; j++) {
            std::string t_id;
            data->getValue(GENERATOR_ID,&t_id,j);
            t_id = p_util.clean2Char(t_id);
            if (tag == t_id) {
              g_id = j;
              break;
            }
          }
        }
        double rval;
        // Add sequence parameters for this generator
        if (p_gen_data[i].use_pos) {
          if (data->getValue(GENERATOR_SEQ_ZRPOS,&rval,g_id)) {
            data->setValue(GENERATOR_SEQ_ZRPOS,p_gen_data[i].zrpos,g_id);
          } else {
            data->addValue(GENERATOR_SEQ_ZRPOS,p_gen_data[i].zrpos,g_id);
          }
          if (data->getValue(GENERATOR_SEQ_ZXPOS,&rval,g_id)) {
            data->setValue(GENERATOR_SEQ_ZXPOS,p_gen_data[i].zxpos,g_id);
          } else {
            data->addValue(GENERATOR_SEQ_ZXPOS,p_gen_data[i].zxpos,g_id);
          }
        }

        if (p_gen_data[i].use_neg) {
          if (data->getValue(GENERATOR_SEQ_ZRNEG,&rval,g_id)) {
            data->setValue(GENERATOR_SEQ_ZRNEG,p_gen_data[i].zrneg,g_id);
          } else {
            data->addValue(GENERATOR_SEQ_ZRNEG,p_gen_data[i].zrneg,g_id);
          }
          if (data->getValue(GENERATOR_SEQ_ZXNEG,&rval,g_id)) {
            data->setValue(GENERATOR_SEQ_ZXNEG,p_gen_data[i].zxneg,g_id);
          } else {
            data->addValue(GENERATOR_SEQ_ZXNEG,p_gen_data[i].zxneg,g_id);
          }
        }

        if (p_gen_data[i].use_zero) {
          if (data->getValue(GENERATOR_SEQ_RZERO,&rval,g_id)) {
            data->setValue(GENERATOR_SEQ_RZERO,p_gen_data[i].zrneg,g_id);
          } else {
            data->addValue(GENERATOR_SEQ_RZERO,p_gen_data[i].zrneg,g_id);
          }
          if (data->getValue(GENERATOR_SEQ_XZERO,&rval,g_id)) {
            data->setValue(GENERATOR_SEQ_XZERO,p_gen_data[i].zxneg,g_id);
          } else {
            data->addValue(GENERATOR_SEQ_XZERO,p_gen_data[i].zxneg,g_id);
          }
        }
      }

      // Store bus and branch parameters in data collection
      nsize = p_bus_data.size();
      buses.clear();
      for (i=0; i<nsize; i++) {
        buses.push_back(p_bus_data[i].bus_id);
      }

      nsize = p_branch_data.size();
      std::vector<std::pair< int,int > > branches;
      for (i=0; i<nsize; i++) {
        branches.push_back(std::pair<int,int>(p_branch_data[i].from_bus,p_branch_data[i].to_bus));
      }
      // Distribute data to correct processors
      std::vector<int> lbranch_idx;
      {
        gridpack::hash_distr::HashDistribution<_network,bus_params,branch_params>
          distr(p_network);
        distr.distributeBusValues(buses,p_bus_data);
        distr.distributeBranchValues(branches,lbranch_idx,p_branch_data);
      }
      // Now match data with corresponding data collection objects
      // Set bus sequence values
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());
        double rval;
        // Add sequence parameters for this bus
        if (p_bus_data[i].use_neg) {
          if (data->getValue(BUS_SEQ_GNEG,&rval)) {
            data->setValue(BUS_SEQ_GNEG,p_bus_data[i].gneg);
          } else {
            data->addValue(BUS_SEQ_GNEG,p_bus_data[i].gneg);
          }
          if (data->getValue(BUS_SEQ_BNEG,&rval)) {
            data->setValue(BUS_SEQ_BNEG,p_bus_data[i].bneg);
          } else {
            data->addValue(BUS_SEQ_BNEG,p_bus_data[i].bneg);
          }
        }

        if (p_bus_data[i].use_zero) {
          if (data->getValue(BUS_SEQ_GZERO,&rval)) {
            data->setValue(BUS_SEQ_GZERO,p_bus_data[i].gzero);
          } else {
            data->addValue(BUS_SEQ_GZERO,p_bus_data[i].gzero);
          }
          if (data->getValue(BUS_SEQ_BZERO,&rval)) {
            data->setValue(BUS_SEQ_BZERO,p_bus_data[i].bzero);
          } else {
            data->addValue(BUS_SEQ_BZERO,p_bus_data[i].bzero);
          }
        }
      }
      // Set branch sequence values
      nsize = lbranch_idx.size();
      int ival;
      for (i=0; i<nsize; i++) {
        int l_idx = lbranch_idx[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBranchData(l_idx).get());
        int nbranch = 0;
        data->getValue(BRANCH_NUM_ELEMENTS, &nbranch);
        // Identify index of branch to which this data applies
        int br_id = -1;
        if (nbranch > 0) {
          // find 2 character tag for branch ID
          std::string tag = p_branch_data[i].branch_id;
          int j;
          for (j=0; j<nbranch; j++) {
            std::string t_id;
            data->getValue(BRANCH_CKT,&t_id,j);
            t_id = p_util.clean2Char(t_id);
            if (tag == t_id) {
              br_id = j;
              break;
            }
          }
        }
        double rval;
        // Add sequence parameters for this branch/transformer
        if (!p_branch_data[i].is_xform) {
          if (data->getValue(BRANCH_SEQ_RLINZ,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_RLINZ,p_branch_data[i].rlinz,br_id);
          } else {
            data->addValue(BRANCH_SEQ_RLINZ,p_branch_data[i].rlinz,br_id);
          }
          if (data->getValue(BRANCH_SEQ_XLINZ,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_XLINZ,p_branch_data[i].xlinz,br_id);
          } else {
            data->addValue(BRANCH_SEQ_XLINZ,p_branch_data[i].xlinz,br_id);
          }
          if (data->getValue(BRANCH_SEQ_BCHZ,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_BCHZ,p_branch_data[i].bchz,br_id);
          } else {
            data->addValue(BRANCH_SEQ_BCHZ,p_branch_data[i].bchz,br_id);
          }

          if (data->getValue(BRANCH_SEQ_GI,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_GI,p_branch_data[i].gi,br_id);
          } else {
            data->addValue(BRANCH_SEQ_GI,p_branch_data[i].gi,br_id);
          }
          if (data->getValue(BRANCH_SEQ_BI,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_BI,p_branch_data[i].bi,br_id);
          } else {
            data->addValue(BRANCH_SEQ_BI,p_branch_data[i].bi,br_id);
          }

          if (data->getValue(BRANCH_SEQ_GJ,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_GJ,p_branch_data[i].gj,br_id);
          } else {
            data->addValue(BRANCH_SEQ_GJ,p_branch_data[i].gj,br_id);
          }
          if (data->getValue(BRANCH_SEQ_BJ,&rval,br_id)) {
            data->setValue(BRANCH_SEQ_BJ,p_branch_data[i].bj,br_id);
          } else {
            data->addValue(BRANCH_SEQ_BJ,p_branch_data[i].bj,br_id);
          }
        } else { 
          if (data->getValue(TRANSFORMER_SEQ_K,&ival,br_id)) {
            data->setValue(TRANSFORMER_SEQ_K,p_branch_data[i].k,br_id);
          } else {
            data->addValue(TRANSFORMER_SEQ_K,p_branch_data[i].k,br_id);
          }
          if (data->getValue(TRANSFORMER_SEQ_CC,&ival,br_id)) {
            data->setValue(TRANSFORMER_SEQ_CC,p_branch_data[i].cc,br_id);
          } else {
            data->addValue(TRANSFORMER_SEQ_CC,p_branch_data[i].cc,br_id);
          }

          if (data->getValue(TRANSFORMER_SEQ_RG,&rval,br_id)) {
            data->setValue(TRANSFORMER_SEQ_RG,p_branch_data[i].rg,br_id);
          } else {
            data->addValue(TRANSFORMER_SEQ_RG,p_branch_data[i].rg,br_id);
          }
          if (data->getValue(TRANSFORMER_SEQ_XG,&rval,br_id)) {
            data->setValue(TRANSFORMER_SEQ_XG,p_branch_data[i].xg,br_id);
          } else {
            data->addValue(TRANSFORMER_SEQ_XG,p_branch_data[i].xg,br_id);
          }

          if (data->getValue(TRANSFORMER_SEQ_R1,&rval,br_id)) {
            data->setValue(TRANSFORMER_SEQ_R1,p_branch_data[i].r1,br_id);
          } else {
            data->addValue(TRANSFORMER_SEQ_R1,p_branch_data[i].r1,br_id);
          }
          if (data->getValue(TRANSFORMER_SEQ_X1,&rval,br_id)) {
            data->setValue(TRANSFORMER_SEQ_X1,p_branch_data[i].x1,br_id);
          } else {
            data->addValue(TRANSFORMER_SEQ_X1,p_branch_data[i].x1,br_id);
          }

          if (!p_branch_data[i].three_winding) {
            if (data->getValue(TRANSFORMER_SEQ_RG2,&rval,br_id)) {
              data->setValue(TRANSFORMER_SEQ_RG2,p_branch_data[i].rg2,br_id);
            } else {
              data->addValue(TRANSFORMER_SEQ_RG2,p_branch_data[i].rg2,br_id);
            }
            if (data->getValue(TRANSFORMER_SEQ_XG2,&rval,br_id)) {
              data->setValue(TRANSFORMER_SEQ_XG2,p_branch_data[i].xg2,br_id);
            } else {
              data->addValue(TRANSFORMER_SEQ_XG2,p_branch_data[i].xg2,br_id);
            }
          } else {

            if (data->getValue(TRANSFORMER_SEQ_R2,&rval,br_id)) {
              data->setValue(TRANSFORMER_SEQ_R2,p_branch_data[i].r2,br_id);
            } else {
              data->addValue(TRANSFORMER_SEQ_R2,p_branch_data[i].r2,br_id);
            }
            if (data->getValue(TRANSFORMER_SEQ_X2,&rval,br_id)) {
              data->setValue(TRANSFORMER_SEQ_X2,p_branch_data[i].x2,br_id);
            } else {
              data->addValue(TRANSFORMER_SEQ_X2,p_branch_data[i].x2,br_id);
            }

            if (data->getValue(TRANSFORMER_SEQ_R3,&rval,br_id)) {
              data->setValue(TRANSFORMER_SEQ_R3,p_branch_data[i].r3,br_id);
            } else {
              data->addValue(TRANSFORMER_SEQ_R3,p_branch_data[i].r3,br_id);
            }
            if (data->getValue(TRANSFORMER_SEQ_X3,&rval,br_id)) {
              data->setValue(TRANSFORMER_SEQ_X3,p_branch_data[i].x3,br_id);
            } else {
              data->addValue(TRANSFORMER_SEQ_X3,p_branch_data[i].x3,br_id);
            }
          }
        }
      }
    }

    /**
     * Test to see if string terminates a section
     * @return: false if first non-blank character is TERM_CHAR
     */
    bool test_end(std::string &str) const
    {
      if (str[0] == TERM_CHAR) {
        return false;
      }
      int len = str.length();
      int i=0;
      while (i<len && str[i] == ' ') {
        i++;
      }
      if (i<len && str[i] != TERM_CHAR) {
        return true;
      } else if (i == len) {
        return true;
      } else if (str[i] == TERM_CHAR) {
        i++;
        if (i>=len || str[i] == ' ' || str[i] == '\\') {
          return false;
        } else {
          return true;
        }
      } else {
        return true;
      }
    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
  private:
    boost::shared_ptr<_network> p_network;

    gridpack::utility::CoarseTimer *p_timer;

    // Vectors of parameters
    std::vector<gen_params> p_gen_data;
    std::vector<bus_params> p_bus_data;
    std::vector<branch_params> p_branch_data;

    // Maps that can be used to fill out parameter arrays
    std::map<std::pair<int,std::string>, int> p_genIDs;
    std::map<int, int> p_busIDs;
    std::map<std::pair<std::pair<int,int>,std::string>, int> p_branchIDs;

    gridpack::utility::StringUtils p_util;

    int p_change_code;

    /**
     * Input stream object
     */
    gridpack::stream::InputStream p_input_stream;
}; /* end of PTI base parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* PSSE_SEQ_PARSER_HPP_ */
