/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: December 30, 2014
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef BASEPTIPARSER_HPP_
#define BASEPTIPARSER_HPP_

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
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "parser_classes/gencls.hpp"
#include "parser_classes/gensal.hpp"
#include "parser_classes/genrou.hpp"
#include "parser_classes/wsieg1.hpp"
#include "parser_classes/exdc1.hpp"
#include "parser_classes/esst1a.hpp"
#include "parser_classes/esst4b.hpp"
#include "parser_classes/ggov1.hpp"
#include "parser_classes/wshygp.hpp"

namespace gridpack {
namespace parser {

template <class _network>
class BasePTIParser : public BaseParser<_network>
{
  public:

    /**
     * Constructor
     */
    explicit BasePTIParser()
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
    }


    /**
     * Destructor
     */
    virtual ~BasePTIParser(){}

    /**
     * Parse a second file after original network has been distributed. This
     * requires the data in the second file to be distributed to all network
     * objects that need the data
     * @param fileName name of file
     */
    void externalParse(const std::string &fileName)
    {
      std::string ext = getExtension(fileName);
      if (ext == "dyr") {
        getDSExternal(fileName);
      } else if (ext == "uc") {
        getUCExternal(fileName);
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
     * This routine opens up a .dyr file with parameters for dynamic
     * simulation. It assumes that a .raw file has already been parsed
     */
    void getDS(const std::string & fileName)
    {
      int t_ds = p_timer->createCategory("Parser:getDS");
      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      if (me == 0) {
        std::ifstream input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          p_timer->stop(t_ds);
          return;
        }
        find_ds_par(input);
        input.close();
      }
      p_timer->stop(t_ds);
#if 0
      int i;
      printf("BUS data size: %d\n",p_busData.size());
      for (i=0; i<p_network->numBuses(); i++) {
        printf("Dumping bus: %d\n",i);
        p_network->getBusData(i)->dump();
      }
#endif
    }

    struct ds_params{
      // Generator parameters
      int bus_id; // ID of bus that owns generator
      char gen_id[3]; // Generator ID
      char gen_model[8];  // Generator model
      double inertia;  // Inertia constant 0
      double damping;  // Damping coefficient
      double reactance; // Transient reactance
      double tdop;
      double tdopp;
      double tqop;
      double tqopp;
      double gn_xd;
      double gn_xq;
      double xdp;
      double xqp;
      double xdpp;
      double gn_xl;
      double gn_s1;
      double s12;
      // Exciter parameters
      char ex_model[8];  // Exciter model
      bool has_exciter;
      double ex_tr;
      double ex_ka;
      double ex_ta;
      double ex_tb;
      double ex_tc;
      double vrmax;
      double vrmin;
      double ex_ke;
      double ex_te;
      double ex_kf;
      double tf1;
      double rswitch;
      double ex_e1;
      double se1;
      double ex_e2;
      double se2;
      double uel;
      double vos;
      double vimax;
      double vimin;
      double tc1;
      double tb1;
      double vamax;
      double vamin;
      double ex_kc;
      double ex_tf;
      double klr;
      double ilr;
      double kpr;
      double kir;
      double kpm;
      double kim;
      double vmmax;
      double vmmin;
      double ex_kg;
      double ex_kp;
      double ex_ki;
      double vbmax;
      double ex_xl;
      double thetap;
      // Governor parameters
      char gov_model[8];  // Governor model
      bool has_governor;
      int jbus;
      int gv_m;
      double gv_k;
      double gv_t1;
      double gv_t2;
      double gv_t3;
      double gv_uo;
      double gv_uc;
      double pmax;
      double pmin;
      double gv_t4;
      double gv_k1;
      double gv_k2;
      double gv_t5;
      double gv_k3;
      double gv_k4;
      double gv_t6;
      double gv_k5;
      double gv_k6;
      double gv_t7;
      double gv_k7;
      double gv_k8;
      double db1;
      double err;
      double db2;
      double gv1;
      double pgv1;
      double gv2;
      double pgv2;
      double gv3;
      double pgv3;
      double gv4;
      double pgv4;
      double gv5;
      double pgv5;
      int iblock;
      double rselect;
      double flagswitch;
      double gv_r;
      double tpelec;
      double maxerr;
      double minerr;
      double kpgov;
      double kigov;
      double kdgov;
      double tdgov;
      double vmax;
      double vmin;
      double tact;
      double kturb;
      double wfnl;
      double gv_tb;
      double gv_tc;
      double teng;
      double tfload;
      double kpload;
      double kiload;
      double ldref;
      double gv_dm;
      double ropen;
      double rclose;
      double kimw;
      double aset;
      double gv_ka;
      double gv_ta;
      double trate;
      double gv_db;
      double tsa;
      double tsb;
      double rup;
      double rdown;
      double gv_td;
      double gv_ki;
      double gv_tf;
      double gv_kd;
      double gv_kp;
      double gv_tt;
      double gv_kg;
      double gv_tp;
      double velopen;
      double velclose;
      double aturb;
      double bturb;
      double tturb;
    };

    /**
     * This routine opens up a .dyr file with parameters for dynamic
     * simulation and distributes the parameters to whatever processor holds the
     * corresponding buses. It assumes that a .raw file has already been parsed
     */
    void getDSExternal(const std::string & fileName)
    {

      //      int t_ds = p_timer->createCategory("Parser:getDS");
      //      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      std::vector<ds_params> ds_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          // p_timer->stop(t_ds);
          return;
        }
        find_ds_vector(input, &ds_data);
        input.close();
      }
      int nsize = ds_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(ds_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,ds_params,ds_params>
        distr(p_network);
      distr.distributeBusValues(buses,ds_data);
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        std::string tag = ds_data[i].gen_id;
        int j;
        for (j=0; j<ngen; j++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,j);
          if (tag == t_id) {
            g_id = j;
            break;
          }
        }
        if (g_id == -1) continue;

        std::string sval;
        double rval;
        int ival;
        bool bval=true;
        if (!strcmp(ds_data[i].gen_model,"GENCLS")) {
          GenclsParser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"GENSAL")) {
          GensalParser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"GENROU")) {
          GenrouParser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"WSIEG1")) {
          Wsieg1Parser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"EXDC1") ||
            !strcmp(ds_data[i].gen_model,"EXDC2")) {
          Exdc1Parser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"ESST1A")) {
          Esst1aParser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"ESST4B")) {
          Esst4bParser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"GGOV1")) {
          Ggov1Parser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        } else if (!strcmp(ds_data[i].gen_model,"WSHYGP")) {
          WshygpParser<ds_params> parser;
          parser.extract(ds_data[i], data, sval, g_id);
        }
      }
      //      p_timer->stop(t_ds);
    }

    struct uc_params{
      int bus_id; // ID of bus that owns generator
      char gen_id[3]; // Generator ID
      int type;
      double init_level; // Initial production level
      double min_gen; // Minimum generation
      double max_gen; // Maximum generation
      double max_oper; // Maximum operating generation
      double min_up;
      double min_down;
      double ramp_up;
      double ramp_down;
      double start_up; // Start up cost
      double const_cost; // Constant cost
      double lin_cost; // Linear cost
      double co_2_cost;
      double init_prd; // Init periods
      double start_cap; // Startup cap
      double shut_cap; // Shutdown cap
    };

    /**
     * This routine opens up a .uc file with parameters for a unit commitment
     * calculationand distributes the parameters to whatever processor holds the
     * corresponding buses. It assumes that a .raw file has already been parsed
     */
    void getUCExternal(const std::string & fileName)
    {

      int me(p_network->communicator().rank());

      std::vector<uc_params> uc_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          return;
        }
        find_uc_vector(input, &uc_data);
        input.close();
      }
      int nsize = uc_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(uc_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,uc_params,uc_params>
        distr(p_network);
      distr.distributeBusValues(buses,uc_data);
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        std::string tag = uc_data[i].gen_id;
        int j;
        for (j=0; j<ngen; j++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,j);
          if (tag == t_id) {
            g_id = j;
            break;
          }
        }
        if (g_id == -1) continue;

        double rval;
        int ival;
        if (!data->getValue("GENERATOR_TYPE",&ival,g_id)) {
          data->addValue("GENERATOR_TYPE", uc_data[i].type, g_id);
        } else {
          data->setValue("GENERATOR_TYPE", uc_data[i].type, g_id);
        }

        if (!data->getValue("GENERATOR_INIT_LEVEL",&rval,g_id)) {
          data->addValue("GENERATOR_INIT_LEVEL", uc_data[i].init_level, g_id);
        } else {
          data->setValue("GENERATOR_INIT_LEVEL", uc_data[i].init_level, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_GEN", uc_data[i].min_gen, g_id);
        } else {
          data->setValue("GENERATOR_MIN_GEN", uc_data[i].min_gen, g_id);
        }

        if (!data->getValue("GENERATOR_MAX_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MAX_GEN", uc_data[i].max_gen, g_id);
        } else {
          data->setValue("GENERATOR_MAX_GEN", uc_data[i].max_gen, g_id);
        }

        if (!data->getValue("GENERATOR_MAX_OP_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MAX_OP_GEN", uc_data[i].max_gen, g_id);
        } else {
          data->setValue("GENERATOR_MAX_OP_GEN", uc_data[i].max_gen, g_id);
        }
        if (!data->getValue("GENERATOR_MIN_UP",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_UP", uc_data[i].min_up, g_id);
        } else {
          data->setValue("GENERATOR_MIN_UP", uc_data[i].min_up, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_DOWN",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_DOWN", uc_data[i].min_down, g_id);
        } else {
          data->setValue("GENERATOR_MIN_DOWN", uc_data[i].min_down, g_id);
        }

        if (!data->getValue("GENERATOR_RAMP_UP",&rval,g_id)) {
          data->addValue("GENERATOR_RAMP_UP", uc_data[i].ramp_up, g_id);
        } else {
          data->setValue("GENERATOR_RAMP_UP", uc_data[i].ramp_up, g_id);
        }

        if (!data->getValue("GENERATOR_RAMP_DOWN",&rval,g_id)) {
          data->addValue("GENERATOR_RAMP_DOWN", uc_data[i].ramp_down, g_id);
        } else {
          data->setValue("GENERATOR_RAMP_DOWN", uc_data[i].ramp_down, g_id);
        }

        if (!data->getValue("GENERATOR_START_UP",&rval,g_id)) {
          data->addValue("GENERATOR_START_UP", uc_data[i].start_up, g_id);
        } else {
          data->setValue("GENERATOR_START_UP", uc_data[i].start_up, g_id);
        }

        if (!data->getValue("GENERATOR_CONST_COST",&rval,g_id)) {
          data->addValue("GENERATOR_CONST_COST", uc_data[i].const_cost, g_id);
        } else {
          data->setValue("GENERATOR_CONST_COST", uc_data[i].const_cost, g_id);
        }
        if (!data->getValue("GENERATOR_LIN_COST",&rval,g_id)) {
          data->addValue("GENERATOR_LIN_COST", uc_data[i].lin_cost, g_id);
        } else {
          data->setValue("GENERATOR_LIN_COST", uc_data[i].lin_cost, g_id);
        }
        if (!data->getValue("GENERATOR_CO_2_COST",&rval,g_id)) {
          data->addValue("GENERATOR_CO_2_COST", uc_data[i].co_2_cost, g_id);
        } else {
          data->setValue("GENERATOR_CO_2_COST", uc_data[i].co_2_cost, g_id);
        }

        if (!data->getValue("GENERATOR_INIT_PRD",&rval,g_id)) {
          data->addValue("GENERATOR_INIT_PRD", uc_data[i].init_prd, g_id);
        } else {
          data->setValue("GENERATOR_INIT_PRD", uc_data[i].init_prd, g_id);
        }

        if (!data->getValue("GENERATOR_START_CAP",&rval,g_id)) {
          data->addValue("GENERATOR_START_CAP", uc_data[i].start_cap, g_id);
        } else {
          data->setValue("GENERATOR_START_CAP", uc_data[i].start_cap, g_id);
        }

        if (!data->getValue("GENERATOR_SHUT_CAP",&rval,g_id)) {
          data->addValue("GENERATOR_SHUT_CAP", uc_data[i].shut_cap, g_id);
        } else {
          data->setValue("GENERATOR_SHUT_CAP", uc_data[i].shut_cap, g_id);
        }
      }
    }

    // Extract extension from file name and convert it to lower case
    std::string getExtension(const std::string file)
    {
      std::string ret;
      std::string line = file;
      int ntok1 = line.find('.',0);
      if (ntok1 == std::string::npos) return
        ret;
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

    void find_ds_par(std::ifstream & input)
    {
      std::string          line;
      gridpack::component::DataCollection *data;
      while(std::getline(input,line)) {
        // Check to see if line is blank
        int idx = line.find_first_not_of(' ');
        if (idx == std::string::npos) continue;

        std::string record = line;
        idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        // GENERATOR_BUSNUMBER               "I"                   integer
        int l_idx, o_idx;
        o_idx = atoi(split_line[0].c_str());
#ifdef OLD_MAP
        std::map<int, int>::iterator it;
#else
        boost::unordered_map<int, int>::iterator it;
#endif
        int nstr = split_line.size();
        it = p_busMap->find(o_idx);
        if (it != p_busMap->end()) {
          l_idx = it->second;
        } else {
          continue;
        }
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        gridpack::utility::StringUtils util;
        std::string tag = util.clean2Char(split_line[2]);
        int i;
        for (i=0; i<ngen; i++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,i);
          if (tag == t_id) {
            g_id = i;
            break;
          }
        }
        if (g_id == -1) continue;

        std::string sval;
        double rval;
        int ival;
        bool bval;

        // GENERATOR_MODEL              "MODEL"                  string
        sval = util.trimQuotes(split_line[1]);
        util.toUpper(sval);

        if (sval == "GENCLS") {
          GenclsParser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "GENSAL") {
          GensalParser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "GENROU") {
          GenrouParser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "WSIEG1") {
          Wsieg1Parser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "EXDC1" || sval == "EXDC2") {
          Exdc1Parser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "ESST1A") {
          Esst1aParser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "ESST4A") {
          Esst4bParser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "GGOV1") {
          Ggov1Parser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        } else if (sval == "WSHYGP") {
          WshygpParser<ds_params> parser;
          parser.parse(split_line, data, sval, g_id);
        }
      }
    }

    void find_ds_vector(std::ifstream & input, std::vector<ds_params> *ds_vector)
    {
      std::string          line;
      ds_vector->clear();
      while(std::getline(input,line)) {
        // Check to see if line is blank
        int idx = line.find_first_not_of(' ');
        if (idx == std::string::npos) continue;

        std::string record = line;
        idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        ds_params data;

        // GENERATOR_BUSNUMBER               "I"                   integer
        int o_idx;
        o_idx = atoi(split_line[0].c_str());
        data.bus_id = o_idx;

        int nstr = split_line.size();

        // Clean up 2 character tag for generator ID
        gridpack::utility::StringUtils util;
        std::string tag = util.clean2Char(split_line[2]);
        strcpy(data.gen_id, tag.c_str());

        std::string sval;
        double rval;
        int ival;

        sval = util.trimQuotes(split_line[1]);
        util.toUpper(sval);

        // GENERATOR_MODEL              "MODEL"                  integer
        strcpy(data.gen_model, sval.c_str());

        if (sval == "GENCLS") {
          GenclsParser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "GENSAL") {
          GensalParser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "GENROU") {
          GenrouParser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "WSIEG1") {
          Wsieg1Parser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "EXDC1" || sval == "EXDC2") {
          Exdc1Parser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "ESST1A") {
          Esst1aParser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "ESST4B") {
          Esst4bParser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "GGOV1") {
          Ggov1Parser<ds_params> parser;
          parser.store(split_line,data);
        } else if (sval == "WSHYGP") {
          WshygpParser<ds_params> parser;
          parser.store(split_line,data);
        }
        ds_vector->push_back(data);
      }
    }

    void find_uc_vector(std::ifstream & input, std::vector<uc_params> *uc_vector)
    {
      std::string          line;
      uc_vector->clear();
      // Ignore first line containing header information
      std::getline(input,line);
      while(std::getline(input,line)) {
        std::vector<std::string>  split_line;
        boost::split(split_line, line, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        uc_params data;

        int nstr = split_line.size();
        if (nstr > 1) {
          data.type = atoi(split_line[1].c_str());
        }
        if (nstr > 2) {
          data.init_level = atof(split_line[2].c_str());
        }
        if (nstr > 3) {
          data.min_gen = atof(split_line[3].c_str());
        }
        if (nstr > 4) {
          data.max_gen = atof(split_line[4].c_str());
        }
        if (nstr > 5) {
          data.max_oper = atof(split_line[5].c_str());
        }
        if (nstr > 6) {
          data.min_up = atof(split_line[6].c_str());
        }
        if (nstr > 7) {
          data.min_down = atof(split_line[7].c_str());
        }
        if (nstr > 8) {
          data.ramp_up = atof(split_line[8].c_str());
        }
        if (nstr > 9) {
          data.ramp_down = atof(split_line[9].c_str());
        }
        if (nstr > 10) {
          data.start_up = atof(split_line[10].c_str());
        }
        if (nstr > 11) {
          data.const_cost = atof(split_line[11].c_str());
        }
        if (nstr > 12) {
          data.lin_cost = atof(split_line[12].c_str());
        }
        if (nstr > 13) {
          data.co_2_cost = atof(split_line[13].c_str());
        }
        if (nstr > 14) {
          data.init_prd = atof(split_line[14].c_str());
        }
        if (nstr > 15) {
          data.start_cap = atof(split_line[15].c_str());
        }
        if (nstr > 16) {
          data.shut_cap = atof(split_line[16].c_str());
        }
        if (nstr > 17) {
          data.bus_id = atoi(split_line[17].c_str());
        }
        if (nstr > 18) {
          // Clean up 2 character tag for generator ID
          gridpack::utility::StringUtils util;
          std::string tag = util.clean2Char(split_line[18]);
          strcpy(data.gen_id, tag.c_str());
        }
        uc_vector->push_back(data);
      }
    }

    /** Store pointer to bus and branch maps
     * @param busMap pointer to map between PTI and local indices
     * @param branchMap pointer to map between PTI bus pairs and local branch
     *        indices
     */
    void setMaps(std::map<int,int> *busMap,
                 std::map<std::pair<int, int>, int> *branchMap)
    {
      p_busMap = busMap;
      p_branchMap = branchMap;
    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
    boost::shared_ptr<_network> p_network;

    gridpack::utility::CoarseTimer *p_timer;

    // Map of PTI indices to index in p_busData
    std::map<int,int> *p_busMap;
    // Map of PTI index pair to index in p_branchData
    std::map<std::pair<int, int>, int> *p_branchMap;

}; /* end of PTI base parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* BASEPTIPARSER_HPP_ */
