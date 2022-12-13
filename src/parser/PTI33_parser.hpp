/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * PTI33parser.hpp
 *
 *  Created on: December 10, 2014
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef PTI33_PARSER_HPP_
#define PTI33_PARSER_HPP_


#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>

#include "gridpack/utilities/exception.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/stream/input_stream.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/base_pti_parser.hpp"
#include "gridpack/parser/block_parsers/case_parser33.hpp"
#include "gridpack/parser/block_parsers/bus_parser33.hpp"
#include "gridpack/parser/block_parsers/load_parser33.hpp"
#include "gridpack/parser/block_parsers/fixed_shunt_parser33.hpp"
#include "gridpack/parser/block_parsers/generator_parser33.hpp"
#include "gridpack/parser/block_parsers/branch_parser33.hpp"
#include "gridpack/parser/block_parsers/transformer_parser33.hpp"
#include "gridpack/parser/block_parsers/area_parser33.hpp"
#include "gridpack/parser/block_parsers/two_term_parser33.hpp"
#include "gridpack/parser/block_parsers/vsc_line_parser33.hpp"
#include "gridpack/parser/block_parsers/imped_corr_parser33.hpp"
#include "gridpack/parser/block_parsers/multi_term_parser33.hpp"
#include "gridpack/parser/block_parsers/multi_section_parser33.hpp"
#include "gridpack/parser/block_parsers/zone_parser33.hpp"
#include "gridpack/parser/block_parsers/interarea_parser33.hpp"
#include "gridpack/parser/block_parsers/owner_parser33.hpp"
#include "gridpack/parser/block_parsers/facts_parser33.hpp"
#include "gridpack/parser/block_parsers/switched_shunt_parser33.hpp"

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt

namespace gridpack {
namespace parser {


template <class _network>
class PTI33_parser : public BasePTIParser<_network>
{
  public:
    /// Constructor 
    /**
     * 
     * @param network network object that will be filled with contents
     * of network configuration file (must be child of network::BaseNetwork<>)
     */
    PTI33_parser(boost::shared_ptr<_network> network)
      : p_network(network), p_maxBusIndex(-1)
    {
      this->setNetwork(network);
      p_network_data = network->getNetworkData();
    }

    /**
     * Destructor
     */
    virtual ~PTI33_parser()
    {
      p_busData.clear();      // unnecessary
      p_branchData.clear();
    }

    /**
     * Parse network configuration file and create network
     * @param fileName name of network file
     */
    void parse(const std::string &fileName)
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
      p_timer->configTimer(false);
      int t_total = p_timer->createCategory("Parser:Total Elapsed Time");
      p_timer->start(t_total);
      gridpack::utility::StringUtils util;
      std::string tmpstr = fileName;
      util.trim(tmpstr);
      std::string ext = this->getExtension(tmpstr);
      if (ext == "raw") {
        openStream(tmpstr);
        getCase();
        this->createNetwork(p_busData,p_branchData);
      } else if (ext == "dyr") {
        this->getDS(tmpstr);
      }
      p_timer->stop(t_total);
      p_timer->configTimer(true);
    }

    /**
     * parse string vector representing a PSS/E RAW or DYR file
     * @param fileVec vector of string representing file
     * @param isRAW flag distinguishing between RAW and DYR files
     */
    void parse(const std::vector<std::string> &fileVec, bool isRAW = false)
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
      p_timer->configTimer(false);
      int t_total = p_timer->createCategory("Parser:Total Elapsed Time");
      p_timer->start(t_total);
      openStream(fileVec);
      if (isRAW) {
        getCase();
        this->createNetwork(p_busData,p_branchData);
      } else {
        this->getDS(fileVec);
      }
      p_timer->stop(t_total);
      p_timer->configTimer(true);
    }

    /**
     * Return values of impedence correction table corresponding to tableID
     * @param tableID ID of correction table
     * @param turns Either turns ration or phase shift
     * @param scale scaling factor for transformer nominal impedence
     */
    void getImpedenceTable(int tableID, std::vector<double> &turns,
        std::vector<double> &scale)
    {
      turns.clear();
      scale.clear();
      if (p_imp_corr_table.find(tableID) != p_imp_corr_table.end()) {
        boost::shared_ptr<gridpack::component::DataCollection> data =
          p_imp_corr_table.find(tableID)->second;
        int i;
        char buf[32];
        for (i=0; i<11; i++) {
          sprintf(buf,"XFMR_CORR_TABLE_T%d",i+1);
          double tval, sval;
          bool ok = data->getValue(buf,&tval);
          sprintf(buf,"XFMR_CORR_TABLE_F%d",i+1);
          ok = ok && data->getValue(buf,&sval);
          if (ok) {
            turns.push_back(tval);
            scale.push_back(sval);
          } else {
            break;
          }
        }
      }
    }

  protected:
    /**
     * Open in input stream object based on accessing a file
     */
    void openStream(const std::string & fileName)
    {
      p_istream.openFile(fileName.c_str());
      if (!p_istream.isOpen()) {
        char buf[512];
        sprintf(buf,"Failed to open network configuration file: %s\n\n",
            fileName.c_str());
        throw gridpack::Exception(buf);
      }
    }

    /**
     * Open in input stream object based on a vector of strings
     */
    void openStream(const std::vector<std::string> & fileVec)
    {
      p_istream.openStringVector(fileVec);
      if (!p_istream.isOpen()) {
        char buf[512];
        sprintf(buf,"Failed to open network configuration stream: %s\n\n");
        throw gridpack::Exception(buf);
      }
    }

    /*
     * A case is the collection of all data associated with a PTI33 file.
     * Each case is a a vector of data_set objects the contain all the data
     * associated with a partition of the PTI file. For example, the bus
     * data in the file constitutes a data_set. Each data_set is a vector of
     * gridpack::component::DataCollection objects. Each of these objects
     * contain a single instance of the data associated with a data_set. For
     * example, each line of the bus partition corresponds to a single
     * DataCollection object.
     */
    void getCase()
    {
      int t_case = p_timer->createCategory("Parser:getCase");
      p_timer->start(t_case);
      p_busData.clear();
      p_branchData.clear();
      p_busMap.clear();

      MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());

      int me(p_network->communicator().rank());

      if (me == 0) {
        gridpack::parser::CaseParser33 case_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        case_parser.parse(p_istream,p_network_data,p_case_sbase,p_case_id);
      } else {
        p_case_sbase = 0.0;
        p_case_id = 0;
      }

      // Transmit CASE_SBASE to all processors
      double sval =  p_case_sbase;
      double rval;
      int nprocs = p_network->communicator().size();
      int ierr = MPI_Allreduce(&sval,&rval,1,MPI_DOUBLE,MPI_SUM,comm);
      p_case_sbase = rval;
      // Transmit CASE_ID to all processors
      int isval, irval;
      isval = p_case_id;
      ierr = MPI_Allreduce(&isval,&irval,1,MPI_INT,MPI_SUM,comm);
      p_case_id = irval;
      this->setCaseID(p_case_id);
      this->setCaseSBase(p_case_sbase);

      if (me == 0) {
        gridpack::parser::BusParser33 bus_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        bus_parser.parse(p_istream,p_busData,p_case_sbase,p_case_id,
            &p_maxBusIndex);
        gridpack::parser::LoadParser33 load_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        load_parser.parse(p_istream,p_busData);
        gridpack::parser::FixedShuntParser33 fixed_shunt_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        fixed_shunt_parser.parse(p_istream,p_busData);
        gridpack::parser::GeneratorParser33 generator_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        generator_parser.parse(p_istream,p_busData);
        gridpack::parser::BranchParser33 branch_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        branch_parser.parse(p_istream,p_branchData);
        gridpack::parser::TransformerParser33 transformer_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        transformer_parser.parse(p_istream,p_busData,p_branchData,p_case_sbase,
            p_maxBusIndex);
        gridpack::parser::AreaParser33 area_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        area_parser.parse(p_istream,p_network_data);
        gridpack::parser::TwoTermParser33 two_term_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        two_term_parser.parse(p_istream);
        gridpack::parser::VSCLineParser33 vsc_line_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        vsc_line_parser.parse(p_istream);
        gridpack::parser::ImpedCorrParser33 imped_corr_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        imped_corr_parser.parse(p_istream,p_imp_corr_table);
        gridpack::parser::MultiTermParser33 multi_term_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        multi_term_parser.parse(p_istream);
        gridpack::parser::MultiSectParser33 multi_section_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        multi_section_parser.parse(p_istream,p_branchData);
        gridpack::parser::ZoneParser33 zone_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        zone_parser.parse(p_istream);
        gridpack::parser::InterAreaParser33 interarea_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        interarea_parser.parse(p_istream);
        gridpack::parser::OwnerParser33 owner_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        owner_parser.parse(p_istream);
        gridpack::parser::FACTSParser33 facts_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        facts_parser.parse(p_istream);
        gridpack::parser::SwitchedShuntParser33 switched_shunt_parser(&p_busMap,
            &p_nameMap, &p_branchMap);
        switched_shunt_parser.parse(p_istream,p_busData);

#if 0
        // debug
        int i;
        printf("BUS data size: %d\n",(int)p_busData.size());
        for (i=0; i<p_busData.size(); i++) {
          printf("Dumping bus: %d\n",i);
          p_busData[i]->dump();
        }
        printf("BRANCH data size: %d\n",(int)p_branchData.size());
        for (i=0; i<p_branchData.size(); i++) {
          printf("Dumping branch: %d\n",i);
          p_branchData[i]->dump();
        }
#endif
        p_istream.close();
      }
      // Distribute impedence correction tables to all processors (if necessary)
      int stables = 0;
      int ntables;
      if (me == 0) {
        stables =  p_imp_corr_table.size();
      }
      ierr = MPI_Allreduce(&stables, &ntables, 1, MPI_INT, MPI_SUM, comm);
#if 0
      if (ntables > 0) {
        // distribute tables
        if (me == 0) {
          int idx;
          for (idx=1; idx<nprocs; idx++) {
            static_cast<boost::mpi::communicator>(
                p_network->communicator()).send(idx,idx,p_imp_corr_table);
          }
        } else {
          p_imp_corr_table.clear();
          static_cast<boost::mpi::communicator>(
              p_network->communicator()).recv(0,me,p_imp_corr_table);
        }
      }
#endif
      p_network->broadcastNetworkData(0);
      p_network_data = p_network->getNetworkData();
      p_timer->stop(t_case);
    }

    // Distribute data uniformly on processors
    void brdcst_data(void)
    {
      int t_brdcst = p_timer->createCategory("Parser:brdcst_data");
      int t_serial = p_timer->createCategory("Parser:data packing and unpacking");
      MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
      int me(p_network->communicator().rank());
      int nprocs(p_network->communicator().size());
      if (nprocs == 1) return;
      p_timer->start(t_brdcst);

      // find number of buses and branches and broadcast this information to
      // all processors
      int sbus, sbranch;
      if (me == 0) {
        sbus = p_busData.size();
        sbranch = p_branchData.size();
      } else {
        sbus = 0;
        sbranch = 0;
      }
      int ierr, nbus, nbranch;
      ierr = MPI_Allreduce(&sbus, &nbus, 1, MPI_INT, MPI_SUM, comm);
      ierr = MPI_Allreduce(&sbranch, &nbranch, 1, MPI_INT, MPI_SUM, comm);
      double rprocs = static_cast<double>(nprocs);
      double rme = static_cast<double>(me);
      int n, i;
      std::vector<gridpack::component::DataCollection>recvV;
      // distribute buses
      if (me == 0) {
        for (n=0; n<nprocs; n++) {
          double rn = static_cast<double>(n);
          int istart = static_cast<int>(static_cast<double>(nbus)*rn/rprocs);
          int iend = static_cast<int>(static_cast<double>(nbus)*(rn+1.0)/rprocs);
          if (n != 0) {
            p_timer->start(t_serial);
            std::vector<gridpack::component::DataCollection> sendV;
            for (i=istart; i<iend; i++) {
              sendV.push_back(*(p_busData[i]));
            }
            p_timer->stop(t_serial);
            static_cast<boost::mpi::communicator>(p_network->communicator()).send(n,n,sendV);
          } else {
            p_timer->start(t_serial);
            for (i=istart; i<iend; i++) {
              recvV.push_back(*(p_busData[i]));
            }
            p_timer->stop(t_serial);
          }
        }
      } else {
        int istart = static_cast<int>(static_cast<double>(nbus)*rme/rprocs);
        int iend = static_cast<int>(static_cast<double>(nbus)*(rme+1.0)/rprocs)-1;
        static_cast<boost::mpi::communicator>(p_network->communicator()).recv(0,me,recvV);
      }
      int nsize = recvV.size();
      p_busData.clear();
      p_timer->start(t_serial);
      for (i=0; i<nsize; i++) {
        boost::shared_ptr<gridpack::component::DataCollection> data(new
            gridpack::component::DataCollection);
        *data = recvV[i];
        p_busData.push_back(data);
      }
      p_timer->stop(t_serial);
      recvV.clear();
      // distribute branches
      if (me == 0) {
        for (n=0; n<nprocs; n++) {
          double rn = static_cast<double>(n);
          int istart = static_cast<int>(static_cast<double>(nbranch)*rn/rprocs);
          int iend = static_cast<int>(static_cast<double>(nbranch)*(rn+1.0)/rprocs);
          if (n != 0) {
            p_timer->start(t_serial);
            std::vector<gridpack::component::DataCollection> sendV;
            for (i=istart; i<iend; i++) {
              sendV.push_back(*(p_branchData[i]));
            }
            p_timer->stop(t_serial);
            static_cast<boost::mpi::communicator>(p_network->communicator()).send(n,n,sendV);
          } else {
            p_timer->start(t_serial);
            for (i=istart; i<iend; i++) {
              recvV.push_back(*(p_branchData[i]));
            }
            p_timer->stop(t_serial);
          }
        }
      } else {
        int istart = static_cast<int>(static_cast<double>(nbranch)*rme/rprocs);
        int iend = static_cast<int>(static_cast<double>(nbranch)*(rme+1.0)/rprocs)-1;
        static_cast<boost::mpi::communicator>(p_network->communicator()).recv(0,me,recvV);
      }
      nsize = recvV.size();
      p_branchData.clear();
      p_timer->start(t_serial);
      for (i=0; i<nsize; i++) {
        boost::shared_ptr<gridpack::component::DataCollection> data(new
            gridpack::component::DataCollection);
        *data = recvV[i];
        p_branchData.push_back(data);
      }
      p_timer->stop(t_serial);
#if 0
      // debug
      printf("p[%d] BUS data size: %d\n",me,p_busData.size());
      for (i=0; i<p_busData.size(); i++) {
        printf("p[%d] Dumping bus: %d\n",me,i);
        p_busData[i]->dump();
      }
      printf("p[%d] BRANCH data size: %d\n",me,p_branchData.size());
      for (i=0; i<p_branchData.size(); i++) {
        printf("p[%d] Dumping branch: %d\n",me,i);
        p_branchData[i]->dump();
      }
#endif
      p_timer->stop(t_brdcst);
    }

  private:

    /*
     * The case_data is the collection of all data points in the case file.
     * Each collection in the case data contains the data associated with a given
     * type. For example, the case is the collection of data describing the
     * current case and the bus data is the collection of data associated with
     * each bus. The type data may consist of zero or more instances of the
     * given type. For example, the bus data may contain several instances of
     * a bus. These type instances are composed of a set of key value pairs.
     * Each column as an associated key and each row is an instance of a given
     * type. When the parser is reading data for a type, the value found in each
     * column associated with the key for that column in a field_data structure.
     *
     * Within the PTI file there are the following group of data sets in order:
     *     case
     *     bus
     *     generator
     *     branch
     *     transformer
     *     dc_line
     *     shunt
     *     impedence corr
     *     multi-terminal
     *     multi-section
     *     zone
     *     inter-area
     *     owner
     *     device driver
     *
     * These data sets are stored in the case data as a collection of
     * data set and each data set is a
     */
    boost::shared_ptr<_network> p_network;

    // Vector of bus data objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
    // Vector of branch data objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_branchData;
    // Vector of Impedence correction tables
    std::map<int,boost::shared_ptr<gridpack::component::DataCollection> >
      p_imp_corr_table;
    // Map of PTI indices to index in p_busData
    std::map<int,int> p_busMap;
    std::map<std::string,int> p_nameMap;
    // Map of PTI index pair to index in p_branchData
    std::map<std::pair<int, int>, int> p_branchMap;

    // Input stream object that obtains next line
    gridpack::stream::InputStream p_istream;

    // Global variables that apply to whole network
    int p_case_id;
    int p_maxBusIndex;
    double p_case_sbase;
    gridpack::utility::CoarseTimer *p_timer;

    /**
     * Data collection object associated with network as a whole
     */
    boost::shared_ptr<gridpack::component::DataCollection> p_network_data;
};

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI33PARSER_HPP_ */
