/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * GOSS_parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef GOSS_PARSER_HPP_
#define GOSS_PARSER_HPP_

#define OLD_MAP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "gridpack/component/base_component.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
/*
 *       <xs:complexType name="transmissionElements">
 *         <xs:sequence>
 *           <xs:element maxOccurs="unbounded" minOccurs="0" name="Line" type="transmissionElementLine"/>
 *           <xs:element maxOccurs="unbounded" minOccurs="0" name="Transformer" type="transmissionElementTransformer"/>
 *         </xs:sequence>
 *       </xs:complexType>
 *
 *       <xs:complexType name="transmissionElementLine">
 *         <xs:complexContent>
 *           <xs:extension base="transmissionElement">
 *             <xs:sequence>
 *               <xs:element name="BRANCH_B" type="xs:double"/>
 *               <xs:element name="BRANCH_SHUNT_ADMTTNC_B1" type="xs:double"/>
 *               <xs:element name="BRANCH_SHUNT_ADMTTNC_B2" type="xs:double"/>
 *               <xs:element name="BRANCH_SHUNT_ADMTTNC_G1" type="xs:double"/>
 *               <xs:element name="BRANCH_SHUNT_ADMTTNC_G2" type="xs:double"/>
 *             </xs:sequence>
 *           </xs:extension>
 *         </xs:complexContent>
 *       </xs:complexType>
 */
namespace gridpack {
namespace parser {

template <class _network>
class GOSS_parser : BaseParser<_network>
{
  public:

    enum XML_TYPE {INTEGER, BOOLEAN, DOUBLE, CHARACTER, STRING, N_TYPES};

    /// Constructor
    explicit GOSS_parser(boost::shared_ptr<_network> network) :
      p_network(network), nBuses(0), nBranches(0)  , totalGenerators(0),
      totalLines(0), totalTransformers(0), totalLoads(0), p_case_sbase(0.0),
      p_case_id(0), p_configExists(false)
    {
      this->setNetwork(network);
      p_comm = network->communicator();
    }


    /**
     * Destructor
     */
    virtual ~GOSS_parser(){}

    /**
     * Parse network configuration file and create network
     * @param fileName name of network file
     */
    virtual void parse(const std::string & xmlFileName)
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
      p_timer->configTimer(false);
      int t_total = p_timer->createCategory("Parser:Total Elapsed Time");
      p_timer->start(t_total);
      // setup parser
      boost::property_tree::ptree xmlTree;

      try {
        setupXMLTree(xmlFileName, xmlTree);
      } catch (boost::property_tree::xml_parser_error & e) {
        std::cout << "ERROR IN PARSER SETUP" << std::endl;
        throw e;
      }

      setTypeAssociations(xmlTree);
      loadCase(xmlTree);
      loadBusData(xmlTree);
      loadBranchData(xmlTree);
      this->createNetwork(p_busCollection, p_branchCollection);
      p_timer->stop(t_total);
      p_timer->configTimer(true);
    }

    int getNBuses()             {return nBuses;};
    int getCaseId()             {return p_case_id;};
    int getNBranches()          {return nBranches;};
    double getCaseSbase()       {return p_case_sbase;};

    void copyDataCollection(std::vector<boost::shared_ptr
        <gridpack::component::DataCollection> > & busCollection,
        std::vector<boost::shared_ptr
        <gridpack::component::DataCollection> > & branchCollection)
    {
      busCollection    = p_busCollection;
      branchCollection = p_branchCollection;
    }

    void test_dumpDataColletionVector(std::vector
        <boost::shared_ptr<gridpack::component::DataCollection> >
        & dataCollectionVector)
    {
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> >::iterator
        dataCollection;
      for (dataCollection = dataCollectionVector.begin();
          dataCollection != dataCollectionVector.end(); ++dataCollection)
      {
        (*dataCollection)->dump();
      }
    }

    void test_dumpDataColletion(gridpack::component::DataCollection & dataCollection)
    {
      dataCollection.dump();
    }

    void test_dumpTypeMap()
    {
      typename std::map<std::string, XML_TYPE>::iterator type;
      for (type = typeMap.begin(); type != typeMap.end(); ++type)
      {
        std::cout << "<" << type->first << "; ";
        if (type->second == BOOLEAN) {
          std::cout << "boolean>" << std::endl;
        } else if (type->second == DOUBLE) {
          std::cout << "double>" << std::endl;
        } else if (type->second == INTEGER) {
          std::cout << "integer>" << std::endl;
        } else if (type->second == STRING) {
          std::cout << "string>" << std::endl;
        } else if (type->second == CHARACTER) {
          std::cout << "character>" << std::endl;
        }
      }
    }

  protected:
    /* ************************************************************************
     **************************************************************************
     ***** PRIVATE SCOPE
     **************************************************************************
     *********************************************************************** */
  private:

    void setupXMLTree(const std::string & xmlFileName,
        boost::property_tree::ptree & xmlTree)
    {
      if (p_comm.rank() != 0) return;
      try {
        read_xml(xmlFileName, xmlTree);
      } catch (boost::property_tree::xml_parser_error  & e) {
        throw e;
      }
    }

    /* ************************************************************************
     **************************************************************************
     ***** Type map setup
     **************************************************************************
     *********************************************************************** */
    void setTypeAssociations(boost::property_tree::ptree &  xmlTree)
    {
      if (p_comm.rank() != 0) return;
      std::string          name("");
      boost::property_tree::ptree schemaStyleAttr;

      try {
        schemaStyleAttr =
          xmlTree.get_child("application.grammars.xs:schema");
      } catch (boost::exception & e){
        std::cout << __FILE__ << ":" << __LINE__ << 
          "\n\tThrowing exception from get_child(\"application.grammars.xs:schema\")"
          << std::endl;
        throw;
      }

      typeMap["mrid"] = STRING;
      typeMap["elementIndex"] = INTEGER;

      BOOST_FOREACH( boost::property_tree::ptree::value_type const& elementAttr,
          schemaStyleAttr)
      {
        try {
          recursiveSetTypeAssociation(elementAttr, name);
        } catch (boost::exception & e) {
          std::cout << __FILE__ << ":" << __LINE__ <<
            "\n\tError in recursive type setting "
            << name << std::endl;
        }
      }
    }

    void recursiveSetTypeAssociation(boost::property_tree::ptree::value_type
        const& elementType, std::string & name)
    {
      XML_TYPE             typeEnum;
      // if the XML type starts with xs:element, then it is component value
      if (elementType.first == "xs:element") {
        // get the elements attribute list
        const boost::property_tree::ptree & attributes
          = elementType.second.get_child("<xmlattr>");

        // search for and extract the name and type values
        BOOST_FOREACH( boost::property_tree::ptree::value_type const& attr,
            attributes)
        {
          if (attr.first == "name") {
            name = attr.second.data()   ;
          } else if (attr.first == "type"){
            std::string               type   = attr.second.data();
            std::vector<std::string>  split_type;
            boost::algorithm::split(split_type, type,
                boost::algorithm::is_any_of(":"), boost::token_compress_on);
            if (split_type[0] == "xs"){
              if (split_type[1] == "int") {
                typeEnum = INTEGER;
              } else if (split_type[1] == "boolean") {
                typeEnum = BOOLEAN;
              } else if (split_type[1] == "double") {
                typeEnum = DOUBLE;
              } else if (split_type[1] == "char") {
                typeEnum = CHARACTER;
              } else if (split_type[1] == "string") {
                typeEnum = STRING;
              } else {
                typeEnum = N_TYPES;
              }
              if (typeEnum < N_TYPES) {
                typeMap[name] = typeEnum;
              }
            }
          }
        }
      } else {
        std::string  element        = elementType.first;
        boost::property_tree::ptree subTree;
        try {
          subTree = elementType.second.get_child("");
        } catch (boost::exception & e) {
          std::cout << __FILE__ << ":" << __LINE__ <<
            "\n\tError in recursive statement " << element
            << std::endl;
          throw;
        }
        BOOST_FOREACH( boost::property_tree::ptree::value_type
            const& subTreeType, subTree)
        {
          recursiveSetTypeAssociation(subTreeType, name);
        }
      }

    }

    /* ************************************************************************
     **************************************************************************
     ***** Load case data
     **************************************************************************
     *********************************************************************** */
    void loadCase(boost::property_tree::ptree & xmlTree)
    {
      if (p_comm.rank() == 0) {
        boost::property_tree::ptree caseXMLTree =
          xmlTree.get_child("GridpackPowergrid");
//        xmlTree.get_child("application.GridpackPowergrid");

        BOOST_FOREACH( boost::property_tree::ptree::value_type
            const& caseAttr, caseXMLTree)
        {
          if (caseAttr.first == "CASE_SBASE")
          {
            p_case_sbase    = atof(caseAttr.second.data().c_str());
          }
          if (caseAttr.first == "CASE_ID")
          {
            p_case_id    = atoi(caseAttr.second.data().c_str());
          }
        }
      }

      // Transmit p_case_sbase and p_case_id to all processors
      double sval =  p_case_sbase;
      double rval;
      MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
      int nprocs = p_network->communicator().size();
      int ierr = MPI_Allreduce(&sval,&rval,1,MPI_DOUBLE,MPI_SUM,comm);
      p_case_sbase = rval;
      this->setCaseSBase(p_case_sbase);
      int isval, irval;
      isval = p_case_id;
      ierr = MPI_Allreduce(&isval,&irval,1,MPI_INT,MPI_SUM,comm);
      p_case_id = irval;
      this->setCaseID(p_case_id);
    }

    /* ************************************************************************
     **************************************************************************
     ***** Load bus data
     **************************************************************************
     *********************************************************************** */

    /* ***********************************************************************
     * Find the "Buses" node in the XML tree and read all of the "Bus" tags
     *********************************************************************** */
    void loadBusData(boost::property_tree::ptree & xmlTree)
    {
      if (p_comm.rank() != 0) return;
      boost::property_tree::ptree busXMLTree =
        xmlTree.get_child("GridpackPowergrid.Buses");
//        xmlTree.get_child("application.GridpackPowergrid.Buses");

      // loop through XML Bus subtree "Buses"
      BOOST_FOREACH( boost::property_tree::ptree::value_type
          const& busTree, busXMLTree)
      {
        readBus(busTree);
      }
#if 0
      boost::shared_ptr<gridpack::component::DataCollection>
          data(new gridpack::component::DataCollection);
      data->addValue("N_GENERATORS", totalGenerators);
      data->addValue("N_LOADS", totalLoads);
      p_busCollection.push_back(data);
#endif

    }

    /* ***********************************************************************
     * Read the "Bus" node of the "Buses" subtree. For each XML tag
     * associated with a "Bus" if the tag is:
     *      "Generators" is the set of generators assocaiated with a given
     *           "Bus." Each generator within a bus is given an index value
     *           corresponding to the order in which the "Generator" is read.
     *           For each "Bus" containing N "Generators," the generator index
     *           values range from 0 to N-1.
     *      "Loads" is the set of loads assocaiated with a given "Bus." Each
     *           load within a bus is given an index value corresponding to the
     *           order in which the "Load" is read. For each "Bus" containing
     *           N "Loads," the load index values range from 0 to N-1.
     *      All other tags are loaded into the "Bus."
     *********************************************************************** */
    void readBus(boost::property_tree::ptree::value_type const & busTree)
    {
      boost::shared_ptr<gridpack::component::DataCollection>
        data(new gridpack::component::DataCollection);

      int                      nLoads     = 0;
      int                      nGenerators = 0;

      BOOST_FOREACH( boost::property_tree::ptree::value_type busAttr,
          busTree.second)
      {
        if (busAttr.first == "Generators")
        {
          BOOST_FOREACH( boost::property_tree::ptree::value_type genSet,
              busAttr.second)
          {
            BOOST_FOREACH( boost::property_tree::ptree::value_type genAttr,
                genSet.second)
            {
              loadCollection(data, genAttr, nGenerators);
            }
            ++nGenerators;
          }

        } else if (busAttr.first == "Loads") {

          BOOST_FOREACH( boost::property_tree::ptree::value_type loadSet,
              busAttr.second)
          {
            BOOST_FOREACH( boost::property_tree::ptree::value_type loadAttr,
                loadSet.second)
            {
              loadCollection(data, loadAttr, nLoads);
            }
            ++nLoads;
          }

        } else {
          loadCollection(data, busAttr);
        }
      }
      data->addValue("GENERATOR_NUMBER", nGenerators);
      data->addValue("LOADS_NUMBER", nLoads);
      p_busCollection.push_back(data);
      ++nBuses;

      totalGenerators += nGenerators;
      totalLoads += nLoads;
    }

    /* ************************************************************************
     **************************************************************************
     ***** Load branch data
     **************************************************************************
     *********************************************************************** */

    /* ***********************************************************************
     * Find the "Buses" node in the XML tree and read all of the "Bus" tags
     *********************************************************************** */
    void loadBranchData(boost::property_tree::ptree & xmlTree)
    {
      if (p_comm.rank() != 0) return;
      boost::property_tree::ptree branchXMLTree =
        xmlTree.get_child("GridpackPowergrid.Branches");
//        xmlTree.get_child("application.GridpackPowergrid.Branches");

      // loop through XML Bus subtree "Buses"
      BOOST_FOREACH( boost::property_tree::ptree::value_type const& branchTree,
          branchXMLTree)
      {
        readBranch(branchTree);
      }
#if 0
      boost::shared_ptr<gridpack::component::DataCollection>
          data(new gridpack::component::DataCollection);
      data->addValue("N_TRANSFORMERS", totalTransformers);
      data->addValue("N_LINES", totalLines);
      p_branchCollection.push_back(data);
#endif
   }

    void readBranch(boost::property_tree::ptree::value_type const & branchTree)
    {
      boost::shared_ptr<gridpack::component::DataCollection>
        data(new gridpack::component::DataCollection);
      int                  nTransformers  = 0;
      int                  nLines         = 0;

      BOOST_FOREACH( boost::property_tree::ptree::value_type branchAttr,
          branchTree.second)
      {

          std::cout << "Branch attr = " << branchAttr.first << std::endl;
          if (branchAttr.first == "TransmissionElements")
          {
              BOOST_FOREACH( boost::property_tree::ptree::value_type transAttr,
                      branchAttr.second)
              {
                  if (transAttr.first == "Line")
                  {
                      BOOST_FOREACH( boost::property_tree::ptree::value_type lineSet,
                              transAttr.second)
                      {
                              loadCollection(data, lineSet, nLines);
                      }
                      ++nLines;
                  } else if (transAttr.first == "Transformer")
                  {
                      BOOST_FOREACH( boost::property_tree::ptree::value_type transformerSet,
                              transAttr.second)
                      {
                          loadCollection(data, transformerSet, nTransformers);
                      }
                      ++nTransformers;
                  }
              } // end search for transmission elements
          } else   // transmission element not found
          loadCollection(data, branchAttr);
      }
      data->addValue("LINE_NUMBER", nLines);
      data->addValue("TRANSFORMER_NUMBER", nTransformers);
      totalLines += nLines;
      totalTransformers += nTransformers;
      p_branchCollection.push_back(data);
      ++nBranches;
    }

    /* ************************************************************************
     **************************************************************************
     ***** Load data into data collection
     **************************************************************************
     *********************************************************************** */

    void loadCollection(boost::shared_ptr<gridpack::component::DataCollection> & data,
        boost::property_tree::ptree::value_type & attr)
    {
      XML_TYPE         type               = typeMap[attr.first];

      switch (type)
      {
        case BOOLEAN:
          {
            bool             value          = true;
            if (attr.second.data() == "false") value = false;
            data->addValue(attr.first.c_str(), value);
            break;
          }
        case INTEGER:
          {
            data->addValue(attr.first.c_str(), atoi(attr.second.data().c_str()));
            break;
          }
        case DOUBLE:
          {
            data->addValue(attr.first.c_str(), atof(attr.second.data().c_str()));
            break;
          }
        case CHARACTER:
          {
            data->addValue(attr.first.c_str(), attr.second.data()[0]);
            break;
          }
        case STRING:
          {
            data->addValue(attr.first.c_str(), attr.second.data().c_str());
            break;
          }
        default:
          {
            break;
          }
      }
    }

    void loadCollection(boost::shared_ptr<gridpack::component::DataCollection> & data,
        boost::property_tree::ptree::value_type & attr, int nElems)
    {
      XML_TYPE         type               = typeMap[attr.first];

      switch (type)
      {
        case BOOLEAN:
          {
            bool             value          = true;
            if (attr.second.data() == "false") value = false;
            data->addValue(attr.first.c_str(), value, nElems);
            break;
          }
        case INTEGER:
          {
            data->addValue(attr.first.c_str(), atoi(attr.second.data().c_str()),
                nElems);
            break;
          }
        case DOUBLE:
          {
            data->addValue(attr.first.c_str(), atof(attr.second.data().c_str()), nElems);
            break;
          }
        case CHARACTER:
          {
            data->addValue(attr.first.c_str(), attr.second.data()[0], nElems);
            break;
          }
        case STRING:
          {
            data->addValue(attr.first.c_str(), attr.second.data().c_str(), nElems);
            break;
          }
        default:
          {
            break;
          }
      }
    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
    boost::shared_ptr<_network> p_network;

    int                      nBuses;
    int                      nBranches;
    int                      totalGenerators;
    int                      totalLines;
    int                      totalLoads;
    int                      totalTransformers;

    // Vector of data collection objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >
      p_busCollection;
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >
      p_branchCollection;
    std::map<std::string, XML_TYPE> typeMap;

    parallel::Communicator   p_comm;
    int                      p_case_id;
    double                   p_case_sbase;
    gridpack::utility::CoarseTimer *p_timer;
    bool                     p_configExists;
}; /* end of GOSS_parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* GOSS_PARSER_HPP_ */
