/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * GOSSparser.hpp
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
class GOSS_parser
{
  public:

    enum XML_TYPE {INTEGER, BOOLEAN, DOUBLE, CHARACTER, STRING, N_TYPES};

    /// Constructor
    explicit GOSS_parser(boost::shared_ptr<_network> network) :nBuses(0),
      p_network(network), nBranches(0), p_case_sbase(0) {};


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
    }

    int getNBuses()             {return nBuses;};
    std::string getCaseId()     {return p_case_id;};
    int getNBranches()          {return nBranches;};
    int getCaseSbase()          {return p_case_sbase;};

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
      std::map<std::string, XML_TYPE>::iterator type;
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
      boost::property_tree::ptree caseXMLTree =
        xmlTree.get_child("application.GridpackPowergrid");

      BOOST_FOREACH( boost::property_tree::ptree::value_type
          const& caseAttr, caseXMLTree)
      {
        if (caseAttr.first == "CASE_SBASE")
        {
          p_case_sbase    = atoi(caseAttr.second.data().c_str());
        }
        if (caseAttr.first == "CASE_ID")
        {
          p_case_id    = caseAttr.second.data();
        }
      }
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
      boost::property_tree::ptree busXMLTree =
        xmlTree.get_child("application.GridpackPowergrid.Buses");

      // loop through XML Bus subtree "Buses"
      BOOST_FOREACH( boost::property_tree::ptree::value_type
          const& busTree, busXMLTree)
      {
        readBus(busTree);
      }
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
      p_busCollection.push_back(data);
      ++nBuses;
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
      boost::property_tree::ptree branchXMLTree =
        xmlTree.get_child("application.GridpackPowergrid.Branches");

      // loop through XML Bus subtree "Buses"
      BOOST_FOREACH( boost::property_tree::ptree::value_type const& branchTree,
          branchXMLTree)
      {
        readBranch(branchTree);
      }
    }

    void readBranch(boost::property_tree::ptree::value_type const & branchTree)
    {
      boost::shared_ptr<gridpack::component::DataCollection>
        data(new gridpack::component::DataCollection);

      BOOST_FOREACH( boost::property_tree::ptree::value_type branchAttr,
          branchTree.second)
      {

        if (branchAttr.first == "TransmissionElements")
        {
          BOOST_FOREACH( boost::property_tree::ptree::value_type lineSet,
              branchAttr.second)
          {
            BOOST_FOREACH( boost::property_tree::ptree::value_type lineAttr,
                lineSet.second)
            {
              loadCollection(data, lineAttr);
            }
          }
        } else {
          loadCollection(data, branchAttr);
        }
      }
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

    void createNetwork(void)
    {
      int t_create = p_timer->createCategory("Parser:createNetwork");
      p_timer->start(t_create);
      int me(p_network->communicator().rank());
      int nprocs(p_network->communicator().size());
      int i;
      // Exchange information on number of buses and branches on each
      // processor
      int sbus[nprocs], sbranch[nprocs];
      int nbus[nprocs], nbranch[nprocs];
      for (i=0; i<nprocs; i++) {
        sbus[i] = 0;
        sbranch[i] = 0;
      }
      sbus[me] = p_busCollection.size();
      sbranch[me] = p_branchCollection.size();
      MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
      int ierr;
      ierr = MPI_Allreduce(sbus,nbus,nprocs,MPI_INT,MPI_SUM,comm);
      ierr = MPI_Allreduce(sbranch,nbranch,nprocs,MPI_INT,MPI_SUM,comm);
      // evaluate offsets for buses and branches
      int offset_bus[nprocs], offset_branch[nprocs];
      offset_bus[0] = 0;
      offset_branch[0] = 0;
      for (i=1; i<nprocs; i++) {
        offset_bus[i] = offset_bus[i-1]+nbus[i-1];
        offset_branch[i] = offset_branch[i-1]+nbranch[i-1];
      }

      int numBus = p_busCollection.size();
      for (i=0; i<numBus; i++) {
        int idx;
        p_busCollection[i]->getValue(BUS_NUMBER,&idx);
        p_network->addBus(idx);
        p_network->setGlobalBusIndex(i,i+offset_bus[me]);
        *(p_network->getBusData(i)) = *(p_busCollection[i]);
        p_network->getBusData(i)->addValue(CASE_ID,p_case_id);
        p_network->getBusData(i)->addValue(CASE_SBASE,p_case_sbase);
      }
      int numBranch = p_branchCollection.size();
      for (i=0; i<numBranch; i++) {
        int idx1, idx2;
        p_branchCollection[i]->getValue(BRANCH_FROMBUS,&idx1);
        p_branchCollection[i]->getValue(BRANCH_TOBUS,&idx2);
        p_network->addBranch(idx1, idx2);
        p_network->setGlobalBranchIndex(i,i+offset_branch[me]);
        int g_idx1, g_idx2;
#ifdef OLD_MAP
        std::map<int, int>::iterator it;
#else
        boost::unordered_map<int, int>::iterator it;
#endif
        it = p_busMap.find(idx1);
        g_idx1 = it->second;
        it = p_busMap.find(idx2);
        g_idx2 = it->second;
        *(p_network->getBranchData(i)) = *(p_branchCollection[i]);
        p_network->getBranchData(i)->addValue(CASE_ID,p_case_id);
        p_network->getBranchData(i)->addValue(CASE_SBASE,p_case_sbase);
      }
      p_configExists = true;
#if 1
      // debug
      printf("Number of buses: %d\n",numBus);
      for (i=0; i<numBus; i++) {
        printf("Dumping bus: %d\n",i);
        p_network->getBusData(i)->dump();
      }
      printf("Number of branches: %d\n",numBranch);
      for (i=0; i<numBranch; i++) {
        printf("Dumping branch: %d\n",i);
        p_network->getBranchData(i)->dump();
      }
#endif
      p_busCollection.clear();
      p_branchCollection.clear();
      p_timer->stop(t_create);
    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
    boost::shared_ptr<_network> p_network;

    int                      nBuses;
    int                      nBranches;

    // Vector of data collection objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >
      p_busCollection;
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> >
      p_branchCollection;
    std::map<std::string, XML_TYPE> typeMap;

    std::string              p_case_id;
    int                      p_case_sbase;
    gridpack::utility::CoarseTimer *p_timer;
}; /* end of GOSS_parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* GOSS_PARSER_HPP_ */
