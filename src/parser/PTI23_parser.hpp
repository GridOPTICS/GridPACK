/*
 * PTI23parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: kglass
 */

#ifndef PTI23_PARSER_HPP_
#define PTI23_PARSER_HPP_

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt

namespace gridpack {
namespace parser {


template <class BUS, class BRANCH>
class PTI23_parser {
public:

    PTI23_parser(boost::shared_ptr<gridpack::network::BaseNetwork<BUS, BRANCH> > network){};
    virtual ~PTI23_parser(){};
	/*
	 * A case is the collection of all data associated with a PTI23 file.
	 * Each case is a a vector of data_set objects the contain all the data
	 * associated with a partition of the PTI file. For example, the bus
	 * data in the file constitutes a data_set. Each data_set is a vector of
	 * gridpack::component::DataCollection objects. Each of these objects
	 * contain a single instance of the data associated with a data_set. For
	 * example, each line of the bus partition corresponds to a single
	 * DataCollection object.
	 */
    std::vector<data_set> * getCase(const std::string & fileName) {

        std::ifstream            input;
        input.open(fileName.c_str());
  	std::vector<data_set> * case_data    = new std::vector<data_set>;
        if (!input.is_open()) {
            throw gridpack::Exception("failed to open case data file");
        }

        find_case(case_data, input);
        find_buses(case_data, input);
        find_loads(case_data, input);
        find_generators(case_data, input);
        find_branches(case_data, input);
        find_transformer(case_data, input);
        find_area(case_data, input);
        find_shunt(case_data, input);
        find_imped_corr(case_data, input);
        find_multi_section(case_data, input);
        find_zone(case_data, input);
        find_interarea(case_data, input);
        find_owner(case_data, input);

        return case_data;
    }
protected:

    void find_case(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                                           case_set;
        std::string                                        line;
        std::vector<gridpack::component::DataCollection>   case_instance;

        gridpack::component::DataCollection                data;

        std::getline(input, line);
        std::vector<std::string>  split_line;

        boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(" "), boost::token_compress_on);

        // CASE_IC             "IC"                   ranged integer
        data.setValue(CASE_ID, atoi(split_line[0].c_str()));
        case_instance.push_back(data);

        // CASE_SBASE          "SBASE"                float
        data.setValue(CASE_SBASE, atof(split_line[1].c_str()));
        case_instance.push_back(data);

/*  These do not appear in the dictionary
        // CASE_RECORD2        "RECORD2"              string
        std::getline(input, line);
        data.setValue(CASE_RECORD2, line.c_str());
        case_instance.push_back(data);

        // CASE_RECORD3        "RECORD3"              string
        std::getline(input, line);
        data.setValue(CASE_RECORD3, line.c_str());
        case_instance.push_back(data);
*/
        case_set.push_back(case_instance);
        case_data->push_back(case_set);

    }

    void find_buses(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        bus_set;
        std::string          line;
        int                  index = 0;
        int                  o_idx;
        std::getline(input, line);


        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   bus_instance;
            gridpack::component::DataCollection     data;

            // BUS_I               "I"                   integer
            o_idx = atoi(split_line[0].c_str());
            data.setValue(BUS_NUMBER, o_idx);
            bus_instance.push_back(data);

            // BUS_NAME             "IDE"                 ranged integer
            data.setValue(BUS_NAME, line.c_str());
            bus_instance.push_back(data);

            // BUS_BASEKV           "BASKV"               float
            data.setValue(BUS_BASEKV, atof(split_line[10].c_str()));
            bus_instance.push_back(data);

            // BUS_TYPE               "I"                   integer
            data.setValue(BUS_TYPE, atoi(split_line[0].c_str()));
            bus_instance.push_back(data);

            // BUS_SHUNT_GL              "GL"                  float
            data.setValue(BUS_SHUNT_GL, atof(split_line[4].c_str()));
            bus_instance.push_back(data);

            // BUS_SHUNT_BL              "BL"                  float
            data.setValue(BUS_SHUNT_BL, atof(split_line[5].c_str()));
            bus_instance.push_back(data);

            // BUS_AREA            "ZONE"                integer
            data.setValue(BUS_AREA, atoi(split_line[11].c_str()));
            bus_instance.push_back(data);

            // BUS_ZONE            "ZONE"                integer
            data.setValue(BUS_ZONE, atoi(split_line[11].c_str()));
            bus_instance.push_back(data);

            // BUS_VOLTAGE_MAG              "PL"                  float
            data.setValue(BUS_VOLTAGE_MAG, atof(split_line[2].c_str()));
            bus_instance.push_back(data);

            // BUS_VOLTAGE_ANG              "QL"                  float
            data.setValue(BUS_VOLTAGE_ANG, atof(split_line[3].c_str()));
            bus_instance.push_back(data);

            // BUS_OWNER              "IA"                  integer
            data.setValue(BUS_OWNER, atoi(split_line[6].c_str()));
            bus_instance.push_back(data);

            network->addBus(o_idx);
            boost::shared_ptr<gridpack::component::DataCollection> netData =
                     network->getBusData(index);
            // This only works because process 0 is only process adding buses
            network->setGlobalBusIndex(index,index);

            *netData = data;
            ++index;

            bus_set.push_back(bus_instance);
            std::getline(input, line);
        }
        case_data->push_back(bus_set);
    }

void find_loads(std::vector<data_set>  * case_data, std::ifstream & input)
{
    data_set                        load_set;
    std::string          line;
    std::getline(input, line); //this should be the first line of the block
    std::getline(input, line);

    while(line[0] != TERM_CHAR) {
        std::vector<std::string>  split_line;
        boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
        std::vector<gridpack::component::DataCollection>   load_instance;
        gridpack::component::DataCollection          data;

        // LOAD_BUSNUMBER               "I"                   integer
        data.setValue(LOAD_BUSNUMBER, atoi(split_line[0].c_str()));
        load_instance.push_back(data);

        // LOAD_ID              "ID"                  integer
        data.setValue(LOAD_ID, atoi(split_line[1].c_str()));
        load_instance.push_back(data);

        // LOAD_STATUS              "ID"                  integer
        data.setValue(LOAD_STATUS, atoi(split_line[1].c_str()));
        load_instance.push_back(data);

        // LOAD_AREA            "ZONE"                integer
        data.setValue(LOAD_AREA, atoi(split_line[11].c_str()));
        load_instance.push_back(data);

        // LOAD_ZONE            "ZONE"                integer
        data.setValue(LOAD_ZONE, atoi(split_line[11].c_str()));
        load_instance.push_back(data);

        // LOAD_PL              "PG"                  float
        data.setValue(LOAD_PL, atof(split_line[2].c_str()));
        load_instance.push_back(data);

        // LOAD_QL              "QG"                  float
        data.setValue(LOAD_QL, atof(split_line[3].c_str()));
        load_instance.push_back(data);

        // LOAD_IP              "QT"                  float
        data.setValue(LOAD_IP, atof(split_line[4].c_str()));
        load_instance.push_back(data);

        // LOAD_IQ              "QB"                  float
        data.setValue(LOAD_IQ, atof(split_line[5].c_str()));
        load_instance.push_back(data);

        // LOAD_YP              "VS"                  float
        data.setValue(LOAD_YP, atof(split_line[6].c_str()));
        load_instance.push_back(data);

        // LOAD_YQ            "IREG"                integer
        data.setValue(LOAD_YQ, atoi(split_line[7].c_str()));
        load_instance.push_back(data);

        // LOAD_OWNER              "IA"                  integer
        data.setValue(LOAD_OWNER, atoi(split_line[6].c_str()));
        load_instance.push_back(data);

        load_set.push_back(load_instance);
        std::getline(input, line);
    }
    case_data->push_back(load_set);
}

void find_generators(std::vector<data_set>  * case_data, std::ifstream & input)
{
    data_set                        generator_set;
    std::string          line;
    std::getline(input, line); //this should be the first line of the block
    std::getline(input, line);

    while(line[0] != TERM_CHAR) {
        std::vector<std::string>  split_line;
        boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
        std::vector<gridpack::component::DataCollection>   gen_instance;
        gridpack::component::DataCollection          data;

        // GENERATOR_BUSNUMBER               "I"                   integer
        data.setValue(GENERATOR_BUSNUMBER, atoi(split_line[0].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_ID              "ID"                  integer
        data.setValue(GENERATOR_ID, atoi(split_line[1].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_PG              "PG"                  float
        data.setValue(GENERATOR_PG, atof(split_line[2].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_QG              "QG"                  float
        data.setValue(GENERATOR_QG, atof(split_line[3].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_QMAX              "QT"                  float
        data.setValue(GENERATOR_QMAX, atof(split_line[4].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_QMIN              "QB"                  float
        data.setValue(GENERATOR_QMIN, atof(split_line[5].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_VS              "VS"                  float
        data.setValue(GENERATOR_VS, atof(split_line[6].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_IREG            "IREG"                integer
        data.setValue(GENERATOR_IREG, atoi(split_line[7].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_MBASE           "MBASE"               float
        data.setValue(GENERATOR_MBASE, atof(split_line[8].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_ZSORCE              "ZR"                  float
        data.setValue(GENERATOR_ZSORCE, atof(split_line[9].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_XTRAN              "ZX"                  float
        data.setValue(GENERATOR_XTRAN, atof(split_line[10].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_GTAP              "RT"                  float
        data.setValue(GENERATOR_GTAP, atof(split_line[11].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_XT              "XT"                  float
        data.setValue(GENERATOR_STAT, atof(split_line[12].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_RMPCT           "RMPCT"               float
        data.setValue(GENERATOR_RMPCT, atof(split_line[15].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_PMAX              "PT"                  float
        data.setValue(GENERATOR_PMAX, atof(split_line[16].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_PMIN              "PB"                  float
        data.setValue(GENERATOR_PMIN, atof(split_line[17].c_str()));
        gen_instance.push_back(data);

        // GENERATOR_OWNER              "IA"                  integer
        data.setValue(GENERATOR_OWNER, atoi(split_line[6].c_str()));
        gen_instance.push_back(data);

        generator_set.push_back(gen_instance);
        std::getline(input, line);
    }
    case_data->push_back(generator_set);
}

    void find_branches(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        branch_set;
        std::string line;
        int  index   = 0;
        int  o_idx1, o_idx2;

        std::getline(input, line); //this should be the first line of the block
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   branch_instance;
            gridpack::component::DataCollection          data;

            // BRANCH_FROMBUS            "I"                   integer
            o_idx1 = atoi(split_line[0].c_str());
            data.setValue(BRANCH_FROMBUS, o_idx1);
            branch_instance.push_back(data);

            // BRANCH_TOBUS            "J"                   integer
            o_idx2 = atoi(split_line[1].c_str());
            data.setValue(BRANCH_TOBUS, o_idx2);
            branch_instance.push_back(data);

            // BRANCH_CKT          "CKT"                 character
            data.setValue(BRANCH_CKT, (split_line[2].c_str()));
            branch_instance.push_back(data);

            // BRANCH_R            "R"                   float
            data.setValue(BRANCH_R, atof(split_line[3].c_str()));
            branch_instance.push_back(data);

            // BRANCH_X            "X"                   float
            data.setValue(BRANCH_X, atof(split_line[4].c_str()));
            branch_instance.push_back(data);

            // BRANCH_B            "B"                   float
            data.setValue(BRANCH_B, atof(split_line[5].c_str()));
            branch_instance.push_back(data);

            // BRANCH_RATING_A        "RATEA"               float
            data.setValue(BRANCH_RATING_A, atof(split_line[6].c_str()));
            branch_instance.push_back(data);

            // BBRANCH_RATING_        "RATEB"               float
            data.setValue(BRANCH_RATING_B, atof(split_line[7].c_str()));
            branch_instance.push_back(data);

            // BRANCH_RATING_C        "RATEC"               float
            data.setValue(BRANCH_RATING_C, atof(split_line[8].c_str()));
            branch_instance.push_back(data);

            // BRANCH_SHUNT_ADMTTNC_G1        "RATIO"               float
            data.setValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            // BRANCH_SHUNT_ADMTTNC_B1        "RATIO"               float
            data.setValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            // BRANCH_SHUNT_ADMTTNC_G2        "RATIO"               float
            data.setValue(BRANCH_SHUNT_ADMTTNC_G2, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            // BRANCH_SHUNT_ADMTTNC_B2        "RATIO"               float
            data.setValue(BRANCH_SHUNT_ADMTTNC_B2, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            // BRANCH_STATUS        "ANGLE"               float
            data.setValue(BRANCH_STATUS, atoi(split_line[10].c_str()));
            branch_instance.push_back(data);

            // BRANCH_LENGTH           "GI"                  float
            data.setValue(BRANCH_LENGTH, atof(split_line[11].c_str()));
            branch_instance.push_back(data);

            // BRANCH_OWNER           "ST"                  integer
            data.setValue(BRANCH_OWNER, atoi(split_line[15].c_str()));
            branch_instance.push_back(data);

            network->addBranch(o_idx1, o_idx2);
            boost::shared_ptr<gridpack::component::DataCollection> netData =
                     network->getBranchData(index);
            // This only works because process 0 is only process adding branches
            network->setGlobalBranchIndex(index,index);
            *netData = data;
            ++index;

            branch_set.push_back(branch_instance);
            std::getline(input, line);
        }
        case_data->push_back(branch_set);
    }

    void find_transformer(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        transformer_set;
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::cout << "transformer block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   transformer_instance;
            gridpack::component::DataCollection          data;

            /*
             * type: integer
             * #define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"
             */
            data.setValue(TRANSFORMER_BUS1, atoi(split_line[0].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"
             */
            data.setValue(TRANSFORMER_BUS2, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_BUS3 "TRANSFORMER_BUS3"
             */
            data.setValue(TRANSFORMER_BUS3, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: string
             * #define TRANSFORMER_CKT "TRANSFORMER_CKT"
             */
            data.setValue(TRANSFORMER_CKT, split_line[2].c_str());
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_CW "TRANSFORMER_CW"
X            */
            data.setValue(TRANSFORMER_CW, atoi(split_line[3].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_CZ "TRANSFORMER_CZ"
             */
            data.setValue(TRANSFORMER_CZ, atoi(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_CM "TRANSFORMER_CM"
             */
            data.setValue(TRANSFORMER_CM, atoi(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"
             */
            data.setValue(TRANSFORMER_MAG1, atof(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"
             */
            data.setValue(TRANSFORMER_MAG2, atof(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_NMETR "TRANSFORMER_NMETR"
             */
            data.setValue(TRANSFORMER_NMETR, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: string
             * #define TRANSFORMER_NAME "TRANSFORMER_NAME"
             */
            data.setValue(TRANSFORMER_NAME, split_line[2].c_str());
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_STATUS "TRANSFORMER_STATUS"
             *
             */
            data.setValue(TRANSFORMER_STATUS, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_OWNER "TRANSFORMER_OWNER"
             */
            data.setValue(TRANSFORMER_OWNER, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"
             */
            data.setValue(TRANSFORMER_R1_2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"
             */
            data.setValue(TRANSFORMER_X1_2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"
             */
            data.setValue(TRANSFORMER_SBASE1_2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_R2_3 "TRANSFORMER_R2_3"
             */
            data.setValue(TRANSFORMER_R2_3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_X2_3 "TRANSFORMER_X2_3"
             */
            data.setValue(TRANSFORMER_X2_3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_SBASE2_3 "TRANSFORMER_SBASE2_3"
             */
            data.setValue(TRANSFORMER_SBASE2_3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_R3_1 "TRANSFORMER_R3_1"
             */
            data.setValue(TRANSFORMER_R3_1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_X3_1 "TRANSFORMER_X3_1"
             */
            data.setValue(TRANSFORMER_X3_1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_SBASE3_1 "TRANSFORMER_SBASE3_1"
             */
            data.setValue(TRANSFORMER_SBASE3_1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMSTAR "TRANSFORMER_VMSTAR"
             */
            data.setValue(TRANSFORMER_VMSTAR, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_ANSTAR "TRANSFORMER_ANSTAR"
             */
            data.setValue(TRANSFORMER_ANSTAR, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_WINDV1 "TRANSFORMER_WINDV1"
             */
            data.setValue(TRANSFORMER_WINDV1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_NOMV1 "TRANSFORMER_NOMV1"
             */
            data.setValue(TRANSFORMER_NOMV1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_ANG1 "TRANSFORMER_ANG1"
             */
            data.setValue(TRANSFORMER_ANG1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATA1 "TRANSFORMER_RATA1"
             */
            data.setValue(TRANSFORMER_RATA1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATB1 "TRANSFORMER_RATB1"
             */
            data.setValue(TRANSFORMER_RATB1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATC1 "TRANSFORMER_RATC1"
             */
            data.setValue(TRANSFORMER_RATC1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_COD1 "TRANSFORMER_COD1"
             */
            data.setValue(TRANSFORMER_COD1, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_CONT1 "TRANSFORMER_CONT1"
             */
            data.setValue(TRANSFORMER_CONT1, atoi(split_line[0].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RMA1 "TRANSFORMER_RMA1"
             */
            data.setValue(TRANSFORMER_RMA1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RMI1 "TRANSFORMER_RMI1"
             */
            data.setValue(TRANSFORMER_RMI1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMA1 "TRANSFORMER_VMA1"
             */
            data.setValue(TRANSFORMER_VMA1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMI1 "TRANSFORMER_VMI1"
             */
            data.setValue(TRANSFORMER_VMI1, atof(split_line[3].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_NTP1 "TRANSFORMER_NTP1"
             */
            data.setValue(TRANSFORMER_NTP1, atoi(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_TAB1 "TRANSFORMER_TAB1"
             */
            data.setValue(TRANSFORMER_TAB1, atoi(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_CR1 "TRANSFORMER_CR1"
             */
            data.setValue(TRANSFORMER_CR1, atof(split_line[5].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_CX1 "TRANSFORMER_CX1"
             */
            data.setValue(TRANSFORMER_CX1, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_WINDV2 "TRANSFORMER_WINDV2"
             */
            data.setValue(TRANSFORMER_WINDV2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_NOMV2 "TRANSFORMER_NOMV2"
             */
            data.setValue(TRANSFORMER_NOMV2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_ANG2 "TRANSFORMER_ANG2"
             */
            data.setValue(TRANSFORMER_ANG2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATA2 "TRANSFORMER_RATA2"
             */
            data.setValue(TRANSFORMER_RATA2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATA2 "TRANSFORMER_RATB2"
             */
            data.setValue(TRANSFORMER_RATB2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATC2 "TRANSFORMER_RATC2"
             */
            data.setValue(TRANSFORMER_RATC2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_COD2 "TRANSFORMER_COD2"
             */
            data.setValue(TRANSFORMER_COD2, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_CONT2 "TRANSFORMER_CONT2"
             */
            data.setValue(TRANSFORMER_CONT2, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RMA2 "TRANSFORMER_RMA2"
             */
            data.setValue(TRANSFORMER_RMA2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RMI2 "TRANSFORMER_RMI2"
             */
            data.setValue(TRANSFORMER_RMI2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMA2 "TRANSFORMER_VMA2"
             */
            data.setValue(TRANSFORMER_VMA2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMI2 "TRANSFORMER_VMI2"
             */
            data.setValue(TRANSFORMER_VMI2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_NTP2 "TRANSFORMER_NTP2"
             */
            data.setValue(TRANSFORMER_NTP2, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_TAB2 "TRANSFORMER_TAB2"
             */
            data.setValue(TRANSFORMER_TAB2, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_CR2 "TRANSFORMER_CR2"
             */
            data.setValue(TRANSFORMER_CR2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_CX2 "TRANSFORMER_CX2"
             */
            data.setValue(TRANSFORMER_CX2, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_WINDV3 "TRANSFORMER_WINDV3"
             */
            data.setValue(TRANSFORMER_WINDV3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_NOMV3 "TRANSFORMER_NOMV3"
             */
            data.setValue(TRANSFORMER_NOMV3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_ANG3 "TRANSFORMER_ANG3"
             */
            data.setValue(TRANSFORMER_ANG3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATA3 "TRANSFORMER_RATA3"
             */
            data.setValue(TRANSFORMER_RATA3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATB3 "TRANSFORMER_RATB3"
             */
            data.setValue(TRANSFORMER_RATB3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RATC3 "TRANSFORMER_RATC3"
             */
            data.setValue(TRANSFORMER_RATC3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_COD3 "TRANSFORMER_COD3"
             */
            data.setValue(TRANSFORMER_COD3, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_CONT3 "TRANSFORMER_CONT3"
             */
            data.setValue(TRANSFORMER_CONT3, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RMA3 "TRANSFORMER_RMA3"
             */
            data.setValue(TRANSFORMER_RMA3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_RMI3 "TRANSFORMER_RMI3"
             */
            data.setValue(TRANSFORMER_RMI3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMA3 "TRANSFORMER_VMA3"
             */
            data.setValue(TRANSFORMER_VMA3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_VMI3 "TRANSFORMER_VMI3"
             */
            data.setValue(TRANSFORMER_VMI3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_NTP3 "TRANSFORMER_NTP3"
             */
            data.setValue(TRANSFORMER_NTP3, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: integer
             * #define TRANSFORMER_TAB3 "TRANSFORMER_TAB3"
             */
            data.setValue(TRANSFORMER_TAB3, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_CR3 "TRANSFORMER_CR3"
             */
            data.setValue(TRANSFORMER_CR3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            /*
             * type: real float
             * #define TRANSFORMER_CX3 "TRANSFORMER_CX3"
             */
            data.setValue(TRANSFORMER_CX3, atof(split_line[1].c_str()));
            transformer_instance.push_back(data);

            transformer_set.push_back(transformer_instance);
            std::getline(input, line);
        }
        case_data->push_back(transformer_set);
    }

    void find_area(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        area_set;
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   area_instance;
            gridpack::component::DataCollection          data;

            // AREAINTG_NUMBER             "I"                    integer
            data.setValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()));
            area_instance.push_back(data);

            // AREAINTG_ISW           "ISW"                  integer
            data.setValue(AREAINTG_ISW, atoi(split_line[1].c_str()));
            area_instance.push_back(data);

            // AREAINTG_PDES          "PDES"                 float
            data.setValue(AREAINTG_PDES, atof(split_line[2].c_str()));
            area_instance.push_back(data);

            // AREAINTG_PTOL          "PTOL"                 float
            data.setValue(AREAINTG_PTOL, atof(split_line[3].c_str()));
            area_instance.push_back(data);

            // AREAINTG_NAME         "ARNAM"                string
            data.setValue(AREAINTG_NAME, split_line[4].c_str());
            area_instance.push_back(data);

            area_set.push_back(area_instance);
            std::getline(input, line);
        }
        case_data->push_back(area_set);
    }


    /*

    */
    void find_shunt(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        shunt_set;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   shunt_instance;
            gridpack::component::DataCollection          data;

            /*
             * type: integer
             * #define SHUNT_BUSNUMBER "SHUNT_BUSNUMBER"
             */
            data.setValue(SHUNT_BUSNUMBER, atoi(split_line[0].c_str()));
            shunt_instance.push_back(data);

            /*

type: integer
#define SHUNT_MODSW "SHUNT_MODSW"
*/
            data.setValue(SHUNT_MODSW, atoi(split_line[0].c_str()));
            shunt_instance.push_back(data);

            /*

type: real float
#define SHUNT_VSWHI "SHUNT_VSWHI"
             *
             */
            data.setValue(SHUNT_VSWHI, atof(split_line[1].c_str()));
            shunt_instance.push_back(data);

            /*

type: real float
#define SHUNT_VSWLO "SHUNT_VSWLO"
             *
             */
            data.setValue(SHUNT_VSWLO, atof(split_line[2].c_str()));
            shunt_instance.push_back(data);

            /*

type: integer
#define SHUNT_SWREM "SHUNT_SWREM"
             */
            data.setValue(SHUNT_SWREM, atoi(split_line[3].c_str()));
            shunt_instance.push_back(data);

            /*

type: real float
#define SHUNT_RMPCT "SHUNT_RMPCT"
             *
             */
            data.setValue(SHUNT_RMPCT, atof(split_line[4].c_str()));
            shunt_instance.push_back(data);

            /*
type: string
#define SHUNT_RMIDNT "SHUNT_RMIDNT"
             *
             */
            data.setValue(SHUNT_RMIDNT, split_line[5].c_str());
            shunt_instance.push_back(data);

            /*
             * type: real float
             * #define SHUNT_BINIT "SHUNT_BINIT"
             */
            data.setValue(SHUNT_BINIT, atof(split_line[5].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N1 "SHUNT_N1"
             *
             */
            data.setValue(SHUNT_N1, atoi(split_line[6].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N2 "SHUNT_N2"
             */
            data.setValue(SHUNT_N2, atoi(split_line[8].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N3 "SHUNT_N3"
             */
            data.setValue(SHUNT_N3, atoi(split_line[10].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N4 "SHUNT_N4"
             */
            data.setValue(SHUNT_N4, atoi(split_line[12].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N5 "SHUNT_N5"
             */
            data.setValue(SHUNT_N5, atoi(split_line[14].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N6 "SHUNT_N6"
             */
            data.setValue(SHUNT_N6, atoi(split_line[16].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N7 "SHUNT_N7"
             */
            data.setValue(SHUNT_N7, atoi(split_line[18].c_str()));
            shunt_instance.push_back(data);

            /*
type: integer
#define SHUNT_N8 "SHUNT_N8"

             */
            data.setValue(SHUNT_N8, atoi(split_line[20].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B1 "SHUNT_B1"
             */
            data.setValue(SHUNT_B1, atof(split_line[7].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B2 "SHUNT_B2"
             */
            data.setValue(SHUNT_B2, atof(split_line[9].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B3 "SHUNT_B3"
             */
            data.setValue(SHUNT_B3, atof(split_line[11].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B4 "SHUNT_B4"
             */
            data.setValue(SHUNT_B4, atof(split_line[13].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B5 "SHUNT_B5"
             */
            data.setValue(SHUNT_B5, atof(split_line[15].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B6 "SHUNT_B6"
             */
            data.setValue(SHUNT_B6, atof(split_line[17].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B7 "SHUNT_B7"
             */
            data.setValue(SHUNT_B7, atof(split_line[19].c_str()));
            shunt_instance.push_back(data);

            /*
type: real float
#define SHUNT_B8 "SHUNT_B8"
             */
            data.setValue(SHUNT_B8, atof(split_line[21].c_str()));
            shunt_instance.push_back(data);

            shunt_set.push_back(shunt_instance);
            std::getline(input, line);
        }
        case_data->push_back(shunt_set);
    }

    void find_imped_corr(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        imped_corr_set;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   imped_corr_instance;
            gridpack::component::DataCollection          data;

            /*
     * type: integer
     * #define XFMR_CORR_TABLE_NUMBER "XFMR_CORR_TABLE_NUMBER"
             */
            data.setValue(XFMR_CORR_TABLE_NUMBER, atoi(split_line[0].c_str()));
            imped_corr_instance.push_back(data);

            /*
     * type: real float
     * #define XFMR_CORR_TABLE_Ti "XFMR_CORR_TABLE_Ti"
             */
            data.setValue(XFMR_CORR_TABLE_Ti, atoi(split_line[0].c_str()));
            imped_corr_instance.push_back(data);

            /*
     * type: real float
     * #define XFMR_CORR_TABLE_Fi "XFMR_CORR_TABLE_Fi"
             */
            data.setValue(XFMR_CORR_TABLE_Fi, atoi(split_line[0].c_str()));
            imped_corr_instance.push_back(data);

            imped_corr_set.push_back(imped_corr_instance);
            std::getline(input, line);
        }
        case_data->push_back(imped_corr_set);
    }

    /*
     *
     */

    void find_multi_section(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        multi_section;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::cout << "multi section block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   multi_section_instance;
            gridpack::component::DataCollection          data;

            /*
     * type: integer
     * #define MULTI_SEC_LINE_FROMBUS "MULTI_SEC_LINE_FROMBUS"

             */
            data.setValue(MULTI_SEC_LINE_FROMBUS, atoi(split_line[0].c_str()));
            multi_section_instance.push_back(data);

            /*
     * type: integer
     * #define MULTI_SEC_LINE_TOBUS "MULTI_SEC_LINE_TOBUS"

             */
            data.setValue(MULTI_SEC_LINE_TOBUS, atoi(split_line[0].c_str()));
            multi_section_instance.push_back(data);

            /*
     * type: string
     * #define MULTI_SEC_LINE_ID "MULTI_SEC_LINE_ID"

             */
            data.setValue(MULTI_SEC_LINE_ID, split_line[0].c_str());
            multi_section_instance.push_back(data);

            /*
     * type: integer
     * #define MULTI_SEC_LINE_DUMi "MULTI_SEC_LINE_DUMi"
             */
            data.setValue(MULTI_SEC_LINE_DUMi, atoi(split_line[0].c_str()));
            multi_section_instance.push_back(data);

            multi_section.push_back(multi_section_instance);
            std::getline(input, line);
        }
        case_data->push_back(multi_section);
    }
    /*
     * ZONE_I          "I"                       integer
     * ZONE_NAME       "NAME"                    string
     */
    void find_zone(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        zone;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   zone_instance;
            gridpack::component::DataCollection          data;

            /*
             *
             */
            data.setValue(ZONE_NUMBER, atoi(split_line[0].c_str()));
            zone_instance.push_back(data);

            // ZONE_NAME       "NAME"                    string
            data.setValue(ZONE_NAME, split_line[1].c_str());
            zone_instance.push_back(data);

            zone.push_back(zone_instance);
            std::getline(input, line);
        }
        case_data->push_back(zone);
    }

    /*
     *
     * type: integer
     * #define INTERAREA_TRANSFER_FROM "INTERAREA_TRANSFER_FROM"

     * type: integer
     * #define INTERAREA_TRANSFER_TO "INTERAREA_TRANSFER_TO"

     * type: character
     * #define INTERAREA_TRANSFER_TRID "INTERAREA_TRANSFER_TRID"

     * type: real float
     * #define INTERAREA_TRANSFER_PTRAN "INTERAREA_TRANSFER_PTRAN"
     */
    void find_interarea(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        inter_area;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   inter_area_instance;
            gridpack::component::DataCollection          data;

            /*
             *
             */
            data.setValue(INTERAREA_TRANSFER_FROM, atoi(split_line[0].c_str()));
            inter_area_instance.push_back(data);

            /*
             *
             */
            data.setValue(INTERAREA_TRANSFER_TO, atoi(split_line[0].c_str()));
            inter_area_instance.push_back(data);

            /*
             *
             */
            data.setValue(INTERAREA_TRANSFER_TRID, split_line[0].c_str()[0]);
            inter_area_instance.push_back(data);

            /*
             *
             */
            data.setValue(INTERAREA_TRANSFER_PTRAN, atof(split_line[0].c_str()));
            inter_area_instance.push_back(data);

            inter_area.push_back(inter_area_instance);
            std::getline(input, line);
        }
        case_data->push_back(inter_area);
    }

    /*
     * type: integer
     * #define OWNER_NUMBER "OWNER_NUMBER"

     * type: integer
     * #define OWNER_NAME "OWNER_NAME"
     */
    void find_owner(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        owner;

        std::string          line;
        std::getline(input, line); //this should be the first line of the block

        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   owner_instance;
            gridpack::component::DataCollection          data;

            data.setValue(OWNER_NUMBER, atoi(split_line[0].c_str()));
            owner_instance.push_back(data);

            data.setValue(OWNER_NAME, split_line[1].c_str());
            owner_instance.push_back(data);

            owner.push_back(owner_instance);
            std::getline(input, line);
        }
        case_data->push_back(owner);
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
    boost::shared_ptr<gridpack::network::BaseNetwork<BUS, BRANCH> >
                                 network;

};

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI23PARSER_HPP_ */
