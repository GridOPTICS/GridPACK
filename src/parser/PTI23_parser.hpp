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

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt
#define     CASE_IC             "IC"
//#define     CASE_SBASE          "SBASE"
#define     CASE_RECORD2        "RECORD2"
#define     CASE_RECORD3        "RECORD3"

#define     BUS_I               "I"
#define     BUS_IDE             "IDE"
#define     BUS_PL              "PL"
#define     BUS_QL              "QL"
#define     BUS_GL              "GL"
#define     BUS_BL              "BL"
#define     BUS_IA              "IA"
#define     BUS_VM              "VM"
#define     BUS_VA              "VA"
//#define     BUS_NAME            "NAME"
#define     BUS_BASKV           "BASKV"
//#define     BUS_ZONE            "ZONE"

#define     GEN_I               "I"
#define     GEN_ID              "ID"
#define     GEN_PG              "PG"
#define     GEN_QG              "QG"
#define     GEN_QT              "QT"
#define     GEN_QB              "QB"
#define     GEN_VS              "VS"
#define     GEN_IREG            "IREG"
#define     GEN_MBASE           "MBASE"
#define     GEN_ZR              "ZR"
#define     GEN_ZX              "ZX"
#define     GEN_RT              "RT"
#define     GEN_XT              "XT"
#define     GEN_GTAP            "GTAP"
#define     GEN_STAT            "STAT"
#define     GEN_RMPCT           "RMPCT"
#define     GEN_PT              "PT"
#define     GEN_PB              "PB"

#define     BRANCH_I            "I"
#define     BRANCH_J            "J"
//#define     BRANCH_CKT          "CKT"
//#define     BRANCH_R            "R"
//#define     BRANCH_X            "X"
//#define     BRANCH_B            "B"
#define     BRANCH_RATEA        "RATEA"
#define     BRANCH_RATEB        "RATEB"
#define     BRANCH_RATEC        "RATEC"
#define     BRANCH_RATIO        "RATIO"
#define     BRANCH_ANGLE        "ANGLE"
#define     BRANCH_GI           "GI"
#define     BRANCH_BI           "BI"
#define     BRANCH_GJ           "GJ"
#define     BRANCH_BJ           "BJ"
#define     BRANCH_ST           "ST"

#define     TRANSF_I           "I"
#define     TRANSF_J           "J"
#define     TRANSF_CKT         "CKT"
#define     TRANSF_ICONT       "ICONT"
#define     TRANSF_RMA         "RMA"
#define     TRANSF_RMI         "RMI"
#define     TRANSF_VMA         "VMA"
#define     TRANSF_VMI         "VMI"
#define     TRANSF_STEP        "STEP"
#define     TRANSF_TABLE       "TABLE"

#define     AREA_I             "I"
#define     AREA_ISW           "ISW"
#define     AREA_PDES          "PDES"
#define     AREA_PTOL          "PTOL"
#define     AREA_ARNAM         "ARNAM"

#define     DL_I               "I"
#define     DL_MDC             "MDC"
#define     DL_RDC             "RDC"
#define     DL_SETVL           "SETVL"
#define     DL_VSCHD           "VSCHD"
#define     DL_VCMOD           "VCMOD"
#define     DL_RCOMP           "RCOMP"
#define     DL_DELTI           "DELTI"
#define     DL_METER           "METER"
#define     DL_IPR             "IPR"
#define     DL_NBR             "NBR"
#define     DL_ALFMAX          "ALFMAX"
#define     DL_ALFMN           "ALFMN"
#define     DL_RCR             "RCR"
#define     DL_XCR             "XCR"
#define     DL_EBASR           "EBASR"
#define     DL_TRR             "TRR"
#define     DL_TAPR            "TAPR"
#define     DL_TPMXR           "TPMXR"
#define     DL_TPMNR           "TPMNR"
#define     DL_TSTPR           "TSTPR"
#define     DL_IPI             "IPI"
#define     DL_NBI             "NBI"
#define     DL_GAMMX           "GAMMX"
#define     DL_GAMMN           "GAMMN"
#define     DL_RCI             "RCI"
#define     DL_XCI             "XCI"
#define     DL_EBASI           "EBASI"
#define     DL_TRI             "TRI"
#define     DL_TAPI            "TAPI"
#define     DL_TPMXI           "TPMXI"
#define     DL_TPMNI           "TPMNI"
#define     DL_TSTPI           "TSTPI"

#define     SHUNT_I            "I"
//#define     SHUNT_MODSW        "MODSW"
//#define     SHUNT_VSWHI        "VSWHI"
#define     SHUNT_VSWLD        "VSWLO"
//#define     SHUNT_SWREM        "SWREM"
#define     SHUNT_VDES         "VDES"
//#define     SHUNT_BINIT        "BINIT"
#define     SHUNT_UNK          "UNK"
//#define     SHUNT_N1           "N1"
//#define     SHUNT_B1           "B1"
//#define     SHUNT_N2           "N2"
//#define     SHUNT_B2           "B2"
//#define     SHUNT_N3           "N3"
//#define     SHUNT_B3           "B3"
//#define     SHUNT_N4           "N4"
//#define     SHUNT_B4           "B4"
//#define     SHUNT_N5           "N5"
//#define     SHUNT_B5           "B5"
//#define     SHUNT_N6           "N6"
//#define     SHUNT_B6           "B6"
//#define     SHUNT_N7           "N7"
//#define     SHUNT_B7           "B7"
//#define     SHUNT_N8           "N8"
//#define     SHUNT_B8           "B8"

#define     IMPED_I            "I"

#define     M_TERM_I           "I"

#define     M_SEC_I            "I"

#define     ZONE_I             "I"
//#define     ZONE_NAME          "NAME"

#define     I_AREA_I           "I"

#define     OWNER_I            "I"
//#define     OWNER_NAME         "NAME"

#include "gridpack/parser/dictionary.hpp"


namespace gridpack {
namespace parser {

class PTI23_parser {
public:

    PTI23_parser(){};
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
        find_generators(case_data, input);
        find_branches(case_data, input);
        find_transformer(case_data, input);
        find_area(case_data, input);
        find_dc_line(case_data, input);
        find_shunt(case_data, input);
        find_imped_corr(case_data, input);
        find_multi_terminal(case_data, input);
        find_multi_section(case_data, input);
        find_zone(case_data, input);
        find_interarea(case_data, input);
        find_owner(case_data, input);
        find_device_data(case_data, input);

        return case_data;
    }
protected:

   /*
#define     CASE_ID             "CASE_ID"
#define     CASE_SBASE          "SBASE"
#define     CASE_RECORD2        "RECORD2"
#define     CASE_RECORD3        "RECORD3"
   */
    void find_case(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                                           case_set;
        std::string                                        line;
        std::vector<gridpack::component::DataCollection>   case_instance;

        gridpack::component::DataCollection                data;

        std::getline(input, line);
        std::vector<std::string>  split_line;

        boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(" "), boost::token_compress_on);

        data.setValue(CASE_ID, atoi(split_line[0].c_str()));
        case_instance.push_back(data);

        data.setValue(CASE_SBASE, atoi(split_line[1].c_str()));
        case_instance.push_back(data);

        std::getline(input, line);
        data.setValue(CASE_RECORD2, line.c_str());
        case_instance.push_back(data);

        std::getline(input, line);
        data.setValue(CASE_RECORD3, line.c_str());
        case_instance.push_back(data);

        case_set.push_back(case_instance);
        case_data->push_back(case_set);

std::cout << "last case line " << line << std::endl;
    }

    /*
#define     BUS_NUMBER          "BUS_NUMBER"
#define     BUS_TYPE            "BUS_TYPE
#define     BUS_PL              "PL"
#define     BUS_QL              "QL"
#define     BUS_SHUNT_GL        "BUS_SHUNT_GL"
#define     BUS_SHUNT_BL        "BUS_SHUNT_BL"
#define     BUS_IA              "IA"
#define     BUS_VM              "VM"
#define     BUS_VA              "VA"
#define     BUS_NAME            "NAME"
#define     BUS_BASKV           "BASKV"
#define     BUS_ZONE            "ZONE"
    */
    void find_buses(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        bus_set;
        std::string          line;

	std::getline(input, line);
std::cout << "bus block " << line << std::endl;


        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   bus_instance;
            gridpack::component::DataCollection     data;

            data.setValue(BUS_NUMBER, atoi(split_line[0].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_TYPE, atoi(split_line[1].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_PL, atof(split_line[2].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_QL, atof(split_line[3].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_SHUNT_GL, atof(split_line[4].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_SHUNT_BL, atof(split_line[5].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_IA, atof(split_line[6].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_VM, atof(split_line[7].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_VA, atof(split_line[8].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_NAME, split_line[9].c_str());
            bus_instance.push_back(data);

            data.setValue(BUS_BASKV, atof(split_line[10].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_ZONE, atoi(split_line[11].c_str()));
            bus_instance.push_back(data);

            bus_set.push_back(bus_instance);
            std::getline(input, line);
        }
        case_data->push_back(bus_set);
    }


    /*

#define     GENERATOR_BUSNBUMBER      "GENERATOR_BUSNBUMBER"
#define     GENERATOR_ID              "GENERATOR_ID"
#define     GENERATOR_PG              "GENERATOR_PG"
#define     GENERATOR_QG              "GENERATOR_QG"
#define     GENERATOR_QT              "GENERATOR_QT"
#define     GEN_QB              "QB"
#define     GENERATOR_VS              "GENERATOR_VS"
#define     GENERATOR_IREG            "GENERATOR_IREG"
#define     GENERATOR_MBASE           "GENERATOR_MBASE"
#define     GEN_ZR              "ZR"
#define     GEN_ZX              "ZX"
#define     GEN_RT              "RT"
#define     GEN_XT              "XT"
#define     GENERATOR_GTAP            "GENERATOR_GTAP"
#define     GENERATOR_STAT            "GENERATOR_STAT"
#define     GENERATOR_RMPCT           "GENERATOR_RMPCT"
#define     GEN_PT              "PT"
#define     GEN_PB              "PB"
    */
    void find_generators(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        generator_set;
        std::string          line;
        std::getline(input, line); //this should be the first line of the block
std::cout << "generator block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   gen_instance;
            gridpack::component::DataCollection          data;

            data.setValue(GENERATOR_BUSNUMBER, atof(split_line[0].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_ID, atof(split_line[1].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_PG, atof(split_line[2].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_QG, atof(split_line[3].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QT, atof(split_line[4].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QB, atof(split_line[5].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_VS, atof(split_line[6].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_IREG, atof(split_line[7].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_MBASE, atof(split_line[8].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ZR, atof(split_line[9].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ZX, atof(split_line[10].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_RT, atof(split_line[11].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_XT, atof(split_line[12].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_GTAP, atof(split_line[13].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_STAT, atof(split_line[14].c_str()));
            gen_instance.push_back(data);

            data.setValue(GENERATOR_RMPCT, atof(split_line[15].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PT, atof(split_line[16].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PB, atof(split_line[17].c_str()));
            gen_instance.push_back(data);

            generator_set.push_back(gen_instance);
            std::getline(input, line);
        }
        case_data->push_back(generator_set);
    }


    /*

#define     BRANCH_FROMBUS            "BRANCH_FROMBUS"
#define     BRANCH_TOBUS              "BRANCH_TOBUS"
#define     BRANCH_CKT                "BRANCH_CKT"
#define     BRANCH_R                  "BRANCH_R"
#define     BRANCH_X                  "BRANCH_X"
#define     BRANCH_B                  "BRANCH_B"
#define     BRANCH_RATING_A           "BRANCH_RATING_A"
#define     BRANCH_RATING_B           "BRANCH_RATING_B"
#define     BRANCH_RATING_C           "BRANCH_RATING_C"
#define     BRANCH_RATIO        "RATIO"
#define     BRANCH_ANGLE        "ANGLE"
#define     BRANCH_SHUNT_ADMTTNC_G1   "BRANCH_SHUNT_ADMTTNC_G1"
#define     BRANCH_SHUNT_ADMTTNC_B1   "BRANCH_SHUNT_ADMTTNC_B1"
#define     BRANCH_SHUNT_ADMTTNC_G2   "BRANCH_SHUNT_ADMTTNC_G2"
#define     BRANCH_SHUNT_ADMTTNC_B2   "BRANCH_SHUNT_ADMTTNC_B2"
#define     BRANCH_STATUS             "BRANCH_STATUS"
    */
    void find_branches(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        branch_set;
        std::string line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "branch block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   branch_instance;
            gridpack::component::DataCollection          data;

            data.setValue(BRANCH_FROMBUS, atoi(split_line[0].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_TOBUS, atoi(split_line[1].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_CKT, atof(split_line[2].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_R, atof(split_line[3].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_X, atof(split_line[4].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_B, atof(split_line[5].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATING_A, atof(split_line[6].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATING_B, atof(split_line[7].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATING_C, atof(split_line[8].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATIO, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_ANGLE, atof(split_line[10].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[11].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[12].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_SHUNT_ADMTTNC_G2, atof(split_line[13].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_SHUNT_ADMTTNC_B2, atof(split_line[14].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_STATUS, atoi(split_line[15].c_str()));
            branch_instance.push_back(data);

            branch_set.push_back(branch_instance);
            std::getline(input, line);
        }
        case_data->push_back(branch_set);
    }


    /*
#define     TRANSFORMER_BUS1       "TRANSFORMER_BUS1"
#define     TRANSFORMER_BUS2       "TRANSFORMER_BUS2"
#define     TRANSFORMER_CKT        "TRANSFORMER_CKT"
#define     TRANSF_ICONT       "ICONT"
#define     TRANSF_RMA         "RMA"
#define     TRANSF_RMI         "RMI"
#define     TRANSF_VMA         "VMA"
#define     TRANSF_VMI         "VMI"
#define     TRANSF_STEP        "STEP"
#define     TRANSF_TABLE       "TABLE"
    */
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

            data.setValue(TRANSFORMER_BUS1, atoi(split_line[0].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSFORMER_BUS2, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSFORMER_CKT, atoi(split_line[2].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_ICONT, atoi(split_line[3].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_RMA, atoi(split_line[4].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_RMI, atoi(split_line[5].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_VMA, atof(split_line[6].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_VMI, atof(split_line[7].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_STEP, atoi(split_line[8].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_TABLE, atoi(split_line[9].c_str()));
            transformer_instance.push_back(data);

            transformer_set.push_back(transformer_instance);
            std::getline(input, line);
        }
        case_data->push_back(transformer_set);
    }


    /*
#define     AREAINTG_NUMBER       "AREAINTG_NUMBER"
#define     AREAINTG_ISW          "AREAINTG_ISW"
#define     AREAINTG_PDES         "AREAINTG_PDES"
#define     AREAINTG_PTOL         "AREAINTG_PTOL"
#define     AREAINTG_NAME         "AREAINTG_NAME"
    */
    void find_area(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        area_set;
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "area  block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   area_instance;
            gridpack::component::DataCollection          data;

            data.setValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()));
            area_instance.push_back(data);

            data.setValue(AREAINTG_ISW, atoi(split_line[1].c_str()));
            area_instance.push_back(data);

            data.setValue(AREAINTG_PDES, atof(split_line[2].c_str()));
            area_instance.push_back(data);

            data.setValue(AREAINTG_PTOL, atof(split_line[3].c_str()));
            area_instance.push_back(data);

            data.setValue(AREAINTG_NAME, split_line[4].c_str());
            area_instance.push_back(data);

            area_set.push_back(area_instance);
            std::getline(input, line);
        }
        case_data->push_back(area_set);
    }


    /*
     *  0 #define     DL_I               "I"
     *  1 #define     DL_MDC             "MDC"
     *  2 #define     DL_RDC             "RDC"
     *  3 #define     DL_SETVL           "SETVL"
     *  4 #define     DL_VSCHD           "VSCHD"
     *  5 #define     DL_VCMOD           "VCMOD"
     *  6 #define     DL_RCOMP           "RCOMP"
     *  7 #define     DL_DELTI           "DELTI"
     *  8 #define     DL_METER           "METER"
     *  9 #define     DL_IPR             "IPR"
     * 10 #define     DL_NBR             "NBR"
     * 11 #define     DL_ALFMAX          "ALFMAX"
     * 12 #define     DL_ALFMN           "ALFMN"
     * 13 #define     DL_RCR             "RCR"
     * 14 #define     DL_XCR             "XCR"
     * 15 #define     DL_EBASR           "EBASR"
     * 16 #define     DL_TRR             "TRR"
     * 17 #define     DL_TAPR            "TAPR"
     * 18 #define     DL_TPMXR           "TPMXR"
     * 19 #define     DL_TPMNR           "TPMNR"
     * 20 #define     DL_TSTPR           "TSTPR"
     * 21 #define     DL_IPI             "IPI"
     * 22 #define     DL_NBI             "NBI"
     * 23 #define     DL_GAMMX           "GAMMX"
     * 24 #define     DL_GAMMN           "GAMMN"
     * 25 #define     DL_RCI             "RCI"
     * 26 #define     DL_XCI             "XCI"
     * 27 #define     DL_EBASI           "EBASI"
     * 28 #define     DL_TRI             "TRI"
     * 29 #define     DL_TAPI            "TAPI"
     * 30 #define     DL_TPMXI           "TPMXI"
     * 31 #define     DL_TPMNI           "TPMNI"
     * 32 #define     DL_TSTPI           "TSTPI"

     */
    void find_dc_line(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        dc_line_set;
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "dc line block " << line << std::endl;

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   dc_line_instance;
            gridpack::component::DataCollection          data;

            data.setValue(DL_I, atoi(split_line[0].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_MDC, atoi(split_line[1].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_RDC, atof(split_line[2].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_SETVL, atof(split_line[3].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_VSCHD, atof(split_line[4].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_VCMOD, atof(split_line[5].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_RCOMP, atof(split_line[6].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_DELTI, atof(split_line[7].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_METER, split_line[8].c_str());
            dc_line_instance.push_back(data);

            getline(input, line);
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

            data.setValue(DL_IPR, atoi(split_line[0].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_NBR, atoi(split_line[1].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_ALFMAX, atof(split_line[2].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_ALFMN, atof(split_line[3].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_RCR, atof(split_line[4].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_XCR, atof(split_line[5].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_EBASR, atof(split_line[6].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TRR, atof(split_line[7].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TAPR, atof(split_line[8].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TPMXR, atof(split_line[9].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TPMNR, atof(split_line[10].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TSTPR, atoi(split_line[11].c_str()));
            dc_line_instance.push_back(data);

            getline(input, line);
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

            dc_line_set.push_back(dc_line_instance);
            std::getline(input, line);

            data.setValue(DL_IPI, atof(split_line[0].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_NBI, atof(split_line[1].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_GAMMX, atoi(split_line[2].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_GAMMN, atof(split_line[3].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_RCI, atof(split_line[4].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_XCI, atoi(split_line[5].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_EBASI, atof(split_line[6].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TRI, atof(split_line[7].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TAPI, atoi(split_line[8].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TPMXI, atof(split_line[9].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TPMNI, atof(split_line[10].c_str()));
            dc_line_instance.push_back(data);

            data.setValue(DL_TSTPI, atoi(split_line[11].c_str()));
            dc_line_instance.push_back(data);
        }
        case_data->push_back(dc_line_set);
    }

    /*
#define     SHUNT_BUSNUMBER "SHUNT_BUSNUMBER"
#define     SHUNT_MODSW     "SHUNT_MODSW"
#define     SHUNT_VSWHI     "SHUNT_VSWHI"
#define     SHUNT_VSWLO     "SHUNT_VSWLO"
#define     SHUNT_SWREM     "SHUNT_SWREM"
#define     SHUNT_VDES      "VDES"
#define     SHUNT_BINIT     "SHUNT_BINIT"
#define     SHUNT_N1        "SHUNT_N1"
#define     SHUNT_B1        "SHUNT_B1"
#define     SHUNT_N2        "SHUNT_N2"
#define     SHUNT_B2        "SHUNT_B2"
#define     SHUNT_N3        "SHUNT_N3"
#define     SHUNT_B3        "SHUNT_B3"
#define     SHUNT_N4        "SHUNT_N4"
#define     SHUNT_B4        "SHUNT_B4"
#define     SHUNT_N5        "SHUNT_N5"
#define     SHUNT_B5        "SHUNT_B5"
#define     SHUNT_N6        "SHUNT_N6"
#define     SHUNT_B6        "SHUNT_B6"
#define     SHUNT_N7        "SHUNT_N7"
#define     SHUNT_B7        "SHUNT_B7"
#define     SHUNT_N8        "SHUNT_N8"
#define     SHUNT_B8        "SHUNT_B8"
    */
    void find_shunt(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        shunt_set;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "shunt block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   shunt_instance;
            gridpack::component::DataCollection          data;

            data.setValue(SHUNT_BUSNUMBER, atoi(split_line[0].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_VSWHI, atof(split_line[1].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_VSWLD, atof(split_line[2].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_SWREM, atof(split_line[3].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_VDES, atof(split_line[4].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_BINIT, atof(split_line[5].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N1, atoi(split_line[6].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B1, atof(split_line[7].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N2, atoi(split_line[8].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B2, atof(split_line[9].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N3, atoi(split_line[10].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B3, atof(split_line[11].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N4, atoi(split_line[12].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B4, atof(split_line[13].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N5, atoi(split_line[14].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B5, atof(split_line[15].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N6, atoi(split_line[16].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B6, atof(split_line[17].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N7, atoi(split_line[18].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B7, atof(split_line[19].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_N8, atoi(split_line[20].c_str()));
            shunt_instance.push_back(data);

            data.setValue(SHUNT_B8, atof(split_line[21].c_str()));
            shunt_instance.push_back(data);

            shunt_set.push_back(shunt_instance);
            std::getline(input, line);
        }
        case_data->push_back(shunt_set);
    }

    /*
#define     IMPED_I         "I"
    */
    void find_imped_corr(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        imped_corr_set;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "impedence block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   imped_corr_instance;
            gridpack::component::DataCollection          data;

            data.setValue(IMPED_I, atoi(split_line[0].c_str()));
            imped_corr_instance.push_back(data);

            imped_corr_set.push_back(imped_corr_instance);
            std::getline(input, line);
        }
        case_data->push_back(imped_corr_set);
    }

/*
#define     M_TERM_I        "I"
*/
    void find_multi_terminal(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        multi_terminal;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "multi terminal block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   multi_terminal_instance;
            gridpack::component::DataCollection          data;

            data.setValue(M_TERM_I, atoi(split_line[0].c_str()));
            multi_terminal_instance.push_back(data);

            multi_terminal.push_back(multi_terminal_instance);
            std::getline(input, line);
        }
        case_data->push_back(multi_terminal);
    }


    /*
#define     M_SEC_I         "I"
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

            data.setValue(M_SEC_I, atoi(split_line[0].c_str()));
            multi_section_instance.push_back(data);

            multi_section.push_back(multi_section_instance);
            std::getline(input, line);
        }
        case_data->push_back(multi_section);
    }

    /*
#define     ZONE_I          "I"
#define     ZONE_NAME       "NAME"
    */
    void find_zone(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        zone;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "zone block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   zone_instance;
            gridpack::component::DataCollection          data;

            data.setValue(ZONE_I, atoi(split_line[0].c_str()));
            zone_instance.push_back(data);

            data.setValue(ZONE_NAME, split_line[1].c_str());
            zone_instance.push_back(data);

            zone.push_back(zone_instance);
            std::getline(input, line);
        }
        case_data->push_back(zone);
    }

    /*
#define     I_AREA_I        "I"
    */
    void find_interarea(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        inter_area;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "inter area block " << line << std::endl;
        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   inter_area_instance;
            gridpack::component::DataCollection          data;

            data.setValue(I_AREA_I, atoi(split_line[0].c_str()));
            inter_area_instance.push_back(data);

            inter_area.push_back(inter_area_instance);
            std::getline(input, line);
        }
        case_data->push_back(inter_area);
    }
    /*
#define     OWNER_I         "I"
#define     OWNER_NAME      "NAME"
    */
    void find_owner(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        owner;

        std::string          line;
        std::getline(input, line); //this should be the first line of the block
std::cout << "owner block " << line << std::endl;

        std::getline(input, line);

        while(line[0] != TERM_CHAR) {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   owner_instance;
            gridpack::component::DataCollection          data;

            data.setValue(OWNER_I, atoi(split_line[0].c_str()));
            owner_instance.push_back(data);

            data.setValue(OWNER_NAME, split_line[1].c_str());
            owner_instance.push_back(data);

            owner.push_back(owner_instance);
            std::getline(input, line);
        }
        case_data->push_back(owner);
    }

    /*
    */
    void find_device_data(std::vector<data_set>  * case_data, std::ifstream & input)
    {
        data_set                        device_set;
       if (!input.eof()) {

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
std::cout << "device  block " << line << std::endl;

        while(line[0] != TERM_CHAR && !input.eof())  {
            std::vector<std::string>  split_line;
            boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
            std::vector<gridpack::component::DataCollection>   device_instance;
            gridpack::component::DataCollection          data;
            device_set.push_back(device_instance);
            std::getline(input, line);
        }
        } else {
            std::vector<gridpack::component::DataCollection>   device_instance;
            gridpack::component::DataCollection          data;
            device_set.push_back(device_instance);

        }
        case_data->push_back(device_set);
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
};

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI23PARSER_HPP_ */
