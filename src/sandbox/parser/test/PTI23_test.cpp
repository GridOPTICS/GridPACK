/**
 * @author William A. Perkins
 * @file   vector_construction_test.cpp
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 * @date   2013-06-17 12:09:38 d3g096
 */

#include <iostream>
#include <fstream>
#include <string>
#include <gridpack/parser/Parser.hpp>
#include <gridpack/parser/PTI23_parser.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#define EPSILON     0.0000001
#define TOLERANCE(x, y, eps) (x - EPSILON > y && x + EPSILON < y)

BOOST_AUTO_TEST_SUITE(Parser)

//____________________________________________________________________________//

void generateCaseData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    std::stringstream   stream;

    stream << index << "1 " << index << ".2 "<< std::endl;
}

void generateBusData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    std::stringstream   stream;

    for (int i = 0; i < 9; i++) {
        // BUS_I           integer
        stream << index << i << "1 ";
        // BUS_NAME        ranged integer
        stream << index << i << "2 ";
        // BUS_BASEKV      float
        stream << index << i << ".3 ";
        // BUS_TYPE        integer
        stream << index << i << "4 ";
        // BUS_SHUNT_GL    float
        stream << index << i << ".05 ";
        // BUS_SHUNT_BL    float
        stream << index << i << ".06 ";
        // BUS_AREA        integer
        stream << index << i << "7 ";
        // BUS_ZONE        integer
        stream << index << i << "8 ";
        // BUS_VOLTAGE_MAG float
        stream << index << i << ".09 ";
        // BUS_VOLTAGE_ANG float
        stream << index << i << ".1 ";
        // BUS_OWNER       integer
        stream << index << i << " " << std::endl;
    }
}


void generateLoadData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{

    for (int i = 0; i < 9; i++) {
        // LOAD_BUSNUMBER   integer
        file << index << i << "1 ";
        // LOAD_ID          integer
        file << index << i << "2 ";
        // LOAD_STATUS      integer
        file << index << i << "3 ";
        // LOAD_AREA        integer
        file << index << i << "4 ";
        // LOAD_ZONE        integer
        file << index << i << "5 ";
        // LOAD_PL          float
        file << index << i << ".06 ";
        // LOAD_QL          float
        file << index << i << ".7 ";
        // LOAD_IP          float
        file << index << i << ".08 ";
        // LOAD_IQ          float
        file << index << i << ".09 ";
        // LOAD_YP          float
        file << index << i << ".1 ";
        // LOAD_YQ          integer
        file << index << i << ".11 ";
        // LOAD_OWNER       integer
        file << index << i << ".12 " << std::endl;
    }
}

void generateGeneratorData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{

    std::stringstream   stream;

    for (int i = 0; i < 9; i++) {
        // GENERATOR_BUSNUMBER      integer
        file << index << i << "1 ";
        // GENERATOR_ID             integer
        file << index << i << "2 ";
        // GENERATOR_PG             float
        file << index << i << ".03 ";
        // GENERATOR_QG             float
        file << index << i << ".04 ";
        // GENERATOR_QMAX           float
        file << index << i << ".05 ";
        // GENERATOR_QMIN           float
        file << index << i << ".06 ";
        // GENERATOR_VS             float
        file << index << i << ".07 ";
        // GENERATOR_IREG           integer
        file << index << i << "8 ";
        // GENERATOR_MBASE          float
        file << index << i << ".09 ";
        // GENERATOR_ZSORCE         float
        file << index << i << ".1 ";
        // GENERATOR_XTRAN          float
        file << index << i << ".11 ";
        // GENERATOR_GTAP           float
        file << index << i << ".12 ";
        // GENERATOR_XT             float
        file << index << i << ".13 ";
        // GENERATOR_RMPCT          float
        file << index << i << ".14 ";
        // GENERATOR_PMAX           float
        file << index << i << ".15 ";
        // GENERATOR_PMIN           float
        file << index << i << ".16 ";
        // GENERATOR_OWNER          integer
        file << index << i << "17 ";
    }
}

void generateBranchData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    // BRANCH_FROMBUS      integer
    file << index << i << "1 ";
    // BRANCH_TOBUS        integer
    file << index << i << "2 ";
    // BRANCH_CKT          character
    file << index << 'c' << " ";
    // BRANCH_R            float
    file << index << i << ".03 ";
    // BRANCH_X            float
    file << index << i << ".04 ";
    // BRANCH_B            float
    file << index << i << ".05 ";
    // BRANCH_RATING_A     float
    file << index << i << ".06 ";
    // BBRANCH_RATING_     float
    file << index << i << ".07 ";
    // BRANCH_RATING_C     float
    file << index << i << ".08 ";
    // BRANCH_SHUNT_ADMTTNC_G1     float
    file << index << i << ".09 ";
    // BRANCH_SHUNT_ADMTTNC_B1     float
    file << index << i << ".1 ";
    // BRANCH_SHUNT_ADMTTNC_G2     float
    file << index << i << ".11 ";
    // BRANCH_SHUNT_ADMTTNC_B2     float
    file << index << i << ".12 ";
    // BRANCH_STATUS       float
    file << index << i << ".13 ";
    // BRANCH_LENGTH       float
    file << index << i << ".14 ";
    // BRANCH_OWNER        integer
    file << index << i << "15 ";
}

void generateTransformerData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    /*
      * type: integer
      * #define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"
      */
    file << index << i << "1 ";
     /*
      * type: integer
      * #define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"
      */
    file << index << i << "2 ";
     /*
      * type: integer
      * #define TRANSFORMER_BUS3 "TRANSFORMER_BUS3"
      */
    file << index << i << "3 ";
     /*
      * type: string
      * #define TRANSFORMER_CKT "TRANSFORMER_CKT"
      */
    file << index << i << "string4 ";
     /*
      * type: integer
      * #define TRANSFORMER_CW "TRANSFORMER_CW"
      */
    file << index << i << "5 ";
     /*
      * type: integer
      * #define TRANSFORMER_CZ "TRANSFORMER_CZ"
      */
    file << index << i << "6 ";

     /*
      * type: integer
      * #define TRANSFORMER_CM "TRANSFORMER_CM"
      */
    file << index << i << "7 ";
     /*
      * type: real float
      * #define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"
      */
    file << index << i << ".08 ";

     /*
      * type: real float
      * #define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"
      */
    file << index << i << ".09 ";

     /*
      * type: integer
      * #define TRANSFORMER_NMETR "TRANSFORMER_NMETR"
      */
    file << index << i << "10 ";

     /*
      * type: string
      * #define TRANSFORMER_NAME "TRANSFORMER_NAME"
      */
    file << index << i << "string11 ";

     /*
      * type: integer
      * #define TRANSFORMER_STATUS "TRANSFORMER_STATUS"
      *
      */
    file << index << i << "12 ";
     /*
      * type: integer
      * #define TRANSFORMER_OWNER "TRANSFORMER_OWNER"
      */
    file << index << i << "13 ";
     /*
      * type: real float
      * #define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"
      */
    file << index << i << ".14 ";
     /*
      * type: real float
      * #define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"
      */
    file << index << i << ".15 ";

     /*
      * type: real float
      * #define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"
      */
    file << index << i << ".16 ";

     /*
      * type: real float
      * #define TRANSFORMER_R2_3 "TRANSFORMER_R2_3"
      */
    file << index << i << ".17 ";

     /*
      * type: real float
      * #define TRANSFORMER_X2_3 "TRANSFORMER_X2_3"
      */
    file << index << i << ".18 ";
     /*
      * type: real float
      * #define TRANSFORMER_SBASE2_3 "TRANSFORMER_SBASE2_3"
      */
    file << index << i << ".19 ";
     /*
      * type: real float
      * #define TRANSFORMER_R3_1 "TRANSFORMER_R3_1"
      */
    file << index << i << ".2 ";
     /*
      * type: real float
      * #define TRANSFORMER_X3_1 "TRANSFORMER_X3_1"
      */
    file << index << i << ".21 ";
     /*
      * type: real float
      * #define TRANSFORMER_SBASE3_1 "TRANSFORMER_SBASE3_1"
      */
    file << index << i << ".22 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMSTAR "TRANSFORMER_VMSTAR"
      */
    file << index << i << ".23 ";
     /*
      * type: real float
      * #define TRANSFORMER_ANSTAR "TRANSFORMER_ANSTAR"
      */
    file << index << i << ".24 ";
     /*
      * type: real float
      * #define TRANSFORMER_WINDV1 "TRANSFORMER_WINDV1"
      */
    file << index << i << ".25 ";
     /*
      * type: real float
      * #define TRANSFORMER_NOMV1 "TRANSFORMER_NOMV1"
      */
    file << index << i << ".26 ";
     /*
      * type: real float
      * #define TRANSFORMER_ANG1 "TRANSFORMER_ANG1"
      */
    file << index << i << ".27 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATA1 "TRANSFORMER_RATA1"
      */
    file << index << i << ".28 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATB1 "TRANSFORMER_RATB1"
      */
    file << index << i << ".29 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATC1 "TRANSFORMER_RATC1"
      */
    file << index << i << ".3 ";
     /*
      * type: integer
      * #define TRANSFORMER_COD1 "TRANSFORMER_COD1"
      */
    file << index << i << "31 ";
     /*
      * type: integer
      * #define TRANSFORMER_CONT1 "TRANSFORMER_CONT1"
      */
    file << index << i << "32 ";
     /*
      * type: real float
      * #define TRANSFORMER_RMA1 "TRANSFORMER_RMA1"
      */
    file << index << i << ".33 ";
     /*
      * type: real float
      * #define TRANSFORMER_RMI1 "TRANSFORMER_RMI1"
      */
    file << index << i << ".34 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMA1 "TRANSFORMER_VMA1"
      */
    file << index << i << ".35 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMI1 "TRANSFORMER_VMI1"
      */
    file << index << i << ".36 ";
     /*
      * type: integer
      * #define TRANSFORMER_NTP1 "TRANSFORMER_NTP1"
      */
    file << index << i << "37 ";
     /*
      * type: integer
      * #define TRANSFORMER_TAB1 "TRANSFORMER_TAB1"
      */
    file << index << i << "38 ";
     /*
      * type: real float
      * #define TRANSFORMER_CR1 "TRANSFORMER_CR1"
      */
    file << index << i << ".39 ";
     /*
      * type: real float
      * #define TRANSFORMER_CX1 "TRANSFORMER_CX1"
      */
    file << index << i << ".4 ";
     /*
      * type: real float
      * #define TRANSFORMER_WINDV2 "TRANSFORMER_WINDV2"
      */
    file << index << i << ".41 ";
     /*
      * type: real float
      * #define TRANSFORMER_NOMV2 "TRANSFORMER_NOMV2"
      */
    file << index << i << ".42 ";
     /*
      * type: real float
      * #define TRANSFORMER_ANG2 "TRANSFORMER_ANG2"
      */
    file << index << i << ".43 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATA2 "TRANSFORMER_RATA2"
      */
    file << index << i << ".44 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATA2 "TRANSFORMER_RATB2"
      */
    file << index << i << ".45 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATC2 "TRANSFORMER_RATC2"
      */
    file << index << i << ".46 ";
     /*
      * type: integer
      * #define TRANSFORMER_COD2 "TRANSFORMER_COD2"
      */
    file << index << i << "47 ";
     /*
      * type: integer
      * #define TRANSFORMER_CONT2 "TRANSFORMER_CONT2"
      */
    file << index << i << "48 ";
     /*
      * type: real float
      * #define TRANSFORMER_RMA2 "TRANSFORMER_RMA2"
      */
    file << index << i << ".49 ";
     /*
      * type: real float
      * #define TRANSFORMER_RMI2 "TRANSFORMER_RMI2"
      */
    file << index << i << ".5 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMA2 "TRANSFORMER_VMA2"
      */
    file << index << i << ".51 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMI2 "TRANSFORMER_VMI2"
      */
    file << index << i << ".52 ";
     /*
      * type: integer
      * #define TRANSFORMER_NTP2 "TRANSFORMER_NTP2"
      */
    file << index << i << "53 ";
     /*
      * type: integer
      * #define TRANSFORMER_TAB2 "TRANSFORMER_TAB2"
      */
    file << index << i << "54 ";
     /*
      * type: real float
      * #define TRANSFORMER_CR2 "TRANSFORMER_CR2"
      */
    file << index << i << ".55 ";
     /*
      * type: real float
      * #define TRANSFORMER_CX2 "TRANSFORMER_CX2"
      */
    file << index << i << ".56 ";
     /*
      * type: real float
      * #define TRANSFORMER_WINDV3 "TRANSFORMER_WINDV3"
      */
    file << index << i << ".57 ";
     /*
      * type: real float
      * #define TRANSFORMER_NOMV3 "TRANSFORMER_NOMV3"
      */
    file << index << i << ".58 ";
     /*
      * type: real float
      * #define TRANSFORMER_ANG3 "TRANSFORMER_ANG3"
      */
    file << index << i << ".59 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATA3 "TRANSFORMER_RATA3"
      */
    file << index << i << ".6 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATB3 "TRANSFORMER_RATB3"
      */
    file << index << i << ".61 ";
     /*
      * type: real float
      * #define TRANSFORMER_RATC3 "TRANSFORMER_RATC3"
      */
    file << index << i << ".62 ";
     /*
      * type: integer
      * #define TRANSFORMER_COD3 "TRANSFORMER_COD3"
      */
    file << index << i << "63 ";
     /*
      * type: integer
      * #define TRANSFORMER_CONT3 "TRANSFORMER_CONT3"
      */
    file << index << i << "64 ";
     /*
      * type: real float
      * #define TRANSFORMER_RMA3 "TRANSFORMER_RMA3"
      */
    file << index << i << ".65 ";
     /*
      * type: real float
      * #define TRANSFORMER_RMI3 "TRANSFORMER_RMI3"
      */
    file << index << i << ".66 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMA3 "TRANSFORMER_VMA3"
      */
    file << index << i << ".67 ";
     /*
      * type: real float
      * #define TRANSFORMER_VMI3 "TRANSFORMER_VMI3"
      */
    file << index << i << ".68 ";
     /*
      * type: integer
      * #define TRANSFORMER_NTP3 "TRANSFORMER_NTP3"
      */
    file << index << i << "69 ";
     /*
      * type: integer
      * #define TRANSFORMER_TAB3 "TRANSFORMER_TAB3"
      */
    file << index << i << "70 ";
     /*
      * type: real float
      * #define TRANSFORMER_CR3 "TRANSFORMER_CR3"
      */
    file << index << i << ".71 ";
     /*
      * type: real float
      * #define TRANSFORMER_CX3 "TRANSFORMER_CX3"
      */
    file << index << i << ".72 ";
}

void generateAreaData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    // AREAINTG_NUMBER             "I"                    integer
    file << index << i << "1 ";
    // AREAINTG_ISW           "ISW"                  integer
    file << index << i << "2 ";
    // AREAINTG_PDES          "PDES"                 float
    file << index << i << ".3 ";
    // AREAINTG_PTOL          "PTOL"                 float
    file << index << i << ".4 ";
    // AREAINTG_NAME         "ARNAM"                string
    file << index << i << "string5 ";
}

void generateShuntData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{

    /*
     * type: integer
     * #define SHUNT_BUSNUMBER "SHUNT_BUSNUMBER"
     */
    file << index << i << "1 ";
    /*
    type: integer
    #define SHUNT_MODSW "SHUNT_MODSW"
     */
    file << index << i << "2 ";
    /*
    type: real float
    #define SHUNT_VSWHI "SHUNT_VSWHI"
     */
    file << index << i << ".3 ";
    /*
        type: real float
        #define SHUNT_VSWLO "SHUNT_VSWLO"
     */
    file << index << i << ".4 ";
    /*
        type: integer
        #define SHUNT_SWREM "SHUNT_SWREM"
     */
    file << index << i << "5 ";
    /*
            type: real float
            #define SHUNT_RMPCT "SHUNT_RMPCT"
     */
    file << index << i << ".6 ";
    /*
        type: string
        #define SHUNT_RMIDNT "SHUNT_RMIDNT"
     */
    file << index << i << "string7 ";
    /*
     * type: real float
     * #define SHUNT_BINIT "SHUNT_BINIT"
     */
    file << index << i << ".8 ";
    /*
        type: integer
        #define SHUNT_N1 "SHUNT_N1"
     */
    file << index << i << "9 ";
    /*
        type: integer
        #define SHUNT_N2 "SHUNT_N2"
     */
    file << index << i << "1 ";
    /*
        type: integer
        #define SHUNT_N3 "SHUNT_N3"
     */
    file << index << i << "11 ";
    /*
        type: integer
        #define SHUNT_N4 "SHUNT_N4"
     */
    file << index << i << "12 ";
    /*
        type: integer
        #define SHUNT_N5 "SHUNT_N5"
     */
    file << index << i << "13 ";
    /*
        type: integer
        #define SHUNT_N6 "SHUNT_N6"
     */
    file << index << i << "14 ";
    /*
        type: integer
        #define SHUNT_N7 "SHUNT_N7"
     */
    file << index << i << "15 ";
    /*
        type: integer
        #define SHUNT_N8 "SHUNT_N8"
     */
    file << index << i << "16 ";
    /*
        type: real float
        #define SHUNT_B1 "SHUNT_B1"
     */
    file << index << i << ".17 ";
    /*
        type: real float
        #define SHUNT_B2 "SHUNT_B2"
     */
    file << index << i << ".18 ";
    /*
        type: real float
        #define SHUNT_B3 "SHUNT_B3"
     */
    file << index << i << ".19 ";
    /*
        type: real float
        #define SHUNT_B4 "SHUNT_B4"
     */
    file << index << i << ".2 ";
    /*
        type: real float
        #define SHUNT_B5 "SHUNT_B5"
     */
    file << index << i << ".21 ";
    /*
        type: real float
        #define SHUNT_B6 "SHUNT_B6"
     */
    file << index << i << ".22 ";
    /*
        type: real float
        #define SHUNT_B7 "SHUNT_B7"
     */
    file << index << i << ".23 ";
    /*
        type: real float
        #define SHUNT_B8 "SHUNT_B8"
     */
    file << index << i << ".24 ";
}

void generateImpedanceCorrData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{

    /*
* type: integer
* #define XFMR_CORR_TABLE_NUMBER "XFMR_CORR_TABLE_NUMBER"
     */
    file << index << i << "1 ";
    /*
* type: real float
* #define XFMR_CORR_TABLE_Ti "XFMR_CORR_TABLE_Ti"
     */
    file << index << i << ".2 ";
    /*
* type: real float
* #define XFMR_CORR_TABLE_Fi "XFMR_CORR_TABLE_Fi"
     */
    file << index << i << ".3 ";
}

void generateMultiSectionData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{

    /*
* type: integer
* #define MULTI_SEC_LINE_FROMBUS "MULTI_SEC_LINE_FROMBUS"

     */
    file << index << i << "1 ";
    /*
* type: integer
* #define MULTI_SEC_LINE_TOBUS "MULTI_SEC_LINE_TOBUS"

     */
    file << index << i << "2 ";
    /*
* type: string
* #define MULTI_SEC_LINE_ID "MULTI_SEC_LINE_ID"

     */
    file << index << i << "string3 ";
    /*
* type: integer
* #define MULTI_SEC_LINE_DUMi "MULTI_SEC_LINE_DUMi"
     */
    file << index << i << "4 ";
}

void generateZoneData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{

    /**
     * Zone Number
     * type: integer
     */
    file << index << i << "1 ";

    /**
     * Zone Name
     * type: string
     */
    file << index << i << "string2 ";
}

void generateInterareaData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    /**
     * From area number of interarea transfer
     * type: integer
    #define INTERAREA_TRANSFER_FROM "INTERAREA_TRANSFER_FROM"
     */
    file << index << i << "1 ";
    /**
     * To area number of interarea transfer
     * type: integer
    #define INTERAREA_TRANSFER_TO "INTERAREA_TRANSFER_TO"
     */
    file << index << i << "2 ";
    /**
     * Single-character (0 through 9 or A through Z) upper case interarea transfer identifier
     * used to distinguish among multiple transfers
     * type: character
    #define INTERAREA_TRANSFER_TRID "INTERAREA_TRANSFER_TRID"
     */
    file << index << "3 ";
    /**
     * MW comprising this transfer
     * type: real float
    #define INTERAREA_TRANSFER_PTRAN "INTERAREA_TRANSFER_PTRAN"
     */
    file << index << i << ".4 ";
}

void generateOwnerData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
    /**
     * Owner number
     * type: integer
    #define OWNER_NUMBER "OWNER_NUMBER"
     */
    file << index << i << "1 ";
    /**
     * Owner name
     * type: integer
    #define OWNER_NAME "OWNER_NAME"
     */
    file << index << i << "2 ";

}

void generateDeviceData(std::ostream & file, std::vector<data_set>  & case_data, int index)
{
}



void generatePTI23File(const char * fileName, )
{
    // open file
    std::vector<data_set>  case_data;
    std::ofstream        file = std::ofstream(fileName);
    int                  index = 0;

    generateCaseData(file, case_data, ++index);
    generateBusData(file, case_data, ++index);
    file << '0' << std::endl;
    generateLoadData(file, case_data, ++index);
    file << '0' << std::endl;
    generateGeneratorData(file, case_data, ++index);
    file << '0' << std::endl;
    generateBranchData(file, case_data, ++index);
    file << '0' << std::endl;
    generateTransformerData(file, case_data, ++index);
    file << '0' << std::endl;
    generateAreaData(file, case_data, ++index);
    file << '0' << std::endl;
    generateShuntData(file, case_data, ++index);
    file << '0' << std::endl;
    generateImpedanceCorrData(file, case_data, ++index);
    file << '0' << std::endl;
    generateMultiSectionData(file, case_data, ++index);
    file << '0' << std::endl;
    generateZoneData(file, case_data, ++index);
    file << '0' << std::endl;
    generateOwnerData(file, case_data, ++index);
    file << '0' << std::endl;
    generateDeviceData(file, case_data, ++index);

    file.close();
}

void checkValue(const char * str, const char * value)
{
    char                * compareValue      = NULL;
    int                   compare           = 0;

    data.getValue(str, compareValue);
    compare = strcmp(str, compareValue);
    std::cout << str << value << compareValue<< std::endl;
    BOOST_CHECK_EQUAL(compare, 0); // values match when compare is zero
}

void checkValue(const char * str, int value)
{
    int                   compareValue      = 0;

    data.getValue(str, compareValue);
    std::cout << str << value << compareValue<< std::endl;
    BOOST_CHECK_EQUAL(value, compareValue);
}

void checkValue(const char * str, double value)
{
    double               compareValue      = 0;

    data.getValue(str, compareValue);
    std::cout << str << value << compareValue<< std::endl;
    BOOST_CHECK_EQUAL(value, compareValue);
}

void checkValue(const char * str, char value)
{
    char               compareValue      = 0;

    data.getValue(str, compareValue);
    std::cout << str << value << compareValue<< std::endl;
    BOOST_CHECK_EQUAL(value, compareValue);
}

void validate_case(std::vector<gridpack::component::DataCollection> & collection)
{
    DataCollection   data    = collection.begin();
    int              iValue      = 0;
    double           dValue      = 0.0;
    char             sValue      = "\0";

    checkValue(CASE_IC,  0);
    checkValue(CASE_SBASE,  100.000);
    checkValue(CASE_RECORD2,  "IEEE 118 bus test case");
    checkValue(CASE_RECORD3,  "IEEE 118 bus test case");
}

void validate_buses(std::vector<gridpack::component::DataCollection> & collection)
{
    int              iValue      = 0;
    double           dValue      = 0.0;
    char             sValue      = "\0";

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        checkValue(BUS_I,  0);
        checkValue(BUS_IDE,  0);
        checkValue(BUS_PL,  0);
        checkValue(BUS_QL,  0);

        checkValue(BUS_GL,  0);
        checkValue(BUS_BL,  0);
        checkValue(BUS_IA,  0);
        checkValue(BUS_VM,  0);

        checkValue(BUS_VA,  0);
        checkValue(BUS_NAME,  0);
        checkValue(BUS_BASKV,  0);
        checkValue(BUS_ZONE,  0);
    }
}


void validate_generator(std::vector<gridpack::component::DataCollection> & collection)
{
    int              iValue      = 0;
    double           dValue      = 0.0;
    char             cValue      = '\0';
    char           * sValue      = "\0";

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        checkValue(GEN_I,  0);
        checkValue(GEN_ID,  0);
        checkValue(GEN_PG,  0);
        checkValue(GEN_QG,  0);

        checkValue(GEN_QT,  0);
        checkValue(GEN_QB,  0);
        checkValue(GEN_VS,  0);
        checkValue(GEN_IREG,  0);

        checkValue(GEN_MBASE,  0);
        checkValue(GEN_ZR,  0);
        checkValue(GEN_ZX,  0);
        checkValue(GEN_RT,  0);

        checkValue(GEN_XT,  0);
        checkValue(GEN_GTAP,  0);
        checkValue(GEN_STAT,  0);
        checkValue(GEN_RMPCT,  0);

        checkValue(GEN_PT,  0);
        checkValue(GEN_PB,  0);
    }
}

void validate_branch(std::vector<gridpack::component::DataCollection> & collection)
{

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        checkValue(BRANCH_I,  0);
        checkValue(BRANCH_J,  0);
        checkValue(BRANCH_CKT,  0);
        checkValue(BRANCH_R,  0);

        checkValue(BRANCH_X,  0);
        checkValue(BRANCH_B,  0);
        checkValue(BRANCH_RATEA,  0);
        checkValue(BRANCH_RATEB,  0);

        checkValue(BRANCH_RATEC,  0);
        checkValue(BRANCH_RATIO,  0);
        checkValue(BRANCH_ANGLE,  0);
        checkValue(BRANCH_GI,  0);

        checkValue(GEN_XT,  0);
        checkValue(GEN_GTAP,  0);
        checkValue(GEN_STAT,  0);
        checkValue(GEN_RMPCT,  0);

        checkValue(BRANCH_BI,  0);
        checkValue(BRANCH_GJ,  0);
        checkValue(BRANCH_BJ,  0);
        checkValue(BRANCH_ST,  0);
    }
}

void validate_tranformer(std::vector<gridpack::component::DataCollection> & collection)
{


    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {
        checkValue(TRANSF_I,  0);
        checkValue(TRANSF_J,  0);
        checkValue(TRANSF_CKT,  0);
        checkValue(TRANSF_VMA,  0);

        checkValue(TRANSF_VMI,  0);
        checkValue(TRANSF_STEP,  0);
        checkValue(TRANSF_TABLE,  0);
    }
}

void validate_area(std::vector<gridpack::component::DataCollection> & data)
{

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        checkValue(AREA_I,  0);
        checkValue(AREA_ISW,  0);
        checkValue(AREA_PDES,  0);
        checkValue(AREA_PTOL,  0);

        checkValue(AREA_ARNAM,  0);
    }
}

void validate_dc_line(std::vector<gridpack::component::DataCollection> & data)
{


    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        checkValue(DL_I,  0);
        checkValue(DL_MDC,  0);
        checkValue(DL_RDC,  0);
        checkValue(DL_SETVL,  0);

        checkValue(DL_VSCHD,  0);
        checkValue(DL_VCMOD,  0);
        checkValue(DL_RCOMP,  0);
        checkValue(DL_DELTI,  0);

        checkValue(DL_METER,  0);
        checkValue(DL_IPR,  0);
        checkValue(DL_NBR,  0);
        checkValue(DL_ALFMAX,  0);

        checkValue(DL_ALFMN,  0);
        checkValue(DL_RCR,  0);
        checkValue(DL_XCR,  0);
        checkValue(DL_EBASR,  0);

        checkValue(DL_TRR,  0);
        checkValue(DL_TAPR,  0);
        checkValue(DL_TPMXR,  0);
        checkValue(DL_TPMNR,  0);

        checkValue(DL_TRR,  0);
        checkValue(DL_TAPR,  0);
        checkValue(DL_TPMXR,  0);
        checkValue(DL_TPMNR,  0);

        checkValue(DL_TSTPR,  0);
        checkValue(DL_IPI,  0);
        checkValue(DL_NBI,  0);
        checkValue(DL_GAMMX,  0);

        checkValue(DL_GAMMN,  0);
        checkValue(DL_RCI,  0);
        checkValue(DL_XCI,  0);
        checkValue(DL_EBASI,  0);


        checkValue(DL_TRI,  0);
        checkValue(DL_TAPI,  0);
        checkValue(DL_TPMXI,  0);
        checkValue(DL_TPMNI,  0);

        checkValue(DL_TSTPI,  0);
    }

}

void validate_shunt(std::vector<gridpack::component::DataCollection> & data)
{

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        data.getValue(SHUNT_I, &iValue);
        std::cout << "SHUNT_I               " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_VSWHI, &dValue);
        std::cout << "SHUNT_VSWHI           " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_VSWLD, &dValue);
        std::cout << "SHUNT_VSWLD           " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_SWREM, &dValue);
        std::cout << "SHUNT_SWREM           " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_VDES, &dValue);
        std::cout << "SHUNT_VDES            " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_BINIT, &dValue);
        std::cout << "SHUNT_BINIT           " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N1, &iValue);
        std::cout << "SHUNT_N1              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B1, &dValue);
        std::cout << "SHUNT_B1              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N2, &iValue);
        std::cout << "SHUNT_N2              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B2, &dValue);
        std::cout << "SHUNT_B2              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N3, &iValue);
        std::cout << "SHUNT_N3              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B3, &dValue);
        std::cout << "SHUNT_B3              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N4, &iValue);
        std::cout << "SHUNT_N4              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B4, &dValue);
        std::cout << "SHUNT_B4              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N5, &iValue);
        std::cout << "SHUNT_N5              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B5, &dValue);
        std::cout << "SHUNT_B5              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N6, &iValue);
        std::cout << "SHUNT_N6              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B6, &dValue);
        std::cout << "SHUNT_B6              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N7, &iValue);
        std::cout << "SHUNT_N7              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B7, &dValue);
        std::cout << "SHUNT_B7              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);

        data.getValue(SHUNT_N8, &iValue);
        std::cout << "SHUNT_N8              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(SHUNT_B8, &dValue);
        std::cout << "SHUNT_B8              " << dValue << std::endl;
        BOOST_CHECK_EQUAL(dValue, i);
    }

}

void validate_imped(std::vector<gridpack::component::DataCollection> & data)
{

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        data.getValue(IMPED_I, &iValue);
        std::cout << "IMPED_I               " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);
    }

    data.setValue(, atoi(split_line[0].c_str()));
}

void validate_terminal(std::vector<gridpack::component::DataCollection> & data)
{

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        data.getValue(M_TERM_I, &iValue);
        std::cout << "M_TERM_I              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);
    }
}

void validate_section(std::vector<gridpack::component::DataCollection> & data)
{
    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        data.getValue(M_SEC_I, &iValue);
        std::cout << "M_SEC_I              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);
    }
}

void validate_zone(std::vector<gridpack::component::DataCollection> & data)
{

    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
            data != collection.end(); ++data, ++i) {

        data.getValue(ZONE_I, &iValue);
        std::cout << "ZONE_I                " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

        data.getValue(ZONE_NAME, &sValue);
        std::cout << "ZONE_NAME             " << sValue << std::endl;
        BOOST_CHECK_EQUAL(sValue, i);
    }

}

void validate_interarea(std::vector<gridpack::component::DataCollection> & data)
{
    for (std::vector<gridpack::component::DataCollection>::iterator data = collection.begin(), int i =0 ;
        data != collection.end(); ++data, ++i) {

        data.getValue(I_AREA_I, &iValue);
        std::cout << "I_AREA_I              " << iValue << std::endl;
        BOOST_CHECK_EQUAL(iValue, i);

   }
}

void validate_owner(std::vector<gridpack::component::DataCollection> & data)
{
}

void validate_device(std::vector<gridpack::component::DataCollection> & data)
{
}


BOOST_AUTO_TEST_CASE(openFailure)
{
    bool                    opened         = true;
    try { 
	std::string         fileName        = "";
        gridpack::parser::Parser<gridpack::parser::PTI23_parser> parser;
	parser.getCaseData(fileName);
    } catch (gridpack::Exception & e) {
        opened     = false;
    }

    BOOST_CHECK_EQUAL(opened, false);
}

BOOST_AUTO_TEST_CASE(openSuccess)
{
    bool                    opened        = true;
    try {
        std::string          fileName      = "../PTI23_seqtest.raw";
        gridpack::parser::Parser<gridpack::parser::PTI23_parser> parser;
        parser.getCaseData(fileName);

    } catch (gridpack::Exception & e) {
        opened     = false;
    }

    BOOST_CHECK_EQUAL(opened, true);
}

BOOST_AUTO_TEST_CASE(readValidFile)
{
    bool                    opened        = false;

    try {
        // loop through each data collection for valid data
        std::string          fileName      = "../PTI23_seqtest.raw";
        gridpack::parser::Parser<gridpack::parser::PTI23_parser> parser;
        data_set      * data = parser.getCaseData(fileName);

        std::vector<gridpack::component::DataCollection>   collection   = (*data).begin();
        validate_case(collection);

        // iterate across data
        collection    = ++(*data);
        validate_buses(collection);

        collection    = ++(*data);
        validate_generator(collection);
        
        collection    = ++(*data);
        validate_branch(collection);
        
        collection    = ++(*data);
        validate_tranformer(collection);
        
        collection    = ++(*data);
        validate_area(collection);
        
        collection    = ++(*data);
        validate_dc_line(collection);
        
        collection    = ++(*data);
        validate_shunt(collection);
        
        collection    = ++(*data);
        validate_imped(collection);
        
        collection    = ++(*data);
        validate_terminal(collection);
        
        collection    = ++(*data);
        validate_section(collection);
        
        collection    = ++(*data);
        validate_zone(collection);
        
        collection    = ++(*data);
        validate_interarea(collection);
        
        collection    = ++(*data);
        validate_owner(collection);
        
        collection    = ++(*data);
        validate_device(collection);
 
    } catch (gridpack::Exception & e) {
        opened     = false;
    }

    BOOST_CHECK_EQUAL(true, true);

}


BOOST_AUTO_TEST_SUITE_END()

// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}

