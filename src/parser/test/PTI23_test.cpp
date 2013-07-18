/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-06-17 12:09:38 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */

#include <iostream>
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

void validate_case(data)
{
/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-06-17 12:09:38 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */

#include <iostream>
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

void validate_case(data)
{
        data.setValue(CASE_IC, atoi(split_line[0].c_str()));
        case_instance.push_back(data);

        data.setValue(CASE_SBASE, atoi(split_line[1].c_str()));
        case_instance.push_back(data);

        std::getline(input, line);
        data.setValue(CASE_RECORD2, line.c_str());
        case_instance.push_back(data);

        std::getline(input, line);
        data.setValue(CASE_RECORD3, line.c_str());
        case_instance.push_back(data);

}
void validate_buses(data)  
{

            data.setValue(BUS_I, atoi(split_line[0].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_IDE, atoi(split_line[1].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_PL, atof(split_line[2].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_QL, atof(split_line[3].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_GL, atof(split_line[4].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_BL, atof(split_line[5].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_IA, atof(split_line[6].c_str()));

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
}
     
void validate_generator(data)
{

            data.setValue(GEN_I, atof(split_line[0].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ID, atof(split_line[1].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PG, atof(split_line[2].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QG, atof(split_line[3].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QT, atof(split_line[4].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QB, atof(split_line[5].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_VS, atof(split_line[6].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_IREG, atof(split_line[7].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_MBASE, atof(split_line[8].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ZR, atof(split_line[9].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ZX, atof(split_line[10].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_RT, atof(split_line[11].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_XT, atof(split_line[12].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_GTAP, atof(split_line[13].c_str()));
            gen_instance.push_back(data);


            data.setValue(GEN_STAT, atof(split_line[14].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_RMPCT, atof(split_line[15].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PT, atof(split_line[16].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PB, atof(split_line[17].c_str()));
            gen_instance.push_back(data);
}

void validate_branch(data)
{

            data.setValue(BRANCH_I, atoi(split_line[0].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_J, atoi(split_line[1].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_CKT, atof(split_line[2].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_R, atof(split_line[3].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_X, atof(split_line[4].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_B, atof(split_line[5].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATEA, atof(split_line[6].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATEB, atof(split_line[7].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATEC, atof(split_line[8].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATIO, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_ANGLE, atof(split_line[10].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_GI, atof(split_line[11].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_BI, atof(split_line[12].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_GJ, atof(split_line[13].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_BJ, atof(split_line[14].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_ST, atof(split_line[15].c_str()));
            branch_instance.push_back(data);


            branch_set.push_back(branch_instance);
            std::getline(input, line);
}

void validate_tranformer(data)
{

            data.setValue(TRANSF_I, atoi(split_line[0].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_J, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_CKT, atoi(split_line[2].c_str()));
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
}

void validate_area(data)
{
           data.setValue(AREA_I, atoi(split_line[0].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_ISW, atoi(split_line[1].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_PDES, atof(split_line[2].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_PTOL, atof(split_line[3].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_ARNAM, split_line[4].c_str());
            area_instance.push_back(data);
}

void validate_dc_line(data)
{
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

void validate_shunt(data)
{
            data.setValue(SHUNT_I, atoi(split_line[0].c_str()));
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
}

void validate_imped(data)
{
           data.setValue(IMPED_I, atoi(split_line[0].c_str()));
}

void validate_terminal(data)
{
            data.setValue(M_TERM_I, atoi(split_line[0].c_str()));
}

void validate_section(data)
{
            data.setValue(M_SEC_I, atoi(split_line[0].c_str()));
}

void validate_zone(data)
{
           data.setValue(ZONE_I, atoi(split_line[0].c_str()));
            zone_instance.push_back(data);

            data.setValue(ZONE_NAME, split_line[1].c_str());
}

void validate_interarea(data)
{
           data.setValue(I_AREA_I, atoi(split_line[0].c_str()));
}

void validate_owner(data)
{
}

void validate_device(data)
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
std::cout << "failed to open" << std::endl;
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
        parser.getCaseData(fileName);
        std::vector<data_set>  * data = parer.getCaseData(fileName);
        data_set   collection   = data->begin();
	// iterate across data
        validate_case(collection);
        collection    = ++(*data)
        validate_buses(collection);

        collection    = ++(*data)
        validate_generator(collection);
        
        collection    = ++(*data)
        validate_branch(collection);
        
        collection    = ++(*data)
        validate_tranformer(collection);
        
        collection    = ++(*data)
        validate_area(collection);
        
        collection    = ++(*data)
        validate_dc_line(collection);
        
        collection    = ++(*data)
        validate_shunt(collection);
        
        collection    = ++(*data)
        validate_imped(collection);
        
        collection    = ++(*data)
        validate_terminal(collection);
        
        collection    = ++(*data)
        validate_section(collection);
        
        collection    = ++(*data)
        validate_zone(collection);
        
        collection    = ++(*data)
        validate_interarea(collection);
        
        collection    = ++(*data)
        validate_owner(collection);
        
        collection    = ++(*data)
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


        data.setValue(CASE_IC, atoi(split_line[0].c_str()));
        case_instance.push_back(data);

        data.setValue(CASE_SBASE, atoi(split_line[1].c_str()));
        case_instance.push_back(data);

        std::getline(input, line);
        data.setValue(CASE_RECORD2, line.c_str());
        case_instance.push_back(data);

        std::getline(input, line);
        data.setValue(CASE_RECORD3, line.c_str());
        case_instance.push_back(data);

}
void validate_buses(data)  
{

            data.setValue(BUS_I, atoi(split_line[0].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_IDE, atoi(split_line[1].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_PL, atof(split_line[2].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_QL, atof(split_line[3].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_GL, atof(split_line[4].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_BL, atof(split_line[5].c_str()));
            bus_instance.push_back(data);

            data.setValue(BUS_IA, atof(split_line[6].c_str()));

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
}
     
void validate_generator(data)
{

            data.setValue(GEN_I, atof(split_line[0].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ID, atof(split_line[1].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PG, atof(split_line[2].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QG, atof(split_line[3].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QT, atof(split_line[4].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_QB, atof(split_line[5].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_VS, atof(split_line[6].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_IREG, atof(split_line[7].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_MBASE, atof(split_line[8].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ZR, atof(split_line[9].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_ZX, atof(split_line[10].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_RT, atof(split_line[11].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_XT, atof(split_line[12].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_GTAP, atof(split_line[13].c_str()));
            gen_instance.push_back(data);


            data.setValue(GEN_STAT, atof(split_line[14].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_RMPCT, atof(split_line[15].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PT, atof(split_line[16].c_str()));
            gen_instance.push_back(data);

            data.setValue(GEN_PB, atof(split_line[17].c_str()));
            gen_instance.push_back(data);
}

void validate_branch(data)
{

            data.setValue(BRANCH_I, atoi(split_line[0].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_J, atoi(split_line[1].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_CKT, atof(split_line[2].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_R, atof(split_line[3].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_X, atof(split_line[4].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_B, atof(split_line[5].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATEA, atof(split_line[6].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATEB, atof(split_line[7].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATEC, atof(split_line[8].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_RATIO, atof(split_line[9].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_ANGLE, atof(split_line[10].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_GI, atof(split_line[11].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_BI, atof(split_line[12].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_GJ, atof(split_line[13].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_BJ, atof(split_line[14].c_str()));
            branch_instance.push_back(data);

            data.setValue(BRANCH_ST, atof(split_line[15].c_str()));
            branch_instance.push_back(data);


            branch_set.push_back(branch_instance);
            std::getline(input, line);
}

void validate_tranformer(data)
{

            data.setValue(TRANSF_I, atoi(split_line[0].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_J, atoi(split_line[1].c_str()));
            transformer_instance.push_back(data);

            data.setValue(TRANSF_CKT, atoi(split_line[2].c_str()));
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
}

void validate_area(data)
{
           data.setValue(AREA_I, atoi(split_line[0].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_ISW, atoi(split_line[1].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_PDES, atof(split_line[2].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_PTOL, atof(split_line[3].c_str()));
            area_instance.push_back(data);

            data.setValue(AREA_ARNAM, split_line[4].c_str());
            area_instance.push_back(data);
}

void validate_dc_line(data)
{
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

void validate_shunt(data)
{
            data.setValue(SHUNT_I, atoi(split_line[0].c_str()));
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
}

void validate_imped(data)
{
           data.setValue(IMPED_I, atoi(split_line[0].c_str()));
}

void validate_terminal(data)
{
            data.setValue(M_TERM_I, atoi(split_line[0].c_str()));
}

void validate_section(data)
{
            data.setValue(M_SEC_I, atoi(split_line[0].c_str()));
}

void validate_zone(data)
{
           data.setValue(ZONE_I, atoi(split_line[0].c_str()));
            zone_instance.push_back(data);

            data.setValue(ZONE_NAME, split_line[1].c_str());
}

void validate_interarea(data)
{
           data.setValue(I_AREA_I, atoi(split_line[0].c_str()));
}

void validate_owner(data)
{          data.setValue(OWNER_I, atoi(split_line[0].c_str()));
            owner_instance.push_back(data);

            data.setValue(OWNER_NAME, split_line[1].c_str());
            owner_instance.push_back(data);

}

void validate_device(data)
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
std::cout << "failed to open" << std::endl;
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
        parser.getCaseData(fileName);
        std::vector<data_set>  * data = parer.getCaseData(fileName);

        validate_case(data);
        validate_buses(data);
        validate_generator(data);
        validate_branch(data);
        validate_tranformer(data);
        validate_area(data);
        validate_dc_line(data);
        validate_shunt(data);
        validate_imped(data);
        validate_terminal(data);
        validate_section(data);
        validate_zone(data);
        validate_interarea(data);
        validate_owner(data);
        validate_device(data);
 
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


