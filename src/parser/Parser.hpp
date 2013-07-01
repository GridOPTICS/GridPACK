/*
 * Parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: kglass
 */

#ifndef PARSER_HPP_
#define PARSER_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gridpack/component/data_collection.hpp>



typedef std::vector<std::string::iterator> string_array;
typedef std::vector<std::vector<gridpack::component::DataCollection> >   data_set;


namespace gridpack {
namespace parser {

template <typename PARSER_TYPE>
class Parser {
public:
    Parser<PARSER_TYPE>(void){};
    std::vector<data_set>  * getCaseData(std::string & fileName)
        {
        // open valid file
        std::ifstream            input;
        std::vector<data_set>  * case_data;
        try {
            input.open(fileName.c_str());
        } catch (std::ios_base::failure & e){
            // let the calling function determine the response
            throw;
        }
        try {
            case_data = parser.getCase(input);
        } catch (std::ios_base::failure & e) {
            // let the calling function determine the response
            throw;
        }
        return case_data;
    }
    virtual ~Parser(){};

protected:
    PARSER_TYPE              parser;
private:

};

} /* namespace parser */
} /* namespace gridpack */
#endif /* PARSER_HPP_ */
