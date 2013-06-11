/*
 * Parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: kglass
 */

#ifndef PARSER_HPP_
#define PARSER_HPP_

#include <iostream>
#include <string>
#include <vector>
#include "gridpack/component/data_collection.hpp"



typedef std::vector<std::string::iterator> string_array;
typedef std::vector<std::vector<gridpack::component::DataCollection> >   data_set;

namespace gridpack {
namespace parser {

class Parser {
public:
    Parser(std::string & file_name)
    {
        // open validated file
        try {
            fh = open_valid_file(file_name);
        } catch () {
            // handle file open error
        }
    }
	virtual ~Parser();

protected:
	std::ifstream                   input;

private:
};

} /* namespace parser */
} /* namespace gridpack */
#endif /* PARSER_HPP_ */
