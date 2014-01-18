/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include "configuration.hpp"
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <fstream>
#include <iostream>
#include <streambuf>

using boost::property_tree::ptree;
using std::string;
using std::ostream;

namespace gridpack {
namespace utility {

class ConfigInternals {
public:
	ConfigInternals() : logging(NULL), path("") { } 
	std::ostream * logging ;
	boost::property_tree::ptree pt;
	string path;
	bool initialize(const std::string &);
} ;


Configuration::Configuration(void)
{
	pimpl = new ConfigInternals;
}


Configuration::~Configuration(void)
{
	delete pimpl;
}

void Configuration::enableLogging(std::ostream * out) {
	pimpl->logging = out;
}

static Configuration * config = NULL; 
Configuration * Configuration::configuration() {
	if(config == NULL) {
		config = new Configuration();
	}
	return config;
}


#ifdef USE_MPI
bool Configuration::open(const std::string & file,MPI_Comm comm) {
#else
bool Configuration::open(const std::string & file) {
#endif

#ifdef USE_MPI
	int rank;
	MPI_Comm_rank(comm,&rank);
	if(rank != 0) {
		return initialize_internal(comm);
	}
#endif 
	std::string str;
	std::ifstream input(file.c_str());
	if(!input.bad()) {
		input.seekg(0, std::ios::end);   
		str.reserve((unsigned) input.tellg());
		input.seekg(0, std::ios::beg);

		str.assign((std::istreambuf_iterator<char>(input)),
					std::istreambuf_iterator<char>());
	}
#ifdef USE_MPI
	int n = str.size();
	MPI_Bcast(&n, 1, MPI_INT, rank, comm);
	MPI_Bcast((void*) str.c_str(), n, MPI_CHAR, rank, comm);
#endif
    // Load the XML file into the property tree. If reading fails
    // (cannot open file, parse error), an exception is thrown.
	try {
		pimpl->initialize(str);
	}
	catch(...) {
		if(pimpl->logging != NULL)
		 (*pimpl->logging) << "Error reading XML file " << file << std::endl;
		return false;
	}
	if(!pimpl->logging && pimpl->pt.get<bool>("Configuration/enableLogging",false))
		pimpl->logging = & std::cout;
	return true;
}

bool ConfigInternals::initialize(const std::string & input) {
	std::istringstream ss(input);
    read_xml(ss, pt);
	return true;
}

#ifdef USE_MPI
bool Configuration::initialize(MPI_Comm comm) {
	std::cout << "warning: Configuration::initialize is deprecated" << std::endl;
	return initialize_internal(comm);
}
bool Configuration::initialize_internal(MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(comm,&rank);
	int n ;
	MPI_Bcast(&n, 1, MPI_INT, 0, comm);
	assert(n > 0 && n < (1<<20)); // sanity check
	char * buffer = new char[n+1];
	MPI_Bcast(buffer, n, MPI_CHAR, 0, comm);
	std::string input(buffer,buffer+n);
	try {
		pimpl->initialize(input);
	}
	catch(...) {
		std::cout << "Configuration::initiaze fails for rank " << rank << std::endl;
		return false;
	}
	return true;
}
#endif

template<typename T>
T get0(ptree & pt, Configuration::KeyType key, const T & default_value) {
	return pt.get<T>(key, default_value);
}
template<typename T>
bool get0_bool(ptree & pt, Configuration::KeyType key, T * output) {

	boost::optional<T> temp = pt.get_optional<T>(key);
	if(temp) { 
		*output = *temp;
		return true;
	}
	return false;
}

bool Configuration::get(Configuration::KeyType key, bool default_value) { return get0(pimpl->pt, key, default_value) ; }
bool Configuration::get(Configuration::KeyType key, bool * output) { return get0_bool(pimpl->pt,key, output); }
int Configuration::get(Configuration::KeyType key, int default_value) { return get0(pimpl->pt, key, default_value) ; }
bool Configuration::get(Configuration::KeyType key, int * output) { return get0_bool(pimpl->pt,key, output); }
double Configuration::get(Configuration::KeyType key, double default_value) { return get0(pimpl->pt, key, default_value) ; }
bool Configuration::get(Configuration::KeyType key, double * output) { return get0_bool(pimpl->pt,key, output); }
std::string Configuration::get(Configuration::KeyType key, const std::string & default_value) {
  std::string ret = get0(pimpl->pt, key, default_value) ;

  // remove leading and trailing white space from string
  ret.replace(0,ret.find_first_not_of(" "), "");
  ret.replace(ret.find_last_not_of(" ")+1, std::string::npos,"");
  return ret;
}
bool Configuration::get(Configuration::KeyType key, std::string * output) { return get0_bool(pimpl->pt,key, output); }

std::vector<double> Configuration::get(Configuration::KeyType key, const std::vector<double> & default_value) { 
	CursorPtr c = getCursor(key);
	if(!c) return default_value;
	std::vector<double> v;
	v.push_back(c->get("x",0.0));
	v.push_back(c->get("y",0.0));
	v.push_back(c->get("z",0.0));
	return v;
}
bool Configuration::get(Configuration::KeyType key, std::vector<double> * output) { 
	CursorPtr c = getCursor(key);
	if(!c) return false;
	output->push_back(c->get("x",0.0));
	output->push_back(c->get("y",0.0));
	output->push_back(c->get("z",0.0));
	return true;
}


Configuration::CursorPtr Configuration::getCursor(Configuration::KeyType key) {
	boost::optional<ptree&> cpt = pimpl->pt.get_child_optional(key);
	if(!cpt) return CursorPtr((Configuration*)NULL);
	Configuration * c = new Configuration;
	c->pimpl->logging = pimpl->logging;
	c->pimpl->pt = *cpt;
	return CursorPtr (c);
}

void Configuration::children(ChildCursors & cs) {
	cs.clear();
	BOOST_FOREACH(ptree::value_type & v, pimpl->pt) {
		boost::shared_ptr<Cursor> c(new Configuration);
		c->pimpl->logging = pimpl->logging;
		c->pimpl->pt = v.second;
		cs.push_back(c);
	}
}

void Configuration::children_with_names(ChildElements & cs) {
	cs.clear();
	BOOST_FOREACH(ptree::value_type & v, pimpl->pt) {
//		boost::shared_ptr<Cursor> c(new Configuration);
		CursorPtr c(new Configuration);
		c->pimpl->logging = pimpl->logging;
		c->pimpl->pt = v.second;
		cs.push_back(ChildElement());
		ChildElement & last = cs.back();
		last.cursor = c;
		last.name = v.first;
	}
}

} // namespace utility
} // namespace gridpack
