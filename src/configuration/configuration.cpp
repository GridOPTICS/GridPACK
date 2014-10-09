/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include "configuration.hpp"
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
// #include <boost/property_tree/json_parser.hpp>
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

static void merge_trees(boost::property_tree::ptree &parent, 
						const boost::property_tree::ptree::path_type &childPath, 
						const boost::property_tree::ptree &child) 
{
	for(ptree::const_iterator it=child.begin();it!=child.end();++it) {
	    ptree::path_type curPath = childPath / ptree::path_type(it->first);
		parent.add_child(curPath,it->second);
//	    merge_trees(parent, curPath, it->second);
  }
}



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

const int TAB_SIZE = 3;
static void indent(std::ostream & out, int indent) {
	while(indent-- > 0)
		out << " " ;
}
static bool dump_xml(boost::property_tree::ptree pt, std::ostream & out, int indent_amount = 0) {
	bool first = true;
 	BOOST_FOREACH(ptree::value_type & v, pt) {
		if(first) {
			first = false;
			if(indent_amount > 0)  out << std::endl;
		}
		indent(out, indent_amount);
		out << "<" << v.first << ">" ;
		bool no_child = dump_xml(v.second, out, indent_amount+TAB_SIZE);
		if(! no_child) {
			indent(out, indent_amount);
		}
		else {
			out << v.second.data();
		}
		out << "</" << v.first << ">" << std::endl;
	}
	return first;
}


#ifdef CONFIGURATION_USE_MPI
bool Configuration::open(const std::string & file,
                         gridpack::parallel::Communicator tcomm) {
#else
bool Configuration::open(const std::string & file) {
#endif

	int rank = 0;
#ifdef CONFIGURATION_USE_MPI
    MPI_Comm comm = static_cast<gridpack::parallel::Communicator>(tcomm);
	MPI_Comm_rank(comm,&rank);
	if(rank != 0) {
		return initialize_internal(comm);
	}
#endif 
	std::string str;
	std::ifstream input(file.c_str());
   int n = -1;
	if(!input.bad()) {
		input.seekg(0, std::ios::end);   
		str.reserve((unsigned) input.tellg());
		input.seekg(0, std::ios::beg);

		str.assign((std::istreambuf_iterator<char>(input)),
					std::istreambuf_iterator<char>());
      n = 0;
   }
#ifdef CONFIGURATION_USE_MPI
   if (n >= 0) {
	  n = str.size();
   }
	MPI_Bcast(&n, 1, MPI_INT, rank, comm);
   if (n > 0) {
	  MPI_Bcast((void*) str.c_str(), n, MPI_CHAR, rank, comm);
   } else {
     std::cout<<"Configure: Unable to open file "<<file<<std::endl;
     return false;
   }
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
	if(!pimpl->logging && pimpl->pt.get<bool>("Configuration.enableLogging",false))
		pimpl->logging = & std::cout;
	if(pimpl->logging != NULL && rank== 0) {
		try {
			dump_xml(pimpl->pt, *pimpl->logging);
		}
		catch(...) {
			 (*pimpl->logging) << "Error writing XML file " << file << std::endl;
		}
    }
	return true;
}

bool ConfigInternals::initialize(const std::string & input) {
	std::istringstream ss(input);
	boost::property_tree::ptree pt0;	
    read_xml(ss, pt0);
	merge_trees(pt,"",pt0);
	return true;
}

#ifdef CONFIGURATION_USE_MPI
bool Configuration::initialize(gridpack::parallel::Communicator tcomm) {
	std::cout << "warning: Configuration::initialize is deprecated" << std::endl;
   MPI_Comm comm = static_cast<gridpack::parallel::Communicator>(tcomm);
	return initialize_internal(comm);
}
bool Configuration::initialize_internal(MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(comm,&rank);
	int n ;
	MPI_Bcast(&n, 1, MPI_INT, 0, comm);
   if (n == 0) return false;
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
bool Configuration::get(Configuration::KeyType key, std::string * output) {
  bool ret = get0_bool(pimpl->pt,key, output);

  // remove leading and trailing white space from string
  output->replace(0,output->find_first_not_of(" "), "");
  output->replace(output->find_last_not_of(" ")+1, std::string::npos,"");
  return ret;
}

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

void Configuration::children(ChildElements & cs) {
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
