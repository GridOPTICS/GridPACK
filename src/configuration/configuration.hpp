#pragma once
#define USE_MPI
#ifndef _configuration_h
// TODO -- any convention on version numbers? 
// TODO -- any coding conventions? 
#define _configuration_h 201307
//#include <boost/mpi/communicator.hpp>
#include <string>
#include <vector>
#include <memory>
#ifdef USE_MPI
#include <mpi.h>
#endif
namespace gridpack {
namespace utility {

class Configuration {
/**
 * We expect execution of a GridPACK application to be determined in part
 * by a single input configuration file (perhaps extracted from a database
 * or web service). We assume that a configuration file is a binding of a
 * hierarchically structured set of keys to a set of values.  The
 * Configuration module is used to access values and abstract the structure
 * of the configuration file which might be concretely represented in a ways
 * including XML or with a custom parser. Different components within
 * GridPACK will share a common configuration file and we expect then the
 * top-level of the key name space to identify the module that is primarily
 * associated with a configuration parameter. This provides extensibility over
 * time and insulates modules from each other’s internals. 
 *
 * Keys are represented by std::string with the character
 * Configuration::KeySeparator,‘.’ (properly escaped in literals) used as
 * a name specifier. So “foo.bar” is a key associated by convention with
 * module “bar”.
 *
 * Values are typed and may be of the following primitive types: bool, int, double
 * or std::string. We also allow a value to be a std::vector<double>. 
 *
 * We use the term cursor to a common prefix on a set of names. Access to a value
 * may be specified by a complete key or by a cursor and the corresponding suffix 
 * relative to that cursor. Thus a cursor may represent the prefix “foo” and
 * relative to that cursor the string “bar” selects the same values as
 * “foo.bar”.  Aside from  factoring out common information for the client,
 * it allows a common configuration substructure do be used by different modules.
 *
 *
 * Sample Usage:
 *
 * string src = "input.xml";
 * Configuration config;
 * config.enable_logging(&std::cout);
 * if(config.open(src)) {
 *    // select a value by path name, with default
 *    string start = config.get("Configuration.Time.Start", "DefaultStart");
 *    cout << "Start " << start << endl;
 *
 *    // select a subtree by path 
 *    Configuration::Cursor * time = config.getCursor("Configuration.Time");
 *
 *    // select a value without specified default
 *    double step = -1.0;
 *    if(time->get("Step", & step)) 
 *      cout << "Step " << step << endl;
 *
 *    std::vector<double> dv ;
 *    if(config.get("Configuration.Option.DefaultVelocity", &dv)) {
 *       cout << "DefaultVelocity: ";
 *       for(double e : dv) 
 *         cout << e << " ";
 *       cout << endl;
 *    }
 * }
 *
 * Sample use of iterating over all children (which might have same XML element name) 
 *
 *
 * 	p = c->getCursor("Configuration.DynamicSimulation.Faults");
 *	Configuration::ChildCursors children;
 *	p->children(children);
 *	int i = 0;
 *	for(auto & c : children) {
 *		string s = c->get("Branch", "No Branch");
 *		cout << i << " " << s << endl;
 *		i += 1;
 *	}
 *
 */

	class ConfigInternals * pimpl;
	bool initialize_internal(MPI_Comm comm); 
public:
	typedef std::string KeyType;
	static const char KeySep = '.';  // inhereted from boost, could change at some cost
	
   /**
    * Simple Constructor
    */
	Configuration(void);

   /**
    * Simple Destructor
    */
	~Configuration(void);

	/**
	 * Access a common instance, shared by all modules in configuration database,
	 */
	static Configuration * configuration();

	/**
    * enable logging for diagnostics and provenence (default is std::cout)
    */
	void enableLogging(std::ostream * = NULL);

	/**
    * read a configuration file. true == success, false == some kind of failure
    */
#ifdef USE_MPI
   /**
    * Open external configuration file on all ranks on communicator MPI_Comm
    * @param file name of external configuration file
    * @param MPI_Comm MPI communicator being used in calculation
    * @return false if there is an error reading XML file
    */
	bool open(std::string file,MPI_Comm);  // on all ranks...
   /**
    * Deprecated method that initializes configuration on all processes except
    * process 0. Can be used in conjunction with "open" call on process 0.
    * @param MPI_Comm MPI communicator being used in calculation
    * @return false if there is an error reading XML file
    */
	bool initialize(MPI_Comm comm);  // deprecated....
#else
   /**
    * Open external configuration file
    * @param file name of external configuration file
    */
	bool open(std::string file);  // rank 0 only
#endif
	/**
	 * For each supported type, there are two variants. One takes a default value 
	 * that is returned if the key is not present in the configuration file,
	 * the other takes a pointer to an output location and returns a boolean. When the boolean is true
	 * the output location is updated with the value
    * @param KeyType data key in key-value pair
    * @param default_value data value in key-value pair
	 */
	bool get(KeyType, bool default_value);
	bool get(KeyType, bool *);
	int get(KeyType, int default_value);
	bool get(KeyType, int *);
	double get(KeyType, double default_value);
	bool get(KeyType, double *);
	std::string get(KeyType, const std::string & default_value);
	// this wrapper makes a common case look a little cleaner 
	std::string get(KeyType key, const char * default_value) { return get(key, std::string(default_value)); }
	bool get(KeyType, std::string *);
	
	// first implementation only works with 3-vectors (which are also individually labeled X,Y,Z)
	std::vector<double> get(KeyType, const std::vector<double> & default_value);
	bool get(KeyType, std::vector<double>*);

	/**
	 * This class represents a prefix of a set of key names.
	 * Conveniently this implementation allows it to be the same class
	 */
	typedef Configuration Cursor ;

   /**
	 * select a prefix, returns NULL of the prefix has no defined keys.
	 * The XML-based implementation on top of boost property maps might
    * have multiple elements with the same name.
	 * This will simply return a pointer to the first, which might not be sufficient.
    * This function will set the cursor so that it matches the XML block pointed
    * to by the KeyType variable. For example, if the key type is set to
    * "Configuration.Powerflow", then subsequent calls to get will look for
    * variables within the block deliminated by
    * <Configuration><Powerflow>...</Powerflow></Configuration>
    * @param KeyType string representing data block to set cursor
    * return cursor pointing to correct data block in configuration file
	 */
	Cursor * getCursor(KeyType);

	/*
  	 * For the root note, return an empty string
	 *   for a cursor, return the key-path as a string that selects this point in the heirarchy
	 */
	const std::string & path();

	/* iterate over children */
	typedef std::vector<std::shared_ptr<Cursor>> ChildCursors;
	void children(ChildCursors &);
};


} // utitilities
} // gridpack
#endif
