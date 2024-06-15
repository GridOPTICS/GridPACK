/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 12, 2024
 */
#ifndef GENPLB_HPP
#define GENPLB_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class GenPlbParser
{
  public:
    /**
     * Constructor
     */
    explicit GenPlbParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~GenPlbParser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
     * @param data_struct data struct object
     * @param data data collection object
     * @param gen_id index of generator
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      std::string stmp;
      // GENERATOR_MODEL              "MODEL"        string
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GENERATOR_MODEL, data_struct.model, g_id);
      }

      // Playback file                string
      if (!data->getValue(GENERATOR_PLAYBACK_FILE,&stmp,g_id)) {
        data->addValue(GENERATOR_PLAYBACK_FILE,
            data_struct.playback_file, g_id);
      } else {
        data->setValue(GENERATOR_PLAYBACK_FILE,
            data_struct.playback_file, g_id);
      }

      // GENERATOR_DAMPING_COEFFICIENT_0             float
      if (!data->getValue(GENERATOR_PLAYBACK_ISCALE,&rval,g_id)) {
        data->addValue(GENERATOR_PLAYBACK_ISCALE,
            data_struct.playback_iscale, g_id);
      } else {
        data->setValue(GENERATOR_PLAYBACK_ISCALE,
            data_struct.playback_iscale, g_id);
      }
    }

    /**
     * Parser list of strings and store results in data collection object
     * @param split_line list of tokens from .dyr file
     * @param data data collection object
     * @param gen_id index of generator
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      int nstr = split_line.size();
      // GENERATOR_MODEL              "MODEL"                  string
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GENERATOR_MODEL, model.c_str(), g_id);
      }

      // PLAYBACK_FILE                           string
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_PLAYBACK_FILE,&stmp,g_id)) {
          data->addValue(GENERATOR_PLAYBACK_FILE,
              split_line[3].c_str(), g_id);
        } else {
          data->setValue(GENERATOR_PLAYBACK_FILE,
              split_line[3].c_str(), g_id);
        }
      } 

      // PLAYBACK_ISCALE                           float
      if (nstr > 4) {
        if (!data->getValue(GENERATOR_PLAYBACK_ISCALE,&rval,g_id)) {
          data->addValue(GENERATOR_PLAYBACK_ISCALE,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_PLAYBACK_ISCALE,
              atof(split_line[4].c_str()), g_id);
        }
      }
    }

    /**
     * Parse list of strings store results in data_struct object
     * @param split_line list of tokens from .dyr file
     * @param data data struct that stores information from file
     */
    void store(std::vector<std::string> &split_line,_data_struct &data)
    {
      // GENERATOR_BUSNUMBER               "I"                   integer
      int o_idx;
      o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      // Clean up 2 character tag for generator ID
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[2]);
      strcpy(data.gen_id, tag.c_str());

      std::string sval;

      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // GENERATOR_MODEL              "MODEL"                 string
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
      // PLAYBACK_FILE                           string
      if (nstr > 3) {
	strcpy(data.playback_file,split_line[3].c_str());
      } 

      // GENERATOR_PLAYBACK_ISCALE                           float
      if (nstr > 4) {
        data.playback_iscale = atof(split_line[4].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
