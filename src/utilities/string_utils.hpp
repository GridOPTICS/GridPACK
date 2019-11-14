//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   string_utils.hpp
 * @author Bruce Palmer
 * @date   2019-11-14 07:01:04 d3g096
 * 
 * @brief  
 * A small utility class that contains various methods for manipulating STL
 * strings in GridPACK applications
 */

#ifndef _string_utils_h_
#define _string_utils_h_

#include <boost/algorithm/string.hpp>
#include <vector>

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class StringUtils
// -------------------------------------------------------------
class StringUtils {
public:

  /**
   * Simple constructor
   */
  StringUtils() {};

  /**
   * Simple destructor
   */
  ~StringUtils() {};

  /**
   * Remove white space from both ends of string
   * @param str string that needs white space removed
   */
  void trim(std::string &str)
  {
    // replace tabs and/or carriage returns with blank space
    int ntok;
    ntok = str.find('\n',0);
    while (ntok != std::string::npos) {
      str[ntok] = ' ';
      ntok++;
      ntok = str.find('\n',ntok);
    }
    ntok = str.find('\t',0);
    while (ntok != std::string::npos) {
      str[ntok] = ' ';
      ntok++;
      ntok = str.find('\t',ntok);
    }
    boost::trim(str);
  };

  /**
   * Convert all characters to upper case
   * @param str characters are converted for string str
   */
  void toUpper(std::string &str)
  {
    boost::to_upper(str);
  }

  /**
   * Convert all characters to lower case
   * @param str characters are converted for string str
   */
  void toLower(std::string &str)
  {
    boost::to_lower(str);
  }

  /**
   * Clean up 2 character tags so that single quotes are removed and single
   * character tags are right-justified. These tags can be delimited by a
   * pair of single quotes, a pair of double quotes, or no quotes
   * @param string original string before reformatting
   * @return 2 character string that is right-justified
   */
  std::string clean2Char(std::string &string)
  {
    std::string tag = string;
    // Find and remove single or double quotes
    int ntok1 = tag.find('\'',0);
    bool sngl_qt = true;
    bool no_qt = false;
    // if no single quote found, then assume double quote or no quote
    if (ntok1 == std::string::npos) {
      ntok1 = tag.find('\"',0);
      // if no double quote found then assume no quote
      if (ntok1 == std::string::npos) {
        ntok1 = tag.find_first_not_of(' ',0);
        no_qt = true;
      } else {
        sngl_qt = false;
      }
    }
    int ntok2;
    if (sngl_qt) {
      ntok1 = tag.find_first_not_of('\'',ntok1);
      ntok2 = tag.find('\'',ntok1);
    } else if (no_qt) {
      ntok2 = tag.find(' ',ntok1);
    } else {
      ntok1 = tag.find_first_not_of('\"',ntok1);
      ntok2 = tag.find('\"',ntok1);
    }
    if (ntok2 == std::string::npos) ntok2 = tag.length();
    std::string clean_tag = tag.substr(ntok1,ntok2-ntok1);
    //get rid of white space
    ntok1 = clean_tag.find_first_not_of(' ',0);
    ntok2 = clean_tag.find(' ',ntok1);
    if (ntok2 == std::string::npos) ntok2 = clean_tag.length();
    tag = clean_tag.substr(ntok1,ntok2-ntok1);
    if (tag.length() == 1) {
      clean_tag = tag;
      clean_tag.append(" ");
    } else {
      clean_tag = tag;
    }
    return clean_tag;
  }

  /**
   * Trim any quotes from the beginning and end of a string
   * @param str string with quotes at beginning and end
   * @return string with quotes and leading and trailing white space removed
   */
  std::string trimQuotes(std::string &str)
  {
    std::string tag = str;
    // Find and remove single or double quotes
    int ntok1 = tag.find('\'',0);
    bool sngl_qt = true;
    bool no_qt = false;
    // if no single quote found, then assume double quote or no quote
    if (ntok1 == std::string::npos) {
      ntok1 = tag.find('\"',0);
      // if no double quote found then assume no quote
      if (ntok1 == std::string::npos) {
        ntok1 = tag.find_first_not_of(' ',0);
        no_qt = true;
      } else {
        sngl_qt = false;
      }
    }
    int ntok2;
    if (sngl_qt) {
      ntok1 = tag.find_first_not_of('\'',ntok1);
      ntok2 = tag.find('\'',ntok1);
    } else if (no_qt) {
      ntok2 = tag.find_last_not_of(' ',ntok1)+1;
    } else {
      ntok1 = tag.find_first_not_of('\"',ntok1);
      ntok2 = tag.find('\"',ntok1);
    }
    if (ntok2 == std::string::npos) ntok2 = tag.length();
    std::string clean_tag = tag.substr(ntok1,ntok2-ntok1);
    //get rid of white space
    trim(clean_tag);
    return clean_tag;
  }

  /**
   * Tokenize a string on blanks and return a vector of strings
   * Blanks within single or double quotes are ignored
   * @param str string of tokens separated by blank characters
   * @return vector containing individual tokens
   */
  std::vector<std::string> blankTokenizer(std::string &str)
  {
    // Replace any tabs with blank spaces
    std::string strcpy = str;
    int slen = strcpy.length();
    int i;
    for (i=0; i<slen; i++) {
      if (strcpy[i] == '\t') strcpy[i] = ' ';
    }
    int ntok1, ntok2;
    std::vector<std::string> ret;
    ntok1 = strcpy.find_first_not_of(' ',0);
    if (strcpy[ntok1] == '\'') {
      ntok2 = strcpy.find('\'',ntok1+1);
      ntok2++;
    } else if (strcpy[ntok1] == '\"') {
      ntok2 = strcpy.find('\"',ntok1+1);
      ntok2++;
    } else if (ntok1 != std::string::npos) {
      ntok2 = strcpy.find(' ',ntok1);
      if (ntok2 == std::string::npos) ntok2 = slen;
    } else {
      return ret;
    }
    ret.push_back(strcpy.substr(ntok1,ntok2-ntok1));
    while (ntok2 < slen-1 && ntok1 != std::string::npos) {
      ntok1 = strcpy.find_first_not_of(' ',ntok2);
      if (strcpy[ntok1] == '\'') {
        ntok2 = strcpy.find('\'',ntok1+1);
        if (ntok2 != std::string::npos) {
          ntok2++;
        } else {
          ntok2 = slen;
        }
      } else if (strcpy[ntok1] == '\"') {
        ntok2 = strcpy.find('\"',ntok1+1);
        if (ntok2 != std::string::npos) {
          ntok2++;
        } else {
          ntok2 = slen;
        }
      } else if (ntok1 != std::string::npos) {
        ntok2 = strcpy.find(' ',ntok1);
        if (ntok2 == std::string::npos) ntok2 = slen;
      } 
      if (ntok2 != std::string::npos) {
        ret.push_back(strcpy.substr(ntok1,ntok2-ntok1));
      }
    }
    return ret;
  }

  /**
   * Tokenize a string on a user specified character string and return a
   * vector of strings. This function does not pay attention to quotes
   * @param str input string
   * @param sep separator string
   * @return vector containing individual tokens
   */
  std::vector<std::string> charTokenizer(std::string &str, const char *sep)
  {
    int slen = str.length();
    int ntok1, ntok2;
    std::vector<std::string> ret;
    ntok1 = str.find_first_not_of(sep,0);
    ntok2 = str.find(sep,ntok1);
    if (ntok2 == std::string::npos) ntok2 = slen;
    if (ntok1<=ntok2) ret.push_back(str.substr(ntok1,ntok2-ntok1));
    while (ntok2 < slen-1 && ntok1 != std::string::npos) {
      ntok1 = str.find_first_not_of(sep,ntok2);
      ntok2 = str.find(sep,ntok1);
      if (ntok2 == std::string::npos) ntok2 = slen;
      if (ntok1<=ntok2) ret.push_back(str.substr(ntok1,ntok2-ntok1));
    }
    return ret;
  }

  /**
   * Convert a string into a bool value. This function will take a variety of
   * strings can generate the corresponding bool. Acceptable values for "true"
   * are "true", "yes", "t", "y", "1". Acceptable values for "false" are
   * "false", "no", "f", "n", "0".
   */
  bool getBool(const char* str) {
    std::string strcpy = str;
    trim(strcpy);
    toLower(strcpy);
    if (strcpy == "true") {
      return true;
    } else if (strcpy == "yes") {
      return true;
    } else if (strcpy == "t") {
      return true;
    } else if (strcpy == "y") {
      return true;
    } else if (strcpy == "1") {
      return true;
    } else if (strcpy == "false") {
      return false;
    } else if (strcpy == "no") {
      return false;
    } else if (strcpy == "f") {
      return false;
    } else if (strcpy == "n") {
      return false;
    } else if (strcpy == "0") {
      return false;
    }
    return false;
  }
  bool getBool(std::string str) {
    return getBool(str.c_str());
  }

};

} // namespace utility
} // namespace gridpack
#endif
