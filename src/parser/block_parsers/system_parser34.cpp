/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * system_parser34.cpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#include "system_parser34.hpp"
#include "gridpack/parser/dictionary.hpp"
#include <cctype>
#include <cstdlib>
#include <sstream>

namespace {

// Trim ASCII whitespace and surrounding quotes in place.
void strip(std::string &s)
{
  while (!s.empty() && (s.back() == ' ' || s.back() == '\t' ||
                        s.back() == '\r' || s.back() == '\n')) {
    s.erase(s.size() - 1);
  }
  size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
  if (i > 0) s.erase(0, i);
  if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') &&
      s.back() == s.front()) {
    s = s.substr(1, s.size() - 2);
  }
}

void upper(std::string &s)
{
  for (size_t i = 0; i < s.size(); ++i) s[i] = std::toupper((unsigned char)s[i]);
}

// Split a System-Wide Data line on commas. Leading record keyword is
// returned in `keyword`; remaining tokens are returned as KEY=VALUE pairs
// in `kv` (KEYs are uppercased, VALUEs are stripped).
void parse_record(const std::string &line, std::string &keyword,
                  std::vector<std::pair<std::string,std::string> > &kv)
{
  keyword.clear();
  kv.clear();
  std::string buf;
  std::vector<std::string> tokens;
  for (size_t i = 0; i < line.size(); ++i) {
    char c = line[i];
    if (c == ',') {
      tokens.push_back(buf);
      buf.clear();
    } else {
      buf.push_back(c);
    }
  }
  tokens.push_back(buf);
  if (tokens.empty()) return;
  std::string first = tokens[0];
  strip(first);
  upper(first);
  keyword = first;
  for (size_t i = 1; i < tokens.size(); ++i) {
    std::string tok = tokens[i];
    size_t eq = tok.find('=');
    if (eq == std::string::npos) continue;
    std::string k = tok.substr(0, eq);
    std::string v = tok.substr(eq + 1);
    strip(k); upper(k);
    strip(v);
    if (!k.empty()) kv.push_back(std::make_pair(k, v));
  }
}

// Locate `key` (uppercased) in `kv`; if found, write its parsed value
// into `*out` and return true. Otherwise return false.
bool kv_get_double(const std::vector<std::pair<std::string,std::string> > &kv,
                   const char *key, double *out)
{
  for (size_t i = 0; i < kv.size(); ++i) {
    if (kv[i].first == key) { *out = std::atof(kv[i].second.c_str()); return true; }
  }
  return false;
}
bool kv_get_int(const std::vector<std::pair<std::string,std::string> > &kv,
                const char *key, int *out)
{
  for (size_t i = 0; i < kv.size(); ++i) {
    if (kv[i].first == key) { *out = std::atoi(kv[i].second.c_str()); return true; }
  }
  return false;
}

}  // namespace

gridpack::parser::SystemParser34::SystemParser34(
    std::map<int,int> *bus_map,
    std::map<std::string,int> *name_map,
    std::map<std::pair<int, int>, int> *branch_map) :
    gridpack::parser::BaseBlockParser(
      bus_map, name_map, branch_map)
{
}


gridpack::parser::SystemParser34::~SystemParser34(void)
{
}

void gridpack::parser::SystemParser34::parse(
    gridpack::stream::InputStream &stream,
    boost::shared_ptr<gridpack::component::DataCollection> network_data)
{
  // PSS/E defaults from the v34 manual.
  double thrshz = 1.0e-4;
  double pqbrak = 0.7;
  int newton_itmxn = 50;
  double newton_toln = 0.1;
  double newton_dvlim = 0.99;
  int solver_flatst = 0;
  int solver_varlim = 99;
  int solver_swshnt = 1;
  int solver_nondiv = 0;
  int solver_actaps = 1;
  int solver_areain = 1;
  int solver_phshft = 1;
  int solver_dctaps = 1;

  std::string line;
  stream.nextLine(line);
  while (test_end(line)) {
    std::string keyword;
    std::vector<std::pair<std::string,std::string> > kv;
    parse_record(line, keyword, kv);
    if (keyword == "GENERAL") {
      kv_get_double(kv, "THRSHZ", &thrshz);
      kv_get_double(kv, "PQBRAK", &pqbrak);
    } else if (keyword == "NEWTON") {
      kv_get_int(kv, "ITMXN", &newton_itmxn);
      kv_get_double(kv, "TOLN", &newton_toln);
      kv_get_double(kv, "DVLIM", &newton_dvlim);
    } else if (keyword == "SOLVER") {
      kv_get_int(kv, "FLATST", &solver_flatst);
      kv_get_int(kv, "VARLIM", &solver_varlim);
      kv_get_int(kv, "SWSHNT", &solver_swshnt);
      kv_get_int(kv, "NONDIV", &solver_nondiv);
      kv_get_int(kv, "ACTAPS", &solver_actaps);
      kv_get_int(kv, "AREAIN", &solver_areain);
      kv_get_int(kv, "PHSHFT", &solver_phshft);
      kv_get_int(kv, "DCTAPS", &solver_dctaps);
    }
    // GAUSS, ADJUST, TYSL, RATING, IMPCOR, OWNERSHIP records are read
    // through silently.
    stream.nextLine(line);
  }

  if (network_data) {
    network_data->addValue(CASE_THRSHZ, thrshz);
    network_data->addValue(CASE_PQBRAK, pqbrak);
    network_data->addValue(CASE_NEWTON_ITMXN, newton_itmxn);
    network_data->addValue(CASE_NEWTON_TOLN, newton_toln);
    network_data->addValue(CASE_NEWTON_DVLIM, newton_dvlim);
    network_data->addValue(CASE_SOLVER_FLATST, solver_flatst);
    network_data->addValue(CASE_SOLVER_VARLIM, solver_varlim);
    network_data->addValue(CASE_SOLVER_SWSHNT, solver_swshnt);
    network_data->addValue(CASE_SOLVER_NONDIV, solver_nondiv);
    network_data->addValue(CASE_SOLVER_ACTAPS, solver_actaps);
    network_data->addValue(CASE_SOLVER_AREAIN, solver_areain);
    network_data->addValue(CASE_SOLVER_PHSHFT, solver_phshft);
    network_data->addValue(CASE_SOLVER_DCTAPS, solver_dctaps);
  }
}

void gridpack::parser::SystemParser34::parse(
    gridpack::stream::InputStream &stream)
{
  boost::shared_ptr<gridpack::component::DataCollection> none;
  parse(stream, none);
}
