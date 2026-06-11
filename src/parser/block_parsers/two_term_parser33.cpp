/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * two_term_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "two_term_parser33.hpp"
#include "gridpack/parser/dictionary.hpp"
#include <cstdlib>
#include <cmath>

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to internal indices
 */
gridpack::parser::TwoTermParser33::TwoTermParser33(
    std::map<int,int> *bus_map,
    std::map<std::string,int> *name_map,
    std::map<std::pair<int, int>, int> *branch_map) :
    gridpack::parser::BaseBlockParser(
      bus_map, name_map, branch_map)
{
}


/**
 * Simple Destructor
 */
gridpack::parser::TwoTermParser33::~TwoTermParser33(void)
{
}

/**
 * parse two term block. Currently does not store data
 * @param stream input stream that feeds lines from RAW file
 */
void gridpack::parser::TwoTermParser33::parse(
    gridpack::stream::InputStream &stream)
{
  // Walk the block correctly: 3 lines per record (header + rect + inv).
  // We don't store anything; this overload exists for callers that don't
  // care about DC modeling.
  std::string line;
  stream.nextLine(line);
  while (test_end(line)) {
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    // Skip rectifier and inverter lines
    stream.nextLine(line);  // rectifier
    stream.nextLine(line);  // inverter
    stream.nextLine(line);  // next header (or end-of-block)
  }
}

void gridpack::parser::TwoTermParser33::parse(
    gridpack::stream::InputStream &stream,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &busData,
    double /*case_sbase*/)
{
  // Optional opt-out: keep the legacy "drop DC" behavior if the user wants it.
  if (std::getenv("GRIDPACK_DISABLE_DC_INJECTION")) {
    parse(stream);
    return;
  }

  bool debug = (std::getenv("GRIDPACK_DEBUG_DC") != NULL);
  // Reactive-load fraction at each terminal. Crude PSS/E-style estimate:
  // converter VAR consumption ~0.5 * P_DC. Override with GRIDPACK_DC_QFRAC.
  double q_frac = 0.5;
  if (const char *e = std::getenv("GRIDPACK_DC_QFRAC")) {
    double v = atof(e);
    if (v >= 0.0) q_frac = v;
  }

  std::string line;
  stream.nextLine(line);  // first line of the block

  int n_records = 0, n_skipped = 0, n_active = 0;
  double total_mw = 0.0;

  while (test_end(line)) {
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    // Header line: 'NAME', MDC, RDC, SETVL, VSCHD, ...
    std::string header_line = line;
    this->cleanComment(header_line);
    std::vector<std::string> hdr = this->splitPSSELine(header_line);

    // Rectifier line:  IPR, NBR, ANMXR, ANMNR, RCR, XCR, EBASR, ...
    stream.nextLine(line);
    if (!test_end(line)) break;
    std::string rect_line = line;
    this->cleanComment(rect_line);
    std::vector<std::string> rect = this->splitPSSELine(rect_line);

    // Inverter line: IPI, NBI, ANMXI, ANMNI, RCI, XCI, EBASI, ...
    stream.nextLine(line);
    if (!test_end(line)) break;
    std::string inv_line = line;
    this->cleanComment(inv_line);
    std::vector<std::string> inv = this->splitPSSELine(inv_line);

    n_records++;

    // Parse fields with bounds checks
    int mdc = (hdr.size() > 1) ? atoi(hdr[1].c_str()) : 0;
    double setvl = (hdr.size() > 3) ? atof(hdr[3].c_str()) : 0.0;
    double vschd = (hdr.size() > 4) ? atof(hdr[4].c_str()) : 0.0;
    int ipr = (rect.size() > 0) ? atoi(rect[0].c_str()) : 0;
    int ipi = (inv.size()  > 0) ? atoi(inv[0].c_str())  : 0;

    if (mdc == 0 || setvl == 0.0 || ipr == 0 || ipi == 0) {
      n_skipped++;
      if (debug) {
        printf("DC_INJECT skip name=%s mdc=%d setvl=%g ipr=%d ipi=%d\n",
               hdr.size()>0 ? hdr[0].c_str() : "?", mdc, setvl, ipr, ipi);
      }
      stream.nextLine(line);
      continue;
    }

    // Convert SETVL to MW. MDC=1 -> SETVL is MW. MDC=2 -> SETVL is kA at
    // rectifier DC bus; approximate P = SETVL * VSCHD (MW = kA * kV).
    double pdc_mw = (mdc == 2) ? std::abs(setvl) * std::abs(vschd) : std::abs(setvl);
    double qdc_mvar = q_frac * pdc_mw;

    // Look up internal bus indices
    std::map<int,int>::iterator it_r = p_busMap->find(ipr);
    std::map<int,int>::iterator it_i = p_busMap->find(ipi);
    if (it_r == p_busMap->end() || it_i == p_busMap->end()) {
      n_skipped++;
      if (debug) {
        printf("DC_INJECT skip-unmapped name=%s ipr=%d ipi=%d\n",
               hdr.size()>0 ? hdr[0].c_str() : "?", ipr, ipi);
      }
      stream.nextLine(line);
      continue;
    }
    int idx_r = it_r->second;
    int idx_i = it_i->second;

    // Append a synthetic LOAD record at each terminal.
    // Rectifier consumes P from AC: positive load.
    // Inverter supplies P to AC: negative load (i.e. injection).
    // Both terminals consume reactive (positive Q load).
    auto add_dc_load = [&](int idx, int bus_num, const std::string &id_tag,
                           double pl_mw, double ql_mvar) {
      gridpack::component::DataCollection *bd = busData[idx].get();
      int nld = 0;
      if (!bd->getValue(LOAD_NUMBER, &nld)) nld = 0;
      bd->addValue(LOAD_BUSNUMBER, bus_num, nld);
      bd->addValue(LOAD_ID, id_tag.c_str(), nld);
      bd->addValue(LOAD_STATUS, 1, nld);
      // PL/QL are in MW/MVAr (LoadParser convention)
      if (nld == 0) {
        bd->addValue(LOAD_PL, pl_mw);
        bd->addValue(LOAD_QL, ql_mvar);
      }
      bd->addValue(LOAD_PL, pl_mw, nld);
      bd->addValue(LOAD_QL, ql_mvar, nld);
      bd->addValue(LOAD_IP, 0.0, nld);
      bd->addValue(LOAD_IQ, 0.0, nld);
      bd->addValue(LOAD_YP, 0.0, nld);
      bd->addValue(LOAD_YQ, 0.0, nld);
      nld++;
      if (!bd->setValue(LOAD_NUMBER, nld)) {
        bd->addValue(LOAD_NUMBER, nld);
      }
    };

    add_dc_load(idx_r, ipr, "DR",  pdc_mw, qdc_mvar);  // rectifier: +P load
    add_dc_load(idx_i, ipi, "DI", -pdc_mw, qdc_mvar);  // inverter:  -P load

    n_active++;
    total_mw += pdc_mw;
    if (debug) {
      printf("DC_INJECT add name=%s rect_bus=%d inv_bus=%d P=%.2f MW Q=%.2f MVAr\n",
             hdr.size()>0 ? hdr[0].c_str() : "?", ipr, ipi, pdc_mw, qdc_mvar);
    }

    stream.nextLine(line);
  }

  if (n_records > 0) {
    printf("Two-terminal DC: %d records, %d active (P-injections synthesized), "
           "%d skipped, total |P| = %.1f MW\n",
           n_records, n_active, n_skipped, total_mw);
  }
}
