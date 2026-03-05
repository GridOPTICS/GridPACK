//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   results_exporter.hpp
 * @date   2026-03-03
 *
 * @brief
 * Data structures and exporter class for writing power flow and
 * contingency analysis results to JSON and CSV formats.
 */

#ifndef _results_exporter_h_
#define _results_exporter_h_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  Result data structures
// -------------------------------------------------------------

struct BusResult {
  int busId;
  int type;           // 1=PQ, 2=PV, 3=Slack
  int area;
  int zone;
  double baseKV;
  double voltage;     // pu
  double angle;       // degrees
  double pInjection;  // MW
  double qInjection;  // MVAr
  double pLoad;       // MW
  double qLoad;       // MVAr
  double pGen;        // MW
  double qGen;        // MVAr
  double shuntMvar;   // MVAr
};

struct BranchResult {
  int fromBus;
  int toBus;
  std::string circuitId;
  double pFrom;           // MW
  double qFrom;           // MVAr
  double pTo;             // MW
  double qTo;             // MVAr
  double pLoss;           // MW
  double qLoss;           // MVAr
  double mvaFrom;         // MVA
  double mvaTo;           // MVA
  double rateA;           // MVA
  double loadingPercent;  // max(|S_from|,|S_to|)/rateA * 100
};

struct GeneratorResult {
  int busId;
  std::string genId;
  double pGen;            // MW
  double qGen;            // MVAr
  double qMax;            // MVAr
  double qMin;            // MVAr
  double voltageSetpoint; // pu
  int status;
};

struct MismatchInfo {
  int maxPBus;
  double maxPMismatch;    // MW
  int maxQBus;
  double maxQMismatch;    // MVAr
};

struct ConvergenceSummary {
  bool converged;
  int iterations;
  double finalTolerance;
  MismatchInfo finalMismatch;
  std::vector<MismatchInfo> perIteration;
};

struct PowerFlowResults {
  ConvergenceSummary convergence;
  std::vector<BusResult> buses;
  std::vector<BranchResult> branches;
  std::vector<GeneratorResult> generators;
};

struct ContingencyResult {
  std::string name;
  std::string type;       // "branch" or "generator"
  PowerFlowResults solution;
  bool hasVoltageViolation;
  bool hasBranchViolation;
};

struct ContingencyAnalysisResults {
  PowerFlowResults baseCase;
  std::vector<ContingencyResult> contingencies;
};

// -------------------------------------------------------------
//  class ResultsExporter
// -------------------------------------------------------------
class ResultsExporter {
public:

  enum Format { JSON, CSV };

  /**
   * Write power flow results to file(s).
   * @param filename base filename (without extension for CSV, with for JSON)
   * @param results power flow results to write
   * @param format output format (JSON or CSV)
   */
  static void writePowerFlow(const std::string& filename,
                             const PowerFlowResults& results,
                             Format format);

  /**
   * Write contingency analysis results to file(s).
   * @param filename base filename (without extension for CSV, with for JSON)
   * @param results contingency analysis results to write
   * @param format output format (JSON or CSV)
   */
  static void writeContingencyAnalysis(const std::string& filename,
                                       const ContingencyAnalysisResults& results,
                                       Format format);

  /**
   * Write the CA JSON header including the base case results and
   * open the contingencies array.
   * @param out output stream
   * @param baseCase base case power flow results
   */
  static void writeCAJSONHeader(std::ostream& out,
                                const PowerFlowResults& baseCase);

  /**
   * Write a single contingency result as a JSON object to a stream.
   * Used by the CA driver for parallel gathering.
   * @param out output stream
   * @param ct contingency result to write
   */
  static void writeContingencyResultJSON(std::ostream& out,
                                         const ContingencyResult& ct);

  /**
   * Write the CA JSON footer (close the contingencies array and root object).
   * @param out output stream
   */
  static void writeCAJSONFooter(std::ostream& out);

  /**
   * Write power flow results for a single case to CSV files.
   * When contingencyName is non-empty, files are opened in append mode
   * and a contingency column is prepended.
   * @param basename base filename (without extension)
   * @param r power flow results
   * @param contingencyName name of contingency (empty for standalone PF)
   */
  static void writePFCSV(const std::string& basename,
                         const PowerFlowResults& r,
                         const std::string& contingencyName = "");

private:

  static std::string escapeJSON(const std::string& s);

  static void writeConvergenceJSON(std::ostream& out,
                                   const ConvergenceSummary& c,
                                   const std::string& indent);

  static void writeBusesJSON(std::ostream& out,
                             const std::vector<BusResult>& buses,
                             const std::string& indent);

  static void writeBranchesJSON(std::ostream& out,
                                const std::vector<BranchResult>& branches,
                                const std::string& indent);

  static void writeGeneratorsJSON(std::ostream& out,
                                  const std::vector<GeneratorResult>& gens,
                                  const std::string& indent);

  static void writePFJSON(std::ofstream& out, const PowerFlowResults& r);

  static void writeCAJSON(std::ofstream& out,
                          const ContingencyAnalysisResults& r);

  static void writeCACSV(const std::string& basename,
                         const ContingencyAnalysisResults& r);
};

} // namespace utility
} // namespace gridpack

#endif // _results_exporter_h_
