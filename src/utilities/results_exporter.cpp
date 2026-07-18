//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
// -------------------------------------------------------------
/**
 * @file   results_exporter.cpp
 * @date   2026-03-03
 *
 * @brief
 * Implementation of ResultsExporter methods for writing power flow and
 * contingency analysis results to JSON and CSV formats.
 */

#include <iomanip>
#include "results_exporter.hpp"

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
// escapeJSON
// -------------------------------------------------------------
std::string ResultsExporter::escapeJSON(const std::string& s)
{
  std::string result;
  result.reserve(s.size() + 8);
  for (size_t i = 0; i < s.size(); ++i) {
    char c = s[i];
    switch (c) {
      case '"':  result += "\\\""; break;
      case '\\': result += "\\\\"; break;
      case '\b': result += "\\b";  break;
      case '\f': result += "\\f";  break;
      case '\n': result += "\\n";  break;
      case '\r': result += "\\r";  break;
      case '\t': result += "\\t";  break;
      default:   result += c;      break;
    }
  }
  return result;
}

// -------------------------------------------------------------
// writeConvergenceJSON
// -------------------------------------------------------------
void ResultsExporter::writeConvergenceJSON(std::ostream& out,
                                           const ConvergenceSummary& c,
                                           const std::string& indent)
{
  out << indent << "\"convergence\": {\n";
  out << indent << "  \"converged\": " << (c.converged ? "true" : "false") << ",\n";
  out << indent << "  \"iterations\": " << c.iterations << ",\n";
  out << indent << "  \"final_tolerance\": " << std::scientific << c.finalTolerance << ",\n";
  out << std::fixed;
  out << indent << "  \"final_mismatch\": {\n";
  out << indent << "    \"max_p_bus\": " << c.finalMismatch.maxPBus << ",\n";
  out << indent << "    \"max_p_mismatch_mw\": " << std::setprecision(4) << c.finalMismatch.maxPMismatch << ",\n";
  out << indent << "    \"max_q_bus\": " << c.finalMismatch.maxQBus << ",\n";
  out << indent << "    \"max_q_mismatch_mvar\": " << std::setprecision(4) << c.finalMismatch.maxQMismatch << "\n";
  out << indent << "  }";

  if (!c.perIteration.empty()) {
    out << ",\n";
    out << indent << "  \"per_iteration\": [\n";
    for (size_t i = 0; i < c.perIteration.size(); ++i) {
      const MismatchInfo& m = c.perIteration[i];
      out << indent << "    {\n";
      out << indent << "      \"max_p_bus\": " << m.maxPBus << ",\n";
      out << indent << "      \"max_p_mismatch_mw\": " << std::setprecision(4) << m.maxPMismatch << ",\n";
      out << indent << "      \"max_q_bus\": " << m.maxQBus << ",\n";
      out << indent << "      \"max_q_mismatch_mvar\": " << std::setprecision(4) << m.maxQMismatch << "\n";
      out << indent << "    }";
      if (i + 1 < c.perIteration.size()) {
        out << ",";
      }
      out << "\n";
    }
    out << indent << "  ]\n";
  } else {
    out << "\n";
  }

  out << indent << "}";
}

// -------------------------------------------------------------
// writeBusesJSON
// -------------------------------------------------------------
void ResultsExporter::writeBusesJSON(std::ostream& out,
                                     const std::vector<BusResult>& buses,
                                     const std::string& indent)
{
  out << std::fixed;
  out << indent << "\"buses\": [\n";
  for (size_t i = 0; i < buses.size(); ++i) {
    const BusResult& b = buses[i];
    out << indent << "  {\n";
    out << indent << "    \"bus_id\": " << b.busId << ",\n";
    out << indent << "    \"type\": " << b.type << ",\n";
    out << indent << "    \"area\": " << b.area << ",\n";
    out << indent << "    \"zone\": " << b.zone << ",\n";
    out << indent << "    \"base_kv\": " << std::setprecision(2) << b.baseKV << ",\n";
    out << indent << "    \"voltage_pu\": " << std::setprecision(6) << b.voltage << ",\n";
    out << indent << "    \"angle_deg\": " << std::setprecision(6) << b.angle << ",\n";
    out << indent << "    \"p_injection_mw\": " << std::setprecision(4) << b.pInjection << ",\n";
    out << indent << "    \"q_injection_mvar\": " << std::setprecision(4) << b.qInjection << ",\n";
    out << indent << "    \"p_load_mw\": " << std::setprecision(4) << b.pLoad << ",\n";
    out << indent << "    \"q_load_mvar\": " << std::setprecision(4) << b.qLoad << ",\n";
    out << indent << "    \"p_gen_mw\": " << std::setprecision(4) << b.pGen << ",\n";
    out << indent << "    \"q_gen_mvar\": " << std::setprecision(4) << b.qGen << ",\n";
    out << indent << "    \"shunt_mvar\": " << std::setprecision(4) << b.shuntMvar << "\n";
    out << indent << "  }";
    if (i + 1 < buses.size()) {
      out << ",";
    }
    out << "\n";
  }
  out << indent << "]";
}

// -------------------------------------------------------------
// writeBranchesJSON
// -------------------------------------------------------------
void ResultsExporter::writeBranchesJSON(std::ostream& out,
                                        const std::vector<BranchResult>& branches,
                                        const std::string& indent)
{
  out << std::fixed;
  out << indent << "\"branches\": [\n";
  for (size_t i = 0; i < branches.size(); ++i) {
    const BranchResult& br = branches[i];
    out << indent << "  {\n";
    out << indent << "    \"from_bus\": " << br.fromBus << ",\n";
    out << indent << "    \"to_bus\": " << br.toBus << ",\n";
    out << indent << "    \"circuit_id\": \"" << escapeJSON(br.circuitId) << "\",\n";
    out << indent << "    \"p_from_mw\": " << std::setprecision(4) << br.pFrom << ",\n";
    out << indent << "    \"q_from_mvar\": " << std::setprecision(4) << br.qFrom << ",\n";
    out << indent << "    \"p_to_mw\": " << std::setprecision(4) << br.pTo << ",\n";
    out << indent << "    \"q_to_mvar\": " << std::setprecision(4) << br.qTo << ",\n";
    out << indent << "    \"p_loss_mw\": " << std::setprecision(4) << br.pLoss << ",\n";
    out << indent << "    \"q_loss_mvar\": " << std::setprecision(4) << br.qLoss << ",\n";
    out << indent << "    \"mva_from\": " << std::setprecision(4) << br.mvaFrom << ",\n";
    out << indent << "    \"mva_to\": " << std::setprecision(4) << br.mvaTo << ",\n";
    out << indent << "    \"rate_a_mva\": " << std::setprecision(4) << br.rateA << ",\n";
    out << indent << "    \"rate_selected_mva\": " << std::setprecision(4) << br.rateSelected << ",\n";
    out << indent << "    \"loading_percent\": " << std::setprecision(2) << br.loadingPercent << "\n";
    out << indent << "  }";
    if (i + 1 < branches.size()) {
      out << ",";
    }
    out << "\n";
  }
  out << indent << "]";
}

// -------------------------------------------------------------
// writeGeneratorsJSON
// -------------------------------------------------------------
void ResultsExporter::writeGeneratorsJSON(std::ostream& out,
                                          const std::vector<GeneratorResult>& gens,
                                          const std::string& indent)
{
  out << std::fixed;
  out << indent << "\"generators\": [\n";
  for (size_t i = 0; i < gens.size(); ++i) {
    const GeneratorResult& g = gens[i];
    out << indent << "  {\n";
    out << indent << "    \"bus_id\": " << g.busId << ",\n";
    out << indent << "    \"gen_id\": \"" << escapeJSON(g.genId) << "\",\n";
    out << indent << "    \"p_gen_mw\": " << std::setprecision(4) << g.pGen << ",\n";
    out << indent << "    \"q_gen_mvar\": " << std::setprecision(4) << g.qGen << ",\n";
    out << indent << "    \"q_max_mvar\": " << std::setprecision(4) << g.qMax << ",\n";
    out << indent << "    \"q_min_mvar\": " << std::setprecision(4) << g.qMin << ",\n";
    out << indent << "    \"voltage_setpoint_pu\": " << std::setprecision(6) << g.voltageSetpoint << ",\n";
    out << indent << "    \"status\": " << g.status << "\n";
    out << indent << "  }";
    if (i + 1 < gens.size()) {
      out << ",";
    }
    out << "\n";
  }
  out << indent << "]";
}

// -------------------------------------------------------------
// writeViolationsJSON
// -------------------------------------------------------------
void ResultsExporter::writeViolationsJSON(std::ostream& out,
    const std::vector<BranchViolation>& brViols,
    const std::vector<VoltageViolation>& vViols,
    const std::string& indent)
{
  out << std::fixed;
  out << indent << "\"violations\": {\n";
  out << indent << "  \"branches\": [";
  for (size_t i = 0; i < brViols.size(); ++i) {
    const BranchViolation& v = brViols[i];
    out << (i == 0 ? "\n" : ",\n");
    out << indent << "    {"
        << "\"from_bus\": " << v.fromBus
        << ", \"to_bus\": " << v.toBus
        << ", \"circuit_id\": \"" << escapeJSON(v.circuitId) << "\""
        << ", \"mva\": "             << std::setprecision(4) << v.mva
        << ", \"rate_mva\": "        << std::setprecision(4) << v.rate
        << ", \"loading_percent\": " << std::setprecision(2) << v.loadingPercent
        << ", \"base_mva\": "        << std::setprecision(4) << v.baseMva
        << ", \"delta_mva\": "       << std::setprecision(4) << v.deltaMva
        << ", \"severity\": \"" << escapeJSON(v.severity) << "\""
        << "}";
  }
  if (!brViols.empty()) out << "\n" << indent << "  ";
  out << "],\n";
  out << indent << "  \"voltages\": [";
  for (size_t i = 0; i < vViols.size(); ++i) {
    const VoltageViolation& v = vViols[i];
    out << (i == 0 ? "\n" : ",\n");
    out << indent << "    {"
        << "\"bus_id\": " << v.busId
        << ", \"v_pu\": "         << std::setprecision(6) << v.vPu
        << ", \"limit_low\": "    << std::setprecision(4) << v.limitLow
        << ", \"limit_high\": "   << std::setprecision(4) << v.limitHigh
        << ", \"deviation_pu\": " << std::setprecision(6) << v.deviationPu
        << ", \"severity\": \"" << escapeJSON(v.severity) << "\""
        << "}";
  }
  if (!vViols.empty()) out << "\n" << indent << "  ";
  out << "]\n";
  out << indent << "}";
}

// -------------------------------------------------------------
// writePFJSON
// -------------------------------------------------------------
void ResultsExporter::writePFJSON(std::ofstream& out,
                                  const PowerFlowResults& r)
{
  out << "{\n";
  writeConvergenceJSON(out, r.convergence, "  ");
  out << ",\n";
  writeBusesJSON(out, r.buses, "  ");
  out << ",\n";
  writeBranchesJSON(out, r.branches, "  ");
  out << ",\n";
  writeGeneratorsJSON(out, r.generators, "  ");
  out << "\n";
  out << "}\n";
}

// -------------------------------------------------------------
// writeCAJSON
// -------------------------------------------------------------
void ResultsExporter::writeCAJSON(std::ofstream& out,
                                  const ContingencyAnalysisResults& r)
{
  out << "{\n";
  out << "  \"base_case\": {\n";
  writeConvergenceJSON(out, r.baseCase.convergence, "    ");
  out << ",\n";
  writeBusesJSON(out, r.baseCase.buses, "    ");
  out << ",\n";
  writeBranchesJSON(out, r.baseCase.branches, "    ");
  out << ",\n";
  writeGeneratorsJSON(out, r.baseCase.generators, "    ");
  out << "\n";
  out << "  },\n";

  out << "  \"contingencies\": [\n";
  for (size_t i = 0; i < r.contingencies.size(); ++i) {
    writeContingencyResultJSON(out, r.contingencies[i]);
    if (i + 1 < r.contingencies.size()) {
      out << ",";
    }
    out << "\n";
  }
  out << "  ]\n";
  out << "}\n";
}

// -------------------------------------------------------------
// writePFCSV
// -------------------------------------------------------------
void ResultsExporter::writePFCSV(const std::string& basename,
                                 const PowerFlowResults& r,
                                 const std::string& contingencyName)
{
  bool isCA = !contingencyName.empty();

  // --- Buses CSV ---
  {
    std::string fname = basename + "_buses.csv";
    std::ofstream out;
    if (isCA) {
      out.open(fname.c_str(), std::ios::app);
    } else {
      out.open(fname.c_str());
    }
    out << std::fixed;

    // Write header only if file is at position 0 (new file) or not appending
    if (!isCA || out.tellp() == 0) {
      if (isCA) {
        out << "contingency,";
      }
      out << "bus_id,type,area,zone,base_kv,"
          << "voltage_pu,angle_deg,"
          << "p_injection_mw,q_injection_mvar,"
          << "p_load_mw,q_load_mvar,"
          << "p_gen_mw,q_gen_mvar,shunt_mvar\n";
    }

    for (size_t i = 0; i < r.buses.size(); ++i) {
      const BusResult& b = r.buses[i];
      if (isCA) {
        out << contingencyName << ",";
      }
      out << b.busId << ","
          << b.type << ","
          << b.area << ","
          << b.zone << ","
          << std::setprecision(2) << b.baseKV << ","
          << std::setprecision(6) << b.voltage << ","
          << std::setprecision(6) << b.angle << ","
          << std::setprecision(4) << b.pInjection << ","
          << std::setprecision(4) << b.qInjection << ","
          << std::setprecision(4) << b.pLoad << ","
          << std::setprecision(4) << b.qLoad << ","
          << std::setprecision(4) << b.pGen << ","
          << std::setprecision(4) << b.qGen << ","
          << std::setprecision(4) << b.shuntMvar << "\n";
    }
    out.close();
  }

  // --- Branches CSV ---
  {
    std::string fname = basename + "_branches.csv";
    std::ofstream out;
    if (isCA) {
      out.open(fname.c_str(), std::ios::app);
    } else {
      out.open(fname.c_str());
    }
    out << std::fixed;

    if (!isCA || out.tellp() == 0) {
      if (isCA) {
        out << "contingency,";
      }
      out << "from_bus,to_bus,circuit_id,"
          << "p_from_mw,q_from_mvar,p_to_mw,q_to_mvar,"
          << "p_loss_mw,q_loss_mvar,"
          << "mva_from,mva_to,rate_a_mva,rate_selected_mva,loading_percent\n";
    }

    for (size_t i = 0; i < r.branches.size(); ++i) {
      const BranchResult& br = r.branches[i];
      if (isCA) {
        out << contingencyName << ",";
      }
      out << br.fromBus << ","
          << br.toBus << ","
          << br.circuitId << ","
          << std::setprecision(4) << br.pFrom << ","
          << std::setprecision(4) << br.qFrom << ","
          << std::setprecision(4) << br.pTo << ","
          << std::setprecision(4) << br.qTo << ","
          << std::setprecision(4) << br.pLoss << ","
          << std::setprecision(4) << br.qLoss << ","
          << std::setprecision(4) << br.mvaFrom << ","
          << std::setprecision(4) << br.mvaTo << ","
          << std::setprecision(4) << br.rateA << ","
          << std::setprecision(4) << br.rateSelected << ","
          << std::setprecision(2) << br.loadingPercent << "\n";
    }
    out.close();
  }

  // --- Generators CSV ---
  {
    std::string fname = basename + "_generators.csv";
    std::ofstream out;
    if (isCA) {
      out.open(fname.c_str(), std::ios::app);
    } else {
      out.open(fname.c_str());
    }
    out << std::fixed;

    if (!isCA || out.tellp() == 0) {
      if (isCA) {
        out << "contingency,";
      }
      out << "bus_id,gen_id,"
          << "p_gen_mw,q_gen_mvar,"
          << "q_max_mvar,q_min_mvar,"
          << "voltage_setpoint_pu,status\n";
    }

    for (size_t i = 0; i < r.generators.size(); ++i) {
      const GeneratorResult& g = r.generators[i];
      if (isCA) {
        out << contingencyName << ",";
      }
      out << g.busId << ","
          << g.genId << ","
          << std::setprecision(4) << g.pGen << ","
          << std::setprecision(4) << g.qGen << ","
          << std::setprecision(4) << g.qMax << ","
          << std::setprecision(4) << g.qMin << ","
          << std::setprecision(6) << g.voltageSetpoint << ","
          << g.status << "\n";
    }
    out.close();
  }

  // --- Convergence CSV ---
  {
    std::string fname = basename + "_convergence.csv";
    std::ofstream out;
    if (isCA) {
      out.open(fname.c_str(), std::ios::app);
    } else {
      out.open(fname.c_str());
    }

    if (!isCA || out.tellp() == 0) {
      if (isCA) {
        out << "contingency,";
      }
      out << "converged,iterations,final_tolerance,"
          << "max_p_bus,max_p_mismatch_mw,"
          << "max_q_bus,max_q_mismatch_mvar\n";
    }

    if (isCA) {
      out << contingencyName << ",";
    }
    out << (r.convergence.converged ? "true" : "false") << ","
        << r.convergence.iterations << ","
        << std::scientific << r.convergence.finalTolerance << ","
        << std::fixed
        << r.convergence.finalMismatch.maxPBus << ","
        << std::setprecision(4) << r.convergence.finalMismatch.maxPMismatch << ","
        << r.convergence.finalMismatch.maxQBus << ","
        << std::setprecision(4) << r.convergence.finalMismatch.maxQMismatch << "\n";
    out.close();
  }
}

// -------------------------------------------------------------
// writeCACSV
// -------------------------------------------------------------
void ResultsExporter::writeCACSV(const std::string& basename,
                                 const ContingencyAnalysisResults& r)
{
  // Write base case first
  writePFCSV(basename, r.baseCase, "base_case");

  // Write each contingency
  for (size_t i = 0; i < r.contingencies.size(); ++i) {
    writePFCSV(basename, r.contingencies[i].solution,
               r.contingencies[i].name);
  }

  // Write summary CSV
  {
    std::string fname = basename + "_summary.csv";
    std::ofstream out(fname.c_str());
    out << "contingency,type,converged,iterations,"
        << "has_voltage_violation,has_branch_violation\n";

    for (size_t i = 0; i < r.contingencies.size(); ++i) {
      const ContingencyResult& ct = r.contingencies[i];
      out << ct.name << ","
          << ct.type << ","
          << (ct.solution.convergence.converged ? "true" : "false") << ","
          << ct.solution.convergence.iterations << ","
          << (ct.hasVoltageViolation ? "true" : "false") << ","
          << (ct.hasBranchViolation ? "true" : "false") << "\n";
    }
    out.close();
  }
}

// -------------------------------------------------------------
// writeCAJSONHeader (public)
// -------------------------------------------------------------
void ResultsExporter::writeCAJSONHeader(std::ostream& out,
    const PowerFlowResults& baseCase)
{
  out << "{\n";
  out << "  \"base_case\": {\n";
  writeConvergenceJSON(out, baseCase.convergence, "    ");
  out << ",\n";
  writeBusesJSON(out, baseCase.buses, "    ");
  out << ",\n";
  writeBranchesJSON(out, baseCase.branches, "    ");
  out << ",\n";
  writeGeneratorsJSON(out, baseCase.generators, "    ");
  out << "\n";
  out << "  },\n";
  out << "  \"contingencies\": [\n";
}

// -------------------------------------------------------------
// writeContingencyResultJSON (public)
// -------------------------------------------------------------
void ResultsExporter::writeContingencyResultJSON(std::ostream& out,
    const ContingencyResult& ct)
{
  out << "    {\n";
  out << "      \"name\": \"" << escapeJSON(ct.name) << "\",\n";
  out << "      \"type\": \"" << escapeJSON(ct.type) << "\",\n";
  out << "      \"has_voltage_violation\": "
      << (ct.hasVoltageViolation ? "true" : "false") << ",\n";
  out << "      \"has_branch_violation\": "
      << (ct.hasBranchViolation ? "true" : "false") << ",\n";
  writeViolationsJSON(out, ct.branchViolations, ct.voltageViolations, "      ");
  out << ",\n";
  out << "      \"solution\": {\n";
  writeConvergenceJSON(out, ct.solution.convergence, "        ");
  out << ",\n";
  writeBusesJSON(out, ct.solution.buses, "        ");
  out << ",\n";
  writeBranchesJSON(out, ct.solution.branches, "        ");
  out << ",\n";
  writeGeneratorsJSON(out, ct.solution.generators, "        ");
  out << "\n";
  out << "      }\n";
  out << "    }";
}

// -------------------------------------------------------------
// writeCAJSONFooter (public)
// -------------------------------------------------------------
void ResultsExporter::writeCAJSONFooter(std::ostream& out)
{
  out << "  ]\n";
  out << "}\n";
}

// -------------------------------------------------------------
// writePowerFlow (public)
// -------------------------------------------------------------
void ResultsExporter::writePowerFlow(const std::string& filename,
                                     const PowerFlowResults& results,
                                     Format format)
{
  if (format == JSON) {
    std::string fname = filename + ".json";
    std::ofstream out(fname.c_str());
    writePFJSON(out, results);
    out.close();
  } else {
    writePFCSV(filename, results);
  }
}

// -------------------------------------------------------------
// writeContingencyAnalysis (public)
// -------------------------------------------------------------
void ResultsExporter::writeContingencyAnalysis(
    const std::string& filename,
    const ContingencyAnalysisResults& results,
    Format format)
{
  if (format == JSON) {
    std::string fname = filename + ".json";
    std::ofstream out(fname.c_str());
    writeCAJSON(out, results);
    out.close();
  } else {
    writeCACSV(filename, results);
  }
}

} // namespace utility
} // namespace gridpack
