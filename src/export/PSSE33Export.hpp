/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * base_export.hpp
 *
 *  This is a utility class that is used to implement the rest of the export
 *  functionality
 *
 *  Created on: April 4, 2019
 *      Author: Bruce Palmer
 */

#ifndef PSSE33EXPORT_HPP_
#define PSSE33EXPORT_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "gridpack/export/data_blocks/export_bus33.hpp"
#include "gridpack/export/data_blocks/export_load33.hpp"
#include "gridpack/export/data_blocks/export_fxshnt33.hpp"
#include "gridpack/export/data_blocks/export_gen33.hpp"
#include "gridpack/export/data_blocks/export_line33.hpp"
#include "gridpack/export/data_blocks/export_xform33.hpp"
#include "gridpack/export/data_blocks/export_area33.hpp"
#include "gridpack/export/data_blocks/export_2term33.hpp"
#include "gridpack/export/data_blocks/export_vscline33.hpp"
#include "gridpack/export/data_blocks/export_icorr33.hpp"
#include "gridpack/export/data_blocks/export_mterm33.hpp"
#include "gridpack/export/data_blocks/export_msect33.hpp"
#include "gridpack/export/data_blocks/export_zone33.hpp"
#include "gridpack/export/data_blocks/export_iarea33.hpp"
#include "gridpack/export/data_blocks/export_owner33.hpp"
#include "gridpack/export/data_blocks/export_facts33.hpp"
#include "gridpack/export/data_blocks/export_swshnt33.hpp"

namespace gridpack {
namespace expnet {

template <class _network>
class PSSE33Export
{
  public:

    /**
     * Constructor
     */
    explicit PSSE33Export(boost::shared_ptr<_network> network) :
      p_network(network), p_comm(network->communicator())
    {
    }

    /**
     * Destructor
     */
    virtual ~PSSE33Export()
    {
    }

    /**
     * Write out data in PSS/E v33 format to a file
     * @param filename name of file that contains output
     */
    void  writeFile(std::string filename) {
      int me = p_comm.rank();
      std::ofstream fout;
      if (me == 0) {
        fout.open(filename.c_str());
      }
      // Write out individual data blocks
      ExportBus33<_network> buses(p_network);
      buses.writeBusBlock(fout);
      ExportLoad33<_network> loads(p_network);
      loads.writeLoadBlock(fout);
      ExportFxShnt33<_network> fxshnts(p_network);
      fxshnts.writeFxShntBlock(fout);
      ExportGen33<_network> generators(p_network);
      generators.writeGenBlock(fout);
      ExportLine33<_network> lines(p_network);
      lines.writeLineBlock(fout);
      ExportXform33<_network> xform(p_network);
      xform.writeXformBlock(fout);
      ExportArea33<_network> area(p_network);
      area.writeAreaBlock(fout);
      Export2Term33<_network> term2(p_network);
      term2.write2TermBlock(fout);
      ExportVSCLine33<_network> vscline(p_network);
      vscline.writeVSCLineBlock(fout);
      ExportImpedCorr33<_network> icorr(p_network);
      icorr.writeImpedCorrBlock(fout);
      ExportMultiTerm33<_network> mterm(p_network);
      mterm.writeMultiTermBlock(fout);
      ExportMultiSect33<_network> msect(p_network);
      msect.writeMultiSectBlock(fout);
      ExportZone33<_network> zone(p_network);
      zone.writeZoneBlock(fout);
      ExportInterArea33<_network> iarea(p_network);
      iarea.writeInterAreaBlock(fout);
      ExportOwner33<_network> owner(p_network);
      owner.writeOwnerBlock(fout);
      ExportFACTS33<_network> facts(p_network);
      facts.writeFACTSBlock(fout);
      ExportSwShnt33<_network> swshnt(p_network);
      swshnt.writeSwShntBlock(fout);
      if (me == 0) {
        // Write closing 'Q'
        fout << "Q" << std::endl;
        fout.close();
      }
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* PSSE33EXPORT_HPP_ */
