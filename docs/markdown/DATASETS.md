## Data Sets

The GridPACK repository contains a large number of data sets that can be used
for testing of GridPACK. These are used in the input files that are typically
included in each application directory when GridPACK is built. A brief
description of these data sets is included here.

## PSS/E raw files

Most of these files are in the version 23 format, but some
files use the version 33 format. Typically, files that don't specify a version
are version 23. All of these files may be used and redistributed freely. Cases
that originate outside GridPACK carry attribution requirements; see
[Third-party data set attributions](#third-party-data-set-attributions) at the
end of this page.

### Simple 9 bus network
* [case9.raw](../../src/applications/data_sets/raw/case9.raw)
* [9b3g.raw](../../src/applications/data_sets/raw/9b3g.raw)
* [IEEE3G9B_V23_bus5smallload.raw](../../src/applications/data_sets/raw/IEEE3G9B_V23_bus5smallload.raw)

### IEEE 14 bus test case
* [IEEE14_ca_mod_rate6.raw](../../src/applications/data_sets/raw/IEEE14_ca_mod_rate6.raw)
* [IEEE14.raw](../../src/applications/data_sets/raw/IEEE14.raw)
* [IEEE14_PTIv33.raw](../../src/applications/data_sets/raw/IEEE14_PTIv33.raw)
* [IEEE14_ca.raw](../../src/applications/data_sets/raw/IEEE14_ca.raw)
* [IEEE14_double.raw](../../src/applications/data_sets/raw/IEEE14_double.raw)
* [IEEE14_kds.raw](../../src/applications/data_sets/raw/IEEE14_kds.raw)

### IEEE 118 bus test case
* [IEEE118.raw](../../src/applications/data_sets/raw/IEEE118.raw)
* [118_PTIv33.raw](../../src/applications/data_sets/raw/118_PTIv33.raw)

### IEEE 145 bus test case
* [IEEE_145bus_v23_PSLF.raw](../../src/applications/data_sets/raw/IEEE_145bus_v23_PSLF.raw)
* [IEEE145.raw](../../src/applications/data_sets/raw/EEE145.raw])

### Larger Networks
* [240busWECC_2018_PSS_fixedshunt.raw](../../src/applications/data_sets/raw/240busWECC_2018_PSS_fixedshunt.raw):
  A reduced model of the WECC network developed at
  [NREL](https://www.nrel.gov/docs/fy21osti/74481.pdf) 
* [kundur-twoarea_v33.raw](../../src/applications/data_sets/raw/kundur-twoarea_v33.raw)
* [300bus_v23_no0imp_pslf.raw](../../src/applications/data_sets/raw/300bus_v23_no0imp_pslf.raw)
* [TAMU_500_rmsmallgen_v23.raw](../../src/applications/data_sets/raw/TAMU_500_rmsmallgen_v23.raw):
  The 500 bus network was created as part of the Texas A &amp; M Electric Grid Datasets
  project and is available
  [here](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/)
* [ACTIVSg500_rmsmallgen_pslfv23.raw](../../src/applications/data_sets/raw/ACTIVSg500_rmsmallgen_pslfv23.raw)
* [Polish_model_v23.raw](../../src/applications/data_sets/raw/Polish_model_v23.raw)
* [Polish_model_v33.raw](../../src/applications/data_sets/raw/Polish_model_v33.raw)
* [bus3000_gen_no0imp_v23_pslf.raw](../../src/applications/data_sets/raw/bus3000_gen_no0imp_v23_pslf.raw)
* [EuropeanOpenModel_v23.raw](../../src/applications/data_sets/raw/EuropeanOpenModel_v23.raw)
* [EuropeanOpenModel_v33.raw](../../src/applications/data_sets/raw/EuropeanOpenModel_v33.raw)
* [case_ACTIVSg2000_1.raw](../../src/applications/data_sets/dsa/case_ACTIVSg2000_1.raw) and
  [texas2000_hr_4882_mod_33.raw](../../src/applications/data_sets/dsa/texas2000_hr_4882_mod_33.raw):
  the Texas A&amp;M 2000 bus synthetic network used by the DSA applications

## PSS/E dyr files
These files are needed to run different dynamic simulation test cases. The dyr
formats are independent of version and any dyr file can be paired with any
version of a raw file.

### 9 Bus Network
Many of these test cases were created to test new models.
* [case9_GENROU.dyr](../../src/applications/data_sets/dyr/case9_GENROU.dyr)
* [case9_GENROU_ESST1A.dyr](../../src/applications/data_sets/dyr/case9_GENROU_ESST1A.dyr)
* [case9_GENROU_ESST1A_WSIEG1.dyr](../../src/applications/data_sets/dyr/case9_GENROU_ESST1A_WSIEG1.dyr])
* [case9_GENROU_EXDC1.dyr](../../src/applications/data_sets/dyr/case9_GENROU_EXDC1.dyr)
* [case9_GENROU_EXDC1_WSIEG1.dyr](../../src/applications/data_sets/dyr/case9_GENROU_EXDC1_WSIEG1.dyr)
* [case9_GENSAL.dyr](../../src/applications/data_sets/dyr/case9_GENSAL.dyr)
* [case9_GENSAL_ESST1A.dyr](../../src/applications/data_sets/dyr/case9_GENSAL_ESST1A.dyr)
* [case9_GENSAL_ESST1A_WSIEG1.dyr](../../src/applications/data_sets/dyr/case9_GENSAL_ESST1A_WSIEG1.dyr)
* [case9_HYGOV.dyr](../../src/applications/data_sets/dyr/case9_HYGOV.dyr)
* [case9_REGCA1_REECA1_REPCA1.dyr](../../src/applications/data_sets/dyr/case9_REGCA1_REECA1_REPCA1.dyr)
* [case9_SEXS.dyr](../../src/applications/data_sets/dyr/case9_SEXS.dyr)
* [case9_GAST.dyr](../../src/applications/data_sets/dyr/case9_GAST.dyr)
* [9b3g_GENSAL.dyr](../../src/applications/data_sets/dyr/9b3g_GENSAL.dyr)
* [9b3g.dyr](../../src/applications/data_sets/dyr/9b3g.dyr)
* [3g9b_gensal_esst1a_wsieg1_1acmotor_bus5smallload.dyr](../../src/applications/data_sets/dyr/3g9b_gensal_esst1a_wsieg1_1acmotor_bus5smallload.dyr)

### IEEE 14 Bus Network
* [IEEE14.dyr](../../src/applications/data_sets/dyr/IEEE14.dyr)
* [IEEE14_classicGen.dyr](../../src/applications/data_sets/dyr/IEEE14_classicGen.dyr)

### Larger Networks
* [IEEE_145b_classical_model.dyr](../../src/applications/data_sets/dyr/IEEE_145b_classical_model.dyr)
* [IEEE145_classicGen.dyr](../../src/applications/data_sets/dyr/IEEE145_classicGen.dyr)
* [kundur-twoarea_4renewable_mech.dyr](../../src/applications/data_sets/dyr/kundur-twoarea_4renewable_mech.dyr)
* [kundur-twoarea_4renewable_regcc1.dyr](../../src/applications/data_sets/dyr/kundur-twoarea_4renewable_regcc1.dyr)
* [kundur-twoarea.dyr](../../src/applications/data_sets/dyr/kundur-twoarea.dyr)
* [tamu_500bus_detail.dyr](../../src/applications/data_sets/dyr/tamu_500bus_detail.dyr)
* [240busWECC_2018_PSS_mod.dyr](../../src/applications/data_sets/dyr/240busWECC_2018_PSS_mod.dyr)
* [300bus_detail_model_cmpld_combine.dyr](../../src/applications/data_sets/dyr/300bus_detail_model_cmpld_combine.dyr)
* [classical_model_3000bus.dyr](../../src/applications/data_sets/dyr/classical_model_3000bus.dyr)

## Third-party data set attributions

The cases below originate outside GridPACK and were modified by the GridPACK
team; they are not the original files. Cite the original source when
publishing results based on them.

* **ACTIVSg500 and ACTIVSg2000 (Texas A&amp;M synthetic grids):** `TAMU_500_*`,
  `ACTIVSg500_*`, `tamu_500bus_detail.dyr`, `case_ACTIVSg2000_1*`,
  `texas2000_*`, `tamu_2000_*`, and the fault, watch-list, and wind-scenario
  files under `data_sets/dsa/`. Created by A.B. Birchfield, T. Xu,
  K.M. Gegner, K.S. Shetye, and T.J. Overbye, Texas A&amp;M University.
  Available for download at <https://electricgrids.engr.tamu.edu/>, free for
  commercial or non-commercial use. Synthetic systems, no CEII. Cite
  A.B. Birchfield et al., "Grid Structural Characteristics as Validation
  Criteria for Synthetic Networks," IEEE Trans. Power Syst., vol. 32, no. 4,
  pp. 3258-3265, 2017, and for the dynamics data T. Xu et al., "Creation of
  Synthetic Electric Grid Models for Transient Stability Studies," IREP
  Symposium, 2017.
* **European Open Model (PEGASE 13659-bus system):**
  `EuropeanOpenModel_v23.raw`, `EuropeanOpenModel_v33.raw`. Copyright (c) 2016
  C. Josz, S. Fliscounakis, J. Maeght, P. Panciatici;
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); distributed with
  MATPOWER as `case13659pegase`. Fictitious data, not for operation or planning
  of the European grid. Cite C. Josz et al., "AC Power Flow Data in MATPOWER
  and QCQP Format: iTesla, RTE Snapshots, and PEGASE,"
  <https://arxiv.org/abs/1603.01533>.
* **Polish model (MATPOWER case3120sp):** `Polish_model_v23.raw`,
  `Polish_model_v33.raw`. Polish 400/220/110 kV network, summer 2008 morning
  peak, provided by Roman Korab, Silesian University of Technology, and
  distributed with his permission through MATPOWER.
* **WECC 240-bus reduced model:** `240busWECC_2018_PSS_*`. Developed at NREL,
  <https://www.nrel.gov/docs/fy21osti/74481.pdf>.
