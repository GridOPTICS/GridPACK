#ifndef __ADMITBUILD_H__
#define __ADMITBUILD_H__

#include <petsc.h>

class AdmitBuild {

public:

  explicit AdmitBuild(const int ome, const int onbus, const int onbrch,
		      const int ons, const int *obus_i, const Vec &obus_gs, const Vec &obus_bs,
		      const int *oline_from, const int *oline_to, const Vec &oline_r,
		      const Vec &oline_x, const Vec &oline_charge,
		      const Vec &oline_ratio, const Vec &oline_shift);
  ~AdmitBuild();

  PetscErrorCode buildMatrix(Mat &ybus);

private:

  const int me;
  const int nbus, nbrch;
  const int ns;
  const int *bus_i;
  Vec bus_gs, bus_bs;
  const int *line_from, *line_to;
  Vec line_r, line_x, line_charge, line_ratio, line_shift;

  AdmitBuild();
  AdmitBuild(const AdmitBuild &othis);
};

inline AdmitBuild::AdmitBuild(const int ome, const int onbus, const int onbrch,
			      const int ons, const int *obus_i, const Vec &obus_gs,
			      const Vec &obus_bs, const int *oline_from, const int *oline_to,
			      const Vec &oline_r, const Vec &oline_x, const Vec &oline_charge,
			      const Vec &oline_ratio, const Vec &oline_shift)
  : me(ome), nbus(onbus), nbrch(onbrch), ns(ons), bus_i(obus_i), bus_gs(obus_gs),
    bus_bs(obus_bs), line_from(oline_from), line_to(oline_to), line_r(oline_r),
    line_x(oline_x), line_charge(oline_charge), line_ratio(oline_ratio),
    line_shift(oline_shift)
{
} // AdmitBuild

inline AdmitBuild::AdmitBuild()
  : me(-1), nbus(0), nbrch(0), ns(0), bus_i(NULL), bus_gs(PETSC_NULL), bus_bs(PETSC_NULL),    
    line_from(NULL), line_to(NULL), line_r(PETSC_NULL), line_x(PETSC_NULL),
    line_charge(PETSC_NULL), line_ratio(PETSC_NULL), line_shift(PETSC_NULL)
{
} // AdmitBuild

inline AdmitBuild::AdmitBuild(const AdmitBuild &othis)
  : me(-1), nbus(0), nbrch(0), ns(0), bus_i(NULL), bus_gs(PETSC_NULL), bus_bs(PETSC_NULL),    
    line_from(NULL), line_to(NULL), line_r(PETSC_NULL), line_x(PETSC_NULL),
    line_charge(PETSC_NULL), line_ratio(PETSC_NULL), line_shift(PETSC_NULL)
{
} // AdmitBuild

inline AdmitBuild::~AdmitBuild()
{
} // ~AdmitBuild

#endif // __ADMITBUILD_H
