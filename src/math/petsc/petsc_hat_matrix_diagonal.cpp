// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_hat_matrix_diagonal.cpp
 * @brief  Diagonal of the weighted least squares hat matrix
 */
// -------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <petscmat.h>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include "hat_matrix_diagonal.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#ifdef GRIDPACK_HAVE_KLU
#include <klu.h>
#endif

namespace gridpack {
namespace math {

#define HAT_CHK(x) do { PetscErrorCode e_ = (x); CHKERRABORT(PETSC_COMM_SELF, e_); } while (0)

// Real H in row storage with the weight of each row
struct HatRows {
  int m, n;
  std::vector<int> ptr, col;
  std::vector<double> val, w;
};

#ifdef GRIDPACK_HAVE_KLU
// Selected inverse of G = H^T W H on the pattern of its factors; false if not usable
static bool
selectedInversionDiagonal(const HatRows& h, std::vector<double>& diag,
                          double& tfactor, std::string& note)
{
  int m = h.m, n = h.n;
  // Column storage of H
  std::vector<int> cp(n+1, 0), ci(h.col.size());
  std::vector<double> cx(h.col.size());
  for (size_t k = 0; k < h.col.size(); k++) cp[h.col[k]+1]++;
  for (int j = 0; j < n; j++) cp[j+1] += cp[j];
  std::vector<int> next(cp.begin(), cp.end()-1);
  for (int i = 0; i < m; i++) {
    for (int k = h.ptr[i]; k < h.ptr[i+1]; k++) {
      int q = next[h.col[k]]++;
      ci[q] = i;
      cx[q] = h.val[k];
    }
  }
  // G column by column, keeping every structural entry
  std::vector<int> gp(n+1, 0), gi, mark(n, -1), list;
  std::vector<double> gx, acc(n, 0.0);
  for (int b = 0; b < n; b++) {
    list.clear();
    for (int q = cp[b]; q < cp[b+1]; q++) {
      int i = ci[q];
      double s = h.w[i]*cx[q];
      for (int k = h.ptr[i]; k < h.ptr[i+1]; k++) {
        int a = h.col[k];
        if (mark[a] != b) {
          mark[a] = b;
          acc[a] = 0.0;
          list.push_back(a);
        }
        acc[a] += s*h.val[k];
      }
    }
    std::sort(list.begin(), list.end());
    for (size_t k = 0; k < list.size(); k++) {
      gi.push_back(list[k]);
      gx.push_back(acc[list[k]]);
    }
    gp[b+1] = static_cast<int>(gi.size());
  }

  // One block, no scaling, pivots kept on the diagonal (G is positive definite)
  klu_common common;
  klu_defaults(&common);
  common.btf = 0;
  common.scale = -1;
  common.tol = 1.0e-12;
  klu_symbolic *symbolic = klu_analyze(n, &gp[0], &gi[0], &common);
  if (symbolic == NULL) {
    note = "KLU analysis failed";
    return false;
  }
  klu_numeric *numeric = klu_factor(&gp[0], &gi[0], &gx[0], symbolic, &common);
  if (numeric == NULL) {
    klu_free_symbolic(&symbolic, &common);
    note = "KLU factorization failed (singular gain matrix)";
    return false;
  }
  klu_sort(symbolic, numeric, &common);
  int lnz = numeric->lnz, unz = numeric->unz;
  std::vector<int> Lp(n+1), Li(lnz), Up(n+1), Ui(unz), P(n), Q(n);
  std::vector<double> Lx(lnz), Ux(unz);
  int ok = klu_extract(numeric, symbolic, &Lp[0], &Li[0], &Lx[0],
      &Up[0], &Ui[0], &Ux[0], NULL, NULL, NULL, &P[0], &Q[0], NULL, NULL, &common);
  klu_free_numeric(&numeric, &common);
  klu_free_symbolic(&symbolic, &common);
  tfactor = MPI_Wtime();
  if (!ok) {
    note = "KLU factors could not be extracted";
    return false;
  }
  for (int k = 0; k < n; k++) {
    if (P[k] != Q[k]) {
      note = "KLU pivoted off the diagonal";
      return false;
    }
  }

  // Strictly lower part of L and the pivots
  std::vector<int> lp(n+1, 0), li;
  std::vector<double> lx, d(n, 0.0);
  li.reserve(lnz);
  lx.reserve(lnz);
  for (int j = 0; j < n; j++) {
    for (int p = Lp[j]; p < Lp[j+1]; p++) {
      if (Li[p] > j) {
        li.push_back(Li[p]);
        lx.push_back(Lx[p]);
      }
    }
    lp[j+1] = static_cast<int>(li.size());
    for (int p = Up[j]; p < Up[j+1]; p++) {
      if (Ui[p] == j) d[j] = Ux[p];
    }
    if (d[j] == 0.0) {
      note = "KLU returned a zero pivot";
      return false;
    }
  }

  // Takahashi recurrence: Z = G^-1 on the pattern of L, last column first
  std::vector<double> zd(n, 0.0), zl(li.size(), 0.0);
  for (int j = n-1; j >= 0; j--) {
    int p0 = lp[j], p1 = lp[j+1];
    for (int pk = p0; pk < p1; pk++) {
      int k = li[pk];
      zl[pk] += lx[pk]*zd[k];
      int q = lp[k], qe = lp[k+1];
      for (int pr = pk+1; pr < p1; pr++) {
        int r = li[pr];
        while (q < qe && li[q] < r) q++;
        if (q == qe || li[q] != r) {
          note = "factor pattern is not symmetric";
          return false;
        }
        zl[pr] += lx[pk]*zl[q];
        zl[pk] += lx[pr]*zl[q];
      }
    }
    double s = 0.0;
    for (int p = p0; p < p1; p++) {
      zl[p] = -zl[p];
      s += lx[p]*zl[p];
    }
    zd[j] = 1.0/d[j] - s;
  }

  // Quadratic form of each row of H with Z
  std::vector<int> pinv(n);
  for (int k = 0; k < n; k++) pinv[P[k]] = k;
  for (int i = 0; i < m; i++) {
    double sum = 0.0;
    for (int ka = h.ptr[i]; ka < h.ptr[i+1]; ka++) {
      int a = pinv[h.col[ka]];
      sum += h.val[ka]*h.val[ka]*zd[a];
      for (int kb = ka+1; kb < h.ptr[i+1]; kb++) {
        int b = pinv[h.col[kb]];
        if (a == b) {
          sum += 2.0*h.val[ka]*h.val[kb]*zd[a];
          continue;
        }
        int lo = std::min(a, b), hi = std::max(a, b);
        std::vector<int>::const_iterator first = li.begin()+lp[lo];
        std::vector<int>::const_iterator last = li.begin()+lp[lo+1];
        std::vector<int>::const_iterator it = std::lower_bound(first, last, hi);
        if (it == last || *it != hi) {
          note = "inverse entry outside the factor pattern";
          return false;
        }
        sum += 2.0*h.val[ka]*h.val[kb]*zl[it-li.begin()];
      }
    }
    diag[i] = h.w[i]*sum;
  }
  return true;
}
#endif

// One solve per row of H with a serial factorization, in blocks
static void
blockedSolvesDiagonal(const HatRows& h, std::vector<double>& diag, double& tfactor)
{
  int m = h.m, n = h.n;
  std::vector<PetscInt> ia(h.ptr.begin(), h.ptr.end()), ja(h.col.begin(), h.col.end());
  std::vector<PetscScalar> va(h.val.begin(), h.val.end()), wva(h.val.size());
  for (int i = 0; i < m; i++) {
    for (PetscInt k = ia[i]; k < ia[i+1]; k++) wva[k] = h.w[i]*va[k];
  }
  Mat Hs, WH, G, F, B, X;
  HAT_CHK(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, m, n, &ia[0], &ja[0], &va[0], &Hs));
  HAT_CHK(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, m, n, &ia[0], &ja[0], &wva[0], &WH));
  HAT_CHK(MatTransposeMatMult(Hs, WH, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &G));
  IS rperm, cperm;
  MatFactorInfo info;
  HAT_CHK(MatGetFactor(G, MATSOLVERPETSC, MAT_FACTOR_LU, &F));
  HAT_CHK(MatGetOrdering(G, MATORDERINGND, &rperm, &cperm));
  HAT_CHK(MatFactorInfoInitialize(&info));
  HAT_CHK(MatLUFactorSymbolic(F, G, rperm, cperm, &info));
  HAT_CHK(MatLUFactorNumeric(F, G, &info));
  tfactor = MPI_Wtime();

  // Columns of each block are rows of H
  int nb = std::min(128, m);
  HAT_CHK(MatCreateSeqDense(PETSC_COMM_SELF, n, nb, NULL, &B));
  HAT_CHK(MatCreateSeqDense(PETSC_COMM_SELF, n, nb, NULL, &X));
  for (int start = 0; start < m; start += nb) {
    int width = std::min(nb, m-start);
    PetscScalar *b;
    HAT_CHK(MatDenseGetArray(B, &b));
    std::fill(b, b+static_cast<size_t>(n)*nb, 0.0);
    for (int c = 0; c < width; c++) {
      int i = start+c;
      for (PetscInt k = ia[i]; k < ia[i+1]; k++) b[static_cast<size_t>(c)*n+ja[k]] = va[k];
    }
    HAT_CHK(MatDenseRestoreArray(B, &b));
    HAT_CHK(MatMatSolve(F, B, X));
    const PetscScalar *x;
    HAT_CHK(MatDenseGetArrayRead(X, &x));
    for (int c = 0; c < width; c++) {
      int i = start+c;
      double kii = 0.0;
      for (PetscInt k = ia[i]; k < ia[i+1]; k++) {
        kii += PetscRealPart(wva[k]*x[static_cast<size_t>(c)*n+ja[k]]);
      }
      diag[i] = kii;
    }
    HAT_CHK(MatDenseRestoreArrayRead(X, &x));
  }
  HAT_CHK(ISDestroy(&rperm));
  HAT_CHK(ISDestroy(&cperm));
  HAT_CHK(MatDestroy(&B));
  HAT_CHK(MatDestroy(&X));
  HAT_CHK(MatDestroy(&F));
  HAT_CHK(MatDestroy(&G));
  HAT_CHK(MatDestroy(&WH));
  HAT_CHK(MatDestroy(&Hs));
}

// -------------------------------------------------------------
// hatMatrixDiagonal
// -------------------------------------------------------------
void
hatMatrixDiagonal(const Matrix& H, const std::vector<double>& localWeights,
                  std::vector<double>& localDiagonal, HatMatrixDiagonalInfo& info,
                  HatMatrixDiagonalMethod method)
{
  boost::mpi::communicator comm =
    static_cast<boost::mpi::communicator>(H.communicator());
  double t0 = MPI_Wtime();
  int m = H.rows();
  int n = H.cols();
  info.rows = m;
  info.cols = n;
  info.method = "";
  info.note = "";
  info.gatherTime = info.factorTime = info.diagonalTime = 0.0;

  // Local rows of H; a complex entry may be stored as a real block
  const Mat *A = PETScMatrix(H);
  PetscInt pm, pn, plo, phi;
  HAT_CHK(MatGetSize(*A, &pm, &pn));
  HAT_CHK(MatGetOwnershipRange(*A, &plo, &phi));
  int es = (m > 0) ? pm/m : 1;
  int lo = plo/es;
  int nlocal = (phi-plo)/es;
  std::vector<int> lcount, lcol;
  std::vector<double> lval;
  for (PetscInt r = plo; r < phi; r += es) {
    PetscInt nc;
    const PetscInt *cols;
    const PetscScalar *vals;
    int count = 0;
    HAT_CHK(MatGetRow(*A, r, &nc, &cols, &vals));
    for (PetscInt k = 0; k < nc; k++) {
      if (cols[k]%es == 0 && PetscRealPart(vals[k]) != 0.0) {
        lcol.push_back(cols[k]/es);
        lval.push_back(PetscRealPart(vals[k]));
        count++;
      }
    }
    HAT_CHK(MatRestoreRow(*A, r, &nc, &cols, &vals));
    lcount.push_back(count);
  }
  std::vector<double> lw(localWeights);
  lw.resize(nlocal, 0.0);
  std::vector<std::vector<int> > gcount, gcol;
  std::vector<std::vector<double> > gval, gw;
  boost::mpi::gather(comm, lcount, gcount, 0);
  boost::mpi::gather(comm, lcol, gcol, 0);
  boost::mpi::gather(comm, lval, gval, 0);
  boost::mpi::gather(comm, lw, gw, 0);

  std::vector<double> diag(m, 0.0);
  if (comm.rank() == 0) {
    // Processes own consecutive rows, so the pieces join in order
    HatRows h;
    h.m = m;
    h.n = n;
    h.ptr.assign(1, 0);
    for (size_t p = 0; p < gcount.size(); p++) {
      for (size_t i = 0; i < gcount[p].size(); i++) {
        h.ptr.push_back(h.ptr.back()+gcount[p][i]);
      }
      h.col.insert(h.col.end(), gcol[p].begin(), gcol[p].end());
      h.val.insert(h.val.end(), gval[p].begin(), gval[p].end());
      h.w.insert(h.w.end(), gw[p].begin(), gw[p].end());
    }
    double t1 = MPI_Wtime();
    double t2 = t1;
    bool done = false;
    if (method == SelectedInversion) {
#ifdef GRIDPACK_HAVE_KLU
      done = selectedInversionDiagonal(h, diag, t2, info.note);
      if (done) info.method = "KLU selected inversion";
#else
      info.note = "GridPACK was built without KLU";
#endif
    }
    if (!done) {
      blockedSolvesDiagonal(h, diag, t2);
      info.method = "blocked solves";
    }
    double t3 = MPI_Wtime();
    info.gatherTime = t1-t0;
    info.factorTime = t2-t1;
    info.diagonalTime = t3-t2;
  }
  boost::mpi::broadcast(comm, diag, 0);
  boost::mpi::broadcast(comm, info.method, 0);
  boost::mpi::broadcast(comm, info.note, 0);
  boost::mpi::broadcast(comm, info.gatherTime, 0);
  boost::mpi::broadcast(comm, info.factorTime, 0);
  boost::mpi::broadcast(comm, info.diagonalTime, 0);
  localDiagonal.assign(diag.begin()+lo, diag.begin()+lo+nlocal);
}

} // namespace math
} // namespace gridpack
