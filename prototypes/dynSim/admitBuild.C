#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <petsc.h>

#include "admitBuild.h"

static PetscErrorCode checkBuses(const int *bus_i, const int nbus, bool &check);
static PetscErrorCode buildYshMatrix(Mat &mysh, const int ns, PetscScalar ysh[],
				     const int nbus, const int nbrch, const int me);
static PetscErrorCode buildConnectMatrix(Mat &cf, Mat &ct, Mat &cft, Mat &ctt,
					 const int nbus, const int nbrch,
					 const int *line_from, const int *line_to,
					 const PetscScalar yff[], const PetscScalar yft[],
					 const PetscScalar ytf[], const PetscScalar ytt[],
					 const int me);
static PetscErrorCode buildDiagMatrices(const Mat &cf, const Mat &ct, const Vec &vyff,
					const Vec &vyft, const Vec &vytf, const Vec &vytt,
					Mat &mcff, Mat &mcft, Mat &mctf, Mat &mctt);
static PetscErrorCode createAndAssembleVectors(Vec &vyff, Vec &vyft, Vec &vytf, Vec &vytt,
					       const PetscScalar yff[], const PetscScalar yft[],
					       const PetscScalar ytf[], const PetscScalar ytt[],
					       const int nbrch, const int nbus);
static PetscErrorCode destroyVectors(Vec &vyff, Vec &vyft, Vec &vytf, Vec &vytt);
static PetscErrorCode matMatSpecialProduct(const Mat &mat, const Vec &vec, Mat &rmat);
static PetscErrorCode vecAllgather(const Vec &vec, Vec &localVec);
static PetscErrorCode nzPrefixSum(const Mat &mat, const int me, const int nproc,
				  std::vector<int> &psum);

PetscErrorCode AdmitBuild::buildMatrix(Mat &ybus)
{
  PetscErrorCode ierr;
  bool check;

  ierr = checkBuses(bus_i, nbus, check);
  CHKERRQ(ierr);

  if (!check)
    SETERRQ(PETSC_COMM_WORLD, 1, "Error: buses must appear in order by bus number");

  PetscScalar *ys = new PetscScalar[nbrch];
  PetscScalar *bc, *dlratio;
  PetscScalar *dlr, *dlx, *dbgs, *dbbs;
  PetscScalar *ytt = new PetscScalar[nbrch];
  PetscScalar *yff = new PetscScalar[nbrch];
  PetscScalar *yft = new PetscScalar[nbrch];
  PetscScalar *ytf = new PetscScalar[nbrch];
  PetscScalar *ysh = new PetscScalar[nbrch];

  ierr = VecGetArray(line_r, &dlr);
  CHKERRQ(ierr);
  ierr = VecGetArray(line_x, &dlx);
  CHKERRQ(ierr);
  ierr = VecGetArray(line_charge, &bc);
  CHKERRQ(ierr);
  ierr = VecGetArray(line_ratio, &dlratio);
  CHKERRQ(ierr);
  ierr = VecGetArray(bus_gs, &dbgs);
  CHKERRQ(ierr);
  ierr = VecGetArray(bus_bs, &dbbs);
  CHKERRQ(ierr);

  for (int i = 0; i < nbrch; i++) {
    ys[i] = 1.0 / (PetscScalar(PetscRealPart(dlr[i]), PetscRealPart(dlx[i])));

    if (dlratio[i] == 0.0)
      dlratio[i] = 1.0;

    ytt[i] = ys[i] + PetscScalar(0.0, PetscRealPart(bc[i] / 2.0));
    yff[i] = ytt[i] / (dlratio[i] * dlratio[i]);
    yft[i] = -(ys[i] / dlratio[i]);
    ytf[i] = yft[i];
    ysh[i] = PetscScalar(PetscRealPart(dbgs[i]), PetscRealPart(dbbs[i]));
  }

  Mat mysh, cf, ct, cft, ctt;

  buildYshMatrix(mysh, ns, ysh, nbus, nbrch, me);
  buildConnectMatrix(cf, ct, cft, ctt, nbus, nbrch, line_from, line_to, yff, yft, ytf, ytt,
		     me);

  ierr = VecRestoreArray(line_r, &dlr);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(line_x, &dlx);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(line_charge, &bc);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(line_ratio, &dlratio);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(bus_gs, &dbgs);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(bus_bs, &dbbs);
  CHKERRQ(ierr);

  Vec vyff, vyft, vytf, vytt;
  Mat mcff, mcft, mctf, mctt;

  createAndAssembleVectors(vyff, vyft, vytf, vytt, yff, yft, ytf, ytt, nbrch, nbus);
  buildDiagMatrices(cf, ct, vyff, vyft, vytf, vytt, mcff, mcft, mctf, mctt);
  destroyVectors(vyff, vyft, vytf, vytt);

  Mat y1, y2, y3, y4;

  ierr = MatMatMult(mcff, cft, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y1);
  CHKERRQ(ierr);
  ierr = MatMatMult(mcft, ctt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y2);
  CHKERRQ(ierr);
  ierr = MatMatMult(mctf, cft, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y3);
  CHKERRQ(ierr);
  ierr = MatMatMult(mctt, ctt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y4);
  CHKERRQ(ierr);

  ierr = MatAXPY(y2, PetscScalar(1.0), y1, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ(ierr);
  ierr = MatAXPY(y2, PetscScalar(1.0), y3, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ(ierr);
  ierr = MatAXPY(y2, PetscScalar(1.0), y4, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ(ierr);

  int y2m, y2n, myshm, myshn;

  ierr = MatGetSize(y2, &y2m, &y2n);
  CHKERRQ(ierr);
  ierr = MatGetSize(mysh, &myshm, &myshn);
  CHKERRQ(ierr);

  if (me == 0) {
    std::cout << "y2: " << y2m << " x " << y2n << std::endl;
    std::cout << "mysh: " << myshm << " x " << myshn << std::endl;
  }

  // ierr = MatAXPY(y2, PetscScalar(1.0), mysh, DIFFERENT_NONZERO_PATTERN);
  // CHKERRQ(ierr);

  ierr = MatDuplicate(y2, MAT_COPY_VALUES, &ybus);
  CHKERRQ(ierr);

  ierr = MatDestroy(&mcff);
  CHKERRQ(ierr);
  ierr = MatDestroy(&mcft);
  CHKERRQ(ierr);
  ierr = MatDestroy(&mctf);
  CHKERRQ(ierr);
  ierr = MatDestroy(&mctt);
  CHKERRQ(ierr);

  ierr = MatDestroy(&y1);
  CHKERRQ(ierr);
  ierr = MatDestroy(&y2);
  CHKERRQ(ierr);
  ierr = MatDestroy(&y3);
  CHKERRQ(ierr);
  ierr = MatDestroy(&y4);
  CHKERRQ(ierr);

  ierr = MatDestroy(&mysh);
  CHKERRQ(ierr);
  ierr = MatDestroy(&cf);
  CHKERRQ(ierr);
  ierr = MatDestroy(&ct);
  CHKERRQ(ierr);
  ierr = MatDestroy(&cft);
  CHKERRQ(ierr);
  ierr = MatDestroy(&ctt);
  CHKERRQ(ierr);

  delete [] ys;
  delete [] ytt;
  delete [] yff;
  delete [] yft;
  delete [] ytf;
  delete [] ysh;

  PetscFunctionReturn(0);
} // buildMatrix

PetscErrorCode checkBuses(const int *bus_i, const int nbus, bool &check)
{
  bool ret = true;

  for (int i = 0; i < nbus; i++) // zero-based array (modified from input data)
    if (bus_i[i] != i) {
      ret = false;
      break;
    }

  check = ret;

  PetscFunctionReturn(0);
} // checkBuses

PetscErrorCode buildYshMatrix(Mat &mysh, const int ns, PetscScalar ysh[],
			      const int nbus, const int nbrch, const int me)
{
  PetscErrorCode ierr;
  int nproc, nrows;

  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

  nrows = (nbus + nproc - 1) / nproc; // ceiling

  ierr = MatCreate(PETSC_COMM_WORLD, &mysh);
  CHKERRQ(ierr);
  ierr = MatSetType(mysh, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(mysh, PETSC_DECIDE, PETSC_DECIDE, nbus, nbus);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(mysh);
  CHKERRQ(ierr);
  ierr = MatSetUp(mysh);
  CHKERRQ(ierr);

  int rpos = me * nrows, nrpos = std::min((me + 1) * nrows, nbus);
  
  for (int i = rpos; i < nrpos; i++) {
    if (i < nbrch) {
      ierr = MatSetValue(mysh, i, i, ysh[i], INSERT_VALUES);
      CHKERRQ(ierr);
    }
    else {
      ierr = MatSetValue(mysh, i, i, PetscScalar(0.0), INSERT_VALUES);
      CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(mysh, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mysh, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
} // buildYshMatrix

PetscErrorCode buildConnectMatrix(Mat &cf, Mat &ct, Mat &cft, Mat &ctt,
				  const int nbus, const int nbrch,
				  const int *line_from, const int *line_to,
				  const PetscScalar yff[], const PetscScalar yft[],
				  const PetscScalar ytf[], const PetscScalar ytt[],
				  const int me)
{
  PetscErrorCode ierr;
  int nproc, nrows;

  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

  nrows = (nbus + nproc - 1) / nproc; // ceiling

  ierr = MatCreate(PETSC_COMM_WORLD, &cf);
  CHKERRQ(ierr);
  ierr = MatSetType(cf, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(cf, PETSC_DECIDE, PETSC_DECIDE, nbus, nbrch);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(cf);
  CHKERRQ(ierr);
  ierr = MatSetUp(cf);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &ct);
  CHKERRQ(ierr);
  ierr = MatSetType(ct, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(ct, PETSC_DECIDE, PETSC_DECIDE, nbus, nbrch);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(ct);
  CHKERRQ(ierr);
  ierr = MatSetUp(ct);
  CHKERRQ(ierr);


  ierr = MatCreate(PETSC_COMM_WORLD, &cft);
  CHKERRQ(ierr);
  ierr = MatSetType(cft, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(cft, PETSC_DECIDE, PETSC_DECIDE, nbrch, nbus);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(cft);
  CHKERRQ(ierr);
  ierr = MatSetUp(cft);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &ctt);
  CHKERRQ(ierr);
  ierr = MatSetType(ctt, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(ctt, PETSC_DECIDE, PETSC_DECIDE, nbrch, nbus);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(ctt);
  CHKERRQ(ierr);
  ierr = MatSetUp(ctt);
  CHKERRQ(ierr);

  int rpos = me * nrows, nrpos = std::min((me + 1) * nrows, nbus);
  unsigned int nz = 0;
  
  for (int i = rpos; i < nrpos; i++) {
    for (int j = 0; j < nbrch; j++) {
      if (line_from[j] == i) { // All arrays are zero-based, including read data
	ierr = MatSetValue(cf, i, j, PetscScalar(1.0), INSERT_VALUES);
	CHKERRQ(ierr);
	nz++;
      }
      if (line_to[j] == i) {
	ierr = MatSetValue(ct, i, j, PetscScalar(1.0), INSERT_VALUES);
	CHKERRQ(ierr);
      }
    }
  }

  nrows = (nbrch + nproc - 1) / nproc; // ceiling

  rpos = me * nrows;
  nrpos = std::min((me + 1) * nrows, nbrch);

  for (int i = rpos; i < nrpos; i++) {
    ierr = MatSetValue(cft, i, line_from[i], PetscScalar(1.0), INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = MatSetValue(ctt, i, line_to[i], PetscScalar(1.0), INSERT_VALUES);
    CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(cf, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(ct, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(cft, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(ctt, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  ierr = MatAssemblyEnd(ct, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(cf, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(cft, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ctt, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // buildConnectMatrix

PetscErrorCode buildDiagMatrices(const Mat &cf, const Mat &ct, const Vec &vyff,
				 const Vec &vyft, const Vec &vytf, const Vec &vytt,
				 Mat &mcff, Mat &mcft, Mat &mctf, Mat &mctt)
{
  PetscErrorCode ierr;

  ierr = MatDuplicate(cf, MAT_DO_NOT_COPY_VALUES, &mcff);
  CHKERRQ(ierr);
  ierr = MatDuplicate(cf, MAT_DO_NOT_COPY_VALUES, &mcft);
  CHKERRQ(ierr);

  ierr = MatDuplicate(ct, MAT_DO_NOT_COPY_VALUES, &mctf);
  CHKERRQ(ierr);
  ierr = MatDuplicate(ct, MAT_DO_NOT_COPY_VALUES, &mctt);
  CHKERRQ(ierr);

  ierr = matMatSpecialProduct(cf, vyff, mcff);
  CHKERRQ(ierr);
  ierr = matMatSpecialProduct(cf, vyft, mcft);
  CHKERRQ(ierr);
  ierr = matMatSpecialProduct(ct, vytf, mctf);
  CHKERRQ(ierr);
  ierr = matMatSpecialProduct(ct, vytt, mctt);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // buildDiagMatrices

PetscErrorCode createAndAssembleVectors(Vec &vyff, Vec &vyft, Vec &vytf, Vec &vytt,
					const PetscScalar yff[], const PetscScalar yft[],
					const PetscScalar ytf[], const PetscScalar ytt[],
					const int nbrch, const int nbus)
{
  PetscErrorCode ierr;
  int *idxs = new int[nbrch];

  for (int i = 0; i < nbrch; i++)
    idxs[i] = i;

  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nbrch, &vyff);
  CHKERRQ(ierr);
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nbrch, &vyft);
  CHKERRQ(ierr);
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nbrch, &vytf);
  CHKERRQ(ierr);
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nbrch, &vytt);
  CHKERRQ(ierr);

  int lo, hi;

  ierr = VecGetOwnershipRange(vyff, &lo, &hi);
  CHKERRQ(ierr);
  ierr = VecSetValues(vyff, hi - lo, &idxs[lo], &yff[lo], INSERT_VALUES);
  CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(vyft, &lo, &hi);
  CHKERRQ(ierr);
  ierr = VecSetValues(vyft, hi - lo, &idxs[lo], &yft[lo], INSERT_VALUES);
  CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(vytf, &lo, &hi);
  CHKERRQ(ierr);
  ierr = VecSetValues(vytf, hi - lo, &idxs[lo], &ytf[lo], INSERT_VALUES);
  CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(vytt, &lo, &hi);
  CHKERRQ(ierr);
  ierr = VecSetValues(vytt, hi - lo, &idxs[lo], &ytt[lo], INSERT_VALUES);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(vyff);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vyft);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vytf);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vytt);
  CHKERRQ(ierr);

  delete [] idxs;

  ierr = VecAssemblyEnd(vyff);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vyft);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vytf);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vytt);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // createAndAssembleVectors

PetscErrorCode destroyVectors(Vec &vyff, Vec &vyft, Vec &vytf, Vec &vytt)
{
  PetscErrorCode ierr;

  ierr = VecDestroy(&vyff);
  CHKERRQ(ierr);
  ierr = VecDestroy(&vyft);
  CHKERRQ(ierr);
  ierr = VecDestroy(&vytf);
  CHKERRQ(ierr);
  ierr = VecDestroy(&vytt);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // destroyVectors

#undef TRACE_LOGS

PetscErrorCode matMatSpecialProduct(const Mat &mat, const Vec &vec, Mat &rmat)
{
  PetscErrorCode ierr;
  int lo, hi, vlo, vhi;
  int me, nproc;

  MPI_Comm_rank(PETSC_COMM_WORLD, &me);
  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

#ifdef TRACE_LOGS
  std::stringstream sstr;
  std::string fname;

  sstr << "trace.dat." << me;

  sstr >> fname;

  std::ofstream ofs(fname.c_str(), std::ios_base::out | std::ios_base::app);
#endif

  ierr = MatGetOwnershipRange(mat, &lo, &hi);
  CHKERRQ(ierr);

  PetscScalar *dlvec;
  VecScatter ctx;
  Vec localVec;

  ierr = VecScatterCreateToAll(vec, &ctx, &localVec);
  CHKERRQ(ierr);

  // vecAllgather(vec, localVec);

  ierr = VecGetOwnershipRange(localVec, &vlo, &vhi);
  CHKERRQ(ierr);

#ifdef TRACE_LOGS
  ofs << "Matrix ownership: " << lo << ", " << hi << std::endl;
  ofs << "localVec indices: (" << vlo << ", " << vhi << ")" << std::endl;
#endif

  ierr = VecScatterBegin(ctx, vec, localVec, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx, vec, localVec, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ(ierr);

#ifdef TRACE_LOGS
  ofs << "Scatter succesful!" << std::endl;
#endif

  ierr = VecGetArray(localVec, &dlvec);
  CHKERRQ(ierr);

#if 0  
  std::vector<int> psum;
  int pos = 0;
  
  nzPrefixSum(mat, me, nproc, psum);
#endif

#ifdef TRACE_LOGS
  for (int i = 0; i < nproc; i++)
    ofs << "Rank: " << i << ", accumulated non-zeroes: " << psum[i] << std::endl;
#endif

  for (int i = lo; i < hi; i++) {
    int ncols;
    const int *cols;
    const PetscScalar *nzdat;

    ierr = MatGetRow(mat, i, &ncols, &cols, &nzdat);
    CHKERRQ(ierr);

    for (int j = 0; j < ncols; j++) {
      if ((cols[j] < vlo) || (cols[j] >= vhi))
	SETERRQ4(PETSC_COMM_WORLD, 1, "Rank: %d, vector index out of range: %d (%d %d)",
		 me, cols[j], vlo, vhi);

      if (nzdat[j] != 0.0) {
	ierr = MatSetValue(rmat, i, cols[j], nzdat[j] * dlvec[cols[j]], INSERT_VALUES);
	CHKERRQ(ierr);
#if 0
	pos++;
#endif
      }
    }

    ierr = MatRestoreRow(mat, i, &ncols, &cols, &nzdat);
    CHKERRQ(ierr);
  }

#ifdef TRACE_LOGS
  ofs << "Computation complete!" << std::endl;
#endif

  ierr = MatAssemblyBegin(rmat, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(localVec, &dlvec);
  CHKERRQ(ierr);

  ierr = VecScatterDestroy(&ctx);
  CHKERRQ(ierr);

  ierr = VecDestroy(&localVec);
  CHKERRQ(ierr);

  ierr = MatAssemblyEnd(rmat, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // matMatSpecialProduct

PetscErrorCode vecAllgather(const Vec &vec, Vec &localVec)
{
  PetscErrorCode ierr;
  int lsz, *lszs, *psum, nproc;

  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
  lszs = new int[nproc];
  psum = new int[nproc];

  ierr = VecGetLocalSize(vec, &lsz);
  CHKERRQ(ierr);
  MPI_Allgather(&lsz, 1, MPI_INT, lszs, 1, MPI_INT, PETSC_COMM_WORLD);

  int sum = 0;

  for (int i = 0; i < nproc; i++) {
    psum[i] = sum;
    sum += lszs[i];
  }

  PetscScalar *ldat = new PetscScalar[sum];
  PetscScalar *dat;
  int *idxs = new int[sum];

  ierr = VecGetArray(vec, &dat);
  CHKERRQ(ierr);

  MPI_Allgatherv(dat, lsz, MPIU_SCALAR, ldat, lszs, psum, MPIU_SCALAR, PETSC_COMM_WORLD);

  ierr = VecRestoreArray(vec, &dat);
  CHKERRQ(ierr);

  delete [] lszs;
  delete [] psum;

  for (int i = 0; i < sum; i++)
    idxs[i] = i;

  ierr = VecCreateSeq(PETSC_COMM_SELF, sum, &localVec);
  CHKERRQ(ierr);

  ierr = VecSetValues(localVec, sum, idxs, ldat, INSERT_VALUES);
  CHKERRQ(ierr);

  delete [] idxs;
  delete [] ldat;  

  PetscFunctionReturn(0);
} // vecAllGather

PetscErrorCode nzPrefixSum(const Mat &mat, const int me, const int nproc,
			   std::vector<int> &psum)
{
  PetscErrorCode ierr;
  int nzcnt = 0;
  int lo, hi;

  ierr = MatGetOwnershipRange(mat, &lo, &hi);
  CHKERRQ(ierr);

  for (int i = lo; i < hi; i++) {
    int ncols;
    const PetscScalar *nzdat;

    ierr = MatGetRow(mat, i, &ncols, PETSC_NULL, &nzdat);
    CHKERRQ(ierr);

    for (int j = 0; j < ncols; j++)
      if (nzdat[j] != 0.0)
	nzcnt++;

    ierr = MatRestoreRow(mat, i, &ncols, PETSC_NULL, &nzdat);
    CHKERRQ(ierr);
  }

  int *nzcnts = new int[nproc];

  MPI_Allgather(&nzcnt, 1, MPI_INT, nzcnts, 1, MPI_INT, PETSC_COMM_WORLD);

  int sum = 0;

  for (int i = 0; i < nproc; i++) {
    psum.push_back(sum);
    sum += nzcnts[i];
  }

  delete [] nzcnts;
  
  PetscFunctionReturn(0);
} // nzPrefixSum
