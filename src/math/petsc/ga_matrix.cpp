/* -------------------------------------------------------------
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 * ------------------------------------------------------------- */
/**
 * @file   ga_matrix.c
 * @author William A. Perkins
 * @date   2019-07-29 12:16:09 d3g096
 * 
 * @brief  
 * 
 * 
 */


#include <iostream>
#include <vector>
#include <boost/assert.hpp>

#include <ga.h>
#include "ga_matrix.hpp"


// -------------------------------------------------------------
// Determine the Global Arrays data type from PETSc configuration
// -------------------------------------------------------------
#if defined(PETSC_USE_REAL_DOUBLE)

#if defined(PETSC_USE_COMPLEX)
typedef DoubleComplex PetscScalarGA;
#define MT_PETSC_SCALAR MT_C_DCPL

#else 
typedef double PetscScalarGA;
#define MT_PETSC_SCALAR MT_C_DBL
#endif

#else

#if defined(PETSC_USE_COMPLEX)
typedef SingleComplex PetscScalarGA;
#define MT_PETSC_SCALAR MT_C_SCPL;
#else
typedef float PetscScalarGA;
#define MT_PETSC_SCALAR MT_C_FLOAT
#endif

#endif

static const PetscScalarGA one =
#if defined(PETSC_USE_COMPLEX)
  { 1.0, 0.0 }
#else
  1.0
#endif
  ;

static const PetscScalarGA zero =
#if defined(PETSC_USE_COMPLEX)
  { 0.0, 0.0 }
#else
  0.0
#endif
  ;
  


// -------------------------------------------------------------
// struct MatGACtx
// -------------------------------------------------------------
struct MatGACtx {
  int gaGroup;                  /**< GA process group handle */
  int ga;                       /**< GA handle */
  int lo[2];                    /**< Lower limits of "local" matrix ownership */
  int hi[2];                    /**< Upper limits of "local" matrix ownership */
};

static PetscErrorCode MatSetOperations_DenseGA(Mat A);

// -------------------------------------------------------------
// CreateMatGA
// -------------------------------------------------------------
static
PetscErrorCode
CreateMatGA(int pgroup, int lrows, int lcols, int grows, int gcols, int *ga)
{
  PetscErrorCode ierr = 0;

  /* Try to honor local ownership request (of rows). */

  int nprocs = GA_Pgroup_nnodes(pgroup);
  int me = GA_Pgroup_nodeid(pgroup);
  std::vector<int> tmapc(nprocs+1);
  std::vector<int> mapc(nprocs+1);
  int i;

  for (i = 0; i < nprocs+1; i++) tmapc[i] = 0;
  tmapc[me] = lrows;
  GA_Pgroup_igop(pgroup, &tmapc[0], nprocs+1, "+");
  mapc[0] = 0;
  for (i = 1; i < nprocs; i++) mapc[i] = mapc[i-1]+tmapc[i-1];
  mapc[nprocs] = 0;

  int dims[2] = {grows, gcols};
  int blocks[2] = { nprocs, 1 };


  *ga = NGA_Create_irreg_config(MT_PETSC_SCALAR, 2, &dims[0],
                                "PETSc Matrix in GA",
                                &blocks[0], &mapc[0],
                                pgroup);

  PetscScalar z(0.0);
  GA_Fill(*ga, &z);
  
  return ierr;
}

// -------------------------------------------------------------
// MPIComm2GApgroup
// -------------------------------------------------------------
static
PetscErrorCode
MPIComm2GApgroup(MPI_Comm comm, int *pGrpHandle)
{
  PetscErrorCode ierr = 0;
  int nproc;
  int me, myGlobalRank;
  int *proclist;
  int p;

  ierr = MPI_Comm_size(comm, &nproc); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &me); CHKERRQ(ierr);
  myGlobalRank = GA_Nodeid();
  ierr = PetscMalloc(nproc*sizeof(int), &proclist); CHKERRQ(ierr);
  for (p = 0; p < nproc; ++p) {
    proclist[p] = 0;
  }
  proclist[me] = myGlobalRank;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &proclist[0], nproc, MPI_INT, MPI_SUM, comm); CHKERRQ(ierr);
  *pGrpHandle = GA_Pgroup_create(&proclist[0], nproc); 
  ierr = PetscFree(proclist); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// Mat2GA
// This should only be called by MatMultbyGA(). It assumes the Mat is
// dense and the GA has the same ownership distribution.
// -------------------------------------------------------------
static
PetscErrorCode
Mat2GA(const Mat& A, int Aga)
{
  PetscErrorCode ierr(0);
  int me(GA_Nodeid());
  int r, rmin, rmax;
  int lo[2], hi[2], ld[2];
  int i, j;
  PetscInt grows, gcols, lrows, lcols;
  PetscScalar *matdata;
  ierr = MatGetSize(A, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A, &rmin, &rmax); CHKERRQ(ierr);

  NGA_Distribution(Aga, me, &lo[0], &hi[0]);
  
  std::vector<PetscInt> cidx(gcols);
  for (j = 0; j < gcols; ++j) cidx[j] = j;
  std::vector<PetscScalar> rdata(gcols);
  PetscScalarGA *gadata;

  for (r = rmin; r < rmax; ++r) {
    ierr = MatGetValues(A, 1, &r, gcols, &cidx[0], &rdata[0]); CHKERRQ(ierr);
    lo[0] = r; hi[0] = r;
    lo[1] = 0; hi[1] = gcols - 1;
    ld[0] = 1;
    NGA_Access(Aga, &lo[0], &hi[0], &gadata, &ld[0]);
    for (j = 0; j < gcols; ++j) {
#if defined(PETSC_USE_COMPLEX)
      gadata[j].real = rdata[j].real();
      gadata[j].imag = rdata[j].imag();
#else
      gadata[j] = rdata[j];
#endif
    }
    NGA_Release_update(Aga, &lo[0], &hi[0]);
  }
  GA_Sync();
  // GA_Print(Aga);
  return ierr;
}

// -------------------------------------------------------------
// GA2Mat
// This should only be called by MatMultbyGA(). It assumes the Mat is
// dense and the GA has the same ownership distribution.
// -------------------------------------------------------------
static
PetscErrorCode
GA2Mat(const int Aga, Mat& A)
{
  PetscErrorCode ierr(0);
  int me(GA_Nodeid());
  int r, rmin, rmax;
  int lo[2], hi[2], ld[2];
  int i, j;
  PetscInt grows, gcols, lrows, lcols;
  PetscScalar *matdata;
  ierr = MatGetSize(A, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A, &rmin, &rmax); CHKERRQ(ierr);

  NGA_Distribution(Aga, me, &lo[0], &hi[0]);
  
  std::vector<PetscInt> cidx(gcols);
  for (j = 0; j < gcols; ++j) cidx[j] = j;
  std::vector<PetscScalar> rdata(gcols);
  PetscScalarGA *gadata;

  for (r = rmin; r < rmax; ++r) {
    lo[0] = r; hi[0] = r;
    lo[1] = 0; hi[1] = gcols - 1;
    ld[0] = 1;
    NGA_Access(Aga, &lo[0], &hi[0], &gadata, &ld[0]);
    for (j = 0; j < gcols; ++j) {
#if defined(PETSC_USE_COMPLEX)
      rdata[j].real(gadata[j].real);
      rdata[j].imag(gadata[j].imag);
#else
      rdata[j] = gadata[j];
#endif
    }
    ierr = MatSetValues(A, 1, &r, gcols, &cidx[0], &rdata[0], INSERT_VALUES); CHKERRQ(ierr);
    NGA_Release(Aga, &lo[0], &hi[0]);
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  return ierr;
}

// -------------------------------------------------------------
// MatMultbyGA
// -------------------------------------------------------------
PetscErrorCode
MatMultbyGA(const Mat& A, const Mat& B, Mat& C)
{
  PetscErrorCode ierr(0);
  MPI_Comm comm;
  int lrows, grows, lcols, gcols;
  int pgp, opgp;
  int Aga, Bga, Cga;
  int m, n, k;

  ierr = PetscObjectGetComm((PetscObject)A, &comm); CHKERRQ(ierr);
  opgp = GA_Pgroup_get_default();

  ierr = MPIComm2GApgroup(comm, &pgp);
  GA_Pgroup_set_default(pgp);

  ierr = MatGetSize(A, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols); CHKERRQ(ierr);
  ierr = CreateMatGA(pgp, lrows, lcols, grows, gcols, &Aga); CHKERRQ(ierr);
  ierr = Mat2GA(A, Aga); CHKERRQ(ierr);
  m = grows;
  k = gcols;

  ierr = MatGetSize(B, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(B, &lrows, &lcols); CHKERRQ(ierr);
  ierr = CreateMatGA(pgp, lrows, lcols, grows, gcols, &Bga); CHKERRQ(ierr);
  ierr = Mat2GA(B, Bga); CHKERRQ(ierr);
  n = gcols;

  ierr = MatGetSize(C, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(C, &lrows, &lcols); CHKERRQ(ierr);
  ierr = CreateMatGA(pgp, lrows, lcols, grows, gcols, &Cga); CHKERRQ(ierr);

  char no('n');
  PetscScalarGA alpha, beta;

#if defined(PETSC_USE_REAL_DOUBLE)
#  if defined(PETSC_USE_COMPLEX)

  alpha.real = 1.0;
  alpha.imag = 0.0;
  beta.real = 0.0;
  beta.imag = 0.0;

  GA_Zgemm(no, no, m, n, k, one, Aga, Bga, beta, Cga);

#  else 

  alpha = 1.0;
  beta = 0.0;

  GA_Dgemm(no, no, m, n, k, one, Aga, Bga, beta, Cga);

#  endif
#else
#  if defined(PETSC_USE_COMPLEX)
#    error "no single precision complex"
#  else

  alpha = 1.0;
  beta = 0.0;
  
  GA_Sgemm(no, no, m, n, k, one, Aga, Bga, beta, Cga);

#  endif
#endif

  ierr = GA2Mat(Cga, C); CHKERRQ(ierr);

  GA_Pgroup_sync(pgp);

  GA_Destroy(Aga);
  GA_Destroy(Bga);
  GA_Destroy(Cga);
  
  GA_Pgroup_set_default(opgp);
  GA_Pgroup_destroy(pgp);
  
  return ierr;
}

// -------------------------------------------------------------
// MatMultbyGA_new
// -------------------------------------------------------------
PetscErrorCode
MatMultbyGA_new(const Mat& A, const Mat& B, Mat& C)
{
  PetscErrorCode ierr(0);
  MPI_Comm comm;
  int nproc;
  PetscInt grows_a, grows_b, grows_c;
  PetscInt lrows_a, lrows_b, lrows_c;
  PetscInt gcols_a, gcols_b, gcols_c;
  PetscInt lcols_a, lcols_b, lcols_c;

  ierr = MatGetSize(A, &grows_a, &gcols_a); CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &lrows_a, &lcols_a); CHKERRQ(ierr);

  ierr = MatGetSize(B, &grows_b, &gcols_b); CHKERRQ(ierr);
  ierr = MatGetLocalSize(B, &lrows_b, &lcols_b); CHKERRQ(ierr);

  grows_c = grows_a;
  lrows_c = lrows_a;
  gcols_c = gcols_b;
  lcols_c = lcols_b;

  ierr = PetscObjectGetComm((PetscObject)A, &comm); CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &nproc); 

  ierr = MatCreate(comm, &C); CHKERRQ(ierr);
  ierr = MatSetSizes(C, lrows_c, lcols_c, grows_c, gcols_c); CHKERRQ(ierr);
  if (nproc == 1) {
    ierr = MatSetType(C, MATSEQDENSE); CHKERRQ(ierr);
    ierr = MatSeqDenseSetPreallocation(C, PETSC_NULL); CHKERRQ(ierr);
  } else {
    ierr = MatSetType(C, MATDENSE); CHKERRQ(ierr);
    ierr = MatMPIDenseSetPreallocation(C, PETSC_NULL); CHKERRQ(ierr);
  }
    
  ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatMultbyGA(A, B, C); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// Vec2GA
// -------------------------------------------------------------
static
PetscErrorCode
Vec2GA(Vec x, int pgroup, int *ga, bool trans = false)
{
  int lrows, rows;
  PetscErrorCode ierr = 0;
  
  ierr = VecGetLocalSize(x, &lrows); CHKERRQ(ierr);
  ierr = VecGetSize(x, &rows); CHKERRQ(ierr);
  
  PetscInt vlo, vhi;
  ierr = VecGetOwnershipRange(x, &vlo, &vhi); CHKERRQ(ierr);
  
  const PetscScalar *v;
  ierr = VecGetArrayRead(x, &v); CHKERRQ(ierr);

  int lo[2] = {0,0}, hi[2] = {0,0}, ld[2] = {1,1};
  if (!trans) {
    ierr = CreateMatGA(pgroup, lrows, 1, rows, 1, ga); CHKERRQ(ierr);
    lo[0] = vlo; 
    hi[0] = vhi-1;
  } else {
    ierr = CreateMatGA(pgroup, 1, lrows, 1, rows, ga); CHKERRQ(ierr);
    lo[1] = vlo; 
    hi[1] = vhi-1;
  }
  // bad cast to get rid of const
  NGA_Put(*ga, lo, hi, (void *)v, ld);
  // GA_Print(*ga);
  ierr = VecRestoreArrayRead(x, &v); CHKERRQ(ierr);

  GA_Pgroup_sync(pgroup);
  return ierr;
}

// -------------------------------------------------------------
// GA2Vec
// -------------------------------------------------------------
static
int
GA2Vec(int ga, Vec x)
{
  int lrows, rows;
  PetscErrorCode ierr = 0;
  
  ierr = VecGetLocalSize(x, &lrows); CHKERRQ(ierr);
  ierr = VecGetSize(x, &rows); CHKERRQ(ierr);
  
  PetscInt vlo, vhi;
  ierr = VecGetOwnershipRange(x, &vlo, &vhi); CHKERRQ(ierr);
  
  PetscScalar *v;
  ierr = VecGetArray(x, &v); CHKERRQ(ierr);
   int lo[2] = {vlo,0}, hi[2] = {vhi-1,0}, ld[2] = {1,1};
  NGA_Get(ga, lo, hi, v, ld);
  ierr = VecRestoreArray(x, &v); CHKERRQ(ierr);
    
  return ierr;
}  

// -------------------------------------------------------------
// MatSetValues_DenseGA
// -------------------------------------------------------------
static 
PetscErrorCode
MatSetValues_DenseGA(Mat mat, 
                   PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], 
                   const PetscScalar v[],InsertMode addv)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;  
  int i, j, idx;
  PetscScalar vij, one(1.0);
  int lo[2], hi[2], ld[2] = {1, 1};
  ierr = MatShellGetContext(mat, (void *)&ctx); CHKERRQ(ierr);

  idx = 0;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j, ++idx) {
      lo[0] = idxm[i];
      hi[0] = idxm[i];
      lo[1] = idxn[j];
      hi[1] = idxn[j];
      vij = v[idx];
      switch (addv) {
      case INSERT_VALUES:
        NGA_Put(ctx->ga, lo, hi, (void *)&vij, ld);
        break;
      case ADD_VALUES:
        NGA_Acc(ctx->ga, lo, hi, (void *)&vij, ld, &one);
        break;
      default:
        BOOST_ASSERT_MSG(false, "Unknown set operation");
      }
    }
  }
  return ierr;
}

// -------------------------------------------------------------
// MatGetValues_DenseGA
// -------------------------------------------------------------
static 
PetscErrorCode
MatGetValues_DenseGA(Mat mat, 
                    PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], 
                    PetscScalar v[])
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;  
  int i, j, idx;
  PetscScalar vij;
  int lo[2], hi[2], ld[2] = {1, 1};
  ierr = MatShellGetContext(mat, (void *)&ctx); CHKERRQ(ierr);

  idx = 0;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j, ++idx) {
      lo[0] = idxm[i];
      hi[0] = idxm[i];
      lo[1] = idxn[j];
      hi[1] = idxn[j];
      NGA_Get(ctx->ga, lo, hi, (void *)&vij, ld);
      v[idx] = vij;
    }
  }
  return ierr;
}

// -------------------------------------------------------------
// MatAssemmblyBegin_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatAssemmblyBegin_DenseGA(Mat mat, MatAssemblyType type)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// MatAssemmblyEnd_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatAssemmblyEnd_DenseGA(Mat mat, MatAssemblyType type)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);
  GA_Pgroup_sync(ctx->gaGroup);
  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)mat,&comm); CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);
  // GA_Print(ctx->ga);
  return ierr;
}


// -------------------------------------------------------------
// MatMult_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatMult_DenseGA(Mat mat, Vec x, Vec y)
{
  // FIXME: I'm assuming the Mat and Vec's are compatible and that's
  // been checked somewhere else. Probably a mistake.

  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);

  PetscInt Arows, Acols;
  ierr = MatGetSize(mat, &Arows, &Acols); CHKERRQ(ierr);

  int g_x, g_y;
  ierr = Vec2GA(x, ctx->gaGroup, &g_x, false); CHKERRQ(ierr);
  ierr = Vec2GA(y, ctx->gaGroup, &g_y, false); CHKERRQ(ierr);

  PetscScalarGA alpha(one), beta(zero);
  int ndim, itype, lo[2] = {0,0}, ahi[2], xhi[2], yhi[2];
  NGA_Inquire(ctx->ga, &itype, &ndim, ahi);
  ahi[0] -= 1; ahi[1] -= 1;
  NGA_Inquire(g_x, &itype, &ndim, xhi);
  xhi[0] -= 1; xhi[1] -= 1;
  NGA_Inquire(g_y, &itype, &ndim, yhi);
  yhi[0] -= 1; yhi[1] -= 1;
  // GA_Print(ctx->ga);
  // GA_Print(g_x);
  NGA_Matmul_patch('N', 'N', &alpha, &beta, 
                   ctx->ga, lo, ahi,
                   g_x, lo, xhi,
                   g_y, lo, yhi);

  GA_Pgroup_sync(ctx->gaGroup);
  // GA_Print(g_y);

  ierr = GA2Vec(g_y, y); CHKERRQ(ierr);

  GA_Destroy(g_y);
  GA_Destroy(g_x);

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)mat,&comm); CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);

  return ierr;
}

// -------------------------------------------------------------
// MatMatMult_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatMatMult_DenseGA(Mat A, Mat B, MatReuse scall, PetscReal fill, Mat *C)
{

  // matrix sizes appear to be checked before here, so we won't do it again

  PetscErrorCode ierr = 0;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)A, &comm); CHKERRQ(ierr);

  MatType atype, btype;
  ierr = MatGetType(A, &atype);
  ierr = MatGetType(B, &btype);

  PetscBool issame;
  PetscStrcmp(atype, btype, &issame);

  Mat Bga, Cga;
  struct MatGACtx *Actx, *Bctx, *Cctx;

  ierr = MatShellGetContext(A, &Actx); CHKERRQ(ierr);

  if (issame) {
    Bga = B;
  } else {
    ierr = MatConvertToDenseGA(B, &Bga); CHKERRQ(ierr);
  }
  ierr = MatShellGetContext(Bga, &Bctx); CHKERRQ(ierr);

  PetscInt lrows, lcols, grows, gcols, junk;

  ierr = MatGetSize(A, &grows, &junk);  CHKERRQ(ierr);
  ierr = MatGetSize(B, &junk, &gcols);  CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &lrows, &junk); CHKERRQ(ierr);
  ierr = MatGetLocalSize(B, &junk, &lcols); CHKERRQ(ierr);

  ierr = MatCreateDenseGA(comm, lrows, lcols, grows, gcols, &Cga); CHKERRQ(ierr);
  ierr = MatShellGetContext(Cga, &Cctx); CHKERRQ(ierr);

  PetscScalarGA alpha(one), beta(zero);
  int ndim, itype, lo[2] = {0,0}, ahi[2], bhi[2], chi[2];
  NGA_Inquire(Actx->ga, &itype, &ndim, ahi);
  ahi[0] -= 1; ahi[1] -= 1;
  NGA_Inquire(Bctx->ga, &itype, &ndim, bhi);
  bhi[0] -= 1; bhi[1] -= 1;
  NGA_Inquire(Cctx->ga, &itype, &ndim, chi);
  chi[0] -= 1; chi[1] -= 1;
  // GA_Print(Actx->ga);
  // GA_Print(Bctx->ga);
  NGA_Matmul_patch('N', 'N', &alpha, &beta, 
                   Actx->ga, lo, ahi,
                   Bctx->ga, lo, bhi,
                   Cctx->ga, lo, chi);
  // GA_Print(Cctx->ga);

  switch (scall) {
  case MAT_REUSE_MATRIX:
    ierr = MatCopy(Cga, *C, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  case MAT_INITIAL_MATRIX:
  default:
    ierr = MatDuplicate(Cga, MAT_COPY_VALUES, C); CHKERRQ(ierr);
    break;
  }

  if (!issame) ierr = MatDestroy(&Bga); CHKERRQ(ierr);
  ierr = MatDestroy(&Cga); CHKERRQ(ierr);
  
  return ierr;
}

// -------------------------------------------------------------
// MatView_DenseGA
// -------------------------------------------------------------


// -------------------------------------------------------------
// MatZeroEntries_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatZeroEntries_DenseGA(Mat A)
{
  PetscErrorCode ierr = 0;

  struct MatGACtx *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  GA_Zero(ctx->ga);

  return ierr;
}


// -------------------------------------------------------------
// MatDestroy_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatDestroy_DenseGA(Mat A)
{
  PetscErrorCode ierr = 0;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)A,&comm); CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);

  struct MatGACtx *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  GA_Pgroup_sync(ctx->gaGroup);
  GA_Destroy(ctx->ga);
  GA_Pgroup_destroy(ctx->gaGroup);
  ierr = PetscFree(ctx);
  return ierr;
}

// -------------------------------------------------------------
// MatTranspose_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatTranspose_DenseGA(Mat mat, MatReuse reuse, Mat *B)
{
  PetscErrorCode ierr = 0;
  
  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)mat, &comm); CHKERRQ(ierr);

  struct MatGACtx *ctx, *newctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);

  PetscInt lrows, grows, lcols, gcols;
  ierr = MatGetSize(mat, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(mat, &lrows, &lcols); CHKERRQ(ierr);

  ierr = PetscMalloc(sizeof(struct MatGACtx), &newctx); CHKERRQ(ierr);
  ierr = MPIComm2GApgroup(comm, &(newctx->gaGroup));

  ierr = CreateMatGA(newctx->gaGroup, lcols, lrows, gcols, grows, &(newctx->ga)); CHKERRQ(ierr);
  GA_Transpose(ctx->ga, newctx->ga);

  ierr = MatCreateShell(comm, lcols, lrows, gcols, grows, newctx, B); CHKERRQ(ierr);
  ierr = MatSetOperations_DenseGA(*B);

  return ierr;
}


// -------------------------------------------------------------
// MatCreateDenseGA
// -------------------------------------------------------------
PetscErrorCode
MatCreateDenseGA(MPI_Comm comm, 
                 PetscInt m,PetscInt n,PetscInt M,PetscInt N, 
                 Mat *A)
{  
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  int lrows, grows, lcols, gcols;

  ierr = PetscMalloc(sizeof(struct MatGACtx), &ctx); CHKERRQ(ierr);
  ierr = MPIComm2GApgroup(comm, &(ctx->gaGroup)); CHKERRQ(ierr);
  
  lrows = m; lcols = n;
  grows = M; gcols = N;

  if (lrows == PETSC_DECIDE || lrows == PETSC_DETERMINE ||
      grows == PETSC_DECIDE || grows == PETSC_DETERMINE) {
    ierr = PetscSplitOwnership(comm, &lrows, &grows); CHKERRQ(ierr);
  }

  if (lcols == PETSC_DECIDE || lcols == PETSC_DETERMINE ||
      gcols == PETSC_DECIDE || gcols == PETSC_DETERMINE) {
    ierr = PetscSplitOwnership(comm, &lcols, &gcols); CHKERRQ(ierr);
  }
  
  ierr = CreateMatGA(ctx->gaGroup, lrows, lcols, grows, gcols, &(ctx->ga)); CHKERRQ(ierr);
  ierr = MatCreateShell(comm, lrows, lcols, grows, gcols, ctx, A); CHKERRQ(ierr);
  ierr = MatSetOperations_DenseGA(*A);

  return ierr;
}

// -------------------------------------------------------------
// MatDuplicate_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatDuplicate_DenseGA(Mat mat, MatDuplicateOption op, Mat *M)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx, *newctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)mat, &comm); CHKERRQ(ierr);

  PetscInt lrows, grows, lcols, gcols;
  ierr = MatGetSize(mat, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(mat, &lrows, &lcols); CHKERRQ(ierr);

  ierr = PetscMalloc(sizeof(struct MatGACtx), &newctx); CHKERRQ(ierr);
  ierr = MPIComm2GApgroup(comm, &(newctx->gaGroup));
  
  ierr = CreateMatGA(newctx->gaGroup, 
                     lrows, lcols, grows, gcols, &(newctx->ga)); CHKERRQ(ierr);
  ierr = MatCreateShell(comm, lrows, lcols, grows, gcols, newctx, M); CHKERRQ(ierr);
  ierr = MatSetOperations_DenseGA(*M);

  PetscScalar z(0.0);
  switch (op) {
  case (MAT_COPY_VALUES):
    GA_Copy(ctx->ga, newctx->ga);
    break;
  default:
    GA_Fill(newctx->ga, &z);
    break;
  }

  GA_Pgroup_sync(newctx->gaGroup);

  return ierr;
}

// -------------------------------------------------------------
// MatConvertToDenseGA
// -------------------------------------------------------------
PetscErrorCode
MatConvertToDenseGA(Mat A, Mat *B)
{
  PetscErrorCode ierr = 0;
  MPI_Comm comm;
  int lrows, grows, lcols, gcols;

  ierr = PetscObjectGetComm((PetscObject)A, &comm); CHKERRQ(ierr);

  ierr = MatGetSize(A, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols); CHKERRQ(ierr);

  ierr = MatCreateDenseGA(comm, lrows, lcols, grows, gcols, B); CHKERRQ(ierr);
  ierr = MatCopy(A, *B, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

  
  return ierr;
}

// -------------------------------------------------------------
// MatConvertGAtoDense
// -------------------------------------------------------------
PetscErrorCode
MatConvertGAToDense(Mat A, Mat *B)
{
  PetscErrorCode ierr = 0;
  MPI_Comm comm;
  int nproc;
  struct MatGACtx *ctx;
  PetscInt lrows, grows, lcols, gcols, lo, hi;

  ierr = PetscObjectGetComm((PetscObject)A, &comm); CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &nproc); 

  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);

  ierr = MatGetSize(A, &grows, &gcols); CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols); CHKERRQ(ierr);

  ierr = MatCreate(comm, B); CHKERRQ(ierr);
  ierr = MatSetSizes(*B, lrows, lcols, grows, gcols); CHKERRQ(ierr);
  if (nproc == 1) {
    ierr = MatSetType(*B, MATSEQDENSE); CHKERRQ(ierr);
    ierr = MatSeqDenseSetPreallocation(*B, PETSC_NULL); CHKERRQ(ierr);
  } else {
    ierr = MatSetType(*B, MATDENSE); CHKERRQ(ierr);
    ierr = MatMPIDenseSetPreallocation(*B, PETSC_NULL); CHKERRQ(ierr);
  }
  ierr = MatGetOwnershipRange(*B, &lo, &hi); CHKERRQ(ierr);

  std::vector<PetscInt> cidx(gcols);
  for (PetscInt c = 0; c < gcols; ++c) {
    cidx[c] = c;
  }
  std::vector<PetscScalar> rowvals(gcols);
  for (PetscInt r = lo; r < hi; ++r) {
    int glo[2] = {r, 0};
    int ghi[2] = {r, gcols - 1};
    int ld[2] = {1,1};
    NGA_Get(ctx->ga, glo, ghi, &rowvals[0], ld);
    ierr = MatSetValues(*B, 1, &r, gcols, &cidx[0], &rowvals[0], INSERT_VALUES); CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(*B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// MatSetOperations_DenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatSetOperations_DenseGA(Mat A)
{
  PetscErrorCode ierr = 0;
  ierr = MatShellSetOperation(A, MATOP_SET_VALUES, (void(*)(void))MatSetValues_DenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_GET_VALUES, (void(*)(void))MatGetValues_DenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_ZERO_ENTRIES, (void(*)(void))MatZeroEntries_DenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))MatMult_DenseGA); CHKERRQ(ierr); 
  ierr = MatShellSetOperation(A, MATOP_MAT_MULT, (void(*)(void))MatMatMult_DenseGA); CHKERRQ(ierr); 
  ierr = MatShellSetOperation(A, MATOP_TRANSPOSE, (void(*)(void))MatTranspose_DenseGA); CHKERRQ(ierr); 
  ierr = MatShellSetOperation(A, MATOP_DUPLICATE, (void(*)(void))MatDuplicate_DenseGA); CHKERRQ(ierr); 
  ierr = MatShellSetOperation(A, MATOP_ASSEMBLY_BEGIN, (void(*)(void))MatAssemmblyBegin_DenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_ASSEMBLY_END, (void(*)(void))MatAssemmblyEnd_DenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_DESTROY,  (void(*)(void))MatDestroy_DenseGA); CHKERRQ(ierr);
  return ierr;
}

