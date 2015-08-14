/* -------------------------------------------------------------
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 * ------------------------------------------------------------- */
/**
 * @file   ga_matrix.c
 * @author William A. Perkins
 * @date   2015-08-06 14:33:06 d3g096
 * 
 * @brief  
 * 
 * 
 */


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
  int tmapc[nprocs+1];
  int mapc[nprocs+1];
  int i;

  for (i = 0; i < nprocs+1; i++) tmapc[i] = 0;
  tmapc[me] = lrows;
  GA_Pgroup_igop(pgroup, tmapc, nprocs+1, "+");
  mapc[0] = 0;
  for (i = 1; i < nprocs; i++) mapc[i] = mapc[i-1]+tmapc[i-1];
  mapc[nprocs] = 0;

  int dims[2] = {grows, gcols};
  int blocks[2] = { nprocs, 1 };

  *ga = GA_Create_handle();
  GA_Set_data(*ga, 2, dims, MT_PETSC_SCALAR);
  GA_Set_irreg_distr(*ga, mapc, blocks);
  GA_Set_pgroup(*ga, pgroup);
  if (!GA_Allocate(*ga)) {
    ierr = 1;
  }
  PetscScalar z(0.0);
  GA_Fill(*ga, &z);
  
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
  
  PetscScalar *v;
  ierr = VecGetArray(x, &v); CHKERRQ(ierr);

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
  NGA_Put(*ga, lo, hi, v, ld);
  // GA_Print(*ga);
  ierr = VecRestoreArray(x, &v); CHKERRQ(ierr);

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
  newctx->gaGroup = ctx->gaGroup;
  ierr = CreateMatGA(newctx->gaGroup, lcols, lrows, gcols, grows, &(newctx->ga)); CHKERRQ(ierr);
  GA_Transpose(ctx->ga, newctx->ga);

  ierr = MatCreateShell(comm, lcols, lrows, gcols, grows, newctx, B); CHKERRQ(ierr);
  ierr = MatSetOperations_DenseGA(*B);

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
  // GA_Pgroup_destroy(ctx->gaGroup);
  ierr = PetscFree(ctx);
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
  newctx->gaGroup = ctx->gaGroup;
  newctx->ga = GA_Duplicate(ctx->ga, "PETSc Dense Matrix");

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
    ierr = PetscSplitOwnership(comm, &lrows, &grows); CHKERRXX(ierr);
  }

  if (lcols == PETSC_DECIDE || lcols == PETSC_DETERMINE ||
      gcols == PETSC_DECIDE || gcols == PETSC_DETERMINE) {
    ierr = PetscSplitOwnership(comm, &lcols, &gcols); CHKERRXX(ierr);
  }
  
  ierr = CreateMatGA(ctx->gaGroup, lrows, lcols, grows, gcols, &(ctx->ga)); CHKERRQ(ierr);
  ierr = MatCreateShell(comm, lrows, lcols, grows, gcols, ctx, A); CHKERRQ(ierr);
  ierr = MatSetOperations_DenseGA(*A);

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

  ierr = MatCreateDenseGA(comm, lrows, lcols, grows, gcols, B); CHKERRXX(ierr);
  ierr = MatCopy(A, *B, SAME_NONZERO_PATTERN); CHKERRXX(ierr);

  
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

  ierr = MatCreateDense(comm, lrows, lcols, grows, gcols, NULL, B); CHKERRQ(ierr);
  
  ierr = MatCreate(comm, B); CHKERRXX(ierr);
  ierr = MatSetSizes(*B, lrows, lcols, grows, gcols); CHKERRXX(ierr);
  if (nproc == 1) {
    ierr = MatSetType(*B, MATSEQDENSE); CHKERRXX(ierr);
    ierr = MatSeqDenseSetPreallocation(*B, PETSC_NULL); CHKERRXX(ierr);
  } else {
    ierr = MatSetType(*B, MATDENSE); CHKERRXX(ierr);
    ierr = MatMPIDenseSetPreallocation(*B, PETSC_NULL); CHKERRXX(ierr);
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

