/* -------------------------------------------------------------
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 * ------------------------------------------------------------- */
/**
 * @file   ga_matrix.c
 * @author William A. Perkins
 * @date   2015-05-01 12:38:40 d3g096
 * 
 * @brief  
 * 
 * 
 */


#include <boost/assert.hpp>

#include <ga.h>
#include "ga_matrix.hpp"


// -------------------------------------------------------------
// Determine the Global Arrays data type from PETSc configuration
// -------------------------------------------------------------
#if defined(PETSC_USE_REAL_DOUBLE)

#if defined(PETSC_USE_COMPLEX)
#define MT_PETSC_SCALAR MT_C_DCPL
#else 
#define MT_PETSC_SCALAR MT_C_DBL
#endif

#else

#if defined(PETSC_USE_COMPLEX)
#define MT_PETSC_SCALAR MT_C_SCPL;
#else
#define MT_PETSC_SCALAR MT_C_FLOAT
#endif

#endif


// -------------------------------------------------------------
// struct MatGACtx
// -------------------------------------------------------------
struct MatGACtx {
  int gaGroup;                  /**< GA process group handle */
  int ga;                       /**< GA handle */
  int lo[2];                    /**< Lower limits of "local" matrix ownership */
  int hi[2];                    /**< Upper limits of "local" matrix ownership */
};

// -------------------------------------------------------------
// CreateMatGA
// -------------------------------------------------------------
static
PetscErrorCode
CreateMatGA(int pgroup, int lrows, int lcols, int grows, int gcols, int *ga)
{
  PetscErrorCode ierr = 0;
  int dims[2] = {grows, gcols};
  int chunk[2] = {1, gcols};

  /* Do not honor local ownership request, initially. */

  *ga = NGA_Create(MT_PETSC_SCALAR, 2, dims, "PETSc Dense Matrix", chunk);
  PetscScalar z(0.0);
  GA_Fill(*ga, &z);
  
  return ierr;
}

// -------------------------------------------------------------
// CreateVecGA
// -------------------------------------------------------------
static
PetscErrorCode
CreateVecGA(int pgroup, int lrows, int grows, int *ga)
{
  PetscErrorCode ierr = 0;
  int dims[1] = {grows};

  /* Do not honor local ownership request, initially. */

  *ga = NGA_Create(MT_PETSC_SCALAR, 1, dims, "PETSc Vector", NULL);
  
  return ierr;
}

// -------------------------------------------------------------
// Vec2GA
// -------------------------------------------------------------
static
PetscErrorCode
Vec2GA(Vec x, int pgroup, int *ga)
{
  int lrows, rows;
  int lo, hi, ld = 1;
  PetscErrorCode ierr = 0;
  PetscScalar *v;
  
  ierr = VecGetLocalSize(x, &lrows); CHKERRQ(ierr);
  ierr = VecGetSize(x, &rows); CHKERRQ(ierr);
  ierr = CreateVecGA(pgroup, lrows, rows, ga); CHKERRQ(ierr);
  
  ierr = VecGetOwnershipRange(x, &lo, &hi); CHKERRQ(ierr);

  ierr = VecGetArray(x, &v); CHKERRQ(ierr);
  NGA_Put(*ga, &lo, &hi, v, &ld);
  ierr = VecRestoreArray(x, &v); CHKERRQ(ierr);
    
  return ierr;
}
  


// -------------------------------------------------------------
// MatSetValuesDenseGA
// -------------------------------------------------------------
static 
PetscErrorCode
MatSetValuesDenseGA(Mat mat, 
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
// MatGetValuesDenseGA
// -------------------------------------------------------------
static 
PetscErrorCode
MatGetValuesDenseGA(Mat mat, 
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
// MatAssemmblyBeginDenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatAssemmblyBeginDenseGA(Mat mat, MatAssemblyType type)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// MatAssemmblyEndDenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatAssemmblyEndDenseGA(Mat mat, MatAssemblyType type)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);
  GA_Pgroup_sync(ctx->gaGroup);
  // GA_Print(ctx->ga);
  return ierr;
}

// -------------------------------------------------------------
// MatMultDenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatMultDenseGA(Mat mat, Vec x, Vec y)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(mat, &ctx); CHKERRQ(ierr);
  

  return ierr;
}

// -------------------------------------------------------------
// MatViewDenseGA
// -------------------------------------------------------------


// -------------------------------------------------------------
// MatDestroyDenseGA
// -------------------------------------------------------------
static
PetscErrorCode
MatDestroyDenseGA(Mat A)
{
  PetscErrorCode ierr = 0;
  struct MatGACtx *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  GA_Destroy(ctx->ga);
  GA_Pgroup_destroy(ctx->gaGroup);
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

  ierr = MatShellSetOperation(*A, MATOP_SET_VALUES, (void(*)(void))MatSetValuesDenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*A, MATOP_GET_VALUES, (void(*)(void))MatGetValuesDenseGA); CHKERRQ(ierr);

  ierr = MatShellSetOperation(*A, MATOP_ASSEMBLY_BEGIN, (void(*)(void))MatAssemmblyBeginDenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*A, MATOP_ASSEMBLY_END, (void(*)(void))MatAssemmblyEndDenseGA); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*A, MATOP_DESTROY,  (void(*)(void))MatDestroyDenseGA); CHKERRQ(ierr);
  
  
  return ierr;
}
