#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ga.h>
//#include <ma.h>
#include <float.h>
#include <time.h>
#include <macdecls.h>
#include <cblas.h>
#include <clapack.h>
//#include <assert.h>
#include <mpi.h>
#define FMT "%lg"

#define CHECK_HANDLE(h)\
    {\
        if ((h)==0) {\
            fprintf(stderr,"found invalid array handle on line %d\n",__LINE__);\
            GA_Error("fatal error",911);\
        }\
    }


#ifdef _DEBUG_ME_
#define CHECK_INDEX(i,j,lo,hi)\
	{\
		if ((i)==-1  || (i)<(lo)[0] || (i)>(hi)[0]) {\
		    fprintf(stderr,"error in indexing i= %d  lo= %d hi = %d  line = %d\n",(i),(lo)[0],(hi)[0],__LINE__);\
		    GA_Error("fatal error",911);\
		}\
		if ((j)==-1  || (j)<(lo)[1] || (j)>(hi)[1]) {\
		    fprintf(stderr,"error in indexing j= %d  lo= %d hi = %d  line = %d\n",(j),(lo)[1],(hi)[1],__LINE__);\
		    GA_Error("fatal error",911);\
		}\
    }
#else
#define CHECK_INDEX(i,j,lo,hi)
#endif

void* replace_malloc ( size_t bytes, int align, char *name )
{
    return malloc ( bytes );
}

void replace_free ( void *ptr )
{
    free ( ptr );
}

void replace_ma()
{
    GA_Register_stack_memory ( replace_malloc, replace_free );
}


static void ga_unit_dmatrix ( int ga_a, double * buffer, int rank )
{
    int64_t lds, i, ixs;
    int64_t lo[2], hi[2];
    size_t asize;
    NGA_Distribution64 ( ga_a, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    asize = lds;
    asize *= hi[0] - lo[0] + 1;
    memset ( buffer, 0, sizeof ( double ) *asize );
    for ( i = lo[0]; i <= hi[0]; ++i ) {
        if ( i < lo[1] || i > hi[1] ) continue;
        CHECK_INDEX ( i, i, lo, hi );
        ixs = ( i - lo[0] ) * lds + i - lo[1];
        buffer[ixs] = 1.;
    }
    NGA_Put64 ( ga_a, lo, hi, buffer, &lds );
}

static inline void * Malloc ( size_t nsz )
{
    void * ptr = malloc ( nsz );
    if ( ptr ) return ptr;
    fprintf ( stderr, "out of memory requested size = %lu\n", nsz );
    GA_Error ( "out of memory", 911 );
    return 0x0;
}

struct Stopwatch_t {
    double acc;
    struct timespec ts, tf;
};

typedef struct Stopwatch_t Stopwatch_t;

static void stopwatch_start ( Stopwatch_t *s )
{
    clock_gettime ( CLOCK_MONOTONIC, & ( s->ts ) );
}

static void stopwatch_stop ( Stopwatch_t *s )
{
    clock_gettime ( CLOCK_MONOTONIC, & ( s->tf ) );
    s->acc += 1.e-9 * ( ( s->tf ).tv_nsec - ( s->ts ).tv_nsec ) + ( ( s->tf ).tv_sec - ( s->ts ).tv_sec );
}

static void stopwatch_clear ( Stopwatch_t *s )
{
    s->acc = 0;
};


static inline double normRand()
{
    double x, y, r;
    const double cut = 1. - DBL_EPSILON;
    do {
        x = ( ( double )  rand() ) / RAND_MAX;
        y = ( ( double )  rand() ) / RAND_MAX;
        y = y + y - 1.;
        r = x * x + y * y;
        if ( r <= cut ) break;
    } while ( 1 );
    x = -log ( r ) / r / log ( 2. );
    x = x + x;
    return ( x * y );
}

inline static void llt_solver ( int rank, 
                                int ga_qmat, double * qmat , double *qinv, int64_t nE )
{
    int err, cid, cnp, zid, pid, kproc;
    int64_t i, lo[2], hi[2], lds, lox[2], hix[2];
    double *qptr;
    cid = GA_Cluster_nodeid();
    cnp = GA_Cluster_nprocs ( cid );
    zid = GA_Cluster_procid ( cid, 0 );
    if ( zid==rank ) {
        lo[0] = lo[1] = 0;
        hi[0] = hi[1] = nE - 1;
        NGA_Get64 ( ga_qmat, lo, hi, qmat, &nE );
        err = clapack_dposv ( CblasRowMajor, CblasUpper, nE, nE, qmat, nE, qinv, nE );
        if ( err ) {
            fprintf (stderr, "matrix inversion failed with error code = %d\n", err );
            GA_Error("failed",911);
        }
        for ( kproc = 0; kproc < cnp; ++kproc ) {
            pid = GA_Cluster_procid ( cid, kproc );
            NGA_Distribution64 ( ga_qmat, pid, lo, hi );
            lds = hi[1] - lo[1] + 1;
            lox[1] = lo[1];
            hix[1] = hi[1];
            for ( i = lo[0]; i <= hi[0]; ++i ) {
                lox[0] = i;
                hix[0] = i;
                qptr =  qinv + ( i - lo[0] ) * lds + lo[1];
                NGA_Put64 ( ga_qmat, lox, hix, qptr , &lds );
            }
        }
    }
    GA_Sync();
}


int ga_block_cyclic_create(int dims[],const char *name,int TYPE)
{
    int err;
    int ga_h=GA_Create_handle();    
    int tile_size[2];
    tile_size[0]=1000;
    tile_size[1]=1000;
    GA_Set_data(ga_h,2,dims,TYPE);
    GA_Set_array_name(ga_h,name);
    GA_Set_block_cyclic_proc_grid(ga_h,pdims,psizes);
    status=GA_Allocate(ga_h);
    if (ga_h && status) return ga_h;
    GA_Error("could not allocate block cylic ga_a",0);
}


int main ( int argc, char **argv )
{
    register int64_t i, j, k, kstep, istep, indx, ix, iz;
    int rank, nproc;
    int64_t nB, nG, nZ, nX, nE, ksize, nsteps, err;
    int64_t chunk[2], dims[2], lo[2], hi[2], lds;
    int64_t xlo[2], xhi[2], hhi[2], hlo[2], vxhi[2], vxlo[2], vhhi[2], vhlo[2], klo[2], khi[2];
    int64_t ldim_x, ldim_h, ldim_vh, ldim_vx;
    int64_t asize, xsize, hsize, vhsize, vxsize;
    double *invHmac, *invXdp, *Pm, *D, *Emag;
    int ga_xpri, ga_xpri_a, ga_recv, ga_recv_f, ga_recv_0, ga_vth,
        ga_merr, ga_hmeas, ga_hmeas_a, ga_ph_t, ga_qmat, ga_zmat, ga_wmat, ga_kmat,
        ga_kalman, ga_xpost, ga_outvh, ga_outvx;
    int local_pid,local_cid,local_np;
    const int y_bus = 1;
    double *xp;
    double *h1;
    double *h2;
    double *dxdt1;
    double *dxdt2;
    const double *delta;
    const double *omega;
    double *d_omega_dt;
    double *d_delta_dt;
    double *xf;
    double *qbuff, *qinv, *wbuff, *qubuff;
    double tfault, tclear, vmag, theta, tr, ti, s, ka, c1, c2;
    double *xbuff, *xfwd, *buff1, *buff2;
    DoubleComplex * vbuff;
    double *zbuff, *hbuff, * zbuff2, *kbuff;
    int bsize = 1 << 31;
    FILE *fp;
    Stopwatch_t timer;
    Stopwatch_t itimer;
    Stopwatch_t oloop_timer;
    Stopwatch_t iloop_timer;
    Stopwatch_t qinv_timer;
    Stopwatch_t mm_timer;
    Stopwatch_t misc_timer;
    Stopwatch_t mio2_timer;
/********** start execution here! */    
    MPI_Init ( &argc, &argv );
    GA_Initialize();
    replace_ma();
    
    rank = GA_Nodeid();
    nproc = GA_Nnodes();
    local_cid=GA_Cluster_nodeid();
    local_np =GA_Cluster_nprocs(local_cid);
    for (i=0;i<local_np;++i) {
        local_pid=GA_Cluster_procid(local_cid,i);
        if (local_pid==rank) {
            local_pid=i;
            break;
        }
    }
    stopwatch_clear ( &timer );
    stopwatch_clear ( &itimer );
    stopwatch_clear ( &iloop_timer );
    stopwatch_clear ( &oloop_timer );
    stopwatch_clear ( &qinv_timer );
    stopwatch_clear ( &mm_timer );
    stopwatch_clear ( &misc_timer );
    stopwatch_clear ( &mio2_timer );
    stopwatch_start ( &timer );
    const double std_dev = 0.0000001;
    const double std_dev_2 = 0.05;
    const double omega0 = 1;
    const double omegaB = 2.*M_PI * 60.0;
    const double time_step = 0.01;
    const double time_step_internal = 0.01;
    const double rfact = 1 / std_dev_2;
    const int64_t ninner = ( int ) rint ( time_step / time_step_internal );
    double mag, ang;
    DoubleComplex z_one, z_zero;
    z_one.real = 1.;
    z_one.imag = 0.;
    z_zero.real = 0.;
    z_zero.imag = 0.;
    fp = fopen ( "data1.txt", "r" );
    if ( feof ( fp ) || !fp ) {
        fprintf ( stderr, "could not open input file for reading\n" );
        GA_Error ( "fatal error", 911 );
    }
    fscanf ( fp, "%ld%lg%lg", &nsteps, &tfault, &tclear );
    fscanf ( fp, "%ld%ld%ld%ld%ld%ld", &nB, &nG, &nZ, &nX, &nE, &ksize );
    if ( !rank ) {
        fprintf ( stderr, "nGen = %ld nBus = %ld  #Ensemble = %ld nX= %ld nZ= %ld\n",
                  nG, nB, nE, nX, nZ );
        fprintf ( stderr, "nsteps = %ld\n", nsteps );
        fprintf ( stderr, "starting\n" );
        fprintf ( stderr, "k size = %ld\n", ksize );
        fprintf ( stderr, "tfault = %lg\n", tfault );
        fprintf ( stderr, "tclear = %lg\n", tclear );
    }
    double oonE = 1 / ( ( double ) nE );
    double nfact = 1 / ( ( double ) ( nE - 1 ) );
    double rnfact = rfact * nfact;
    const int64_t kfault = ( int64_t ) rint ( tfault / time_step );
    const int64_t kclear = ( int64_t ) rint ( tclear / time_step );
    /*****************************************************************
     * create global arrays for data
     ****************************************************************/
    dims[0] = nB;
    dims[1] = nE;
    chunk[0] = -1;
    chunk[1] = nE;
    ga_outvh = NGA_Create64 ( C_DCPL, 2, dims, "ga_outvh", chunk );
    CHECK_HANDLE ( ga_outvh );
    dims[0] = nG;
    dims[1] = nE;
    chunk[0] = -1;
    chunk[1] = nE;
    ga_vth = NGA_Create64 ( C_DCPL, 2, dims, "ga_vth", chunk );
    ga_outvx = GA_Duplicate ( ga_vth, "ga_outvx" );
    CHECK_HANDLE ( ga_vth );
    CHECK_HANDLE ( ga_outvx );
    dims[0] = nB;
    dims[1] = nG;
    chunk[0] = -1;
    chunk[1] = -1;
    ga_recv_f = NGA_Create64 ( C_DCPL, 2, dims, "ga_recv_f", chunk );
    ga_recv_0 = GA_Duplicate ( ga_recv_f, "ga_recv_0" );
    CHECK_HANDLE ( ga_recv_f );
    CHECK_HANDLE ( ga_recv_0 );
    int64_t *map = ( int64_t* ) Malloc ( ( nproc + 1 ) * sizeof ( int64_t ) );
    int64_t *nblock = ( int64_t* ) Malloc ( 2 * sizeof ( int64_t ) );
    nblock[0] = nproc;
    nblock[1] = 1;
    for ( j = 0; j < nproc; ++j ) {
        NGA_Distribution64 ( ga_vth, j, vxlo, vxhi );
        map[j] = vxlo[0] + vxlo[0];
    }
    map[nproc] = 0;
    dims[0] = nX;
    dims[1] = nE;
    ga_xpri = NGA_Create_irreg64 ( C_DBL, 2, dims, "ga_xpri", nblock, map );
    CHECK_HANDLE ( ga_xpri );
    ga_xpri_a = GA_Duplicate ( ga_xpri, "ga_xpri_a" );
    CHECK_HANDLE ( ga_xpri_a );
    GA_Sync();
    nblock[0] = nproc;
    nblock[1] = 1;
    for ( j = 0; j < nproc; ++j ) {
        NGA_Distribution64 ( ga_outvh, j, vhlo, vhhi );
        map[j] = vhlo[0] + vhlo[0];
    }
    map[nproc] = 0;
    dims[0] = nZ;
    dims[1] = nE;
    ga_hmeas = NGA_Create_irreg64 ( C_DBL, 2, dims, "ga_hmeas", nblock, map );
    CHECK_HANDLE ( ga_hmeas );
    ga_hmeas_a = NGA_Create_irreg64 ( C_DBL, 2, dims, "ga_hmeas_a", nblock, map );
    ga_merr = NGA_Create_irreg64 ( C_DBL, 2, dims, "ga_merr", nblock, map );
    CHECK_HANDLE ( ga_hmeas_a );
    CHECK_HANDLE ( ga_merr );
    chunk[0] = -1;
    chunk[1] = -1;
    dims[0] = nX;
    dims[1] = nZ;
    ga_kmat = NGA_Create64 ( C_DBL, 2, dims, "ga_kmat", chunk );
    dims[0] = nE;
    dims[1] = nE;
    ga_qmat = NGA_Create64 ( C_DBL, 2, dims, "ga_qmat", chunk );
    dims[0] = nZ;
    dims[1] = nE;
    ga_zmat = NGA_Create64 ( C_DBL, 2, dims, "ga_zmat", chunk );
    dims[0] = nZ;
    dims[1] = nZ;
    ga_wmat = NGA_Create64 ( C_DBL, 2, dims, "ga_wmat", chunk );
    dims[0] = nX;
    dims[1] = nZ;
    ga_ph_t = NGA_Create64 ( C_DBL, 2, dims, "ga_ph_t", chunk );
    CHECK_HANDLE ( ga_kmat );
    CHECK_HANDLE ( ga_qmat );
    CHECK_HANDLE ( ga_zmat );
    CHECK_HANDLE ( ga_wmat );
    CHECK_HANDLE ( ga_ph_t );
    GA_Sync();
    nblock[0] = 1;
    nblock[1] = nproc;
    map[0] = 0;
    for ( j = 0; j < nproc; ++j ) {
        NGA_Distribution64 ( ga_hmeas, j, hlo, hhi );
        map[j+1] = hlo[0];
    }
    dims[0] = nsteps;
    dims[1] = nZ;
    ga_kalman = NGA_Create_irreg64 ( C_DBL, 2, dims, "ga_kalman", nblock, map );
    CHECK_HANDLE ( ga_kalman );
    GA_Sync();
    nblock[0] = 1;
    nblock[1] = nproc;
    map[0] = 0;
    for ( j = 0; j < nproc; ++j ) {
        NGA_Distribution64 ( ga_xpri, j, xlo, xhi );
        map[j+1] = xlo[0];
    }
    dims[0] = nsteps;
    dims[1] = nX;
    ga_xpost = NGA_Create_irreg64 ( C_DBL, 2, dims, "ga_xpost", nblock, map );
    CHECK_HANDLE ( ga_xpost );
#ifdef _DEBUG_ME_
    /**************
     * write out distributions for check
     *************/
    NGA_Distribution64 ( ga_outvh, rank, lo, hi );
    fprintf ( stderr, "ga_outvh %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_outvx, rank, lo, hi );
    fprintf ( stderr, "ga_outvx %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_xpri, rank, lo, hi );
    fprintf ( stderr, "ga_xpri %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_hmeas, rank, lo, hi );
    fprintf ( stderr, "ga_hmeas %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_kalman, rank, lo, hi );
    fprintf ( stderr, "ga_kalman %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_xpost, rank, lo, hi );
    fprintf ( stderr, "ga_xpost %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_qmat, rank, lo, hi );
    fprintf ( stderr, "ga_qmat %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_recv_0, rank, lo, hi );
    fprintf ( stderr, "ga_recv_0 %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    NGA_Distribution64 ( ga_recv_f, rank, lo, hi );
    fprintf ( stderr, "ga_recv_f %ld (%ld,%ld)(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
    fprintf ( stderr, "ga dict\n" );
    fprintf ( stderr, " ga_outvh %d\n", ga_outvh );
    fprintf ( stderr, " ga_outvx %d\n", ga_outvx );
    fprintf ( stderr, " ga_vth %d\n", ga_vth );
    fprintf ( stderr, " ga_xpri %d\n", ga_xpri );
    fprintf ( stderr, " ga_xpri_a %d\n", ga_xpri_a );
    fprintf ( stderr, " ga_hmeas  %d\n", ga_hmeas );
    fprintf ( stderr, " ga_hmeas_a %d\n", ga_hmeas_a );
    fprintf ( stderr, " ga_merr %d\n", ga_merr );
    fprintf ( stderr, " ga_xpost %d\n", ga_xpost );
    fprintf ( stderr, " ga_kalman %d\n", ga_kalman );
    fprintf ( stderr, " ga_qmat %d\n", ga_qmat );
    fprintf ( stderr, " ga_ph_t %d\n", ga_ph_t );
    fprintf ( stderr, " ga_wmat %d\n", ga_wmat );
    fprintf ( stderr, " ga_zmat %d\n", ga_zmat );
    fprintf ( stderr, " ga_kmat %d\n", ga_kmat );
#endif
    /*******************
     * create buffers for local data
     ******************/
    NGA_Distribution64 ( ga_xpri, rank, xlo, xhi );
    ldim_x = xhi[1] - xlo[1] + 1L;
    asize = ldim_x;
    asize *= ( xhi[0] - xlo[0] + 1L );
    xbuff = ( double* ) Malloc ( sizeof ( double ) * asize );
    xfwd = ( double* ) Malloc ( sizeof ( double ) * asize );
    buff1 = ( double* ) Malloc ( sizeof ( double ) * asize );
    buff2 = ( double* ) Malloc ( sizeof ( double ) * asize );
    xsize = asize;
    NGA_Distribution64 ( ga_recv_0, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    asize = lds;
    asize = asize * ( hi[0] - lo[0] + 1 );
    NGA_Distribution64 ( ga_outvh, rank, vhlo, vhhi );
    ldim_vh = vhhi[1] - vhlo[1] + 1;
    vhsize = ldim_vh;
    vhsize *= ( vhhi[0] - vhlo[0] + 1 );
    asize = ( asize > vhsize ) ? asize : vhsize;
    NGA_Distribution64 ( ga_vth, rank, vxlo, vxhi );
    ldim_vx = vxhi[1] - vxlo[1] + 1;
    vxsize = ldim_vx;
    vxsize *= ( vxhi[0] - vxlo[0] + 1 );
    asize = ( asize > vxsize ) ? asize : vxsize;
    vbuff = ( DoubleComplex* ) Malloc ( sizeof ( DoubleComplex ) * asize );
    NGA_Distribution64 ( ga_hmeas, rank, hlo, hhi );
    ldim_h = hhi[1] - hlo[1] + 1;
    hsize = ldim_h;
    hsize = hsize * ( hhi[0] - hlo[0] + 1 );
    zbuff = ( double* ) Malloc ( sizeof ( double ) * hsize );
    zbuff2 = ( double* ) Malloc ( sizeof ( double ) * hsize );
    hbuff = ( double* ) Malloc ( sizeof ( double ) * hsize );
    NGA_Distribution64 ( ga_wmat, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    asize = lds;
    asize *= hi[0] - lo[0] + 1;
    wbuff = ( double* ) Malloc ( sizeof ( double ) * asize );
    NGA_Distribution64 ( ga_qmat, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    asize = lds;
    asize *= hi[0] - lo[0] + 1;
    qubuff = ( double* ) Malloc ( sizeof ( double ) * asize );
    if ( ( xlo[0] != ( 2 * vxlo[0] ) ) || ( xlo[1] != vxlo[1] ) || ( xhi[1] != vxhi[1] ) ) {
        fprintf ( stderr, "v - x distribution is messed up!\n" );
        fprintf ( stderr, "rank = %d  xdims = (%ld,%ld),(%ld,%ld)\n", rank, xlo[0], xhi[0], xlo[1], xhi[1] );
        fprintf ( stderr, "rank = %d  vdims = (%ld,%ld),(%ld,%ld)\n", rank, vxlo[0], vxhi[0], vxlo[1], vxhi[1] );
        fprintf ( stderr, "rank = %d  =xdims = (%ld,%ld),(%ld,%ld)\n", rank, 2 * vxlo[0], 2 * vxhi[0], vxlo[1], vxhi[1] );
        GA_Error ( "fatal error", 911 );
    }
    if ( ( hlo[0] != ( 2 * vhlo[0] ) ) || ( hlo[1] != vhlo[1] ) || ( hhi[1] != vhhi[1] ) ) {
        fprintf ( stderr, "v - z distribution is messed up!\n" );
        fprintf ( stderr, "rank = %d  zdims = (%ld,%ld),(%ld,%ld)\n", rank, hlo[0], hhi[0], hlo[1], hhi[1] );
        fprintf ( stderr, "rank = %d  hvdims = (%ld,%ld),(%ld,%ld)\n", rank, vhlo[0], vhhi[0], vhlo[1], vhhi[1] );
        fprintf ( stderr, "rank = %d  =zdims = (%ld,%ld),(%ld,%ld)\n", rank, 2 * vhlo[0], 2 * vhhi[0], vhlo[1], vhhi[1] );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_kalman, rank, klo, khi );
    if ( ( hlo[0] !=  klo[1] ) || ( hhi[0] != khi[1] ) ) {
        fprintf ( stderr, "kalman - z distribution is messed up!\n" );
        fprintf ( stderr, "rank = %d  zdims = (%ld,%ld),(%ld,%ld)\n", rank, hlo[0], hhi[0], hlo[1], hhi[1] );
        fprintf ( stderr, "rank = %d  kdims = (%ld,%ld),(%ld,%ld)\n", rank, klo[0], khi[0], klo[1], khi[1] );
        fprintf ( stderr, "rank = %d  =zdims = (%ld,%ld),(%ld,%ld)\n", rank, 0L , nsteps, hlo[0], hhi[0] );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_xpost, rank, lo, hi );
    if ( ( xlo[0] !=  lo[1] ) || ( xhi[0] != hi[1] ) ) {
        fprintf ( stderr, "xpost - x distribution is messed up!\n" );
        fprintf ( stderr, "rank = %d  xdims = (%ld,%ld),(%ld,%ld)\n", rank, xlo[0], xhi[0], xlo[1], xhi[1] );
        fprintf ( stderr, "rank = %d  pdims = (%ld,%ld),(%ld,%ld)\n", rank, lo[0], hi[0], lo[1], hi[1] );
        fprintf ( stderr, "rank = %d  =dims = (%ld,%ld),(%ld,%ld)\n", rank, 0L , nsteps, xlo[0], xhi[0] );
        GA_Error ( "fatal error", 911 );
    }
    /***************************
     * initialize parameter and buffer arrays
     **************************/
    Pm = ( double* ) Malloc ( sizeof ( double ) * nG );
    D = ( double* ) Malloc ( sizeof ( double ) * nG );
    Emag = ( double* ) Malloc ( sizeof ( double ) * nG );
    invHmac = ( double* ) Malloc ( sizeof ( double ) * nG );
    invXdp = ( double* ) Malloc ( sizeof ( double ) * nG );
//    if (!local_pid) {
        qinv  = ( double* ) Malloc ( sizeof ( double ) * nE * nE );
        qbuff = ( double* ) Malloc ( sizeof ( double ) * nE * nE );
//    }
    /*****************************
     * some checks for correctness
     *
     ****************************/
    NGA_Distribution64 ( ga_outvh, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "outvh distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_outvx, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "outvx distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_xpri, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "xpri distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_xpri_a, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "xpri_a distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_hmeas, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "hmeas distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_hmeas_a, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "hmeas_a distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    NGA_Distribution64 ( ga_merr, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    if ( lds != nE ) {
        fprintf ( stderr, "merr distribution wrong!\n" );
        GA_Error ( "fatal error", 911 );
    }
    /*****************************
     * read in data
     ****************************/
    if ( !rank ) fprintf ( stderr, "starting read of file!\n" );
    NGA_Distribution64 ( ga_recv_0, rank, lo, hi );
    lds = hi[1] - lo[1] + 1;
    for ( j = 0; j < nG; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( i < lo[0] ||  i > hi[0] ) continue;
            if ( j < lo[1] ||  j > hi[1] ) continue;
            CHECK_INDEX ( i, j, lo, hi );
            indx = ( i - lo[0] ) * lds + j - lo[1];
            vbuff[indx].real = s;
        }
    }
    for ( j = 0; j < nG; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( i < lo[0] ||  i > hi[0] ) continue;
            if ( j < lo[1] ||  j > hi[1] ) continue;
            CHECK_INDEX ( i, j, lo, hi );
            indx = ( i - lo[0] ) * lds + j - lo[1];
            vbuff[indx].imag = s;
        }
    }
    k = 0;
    for ( i = lo[0]; i <= hi[0]; ++i ) {
        for ( j = lo[1]; j <= hi[1]; ++j ) {
            double tr = vbuff[k].real;
            double ti = vbuff[k].imag;
            vbuff[k].real = tr * cos ( ti );
            vbuff[k].imag = tr * sin ( ti );
            ++k;
        }
    }
    NGA_Put64 ( ga_recv_0, lo, hi, vbuff, &lds );
    GA_Sync();
    for ( j = 0; j < nG; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( i < lo[0] ||  i > hi[0] ) continue;
            if ( j < lo[1] ||  j > hi[1] ) continue;
            CHECK_INDEX ( i, j, lo, hi );
            indx = ( i - lo[0] ) * lds + j - lo[1];
            vbuff[indx].real = s;
        }
    }
    for ( j = 0; j < nG; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( i < lo[0] ||  i > hi[0] ) continue;
            if ( j < lo[1] ||  j > hi[1] ) continue;
            CHECK_INDEX ( i, j, lo, hi );
            indx = ( i - lo[0] ) * lds + j - lo[1];
            vbuff[indx].imag = s;
        }
    }
    k = 0;
    for ( i = lo[0]; i <= hi[0]; ++i ) {
        for ( j = lo[1]; j <= hi[1]; ++j ) {
            CHECK_INDEX ( i, j, lo, hi );
            tr = vbuff[k].real;
            ti = vbuff[k].imag;
            vbuff[k].real = tr * cos ( ti );
            vbuff[k].imag = tr * sin ( ti );
            ++k;
        }
    }
    NGA_Put64 ( ga_recv_f, lo, hi, vbuff, &lds );
    GA_Sync();
    NGA_Distribution64 ( ga_xpri, rank, xlo, xhi );
    ldim_x = xhi[1] - xlo[1] + 1;
    for ( j = 0; j < nE; ++j ) {
        for ( i = 0; i < nX; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( i < xlo[0] ||  i > xhi[0] ) continue;
            if ( j < xlo[1] ||  j > xhi[1] ) continue;
            CHECK_INDEX ( i, j, xlo, xhi );
            indx = ( i - xlo[0] ) * ldim_x + j - xlo[1];
            xbuff[indx] = s;
        }
    }
    NGA_Put64 ( ga_xpri, xlo, xhi, xbuff, &ldim_x );
    GA_Sync();
    for ( i = 0; i < nG; ++i ) {
        fscanf ( fp, FMT, &s );
        invHmac[i] = 1 / s;
    }
    for ( i = 0; i < nG; ++i ) {
        fscanf ( fp, FMT, Pm + i );
    }
    for ( i = 0; i < nG; ++i ) {
        fscanf ( fp, FMT, D + i );
    }
    for ( i = 0; i < nG; ++i ) {
        fscanf ( fp, FMT, Emag + i );
    }
    for ( i = 0; i < nG; ++i ) {
        fscanf ( fp, FMT, &s );
        invXdp[i] = 1 / s;
    }
    NGA_Distribution64 ( ga_kalman, rank, klo, khi );
    lds = khi[1] - klo[1] + 1;
    asize = lds * ( khi[0] - klo[0] + 1 );
    kbuff = ( double* ) Malloc ( sizeof ( double ) * asize );
    for ( j = 0; j < nsteps; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( j < klo[0] || j > khi[0] ) continue;
            if ( i < klo[1] || i > khi[1] ) continue;
            CHECK_INDEX ( j, i, klo, khi );
            indx = ( j - klo[0] ) * lds + i - klo[1];
            kbuff[indx] = s;
        }
    }
    for ( j = nsteps; j < ksize; ++j ) {
        for ( i = 0; i < nB; ++i ) fscanf ( fp, FMT, &s );
    }
    for ( j = 0; j < nsteps; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            if ( j < klo[0] || j > khi[0] ) continue;
            if ( ( i + nB ) < klo[1] || ( i + nB ) > khi[1] ) continue;
            CHECK_INDEX ( j, ( i + nB ), klo, khi );
            indx = ( j - klo[0] ) * lds + i + nB - klo[1];
            kbuff[indx] = s;
        }
    }
    fclose ( fp );
    NGA_Put64 ( ga_kalman, klo, khi, kbuff, &lds );
    GA_Sync();
    if ( !rank ) fprintf ( stderr, "ended reading of file!\n" );
///////// find averages for initial state
    k = 0;
    for ( i = xlo[0]; i <= xhi[0]; ++i ) {
        ix = i - xlo[0];
        xp = xbuff + ix * ldim_x;
        s = 0;
        for ( j = 0; j < nE; ++j ) {
            CHECK_INDEX ( i, j, xlo, xhi );
            s += xp[j];
        }
        s *= oonE;
        buff1[ix] = s;
        ++k;
    }
    lo[0] = 0;
    hi[0] = 0;
    lo[1] = xlo[0];
    hi[1] = xhi[0];
    lds = hi[1] - lo[1] + 1;
    NGA_Put64 ( ga_xpost, lo, hi, buff1, &lds );
    GA_Sync();
    if ( !rank ) {
        fprintf ( stderr, "starting run\n" );
        fprintf ( stderr, "kfault = %ld kclear = %ld nsteps = %ld\n", kfault, kclear, nsteps );
        fprintf ( stderr, "# inner steps = %ld\n", ninner );
        fprintf ( stderr, "# nodes       = %d\n", nproc );
    }
    stopwatch_start ( &itimer );
    NGA_Distribution64 ( ga_kalman, rank, klo, khi );
    for ( kstep = 1; kstep < nsteps; ++kstep ) {
        stopwatch_start ( &iloop_timer );
        if ( ( kstep >= kfault ) && ( kstep < kclear ) ) {
            if ( y_bus ) {
                ga_recv = ga_recv_f;
            } else {
                ga_recv = ga_recv_0;
            }
        }  else {
            ga_recv = ga_recv_0;
        }
        NGA_Distribution64 ( ga_kalman, rank, klo, khi );
        lo[0] = kstep;
        hi[0] = kstep;
        lo[1] = klo[1];
        hi[1] = khi[1];
        lds = khi[1] - klo[1] + 1;
        NGA_Get64 ( ga_kalman, lo, hi, kbuff, &lds );
        NGA_Distribution64 ( ga_xpri, rank, xlo, xhi );
        ldim_x = xhi[1] - xlo[1] + 1;
        NGA_Get64 ( ga_xpri, xlo, xhi, xbuff, &nE  );
        GA_Sync();
        for ( istep = 0; istep < ninner; ++istep ) {
            NGA_Distribution64 ( ga_vth, rank, vxlo, vxhi );
            ldim_vx = vxhi[1] - vxlo[1] + 1;
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                xp = xbuff + ix * nE;
                for ( j = 0; j < nE; ++j, ++k ) {
                    tr = Emag[i];
                    ti = xp[j];
                    vbuff[k].real = tr * cos ( ti );
                    vbuff[k].imag = tr * sin ( ti );
                }
            }
            NGA_Put64 ( ga_vth, vxlo, vxhi, vbuff, &nE );
            GA_Sync();
            GA_Zgemm64 ( 'n', 'n', nG, nE, nG, z_one, ga_recv, ga_vth, z_zero, ga_outvx );
            NGA_Get64 ( ga_outvx, vxlo, vxhi, vbuff, &nE );
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                delta = xbuff + ix * nE;
                omega = xbuff + ix * nE + nE;
                d_delta_dt = buff1 + ix * nE;
                d_omega_dt = buff1 + ix * nE + nE;
                for ( j = 0; j < nE; ++j, ++k ) {
                    tr = vbuff[k].real;
                    ti = vbuff[k].imag;
                    theta = atan2 ( ti, tr );
                    vmag = sqrt ( tr * tr + ti * ti );
                    c1 = omega[j] - omega0;
                    d_delta_dt[j] = omegaB * c1;
                    c2 = Pm[i] - D[i] * c1 - invXdp[i] * vmag * Emag[i] * sin ( delta[j] - theta );
                    d_omega_dt[j] = omega0 * invHmac[i] * 0.5 * c2;
                }
            }
            for ( i = xlo[0]; i <= xhi[0]; ++i ) {
                ix = i - xlo[0];
                xp = xbuff + ix * nE;
                xf = xfwd + ix * nE;
                dxdt1 = buff1 + ix * nE;
                for ( j = 0; j < nE; ++j ) {
                    xf[j] = xp[j] + ( dxdt1[j] * time_step_internal );
                }
            }
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                xp = xfwd + ix * nE;
                for ( j = 0; j < nE; ++j, ++k ) {
                    tr = Emag[i];
                    ti = xp[j];
                    vbuff[k].real = tr * cos ( ti );
                    vbuff[k].imag = tr * sin ( ti );
                }
            }
            NGA_Put64 ( ga_vth, vxlo, vxhi, vbuff, &nE );
            GA_Sync();
            GA_Zgemm64 ( 'n', 'n', nG, nE, nG, z_one, ga_recv, ga_vth, z_zero, ga_outvx );
            NGA_Get64 ( ga_outvx, vxlo, vxhi, vbuff, &nE );
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                delta = xfwd + ix * nE;
                omega = xfwd + ix * nE + nE;
                d_delta_dt = buff2 + ix * nE;
                d_omega_dt = buff2 + ix * nE + nE;
                for ( j = 0; j < nE; ++j, ++k ) {
                    tr = vbuff[k].real;
                    ti = vbuff[k].imag;
                    theta = atan2 ( ti, tr );
                    vmag = sqrt ( tr * tr + ti * ti );
                    c1 = omega[j] - omega0;
                    d_delta_dt[j] = omegaB * c1;
                    c2 = Pm[i] - D[i] * c1 - invXdp[i] * vmag * Emag[i] * sin ( delta[j] - theta );
                    d_omega_dt[j] = omega0 * invHmac[i] * 0.5 * c2;
                }
            }
            for ( i = xlo[0]; i <= xhi[0]; ++i ) {
                ix = i - xlo[0];
                xp = xbuff + ix * nE;
                dxdt1 = buff1 + ix * nE;
                dxdt2 = buff2 + ix * nE;
                for ( j = 0; j < nE; ++j ) {
                    xp[j] += ( dxdt1[j] + dxdt2[j] ) * time_step_internal * 0.5;
                }
            }
        }
        GA_Sync();
        stopwatch_stop ( &iloop_timer );
        stopwatch_start ( &oloop_timer );
        stopwatch_start ( &misc_timer );
        stopwatch_start ( &mio2_timer );
        NGA_Put64 ( ga_xpri, xlo, xhi, xbuff, &nE );
        for ( i = xlo[0]; i <= xhi[0]; ++i ) {
            ix = i - xlo[0];
            const double *xpr = xbuff + ix * nE;
            double *xpa = buff1 + ix * nE;
            s = 0;
            for ( j = 0; j < nE; ++j ) {
                CHECK_INDEX ( i, j, xlo, xhi );
                s += xpr[j];
            }
            s *= oonE;
            for ( j = 0; j < nE; ++j ) {
                CHECK_INDEX ( i, j, xlo, xhi );
                xpa[j] = xpr[j] - s;
            }
        }
        NGA_Put64 ( ga_xpri_a, xlo, xhi, buff1, &nE );
        k = 0;
        for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
            xp = xbuff + ( i + i - xlo[0] ) * nE;
            for ( j = 0; j < nE; ++j, ++k ) {
                CHECK_INDEX ( i, j, vxlo, vxhi );
                CHECK_INDEX ( ( i + i ), j, xlo, xhi );
                mag = Emag[i];
                ang = xp[j];
                vbuff[k].real = mag * cos ( ang );
                vbuff[k].imag = mag * sin ( ang );
            }
        }
        stopwatch_start ( &mio2_timer );
        NGA_Put64 ( ga_vth, vxlo, vxhi, vbuff, &nE );
        GA_Sync();
        GA_Zgemm64 ( 'n', 'n', nB, nE, nG, z_one, ga_recv, ga_vth, z_zero, ga_outvh );
        NGA_Distribution64 ( ga_outvh, rank, vhlo, vhhi );
        ldim_vh = vhhi[1] - vhlo[1] + 1;
        NGA_Get64 ( ga_outvh, vhlo, vhhi, vbuff, &nE );
        k = 0;
        for ( i = vhlo[0]; i <= vhhi[0]; i ++ ) {
            iz = i + i - hlo[0];
            h1 = hbuff + iz * nE;
            h2 = h1 + nE;
            for ( j = 0; j < nE; ++j, ++k ) {
                CHECK_INDEX ( i, j, vhlo, vhhi );
                CHECK_INDEX ( ( i + i ), j, hlo, hhi );
                tr = vbuff[k].real;
                ti = vbuff[k].imag;
                h1[j] = sqrt ( tr * tr + ti * ti );
                h2[j] = atan2 ( ti, tr );
            }
        }
        NGA_Put64 ( ga_hmeas, hlo, hhi, hbuff, &nE );
        for ( i = hlo[0]; i <= hhi[0]; ++i ) {
            iz = i - hlo[0];
            h1 = hbuff + iz * nE;
            h2 = zbuff + iz * nE;
            s = 0;
            for ( j = 0; j < nE; ++j ) {
                CHECK_INDEX ( i, j, hlo, hhi );
                s += h1[j];
            }
            s *= oonE;
            for ( j = 0; j < nE; ++j ) {
                CHECK_INDEX ( i, j, hlo, hhi );
                h2[j] = h1[j] - s;
            }
        }
        NGA_Put64 ( ga_hmeas_a, hlo, hhi, zbuff, &nE );
        for ( i = hlo[0]; i <= hhi[0]; ++i ) {
            iz = i - hlo[0];
            ka = kbuff[iz];
            h1 = hbuff + iz * nE;
            h2 = zbuff2 + iz * nE;
            for ( j = 0; j < nE; ++j ) {
                CHECK_INDEX ( i, j, hlo, hhi )
                h2[j] = ka * ( 1 + 0.0 * normRand() ) - h1[j];
            }
        }
        NGA_Put64 ( ga_merr, hlo, hhi, zbuff2, &nE );
        stopwatch_stop ( &mio2_timer );
        GA_Sync();
        ga_unit_dmatrix ( ga_qmat, qubuff, rank );
        ga_unit_dmatrix ( ga_wmat, wbuff, rank );
        GA_Sync();
        stopwatch_stop ( &misc_timer );
        stopwatch_start ( &mm_timer );
        GA_Dgemm64 ( 'n', 't', nX, nZ, nE, nfact, ga_xpri_a, ga_hmeas_a, 0.0, ga_ph_t );
        GA_Dgemm64 ( 't', 'n', nE, nE, nZ, rnfact, ga_hmeas_a, ga_hmeas_a, 1.0, ga_qmat );
        stopwatch_stop ( &mm_timer );
        stopwatch_start ( &qinv_timer );
        /**
                GA_Llt_solve ( ga_qmat, ga_qinv );
        **/
        llt_solver ( rank, ga_qmat, qbuff , qinv, nE );
        stopwatch_stop ( &qinv_timer );
        stopwatch_start ( &mm_timer );
        GA_Dgemm64 ( 'n', 'n', nZ, nE, nE, rnfact, ga_hmeas_a, ga_qmat, 0.0, ga_zmat );
        GA_Dgemm64 ( 'n', 't', nZ, nZ, nE, -1.0, ga_zmat, ga_hmeas_a, 1.0, ga_wmat );
        GA_Dgemm64 ( 'n', 'n', nX, nZ, nZ, rfact, ga_ph_t, ga_wmat, 0.0, ga_kmat );
        GA_Dgemm64 ( 'n', 'n', nX, nE, nZ, 1.0, ga_kmat, ga_merr, 1.0, ga_xpri );
        stopwatch_stop ( &mm_timer );
        stopwatch_start ( &misc_timer );
        NGA_Get64 ( ga_xpri, xlo, xhi, xbuff, &nE );
        lo[0] = kstep;
        hi[0] = kstep;
        lo[1] = xlo[0];
        hi[1] = xhi[0];
        lds = hi[1] - lo[1] + 1;
        for ( i = xlo[0]; i <= xhi[0]; ++i ) {
            k = i - xlo[0];
            xp = xbuff + k * nE;
            s = 0;
            for ( j = 0; j < nE; ++j ) s += xp[j];
            s *= oonE;
            buff1[k] = s;
        }
        NGA_Put64 ( ga_xpost, lo, hi, buff1, &lds );
        GA_Sync();
        stopwatch_stop ( &misc_timer );
        stopwatch_stop ( &oloop_timer );
    }
    stopwatch_stop ( &itimer );
    GA_Print ( ga_xpost );
    stopwatch_stop ( &timer );
    if ( rank == 0 ) {
        fprintf ( stderr, "tot  time = %lg\n", timer.acc );
        fprintf ( stderr, "init      = %lg\n", timer.acc - itimer.acc );
        fprintf ( stderr, "loop time = %lg\n", itimer.acc );
        fprintf ( stderr, "    outer = %lg\n", oloop_timer.acc );
        fprintf ( stderr, "    inner = %lg\n", iloop_timer.acc );
        fprintf ( stderr, "inverse   = %lg\n", qinv_timer.acc );
        fprintf ( stderr, "mmult     = %lg\n", mm_timer.acc );
        fprintf ( stderr, "misc outer= %lg\n", misc_timer.acc );
        fprintf ( stderr, "mio2 outer= %lg\n", mio2_timer.acc );
    }
    GA_Sync();
    GA_Destroy ( ga_xpost );
    GA_Destroy ( ga_kalman );
    GA_Destroy ( ga_ph_t );
    GA_Destroy ( ga_wmat );
    GA_Destroy ( ga_zmat );
    GA_Destroy ( ga_qmat );
    GA_Destroy ( ga_kmat );
    GA_Destroy ( ga_recv_0 );
    GA_Destroy ( ga_recv_f );
    GA_Destroy ( ga_vth );
    GA_Destroy ( ga_outvh );
    GA_Destroy ( ga_outvx );
    GA_Destroy ( ga_merr );
    GA_Destroy ( ga_hmeas_a );
    GA_Destroy ( ga_hmeas );
    GA_Destroy ( ga_xpri_a );
    GA_Destroy ( ga_xpri );
    GA_Sync();
    fprintf ( stderr, "rank %d done x3\n", rank );
    GA_Terminate();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
