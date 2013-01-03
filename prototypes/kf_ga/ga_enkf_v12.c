#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ga.h>
#include <float.h>
#include <time.h>
#include <macdecls.h>
#include <cblas.h>
#include <clapack.h>
#include <mpi.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>
#include <sched.h>
#include <unistd.h>
#include <sys/types.h>


#define FMT "%lg"
#define ONE 1.0
#define ZERO 0.0
#define HALF 0.5

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
    void * ptr  = malloc(bytes);
    if ( ptr != NULL ) return ptr;
    fprintf ( stderr, "could not allocate %lu bytes for %s\n", (unsigned long)bytes, name );
    exit ( 1 );
}

void replace_free ( void *ptr )
{
    free ( ptr );
}

void replace_ma()
{
    GA_Register_stack_memory ( replace_malloc, replace_free );
}



size_t GAXsizeof ( int type )
{
    switch ( type ) {
    case C_DBL  :
        return ( sizeof ( double ) );
    case C_INT  :
        return ( sizeof ( int ) );
    case C_SCPL :
        return ( sizeof ( SingleComplex ) );
    case C_DCPL :
        return ( sizeof ( DoubleComplex ) );
    case C_FLOAT :
        return ( sizeof ( float ) );
    case C_LONG :
        return ( sizeof ( long ) );
    case C_LONGLONG :
        return ( sizeof ( long long ) );
    default   :
        return 0;
    }
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
        buffer[ixs] = ONE;
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

int str2int(const char *str)
{
    long x;
    char *endptr;

    errno = 0;    /* To distinguish success/failure after call */
    x = strtol(str, &endptr, 10);

    /* Check for various possible errors */

    if ((errno == ERANGE && (x == LONG_MAX || x == LONG_MIN))
            || (errno != 0 && x == 0)) {
        fprintf(stderr,"error parsing %s to an integer\n",str);
        perror("strtol");
        exit(EXIT_FAILURE);
    }

    if (endptr == str) {
        fprintf(stderr, "No digits were found in %s\n",str);
        exit(EXIT_FAILURE);
    }

    /* If we got here, strtol() successfully parsed a number */

//    printf("strtol() returned %ld\n", x);

    if (*endptr != '\0')        /* Not necessarily an error... */
        printf("Further characters after number: %s\n", endptr);

    if (x > INT_MAX || x < INT_MIN) {
        fprintf(stderr,"value %s overflows int in str2int\n",str);
        exit(EXIT_FAILURE);
    }
    return (int)x;
}

void parse_args(int argc,char **argv, int *ns_max, char **input_name)
{
    int ich;
    char ch;
    *ns_max= 0;
    *input_name = 0x0;
    while (1) {
        ich=getopt(argc,argv,"n:i:");
        if (ich==-1) break;
        ch=(int)ich;
        switch (ch) {
        case 'n':
            *ns_max = str2int(optarg);
            break;
        case 'i':
            *input_name = strdup(optarg);
            break;
        default:
            fprintf(stderr,"error: no such input option -%d\n",optopt);
            exit(EXIT_FAILURE);
        }
    }
}


int expandedGA ( int ga_in, int rank,
                 int type, int64_t dims[], char *name )
{
    int j;
    int ga_out;
    int64_t bsz,xsz,slo[2], shi[2], xlo[2], xhi[2], *map, nb[2], asize;
    int nblock[2],knt=0;

    GA_Nblock ( ga_in, nblock );
    asize = ( nblock[0] + nblock[1] +1 ) * sizeof ( int64_t );
    map = ( int64_t* ) replace_malloc ( asize, 8, "xmap" );
    bsz=(dims[0]/2)/nblock[0];
    xsz=(dims[0]/2)%nblock[0];    
    for ( j = 0; j < nblock[0]; ++j ) {
        if (j<xsz) {
            slo[0]=(bsz+1)*j;
        }else{
            slo[0]=(bsz)*j+xsz;
        }
        map[knt] = slo[0] + slo[0];
        ++knt;
    }
    bsz=dims[1]/nblock[1];
    xsz=dims[1]%nblock[1];    
    for ( j = 0; j < nblock[1]; ++j ) {
        if (j < xsz ) {
        map[knt] = (bsz + 1 ) * j; 
        }else{
        map[knt] = (bsz  * j) + xsz;
        }
        ++knt;
    }
    nb[0] = nblock[0];
    nb[1] = nblock[1];
    ga_out = NGA_Create_irreg64 ( type, 2, dims, name, nb, map );
    CHECK_HANDLE ( ga_out );
    free ( map );
    NGA_Distribution64(ga_in,rank,slo,shi);
    NGA_Distribution64(ga_out,rank,xlo,xhi);
    if ( ( xlo[0] != ( 2 * slo[0] ) ) || ( xlo[1] != slo[1] ) || ( xhi[1] != shi[1] ) ) {
        fprintf ( stderr, "v - x distribution is messed up for %s!\n",name );
        fprintf ( stderr, "rank = %d  xdims = (%ld,%ld),(%ld,%ld)\n", rank, xlo[0], xhi[0], xlo[1], xhi[1] );
        fprintf ( stderr, "rank = %d  vdims = (%ld,%ld),(%ld,%ld)\n", rank, slo[0], shi[0], slo[1], shi[1] );
        fprintf ( stderr, "rank = %d  =xdims = (%ld,%ld),(%ld,%ld)\n", rank, 2 * slo[0], 2 * shi[0], slo[1], shi[1] );
        GA_Error ( "fatal error", 911 );
    }
    return ga_out;
}

int expandedGA2 ( int ga_in, int rank,
                 int type, int64_t dims[], char *name )
{
    int j;
    int ga_out;
    int64_t bsz,xsz,slo[2], shi[2], xlo[2], xhi[2], *map, nb[2], asize;
    int nblock[2],knt=0;

    GA_Nblock ( ga_in, nblock );
    asize = ( nblock[0] + 2 ) * sizeof ( int64_t );
    map = ( int64_t* ) replace_malloc ( asize, 8, "xmap" );
    map[0]=0;
    bsz=(dims[0])/nblock[0];
    xsz=(dims[0])%nblock[0];    
    for ( j = 0; j < nblock[0]; ++j ) {
        if (j<xsz) {
            slo[0]=(bsz+1)*j;
        }else{
            slo[0]=(bsz*j)+xsz;
        }
        ++knt;
        map[knt] = slo[0];
    }
    nb[0] = 1;
    nb[1] = nblock[0];
    ga_out = NGA_Create_irreg64 ( type, 2, dims, name, nb, map );
    CHECK_HANDLE ( ga_out );
    free ( map );
    NGA_Distribution64(ga_in,rank,slo,shi);
    NGA_Distribution64(ga_out,rank,xlo,xhi);
    if ( (slo[0]!=xlo[1]) ||  (shi[0]!=xlo[1]) ) {
        fprintf ( stderr, "v - x distribution is messed up for %s!\n",name );
        fprintf ( stderr, "rank = %d  in dims = (%ld,%ld),(%ld,%ld)\n", rank, xlo[0], xhi[0], xlo[1], xhi[1] );
        fprintf ( stderr, "rank = %d out dims = (%ld,%ld),(%ld,%ld)\n", rank, slo[0], shi[0], slo[1], shi[1] );
        GA_Error ( "fatal error", 911 );
    }
    return ga_out;
}

static inline void * create_local_buffer ( int ga_in, int type, int ndims, int64_t dims[] )
{
    int k;
    int nblock[10];
    int64_t bsz, xsz;
    size_t asize;
    GA_Nblock ( ga_in, nblock );
    asize = 1;
    for ( k = 0; k < ndims; ++k ) {
        bsz = dims[0] / ( nblock[0] );
        xsz = dims[0] % nblock[0];
        if ( xsz ) asize *= bsz + 1;
        else asize *= bsz;
    }
    asize *= GAXsizeof ( type );
    return replace_malloc ( asize, 8, GA_Inquire_name ( ga_in ) );
}

int get_row_pgroup ( int ga_in, int rank, int nproc )
{
    int k, knt, nblock[2], pgrp;
    int proclist[nproc+1];
    int64_t xlo[2], xhi[2], slo[2], shi[2];
    knt = 0;
    GA_Nblock ( ga_in, nblock );
    NGA_Distribution64 ( ga_in, rank, xlo, xhi );
    for ( k = 0; k < nproc; ++k ) {
        NGA_Distribution64 ( ga_in, k, slo, shi );
        if ( xlo[0] == slo[0] && xhi[0] == shi[0] ) {
            proclist[knt] = k;
            ++knt;
        }
        if ( ( shi[0] > xlo[0] ) && ( slo[0] < xlo[0] ) ) {
            fprintf ( stderr, "overlapping nonequal data partitions!\n" );
            fprintf ( stderr, "x -> ( %ld, %ld) and (%ld,%ld)\n",
                      slo[0], shi[0], xlo[0], xhi[0] );
            exit ( 1 );
        }
        if ( ( slo[0] < xhi[0] ) && ( shi[0] > xhi[0] ) ) {
            fprintf ( stderr, "overlapping nonequal data partitions!\n" );
            fprintf ( stderr, "x -> ( %ld, %ld) and (%ld,%ld)\n",
                      slo[0], shi[0], xlo[0], xhi[0] );
            exit ( 1 );
        }
    }
    pgrp = GA_Pgroup_create ( proclist, knt );
    return pgrp;
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
    s->acc = ZERO;
};


static inline double normRand()
{
    double x, y, r;
    const double cut = ONE - DBL_EPSILON;
    do {
        x = ( ( double )  rand() ) / RAND_MAX;
        y = ( ( double )  rand() ) / RAND_MAX;
        y = y + y - ONE;
        r = x * x + y * y;
        if ( r <= cut ) break;
    } while ( 1 );
    x = -log ( r ) / r / log ( 2 );
    x = x + x;
    return ( x * y );
}


void store_xpost( double *buffer, double *xbuff, int ga_xpost, int ga_xpri,
int xrow_pgrp, int xrow_rank, int rank,
int kstep,int64_t nE)
{
    int64_t xlo[2],xhi[2],ldim_x,sdim_x,i,j,k,plo[2],phi[2];
    register double s;
    const double oonE=1.0/((double)nE);
        
    NGA_Distribution64(ga_xpri,rank,xlo,xhi);
    ldim_x = xhi[1]-xlo[1]+1;
    sdim_x = xhi[0]-xlo[0]+1;
    NGA_Get64(ga_xpri,xlo,xhi,xbuff,&ldim_x);
    memset(buffer,0,sizeof(double)*sdim_x);
    k=0;
    for (i=xlo[0];i<=xhi[0];++i) {
        s=0;
        for (j=xlo[1];j<=xhi[1];++j) {
            s+=xbuff[k];
            ++k;
        }
        buffer[(i-xlo[0])]=s*oonE;
    }
    GA_Pgroup_dgop(xrow_pgrp,buffer,sdim_x,"+");
    plo[0]=kstep;
    phi[0]=kstep;
    plo[1]=xlo[0];
    phi[1]=xhi[0];
    if (xrow_rank==0) {
        NGA_Put64(ga_xpost,plo,phi,buffer,&sdim_x);
    }
    GA_Sync();
}

void get_variant( double *buffer, double *abuff, double *cbuff, int ga_a, int ga_c,
int row_pgrp, int row_rank, int rank,int64_t nE)
{
    int64_t alo[2],ahi[2],ldim,sdim,i,j,k;
    register double s;
    const double oonE=1.0/((double)nE);    
    NGA_Distribution64(ga_a,rank,alo,ahi);
    ldim = ahi[1]-alo[1]+1;
    sdim = ahi[0]-alo[0]+1;
    NGA_Get64(ga_a,alo,ahi,abuff,&ldim);
    k=0;
    for (i=alo[0];i<=ahi[0];++i) {
        s=0;
        for (j=alo[1];j<=ahi[1];++j) {
            s+=abuff[k];
            ++k;
        }
        buffer[(i-alo[0])]=s*oonE;
    }
    GA_Pgroup_dgop(row_pgrp,buffer,sdim,"+");
    k=0;
    for (i=alo[0];i<=ahi[0];++i) {
        s=buffer[(i-alo[0])];
        for (j=alo[1];j<=ahi[1];++j) {
            cbuff[k]=abuff[k]-s;
            ++k;
        }
    }    
    if (row_rank==0) {
        NGA_Put64(ga_c,alo,ahi,cbuff,&ldim);
    }
    GA_Sync();
}


void check_cpu_bind_pid ( int num_thds, int num_cores, int thd_id )
{
    int k, j, e;
    cpu_set_t set;
    pid_t pid= getpid();
    sched_getaffinity( pid, sizeof ( set ), &set );
    e = -1;
    for ( j = 0; j < num_cores; ++j ) {
        e = CPU_ISSET ( j, &set );
        if ( e ) {
            flockfile(stderr);
            fprintf ( stderr, "thread %d is bound to core %d\n", thd_id, j );
            funlockfile(stderr);
            break;
        }
    }
    if ( e == -1 ) {
        flockfile(stderr);
        fprintf ( stderr, "thread %d is not bound to any core!\n", thd_id );
        funlockfile(stderr);
    }
    return;
}


void cpu_bind_pid( int num_thds, int num_cores, int thd_id, int strat )
{
    cpu_set_t set;
    const size_t cpu_set_size = sizeof ( cpu_set_t );
    int kcore, k, st_core;
    pid_t pid = getpid();
//    if (thd_id==0) {
//        flockfile(stderr);
//        fprintf ( stderr, "# threads = %d  # cores = %d  strategy = %d\n", num_thds, num_cores, strat );
//        funlockfile(stderr);
//    }
    kcore=0xFFFF;
    switch ( strat ) {
    case 0:
        kcore = (thd_id % num_cores);
        break;
    case 1:
        st_core = num_cores - 1;
        kcore = num_cores - 1;
        for ( k = 0; k < thd_id; ++k ) {
            kcore -= 1;
            if ( kcore < 0 ) {
                kcore = st_core;
            }
        }
        break;
    case 2:
        st_core = num_cores - 2;
        kcore = num_cores - 1;
        for ( k = 0; k < thd_id; ++k ) {
            kcore -= 2;
            if ( kcore < 0 ) {
                kcore = st_core;
                st_core= (st_core==(num_cores-1))?(num_cores-2):(num_cores-1);  
            }
        }
        break;
    default:
        st_core = 1;
        kcore = 0;
        for ( k = 0; k < thd_id; ++k ) {
                kcore += 2;
                if ( kcore >= num_cores ) {
                    kcore = st_core;
                    st_core = ( st_core == 0 ) ? 1 : 0;
                }
        }
        break;
    }
    CPU_ZERO ( &set );
    CPU_SET ( kcore, &set );
    sched_setaffinity ( pid, cpu_set_size, &set );
    return;
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
    double *rowBuffer, *kalman;
    int ga_xpri, ga_xpri_a, ga_recv, ga_recv_f, ga_recv_0, ga_vth,
        ga_merr, ga_hmeas, ga_hmeas_a, ga_ph_t, ga_qmat, ga_zmat, ga_wmat, ga_kmat,
        ga_kalman, ga_xpost, ga_outvh, ga_outvx;
    int local_pid, local_cid, local_np;
    int xrow_rank;
    int zrow_rank;
    int xrow_pgrp;
    int zrow_pgrp;
    int nblocks[2];
    int proc_list[32];
    int phandle;
    int ga_qmatL;
    int ga_qinvL;
    long heap_size=1048576L*1024L;
    long stack_size=1048576L*512L;
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
    double *zbuff, *hbuff, * zbuff2;
    int64_t nXrow,nZrow;
    FILE *fp;
    Stopwatch_t timer;
    Stopwatch_t itimer;
    Stopwatch_t oloop_timer;
    Stopwatch_t iloop_timer;
    Stopwatch_t qinv_timer;
    Stopwatch_t mm1_timer,mm2_timer,mm3_timer,umatrix_timer;
    Stopwatch_t mm4_timer;
    Stopwatch_t mio2_timer;
    int nsteps_max;
    char * input_name;
    
/**********************************************************************/
/********** start execution here! */
/**********************************************************************/
    MPI_Init ( &argc, &argv );
    GA_Initialize();
    replace_ma();
    MA_init(MT_DBL,stack_size,heap_size);
    rank = GA_Nodeid();
    nproc = GA_Nnodes();
    local_cid = GA_Cluster_nodeid();
    local_np = GA_Cluster_nprocs ( local_cid );
    local_pid = -1;
    for (j=0;j<local_np;++j) {
        if (GA_Cluster_procid(local_cid,j)==rank) {
            local_pid = j;
        }
    }
    if (local_pid==-1) {
        fprintf(stderr,"local pid = -1 , problem!\n");
        exit(EXIT_FAILURE);
    }
    parse_args(argc,argv,&nsteps_max,&input_name);
    cpu_bind_pid(local_np,32,local_pid,2);
    if (!local_pid) {
        fprintf(stderr,"there are %d procs on this node %d %d \n",local_np,local_cid,rank);
    }

    for (j=0;j<local_np;++j) proc_list[j]= GA_Cluster_procid(local_cid,j);
    phandle=GA_Pgroup_create(proc_list,local_np);

    stopwatch_clear ( &timer );
    stopwatch_clear ( &itimer );
    stopwatch_clear ( &iloop_timer );
    stopwatch_clear ( &oloop_timer );
    stopwatch_clear ( &qinv_timer );
    stopwatch_clear ( &mm1_timer );
    stopwatch_clear ( &mm2_timer );
    stopwatch_clear ( &mm3_timer );
    stopwatch_clear ( &umatrix_timer );
    stopwatch_clear ( &mm4_timer );
    stopwatch_clear ( &mio2_timer );
    stopwatch_start ( &timer );
    const double std_dev = 0.001;
    const double std_dev_2 = 0.05;
    const double omega0 = ONE;
    const double omegaB = 2.*M_PI * 60.0;
    const double time_step = 0.01;
    const double time_step_internal = 0.01;
    const double rfact = ONE / std_dev_2;
    const int64_t ninner = ( int64_t ) rint ( time_step / time_step_internal );
    double mag, ang;
    DoubleComplex z_one, z_zero;
    z_one.real = ONE;
    z_one.imag = ZERO;
    z_zero.real = ZERO;
    z_zero.imag = ZERO;
    if (input_name==0x0) {
        fp = fopen ( "data1.txt", "r" );
        if ( feof ( fp ) || !fp ) {
        fprintf ( stderr, "could not open input file for reading\n" );
        GA_Error ( "fatal error", 911 );
        }
        if (rank==0) fprintf(stderr,"input file data1.txt (default)\n");
    }else{
        fp = fopen ( input_name, "r" );
        if ( feof ( fp ) || !fp ) {
        fprintf ( stderr, "could not open input file %s for reading\n",input_name );
        GA_Error ( "fatal error", 911 );
        }
        if (rank==0) fprintf(stderr,"input file %s \n",input_name);
    }
    fscanf ( fp, "%ld%lg%lg", &nsteps, &tfault, &tclear );
    if (nsteps > nsteps_max) nsteps=nsteps_max;
    fscanf ( fp, "%ld%ld%ld%ld%ld%ld", &nB, &nG, &nZ, &nX, &nE, &ksize );
    if ( !rank ) {
        fprintf ( stderr, "nGen = %ld nBus = %ld  #Ensemble = %ld nX= %ld nZ= %ld\n",
                  nG, nB, nE, nX, nZ );
        fprintf ( stderr, "nsteps = %ld\n", nsteps );
        fprintf ( stderr, "nsteps max = %ld\n",(long)nsteps_max);
        fprintf ( stderr, "starting\n" );
        fprintf ( stderr, "k size = %ld\n", (long)ksize );
        fprintf ( stderr, "tfault = %lg\n", tfault );
        fprintf ( stderr, "tclear = %lg\n", tclear );
    }
    double oonE = ONE / ( ( double ) nE );
    double nfact = ONE / ( ( double ) ( nE - 1 ) );
    double rnfact = rfact * nfact;
    const int64_t kfault = ( int64_t ) rint ( tfault / time_step );
    const int64_t kclear = ( int64_t ) rint ( tclear / time_step );
    /*****************************************************************
     * create global arrays for data
     ****************************************************************/
    dims[0] = nB;
    dims[1] = nE;
    chunk[0] = -1;
    chunk[1] = -1;
    ga_outvh = NGA_Create64 ( C_DCPL, 2, dims, "ga_outvh", chunk );
    CHECK_HANDLE ( ga_outvh );
    dims[0] = nG;
    dims[1] = nE;
    ga_vth = NGA_Create64 ( C_DCPL, 2, dims, "ga_vth", chunk );
    ga_outvx = GA_Duplicate ( ga_vth, "ga_outvx" );
    CHECK_HANDLE ( ga_vth );
    CHECK_HANDLE ( ga_outvx );
    dims[0] = nB;
    dims[1] = nG;
    ga_recv_f = NGA_Create64 ( C_DCPL, 2, dims, "ga_recv_f", chunk );
    ga_recv_0 = GA_Duplicate ( ga_recv_f, "ga_recv_0" );
    CHECK_HANDLE ( ga_recv_f );
    CHECK_HANDLE ( ga_recv_0 );
    dims[0] = nX;
    dims[1] = nE;
    ga_xpri = expandedGA ( ga_outvx, rank, C_DBL, dims, "ga_xpri" );
    ga_xpri_a = GA_Duplicate ( ga_xpri, "ga_xpri_a" );
    GA_Sync();
    dims[0] = nZ;
    dims[1] = nE;
    ga_hmeas = expandedGA ( ga_outvh, rank, C_DBL, dims, "ga_hmeas" );
    ga_hmeas_a = GA_Duplicate ( ga_hmeas, "ga_hmeas_a" );
    ga_merr = GA_Duplicate ( ga_hmeas, "ga_merr" );
 
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
    dims[0] = nE;
    dims[1] = nE;
    ga_qmatL= NGA_Create_config64( C_DBL, 2, dims, "ga_qmatL", chunk, phandle );
    dims[0] = nE;
    dims[1] = nE;
    ga_qinvL= GA_Duplicate( ga_qmatL, "ga_qinvL");
    
    CHECK_HANDLE ( ga_kmat );
    CHECK_HANDLE ( ga_qmat );
    CHECK_HANDLE ( ga_zmat );
    CHECK_HANDLE ( ga_wmat );
    CHECK_HANDLE ( ga_ph_t );
    CHECK_HANDLE ( ga_qmatL );
    CHECK_HANDLE ( ga_qinvL );
    dims[0]=nsteps;
    dims[1]=nX;
    chunk[0]=-1;
    chunk[1]=-1;
    ga_xpost= NGA_Create64( C_DBL, 2, dims, "ga_xpost", chunk);
    CHECK_HANDLE(ga_xpost);

    /********************
     * write out distributions
     *******************/

    if (rank==0) {
      GA_Nblock(ga_xpri,nblocks);
      fprintf(stderr,"ga_xpri blocks %d %d \n",nblocks[0],nblocks[1]);
      GA_Nblock(ga_hmeas,nblocks);
      fprintf(stderr,"ga_hmeas blocks %d %d \n",nblocks[0],nblocks[1]);
      GA_Nblock(ga_qmat,nblocks);
      fprintf(stderr,"ga_qmat blocks %d %d \n",nblocks[0],nblocks[1]);
      GA_Nblock(ga_zmat,nblocks);
      fprintf(stderr,"ga_zmat blocks %d %d \n",nblocks[0],nblocks[1]);
      GA_Nblock(ga_wmat,nblocks);
      fprintf(stderr,"ga_wmat blocks %d %d \n",nblocks[0],nblocks[1]);
      GA_Nblock(ga_ph_t,nblocks);
      fprintf(stderr,"ga_ph_t blocks %d %d \n",nblocks[0],nblocks[1]);
    } 
    GA_Sync();     
    /*******************
     * create buffers for local data
     ******************/
    NGA_Distribution64 ( ga_xpri, rank, xlo, xhi );
    nXrow = xhi[0] - xlo[0] + 1L;
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
    ldim_h = hhi[1] - hlo[1] + 1L;
    nZrow = hhi[0]-hlo[0]+1L;
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
    /***************************
     * initialize parameter and buffer arrays
     **************************/
    Pm = ( double* ) Malloc ( sizeof ( double ) * nG );
    D = ( double* ) Malloc ( sizeof ( double ) * nG );
    Emag = ( double* ) Malloc ( sizeof ( double ) * nG );
    invHmac = ( double* ) Malloc ( sizeof ( double ) * nG );
    invXdp = ( double* ) Malloc ( sizeof ( double ) * nG );
    kalman = ( double* ) Malloc ( (sizeof ( double ) * nB * 2 *nsteps) );
    qinv = ( double* ) Malloc ( sizeof ( double ) * nE * nE );
    qbuff = ( double* ) Malloc ( sizeof ( double ) * nE * nE );
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
    for ( j = 0; j < nsteps; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            indx = j * nZ + i;
            kalman[indx] = s;
        }
    }
    for ( j = nsteps; j < ksize; ++j ) {
        for ( i = 0; i < nB; ++i ) fscanf ( fp, FMT, &s );
    }
    for ( j = 0; j < nsteps; ++j ) {
        for ( i = 0; i < nB; ++i ) {
            fscanf ( fp, FMT, &s );
            indx = j * nZ + i + nB; 
            kalman[indx] = s;
        }
    }
    fclose ( fp );
    if ( !rank ) fprintf ( stderr, "ended reading of file!\n" );
///////// find averages for initial state
    asize = ( nZrow > nXrow ) ? nZrow : nXrow;
    rowBuffer = ( double* ) replace_malloc ( ( sizeof ( double ) * asize ), 8, "rowBuffer" );
    xrow_pgrp = get_row_pgroup ( ga_xpri, rank, nproc );
    zrow_pgrp = get_row_pgroup ( ga_hmeas, rank, nproc );
    xrow_rank = GA_Pgroup_nodeid(xrow_pgrp);
    zrow_rank = GA_Pgroup_nodeid(zrow_pgrp);
    store_xpost(rowBuffer,xbuff,ga_xpost,ga_xpri,xrow_pgrp,xrow_rank,rank,0,nE);
    if ( !rank ) {
        fprintf ( stderr, "starting run\n" );
        fprintf ( stderr, "kfault = %ld kclear = %ld nsteps = %ld\n", kfault, kclear, nsteps );
        fprintf ( stderr, "# inner steps = %ld\n", ninner );
        fprintf ( stderr, "# nodes       = %d\n", nproc );
    }
    stopwatch_start ( &itimer );
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
        NGA_Distribution64 ( ga_xpri, rank, xlo, xhi );
        ldim_x = xhi[1] - xlo[1] + 1;
        NGA_Get64 ( ga_xpri, xlo, xhi, xbuff, &ldim_x  );
        GA_Sync();
        for ( istep = 0; istep < ninner; ++istep ) {
            NGA_Distribution64 ( ga_vth, rank, vxlo, vxhi );
            ldim_vx = vxhi[1] - vxlo[1] + 1;
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                xp = xbuff + ix * ldim_x;
                for ( j = 0; j < ldim_x; ++j, ++k ) {
                    vbuff[k].real = Emag[i] * cos ( xp[j] );
                    vbuff[k].imag = Emag[i] * sin ( xp[j] );
                }
            }
            NGA_Put64 ( ga_vth, vxlo, vxhi, vbuff, &ldim_vx );
            GA_Sync();
            GA_Zgemm64 ( 'n', 'n', nG, nE, nG, z_one, ga_recv, ga_vth, z_zero, ga_outvx );
            NGA_Get64 ( ga_outvx, vxlo, vxhi, vbuff, &ldim_vx );
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                delta = xbuff + ix * ldim_x;
                omega = delta + ldim_x;
                d_delta_dt = buff1 + ix * ldim_x;
                d_omega_dt = d_delta_dt + ldim_x;
                for ( j = 0; j < ldim_x; ++j, ++k ) {
                    tr = vbuff[k].real;
                    ti = vbuff[k].imag;
                    theta = atan2 ( ti, tr );
                    vmag = sqrt ( tr * tr + ti * ti );
                    c1 = omega[j] - omega0;
                    d_delta_dt[j] = omegaB * c1;
                    c2 = Pm[i] - D[i] * c1 - invXdp[i] * vmag * Emag[i] * sin ( delta[j] - theta );
                    d_omega_dt[j] = omega0 * invHmac[i] * HALF * c2;
                }
            }
            for ( i = xlo[0]; i <= xhi[0]; ++i ) {
                ix = i - xlo[0];
                xp = xbuff + ix * ldim_x;
                xf = xfwd + ix * ldim_x;
                dxdt1 = buff1 + ix * ldim_x;
                for ( j = 0; j < ldim_x; ++j ) {
                    xf[j] = xp[j] + ( dxdt1[j] * time_step_internal );
                }
            }
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                xp = xfwd + ix * ldim_x;
                for ( j = 0; j < ldim_x; ++j, ++k ) {
                    tr = Emag[i];
                    ti = xp[j];
                    vbuff[k].real = tr * cos ( ti );
                    vbuff[k].imag = tr * sin ( ti );
                }
            }
            NGA_Put64 ( ga_vth, vxlo, vxhi, vbuff, &ldim_x );
            GA_Sync();
            GA_Zgemm64 ( 'n', 'n', nG, nE, nG, z_one, ga_recv, ga_vth, z_zero, ga_outvx );
            NGA_Get64 ( ga_outvx, vxlo, vxhi, vbuff, &ldim_x );
            k = 0;
            for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
                ix = i + i - xlo[0];
                delta = xfwd + ix * ldim_x;
                omega = delta + ldim_x;
                d_delta_dt = buff2 + ix * ldim_x;
                d_omega_dt = d_delta_dt + ldim_x;
                for ( j = 0; j < ldim_x; ++j, ++k ) {
                    tr = vbuff[k].real;
                    ti = vbuff[k].imag;
                    theta = atan2 ( ti, tr );
                    vmag = sqrt ( tr * tr + ti * ti );
                    c1 = omega[j] - omega0;
                    d_delta_dt[j] = omegaB * c1;
                    c2 = Pm[i] - D[i] * c1 - invXdp[i] * vmag * Emag[i] * sin ( delta[j] - theta );
                    d_omega_dt[j] = omega0 * invHmac[i] * HALF * c2;
                }
            }
            for ( i = xlo[0]; i <= xhi[0]; ++i ) {
                ix = i - xlo[0];
                xp = xbuff + ix * ldim_x;
                dxdt1 = buff1 + ix * ldim_x;
                dxdt2 = buff2 + ix * ldim_x;
                for ( j = 0; j < ldim_x; ++j ) {
                    xp[j] += ( dxdt1[j] + dxdt2[j] ) * time_step_internal * HALF;
                }
            }
        }
/*** end of inner loop, form correction */        
        GA_Sync();
        stopwatch_stop ( &iloop_timer );
        stopwatch_start ( &oloop_timer );
        stopwatch_start ( &mio2_timer );
        NGA_Put64 ( ga_xpri, xlo, xhi, xbuff, &ldim_x );
        get_variant( rowBuffer, xbuff, buff1, ga_xpri, ga_xpri_a, xrow_pgrp, xrow_rank, rank, nE); 
        k = 0;
        for ( i = vxlo[0]; i <= vxhi[0]; ++i ) {
            ix = i + i - xlo[0];
            xp = xbuff + ix * ldim_x;
            for ( j = 0; j < ldim_x; ++j, ++k ) {
                CHECK_INDEX ( i, j, vxlo, vxhi );
                CHECK_INDEX ( ( i + i ), j, xlo, xhi);
                mag = Emag[i];
                ang = xp[j];
                vbuff[k].real = mag * cos ( ang );
                vbuff[k].imag = mag * sin ( ang );
            }
        }
        stopwatch_start ( &mio2_timer );
        NGA_Put64 ( ga_vth, vxlo, vxhi, vbuff, &ldim_x );
        GA_Sync();
        stopwatch_stop(&mio2_timer);
        stopwatch_start(&mm1_timer);
        GA_Zgemm64 ( 'n', 'n', nB, nE, nG, z_one, ga_recv, ga_vth, z_zero, ga_outvh );
        stopwatch_stop(&mm1_timer);
        stopwatch_start(&mio2_timer);
        NGA_Distribution64 ( ga_outvh, rank, vhlo, vhhi );
        ldim_vh = vhhi[1] - vhlo[1] + 1;
        NGA_Get64 ( ga_outvh, vhlo, vhhi, vbuff, &ldim_vh );
        k = 0;
        for ( i = vhlo[0]; i <= vhhi[0]; i ++ ) {
            iz = i + i - hlo[0];
            h1 = hbuff + iz * ldim_vh;
            h2 = h1 + ldim_vh;
            for ( j = 0; j < ldim_vh; ++j, ++k ) {
                CHECK_INDEX ( i, j, vhlo, vhhi );
                CHECK_INDEX ( ( i + i ), j, hlo, hhi );
                tr = vbuff[k].real;
                ti = vbuff[k].imag;
                h1[j] = sqrt ( ( tr * tr + ti * ti ) );
                h2[j] = atan2 ( ti, tr );
            }
        }
        NGA_Put64 ( ga_hmeas, hlo, hhi, hbuff, &ldim_h );
        ldim_h = hhi[1] - hlo[1] + 1;
        get_variant(rowBuffer,hbuff,zbuff,ga_hmeas,ga_hmeas_a,zrow_pgrp,zrow_rank,rank,nE);
        for ( i = hlo[0]; i <= hhi[0]; ++i ) {
            iz = i - hlo[0];
            ka = kalman [(kstep * nZ + i)]; 
            h1 = hbuff + iz * ldim_h;
            h2 = zbuff2 + iz * ldim_h;
            for ( j = 0; j < ldim_h; ++j ) {
                CHECK_INDEX ( i, j, hlo, hhi )
                h2[j] = ka * ( 1 + 0.0 * normRand() ) - h1[j];
            }
        }
        NGA_Put64 ( ga_merr, hlo, hhi, zbuff2, &ldim_h );
        stopwatch_stop ( &mio2_timer );
        GA_Sync();
        stopwatch_start(&umatrix_timer);
        ga_unit_dmatrix ( ga_qmat, qubuff, rank );
        ga_unit_dmatrix ( ga_wmat, wbuff, rank );
        GA_Sync();
        stopwatch_stop(&umatrix_timer);
        stopwatch_start ( &mm2_timer );
        GA_Dgemm64 ( 'n', 't', nX, nZ, nE, nfact, ga_xpri_a, ga_hmeas_a, 0.0, ga_ph_t );
        GA_Dgemm64 ( 't', 'n', nE, nE, nZ, rnfact, ga_hmeas_a, ga_hmeas_a, ONE, ga_qmat );
        stopwatch_stop ( &mm2_timer );
        stopwatch_start ( &qinv_timer );
        if (local_pid == 0) {
/**
 *   Take data out of global GA ga_qmat and put it into local one ga_qmaL
 */
          lo[0]=0;
          lo[1]=0;
          hi[0]=nE-1;
          hi[1]=nE-1;
          NGA_Get64(ga_qmat,lo,hi,qbuff,&nE);
          NGA_Put64(ga_qmatL,lo,hi,qbuff,&nE);
        }
        GA_Pgroup_sync(phandle);
        GA_Lu_solve('u',ga_qmatL,ga_qinvL);
        if (local_pid==0) {
          NGA_Get64(ga_qinvL,lo,hi,qbuff,&nE);
          NGA_Put64(ga_qmat,lo,hi,qbuff,&nE);
        }
        GA_Pgroup_sync(phandle);
        /**
                GA_Llt_solve ( ga_qmat, ga_qinv );
        **/
        stopwatch_stop ( &qinv_timer );
        stopwatch_start ( &mm3_timer );
        GA_Dgemm64 ( 'n', 'n', nZ, nE, nE, rnfact, ga_hmeas_a, ga_qmat, 0.0, ga_zmat );
        GA_Dgemm64 ( 'n', 't', nZ, nZ, nE, -ONE, ga_zmat, ga_hmeas_a, ONE, ga_wmat );
        stopwatch_stop(&mm3_timer);
        stopwatch_start(&mm4_timer);
        GA_Dgemm64 ( 'n', 'n', nX, nZ, nZ, rfact, ga_ph_t, ga_wmat, 0.0, ga_kmat );
        GA_Dgemm64 ( 'n', 'n', nX, nE, nZ, ONE, ga_kmat, ga_merr, ONE, ga_xpri );
        stopwatch_stop ( &mm4_timer );
        NGA_Get64 ( ga_xpri, xlo, xhi, xbuff, &ldim_x );
        store_xpost(rowBuffer,xbuff,ga_xpost,ga_xpri,xrow_pgrp,xrow_rank,rank,kstep,nE);
        stopwatch_stop ( &oloop_timer );
    }
    stopwatch_stop ( &itimer );
    GA_Sync();
    GA_Print ( ga_xpost );
    stopwatch_stop ( &timer );
    if ( rank == 0 ) {
        fprintf ( stderr, "tot  time = %lg\n", timer.acc );
        fprintf ( stderr, "init      = %lg\n", timer.acc - itimer.acc );
        fprintf ( stderr, "loop time = %lg\n", itimer.acc );
        fprintf ( stderr, "    outer = %lg\n", oloop_timer.acc );
        fprintf ( stderr, "    inner = %lg\n", iloop_timer.acc );
        fprintf ( stderr, "inverse   = %lg\n", qinv_timer.acc );
        fprintf ( stderr, "mmult1    = %lg\n", mm1_timer.acc );
        fprintf ( stderr, "mmult2    = %lg\n", mm2_timer.acc );
        fprintf ( stderr, "mmult3    = %lg\n", mm3_timer.acc );
        fprintf ( stderr, "mmult4    = %lg\n", mm4_timer.acc );
        fprintf ( stderr, "umatrix   = %lg\n", umatrix_timer.acc );
    }
    GA_Sync();
    GA_Pgroup_destroy ( zrow_pgrp);
    GA_Pgroup_destroy ( xrow_pgrp); 
    GA_Destroy ( ga_xpost );
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
//    fprintf ( stderr, "rank %d done x3\n", rank );
    GA_Terminate();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
