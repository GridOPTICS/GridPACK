#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cfloat>
#include <ctime>
extern "C" {
#include <cblas.h>
#include <clapack.h>
}
using namespace std;

#ifdef _USE_FLOAT_
#define real_t float
#define EPS FLT_EPSILON
#else
#define real_t double
#define EPS DBL_EPSILON
#endif

struct Complex_t {
    real_t real,imag;
};

typedef struct Complex_t Complex_t;

struct pstopwatch_t {
    double acc;
    timespec ts, tf;

    pstopwatch_t() {
        acc = 0;
    };

    inline void start() {
        clock_gettime ( CLOCK_MONOTONIC, &ts );
    }
    inline void stop ( ) {
        clock_gettime ( CLOCK_MONOTONIC, &tf );
        acc += ( tf.tv_sec - ts.tv_sec ) + 1.e-9 * ( tf.tv_nsec - ts.tv_nsec );
        return;
    }
    double elapsed_time ( ) {
        return acc;
    };

    void pstopwatch_clear () {
        acc = 0;
    };

};


FILE * openFile ( const char *name, const char *mode )
{
    FILE *fp = fopen ( name, mode );
    if ( fp ) return fp;
    cerr << " could not open file " << name << " in mode " << mode << endl;
    exit ( EXIT_FAILURE );
}

inline real_t randXNormal ( )
{
    real_t x, y, r;
    const real_t cut = real_t(1) - EPS ;
    do {
        x = real_t ( rand() ) / RAND_MAX;
        y = real_t ( rand() ) / RAND_MAX;
        y = y + y - 1;
        r = x * x + y * y;
        if ( r <= cut ) break;
    } while ( 1 );
#ifdef _USE_FLOAT_
    x = -log2f ( r ) / r;
#else
    x = -log2 ( r ) / r;
#endif
    x = x + x;
    return ( x * y );
}


void getData ( FILE *fp,
               Complex_t *RecV_0,
               Complex_t *RecV_f,
               real_t *Kalman_input,
               real_t *invH_mac,
               real_t *Pm,
               real_t *D,
               real_t *Emag,
               real_t *invXdp,
               real_t *xpost,
               const int nB, const int nG, const int nX,
               const int nE,
               const int ksize,
               const int nsteps
             )
{
    int i, j;
    real_t din,tr;
#ifdef _USE_FLOAT_
#define FMT "%g"
#else
#define FMT "%lg"
#endif


    if ( !fp ) {
        fprintf ( stderr, "File pointer to data is null" );
        exit ( EXIT_FAILURE );
    }
    if ( feof ( fp ) ) {
        fprintf ( stderr, "data file  is empty" );
        exit ( EXIT_FAILURE );
    }
    // Reading data from files
    for ( j = 0; j < nG; j++ ) {
        for ( i = 0; i < nB; i++ ) {
            fscanf ( fp, FMT, &din);
            RecV_0[i*nG+j].real = din;
        }
    }
    for ( j = 0; j < nG; j++ ) {
        for ( i = 0; i < nB; i++ ) {
            fscanf ( fp, FMT, &din);
            tr= RecV_0[i*nG+j].real;
            RecV_0[i*nG+j].real = tr * cos(din);
            RecV_0[i*nG+j].imag = tr * sin(din);
        }
    }
    for ( j = 0; j < nG; j++ ) {
        for ( i = 0; i < nB; i++ ) {
            fscanf ( fp, FMT, &din);
            RecV_f[i*nG+j].real = din;
        }
    }
    for ( j = 0; j < nG; j++ ) {
        for ( i = 0; i < nB; i++ ) {
            fscanf ( fp, FMT, &din);
            tr= RecV_f[i*nG+j].real;
            RecV_f[i*nG+j].real = tr * cos(din);
            RecV_f[i*nG+j].imag = tr * sin(din);
        }
    }
    for ( j = 0; j < nE; j++ ) {
        for ( i = 0; i < nX; i++ ) {
            fscanf ( fp, FMT, ( xpost + nE * i + j ) );
        }
    }
    for ( i = 0; i < nG; i++ ) {
        fscanf ( fp, FMT, &din );
        invH_mac[i] = 1. / din;
    }
    for ( i = 0; i < nG; i++ ) {
        fscanf ( fp, FMT, ( Pm + i ) );
    }
    for ( i = 0; i < nG; i++ ) {
        fscanf ( fp, FMT, ( D + i ) );
    }
    for ( i = 0; i < nG; i++ ) {
        fscanf ( fp, FMT, ( Emag + i ) );
    }
    for ( i = 0; i < nG; i++ ) {
        fscanf ( fp, FMT, &din );
        invXdp[i] = 1. / din;
    }
    for ( j = 0; j < nsteps; j++ ) {
        for ( i = 0; i < nB; i++ ) {
            fscanf ( fp, FMT, ( Kalman_input + nsteps * i + j ) );
        }
    }
    for ( j = nsteps; j < ksize ; ++ j) {
        for (i=0; i<nB; ++i) {
            fscanf(fp,FMT,&din);
        }
    }
    for ( j = 0; j < nsteps; j++ ) {
        for ( i = 0; i < nB; i++ ) {
            fscanf ( fp, FMT, ( Kalman_input + nsteps * ( i + nB ) + j ) );
        }
    }
}

int main()
{
    int i, j, offset;
    int nsteps, nB, nE, nX, nG, nZ, ksize;
    real_t tfault, tclear;
    Complex_t *RecV_0,*RecV_f,*Vin,*Vout,*RecV;
    const real_t *delta,*omega;
    real_t *xpri, *xfwd, *xpost_all, *xpri_a,
    *dx_dt_e1, *dx_dt_e2, *invH_mac, *D, *Pm, *Emag, *invXdp, *h_meas, *h_meas_a,
    *PH_T, *Kgain, *Kalman_input,
    *meas_err, *qinv, *qmat, *zmat, *wmat;

    FILE *fp;
    pstopwatch_t timer1;
    pstopwatch_t timer2;
    pstopwatch_t timer3;
    pstopwatch_t timer4;
    pstopwatch_t timer5;
    pstopwatch_t timer6;
    pstopwatch_t itimer;
    itimer.start();
    const real_t w0 = 1.;
    const real_t wB = 2.*M_PI * 60.0;
    const real_t tstep_internal = 0.01;
    const real_t tstep = 0.01;
    const int Y_bus = 1;
    const real_t std_dev = 0.1;
    const real_t std_dev_2 = 0.1;
    fp = openFile ( "data1.txt", "r" );
#ifdef _USE_FLOAT_
    fscanf ( fp, "%d%g%g", &nsteps, &tfault, &tclear );
#else
    fscanf ( fp, "%d%lg%lg", &nsteps, &tfault, &tclear );
#endif
    fscanf ( fp, "%d%d%d%d%d%d", &nB, &nG, &nZ, &nX, &nE, &ksize );
    cerr << "nsteps = " << nsteps << " tfault = " << tfault << " tclear " << tclear << endl;
    cerr << " nB = " << nB << " nG = " << nG << " nX = " << nX << " nE = " << nE << " ksize = " << ksize << endl;
    const real_t oonE = 1. / real_t ( nE );
    const real_t tfact = 1. / real_t ( nE - 1 );
    const real_t rfact = 1./ std_dev_2;
    const real_t rtfact = tfact * rfact;
#ifdef _USE_FLOAT_
    const int n_inner = static_cast<int>( rintf ( tstep / tstep_internal ) );
    const int kfault = static_cast<int> ( rintf ( tfault / tstep ) );
    const int kclear = static_cast<int> ( rintf ( tclear / tstep ) );
#else
    const int n_inner = static_cast<int>( rint ( tstep / tstep_internal ) );
    const int kfault = static_cast<int> ( rint ( tfault / tstep ) );
    const int kclear = static_cast<int> ( rint ( tclear / tstep ) );
#endif
    cerr << " n_inner = " << n_inner << endl;
    cerr << " sizeof(real_t) = " << sizeof(real_t) << endl;
    try {
        RecV_0 = new Complex_t[nB*nG];
        RecV_f = new Complex_t[nB*nG];
        Vin = new Complex_t[nG*nE];
        Vout = new Complex_t[nB*nE];
        xpri = new real_t[nX*nE];
        xpri_a = new real_t[nX*nE];
        xpost_all = new real_t[nX*nsteps];
        xfwd = new real_t[nX*nE];
        dx_dt_e1 = new real_t[nX*nE];
        dx_dt_e2 = new real_t[nX*nE];
        invH_mac = new real_t[nG];
        Pm = new real_t[nG];
        D = new real_t[nG];
        Emag = new real_t[nG];
        invXdp = new real_t[nG];
        h_meas = new real_t[nZ*nE];
        h_meas_a = new real_t[nZ*nE];
        meas_err = new real_t[nZ*nE];
        Kalman_input = new real_t[nsteps*nZ];
        qinv=new real_t[nE*nE];
        qmat=new real_t[nE*nE];
        zmat=new real_t[nZ*nE];
        wmat=new real_t[nZ*nZ];
        PH_T = new real_t[ nX * nZ];
        Kgain = new real_t[nX*nZ];
    } catch ( ... ) {
        cerr << " ran out of memory!\n";
        exit ( EXIT_FAILURE );
    }
    getData ( fp, RecV_0, RecV_f, Kalman_input,
              invH_mac, Pm, D, Emag, invXdp, xpri, nB, nG, nX, nE, ksize, nsteps );
    fclose ( fp );
/// record initial state in xpost_all
    for ( i = 0; i < nX; ++i ) {
        offset = i * nE;
        real_t s = 0;
        for ( j = 0; j < nE; ++j ) s += xpri[offset+j];
        s *= oonE;
        xpost_all[nsteps*i] = s;
    }
    itimer.stop();
    timer1.start();
    int ii, kstep;
    const Complex_t zzero={ 0, 0 };
    const Complex_t zone ={ real_t(1), 0 };
    real_t s,cr,ci,c1,c2,vmag,theta;
    real_t sre, sim;
    real_t *h1, *h2, *dxdt1, *dxdt2;
    Complex_t *vptr;
#pragma omp parallel default(shared) private(i,j,offset,delta,omega,dxdt1,dxdt2,cr,ci,c1,c2,\
    vmag,theta,h1,h2,vptr,s,sre,sim)
    {
        for ( kstep = 1; kstep < nsteps; ++kstep ) {
            if ( ( kstep > kfault ) && ( kstep <= kclear ) ) {
                if ( Y_bus == 0 ) {
                    RecV = RecV_0;
                } else {
                    RecV = RecV_f;
                }
            } else {
                RecV = RecV_0;
            }
            timer2.start();
            for ( ii = 0; ii < n_inner; ++ii ) {
                /// generate V and theta
                // transpose E_gen
#pragma omp for
                for ( i = 0; i < nG; ++i ) {
                    offset = i * nE;
                    delta = xpri + offset + offset;
                    for ( j = 0; j < nE; ++j ) {
                        Vin[ offset + j ].real = Emag[i] * cos (delta[j] );
                        Vin[ offset + j ].imag = Emag[i] * sin (delta[j] );
                    }
                }
#pragma omp single
                {
#ifdef _USE_FLOAT_
                    cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                    nG,nE,nG,&zone,RecV,nG,Vin,nE,&zzero,Vout,nE);
#else
                    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                                nG,nE,nG,&zone,RecV,nG,Vin,nE,&zzero,Vout,nE);
#endif
                }
                // find dx_dt_e1
#pragma omp for
                for ( i = 0; i < nG; ++i ) {
                    offset = i * nE;
                    dxdt1 = dx_dt_e1 + offset + offset;
                    dxdt2 = dxdt1 + nE;
                    delta = xpri + offset + offset;
                    omega = delta + nE;
                    for ( j = 0; j < nE; ++j ) {
                        cr = Vout[offset+j].real;
                        ci = Vout[offset+j].imag;
                        vmag= sqrt(cr*cr+ci*ci);
                        theta = atan2(ci,cr);
                        c1 = omega[j] - w0;
                        dxdt1[j] = wB * c1;
                        c2 = Pm[i] - D[i] * c1 - Emag[i] * vmag * invXdp[i] * sin (delta[j]-theta);
                        dxdt2[j] = w0 * invH_mac[i] * 0.5f * c2;
                    }
                }
#pragma omp for
                for ( i = 0; i < nX; ++i ) {
                    offset = i * nE;
                    for ( j = 0; j < nE; ++j ) {
                        xfwd[offset+j] = xpri[offset+j] + dx_dt_e1[offset+j] * tstep_internal;
                    }
                }
                // Repeat for x(t+dt) step
#pragma omp for
                for ( i = 0; i < nG; ++i ) {
                    offset = i * nE;
                    delta = xfwd + offset + offset;
                    for ( j = 0; j < nE; ++j ) {
                        Vin[ offset + j ].real = Emag[i] * cos (delta[j] );
                        Vin[ offset + j ].imag = Emag[i] * sin (delta[j] );
                    }
                }
#pragma omp single
                {
#ifdef _USE_FLOAT_
                    cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                    nG,nE,nG,&zone,RecV,nG,Vin,nE,&zzero,Vout,nE);
#else
                    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                                nG,nE,nG,&zone,RecV,nG,Vin,nE,&zzero,Vout,nE);
#endif
                }
                // find dx_dt_e1
#pragma omp for
                for ( i = 0; i < nG; ++i ) {
                    offset = i * nE;
                    dxdt1 = dx_dt_e2 + offset + offset;
                    dxdt2 = dxdt1 + nE;
                    delta = xfwd + offset + offset;
                    omega = delta + nE;
                    for ( j = 0; j < nE; ++j ) {
                        cr = Vout[offset+j].real;
                        ci = Vout[offset+j].imag;
                        vmag= sqrt(cr*cr+ci*ci);
                        theta = atan2(ci,cr);
                        c1 = omega[j] - w0;
                        dxdt1[j] = wB * c1;
                        c2 = Pm[i] - D[i] * c1 - Emag[i] * vmag * invXdp[i] * sin (delta[j]-theta);
                        dxdt2[j] = w0 * invH_mac[i] * 0.5 * c2;
                    }
                }
#pragma omp for
                for ( i = 0; i < nX; ++i ) {
                    offset = i * nE;
                    for ( j = 0; j < nE; ++j ) {
                        xpri[offset+j] += 0.5f * ( dx_dt_e1[offset+j] + dx_dt_e2[offset+j] ) * tstep_internal;
                    }
                }
            }
            /// generate V and theta
            // transpose E_gen
#pragma omp for
            for ( i = 0; i < nG; ++i ) {
                offset = i * nE;
                delta = xfwd + offset + offset;
                for ( j = 0; j < nE; ++j ) {
                    Vin[ offset + j ].real = Emag[i] * cos (delta[j] );
                    Vin[ offset + j ].imag = Emag[i] * sin (delta[j] );
                }
            }
#pragma omp single
            {
#ifdef _USE_FLOAT_
                cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                nB,nE,nG,&zone,RecV,nG,Vin,nE,&zzero,Vout,nE);
#else
                cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            nB,nE,nG,&zone,RecV,nG,Vin,nE,&zzero,Vout,nE);
#endif
            }
            /// form product E_gen * RecV -> V and theta
#pragma omp for
            for ( i = 0; i < nB; ++i ) {
                offset = i * nE;
                h1 = h_meas + offset + offset;
                h2 = h1 + nB * nE;
                vptr = Vout + offset;
                for ( j = 0; j < nE; ++j ) {
                    cr = vptr[j].real;
                    ci=  vptr[j].imag;
                    h1[j] = sqrt ( ( cr * cr + ci * ci ) );
                    h2[j] = atan2 ( ci, cr );
                }
            }
            timer6.stop();
            // find mean of X and form xpri_a
#pragma omp for
            for ( i = 0; i < nX; ++i ) {
                offset = i * nE;
                s = 0;
                for ( j = 0; j < nE; ++j ) s += xpri[offset+j];
                s *= oonE;
                for ( j = 0; j < nE; ++j ) {
                    xpri_a[offset+j] = xpri[offset+j] - s;
                }
            }
            // find mean of X and form xpri_a
#pragma omp for
            for ( i = 0; i < nZ; ++i ) {
                offset = i * nE;
                s = 0;
                for ( j = 0; j < nE; ++j ) s += h_meas[offset+j];
                s *= oonE;
                for ( j = 0; j < nE; ++j ) {
                    h_meas_a[offset+j] = h_meas[offset+j] - s;
                }
            }

            /// generate a noisy measurement and subtract predicted measurement
#pragma omp for
            for ( i = 0; i < nZ; ++i ) {
                c2 = Kalman_input[ nsteps*i + kstep];
                for ( j = 0; j < nE; ++j ) {
                    c1 = randXNormal ( );
                    meas_err[i*nE+j] = c2 * ( 1 + std_dev * c1 ) - h_meas[i*nE+j];
                }
            }
#pragma omp for
            for (i=0; i<nZ; ++i) {
                for (j=0; j<nZ; ++j) wmat[i*nE+j]=0;
                wmat[i*nE+i]=real_t(1);
            }
#pragma omp for
            for (i=0; i<nE; ++i) {
                for (j=0; j<nE; ++j) qmat[i*nE+j]=0;
                qmat[i*nE+i]=real_t(1);
            }


#pragma omp single
            {
#ifdef _USE_FLOAT_
                cblas_sgemm ( CblasRowMajor, CblasNoTrans, CblasTrans, nX, nZ, nE,
                tfact, xpri_a, nE, h_meas_a, nE, 0., PH_T, nZ );


                cblas_sgemm( CblasRowMajor, CblasTrans, CblasNoTrans, nE, nE, nZ, rtfact,
                             h_meas_a, nE, h_meas_a, nE, 1, qmat, nE);

                // find the inverse of 1 + H_a' Rinv H_a
                int err = clapack_sposv ( CblasRowMajor, CblasUpper, nE, nE, qmat, nE, qinv, nE );
                if ( err ) {
                    cerr << " could not find inverse for run " << kstep << " err = " << err << endl;
                    exit ( EXIT_FAILURE );
                }

                cblas_sgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, nZ, nE, nE, rtfact,
                             h_meas_a, nE, qinv, nE, 0, zmat, nE);


                cblas_sgemm( CblasRowMajor, CblasNoTrans, CblasTrans, nZ, nZ, nE, -1,
                             zmat, nE, h_meas_a, nE, 1, wmat , nZ);

                cblas_sgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, nX, nZ, nZ, rfact,
                             PH_T, nZ, wmat, nZ, 0, Kgain, nZ);

                cblas_sgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, nX, nE, nZ, 1,
                             Kgain, nZ, meas_err, nE, 1, xpri, nE);
#else
                cblas_dgemm ( CblasRowMajor, CblasNoTrans, CblasTrans, nX, nZ, nE,
                              tfact, xpri_a, nE, h_meas_a, nE, 0., PH_T, nZ );

                cblas_dgemm( CblasRowMajor, CblasTrans, CblasNoTrans, nE, nE, nZ, rtfact,
                             h_meas_a, nE, h_meas_a, nE, 1, qmat, nE);

                // find the inverse of 1 + H_a' Rinv H_a
                int err = clapack_dposv ( CblasRowMajor, CblasUpper, nE, nE, qmat, nE, qinv, nE );
                if ( err ) {
                    cerr << " could not find inverse for run " << kstep << " err = " << err << endl;
                    exit ( EXIT_FAILURE );
                }

                cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, nZ, nE, nE, rtfact,
                             h_meas_a, nE, qinv, nE, 0, zmat, nE);


                cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasTrans, nZ, nZ, nE, -1,
                             zmat, nE, h_meas_a, nE, 1, wmat , nZ);

                cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, nX, nZ, nZ, rfact,
                             PH_T, nZ, wmat, nZ, 0, Kgain, nZ);

                cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, nX, nE, nZ, 1,
                             Kgain, nZ, meas_err, nE, 1, xpri, nE);

#endif
            }
            /// store mean values
#pragma omp for 
            for ( i = 0; i < nX; ++i ) {
                offset = i * nE;
                s = 0;
                for ( j = 0; j < nE; ++j ) s += xpri[offset+j];
                s *= oonE;
                xpost_all[nsteps*i + kstep  ] = s;
            }
        }
    }
    timer1.stop();
    double t = timer1.elapsed_time();
    fprintf ( stderr, "loop time = %lg\n", t );
    fprintf ( stderr, "init time = %lg\n", itimer.elapsed_time() );
    fprintf ( stderr, "tot  time = %lg\n", ( itimer.elapsed_time() + t ) );
#ifdef _USE_FLOAT_
    for ( i = 0; i < nX ; ++i ) {
        for ( j = 0; j < nsteps; ++j ) fprintf ( stdout, "%20.6f\t", xpost_all[i*nsteps+j] );
        fprintf ( stdout, "\n" );
    }
#else
    for ( i = 0; i < nX ; ++i ) {
        for ( j = 0; j < nsteps; ++j ) fprintf ( stdout, "%20.12lf\t", xpost_all[i*nsteps+j] );
        fprintf ( stdout, "\n" );
    }
#endif
    return EXIT_SUCCESS;
}




