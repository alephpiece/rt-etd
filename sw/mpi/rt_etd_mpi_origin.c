#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <simd.h>
#include <mpi.h>

#define KB 8.625e-5
#define T 823
#define PI 3.1415926536
#define TIME 1.0e5

void *alloc_aligned(size_t alignment, size_t size) {
    void *ptr;
    posix_memalign(&ptr, alignment, size);
    memset(ptr, 0, size);
    return ptr;
}


void computeF(double *fx,
              size_t n,
              double x0, double k,
              double *alpha, double *beta,
              double *x) {
    size_t i;
    //for (i = 0; i < n; i++) {
    //    fx[i] = x0 * (beta[i] * x[i]) - x[i+1] * (x0 * beta[i+1] - k) + alpha[i+1] * x[i+2];
    //}
    size_t nblocks = (n/4)*4;
    doublev4 vfx, vx, vx1, vx2, vbeta, vbeta1, valpha1;
    doublev4 vx0 = x0;
    doublev4 vk = k;
    for (i = 0; i < nblocks; i+=4) {
        simd_load(vx, &(x[i]));
        simd_load(vbeta, &(beta[i]));
        simd_loadu(vx1, &(x[i+1]));
        simd_loadu(vbeta1, &(beta[i+1]));
        simd_loadu(vx2, &(x[i+2]));
        simd_loadu(valpha1, &(alpha[i+1]));

        vfx = simd_vmad(valpha1, vx2, simd_vmss(vx0, vbeta * vx, vx1 * simd_vmss(vx0, vbeta1, vk)));

        simd_storeu(vfx, &(fx[i]));
    }

    for (i = nblocks; i < n; i++) {
        fx[i] = x0 * (beta[i] * x[i]) - x[i+1] * (x0 * beta[i+1] - k) + alpha[i+1] * x[i+2];
    }
}


void computeSumAlphaBeta(double *sum_alpha, double *sum_beta,
                         double *alpha, double *beta,
                         double *x, const double x0,
                         const size_t N,
                         const int rank) {
    size_t i;
    for (i = 0; i < N; i++) {
        *sum_alpha += alpha[i] * x[i+1];
        *sum_beta += beta[i+1] * x[i+1] * x0;
    }
    //size_t nblocks = (N/4)*4;
    //doublev4 vsum_a = 0., vsum_b = 0.;
    //doublev4 valpha, vbeta1, vx1;
    //for (i = 0; i < nblocks; i+=4) {
    //    simd_load(valpha, &(alpha[i]));
    //    simd_loadu(vx1, &(x[i+1]));
    //    vsum_a = simd_vmad(valpha, vx1, vsum_a);
    //}
    //*sum_alpha += (double)(simd_vextf0(vsum_a) + simd_vextf1(vsum_a)
    //                + simd_vextf2(vsum_a) + simd_vextf3(vsum_a));

    //for (i = 0; i < nblocks; i+=4) {
    //    simd_loadu(vbeta1, &(beta[i+1]));
    //    simd_loadu(vx1, &(x[i+1]));
    //    vsum_b = simd_vmad(vbeta1, vx1, vsum_b);
    //}
    //*sum_beta += (double)(simd_vextf0(vsum_b) + simd_vextf1(vsum_b)
    //                + simd_vextf2(vsum_b) + simd_vextf3(vsum_b));

    //for (i = nblocks; i < N; i++) {
    //    *sum_alpha += alpha[i] * x[i+1];
    //    *sum_beta += beta[i+1] * x[i+1];
    //}
    //*sum_beta *= x0;

    if (rank == 0) {
        *sum_alpha -= alpha[0] * x[1];
        *sum_beta -= beta[1] * x[1] * x0;
    }
}


size_t getBlockSize(size_t N, int rank, int nprocs) {
    size_t quo = N / nprocs;
    size_t res = N % nprocs;
    int first = (quo * rank) + (res * rank / nprocs);   // int( r*n/p )
    int last = (quo * (rank + 1)) + (res * (rank + 1) / nprocs) - 1;   // int( (r+1)*n/p ) - 1
    
    return last - first + 1;
}


void mpiSendRecvHalos(double *x, size_t n, int left, int right,
                      MPI_Comm comm, MPI_Status *status) {
    MPI_Sendrecv(&x[n],1,MPI_DOUBLE,right,0,&x[0],1,MPI_DOUBLE,left,0,comm,status);
    MPI_Sendrecv(&x[1],1,MPI_DOUBLE,left,1,&x[n+1],1,MPI_DOUBLE,right,1,comm,status);
}


int main(int argc,char *argv[]) {
    MPI_Init(&argc,&argv);

    // MPI initialization
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;

    int myid = 0, numprocs = 1;
    MPI_Comm_rank(comm,&myid);
    MPI_Comm_size(comm,&numprocs);

    // MPI rank neighbors
    int left = myid - 1;
    int right = myid + 1;

    if (myid == 0) {
        left = MPI_PROC_NULL;
    }

    if (myid == numprocs-1) {
        right = MPI_PROC_NULL;
    }

    // Start the timer
    clock_t start_time = clock();

    if (!myid) {
        printf("[Rank %d] Reading arguments...\n", myid);
    }

    // Get inputs
    size_t N = 64;
    if (argc > 1) {
        N = atoi(argv[1]);
    }


    // Decomposition
    size_t n = getBlockSize(N, myid, numprocs);

    size_t checksum;
    MPI_Reduce(&n, &checksum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!myid) {
        printf("[Rank %d] Distributed N = %lu to %d processes\n", myid, N, numprocs);

        if (checksum != N) {
            printf("[Rank %d] Summation of n over processes = %lu, which is not %lu\n", myid, checksum, N);
            MPI_Barrier(comm);
            MPI_Abort(comm, 1);
        }
    }


    // Does this necessary?
    // We could compute F(0) twice for rank 0
    //if (myid == 0)
    //    eqn_start=2;

    // Does this necessary?
    // We could ensure each array involved has length n+1 with an ending 0 value.
    //if (myid==numprocs-1)
    //    eqn_end=n-1;

    if (!myid) {
        printf("[Rank %d] Initializing constants...\n", myid);
    }

    size_t i;
    //double sumalpha=0.0,sumalpha_a=0.0,sumalpha_b=0.0,sumalpha_d=0.0;
    //double sumbeta=0.0,sumbeta_a=0.0,sumbeta_b=0.0,sumbeta_d=0.0;
    double AlphaBeta[8]= {0.0}; //sumalpha,sumalpha_a,sumalpha_b,sumalpha_d,sumbeta,sumbeta_a,sumbeta_b,sumbeta_d
    double tot_AB[8]= {0.0};
    double t = 0.0;
    const double dt = 1.0e-1;
    const double k = 0.0;
    double Alpha0,Beta0;

    double Vat = 1.205e-29;								//m^3¡£1.205e-29m^3=1.205e-23cm^3
    double EFvac = 1.77;								//¿ÕÎ»µÄÐÎ³ÉÄÜeV
    double EMvac = 1.1;                                //¿ÕÎ»µÄÇ¨ÒÆÄÜeV
    double Dvac= 1.0e-6*exp(-EMvac/(KB*T));    //m^2/s¡£1.0e-6*exp(-EMvac/(KB*T))m^2/s ¿ÕÎ»µÄÀ©É¢ÏµÊý,m^2/s=1.0e+4nm^2/s
    double Gama = 6.25e18;						//eV/m^2¡£6.25e+18eV/m^2,eV/m^2=1.0e-4eV/nm^2
    double Cinit = 1.0e-7;                              //atom-1 ³õÊ¼Å¨¶È
    Beta0 = Alpha0 = pow((double)(48*PI*PI/Vat/Vat),(double)1/3)*Dvac;

    if (!myid) {
        printf("[Rank %d] Allocating memory...\n", myid);
    }

    double *C = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *A = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *B = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *D = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *Fa = (double *)alloc_aligned(32, (n+3)*sizeof(double));    // F[i] <- F(i'-1), e.g., F[1] is actually F(0)
    double *Fb = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *Fc = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *Fd = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *Alpha = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *Beta = (double *)alloc_aligned(32, (n+3)*sizeof(double));
    double *L3 = (double *)alloc_aligned(32, n*sizeof(double));
    double *phi1 = (double *)alloc_aligned(32, n*sizeof(double));
    double *phi2 = (double *)alloc_aligned(32, n*sizeof(double));
    double *phi3 = (double *)alloc_aligned(32, n*sizeof(double));
    double *poly1 = (double *)alloc_aligned(32, n*sizeof(double));
    double *poly2 = (double *)alloc_aligned(32, n*sizeof(double));
    double *poly3 = (double *)alloc_aligned(32, n*sizeof(double));

    if (!myid) {
        printf("[Rank %d] Initializing vectors...\n", myid);
    }

    //double acons = 6.02e23;
    //double at2cb = 8.385e22;
    //double cb2at = 1./at2cb;

    for (i=0; i<n; i++) {
        double r = pow((double)(3*(myid*n+i+1)*Vat/(4*PI)),(double)1/3);			//nm
        double EBvac = EFvac-2*Gama*Vat/r;								//eV¿ÕÎ»È±ÏÝ½áºÏÄÜ
        Alpha[i] = Alpha0*pow((double)(myid*n+i+1),(double)1/3)*exp(-EBvac/(KB*T));
        Beta[i+1] = Beta0*pow((double)(myid*n+i+1),(double)1/3);				//¹«Ê½2

        L3[i] = -pow((double)(Alpha[i]+k),(double)-3);
        phi1[i] = exp(-dt*(Alpha[i]+k)/2);
        phi2[i] = (1-exp(-dt*(Alpha[i]+k)/2))/(Alpha[i]+k);
        phi3[i] = exp(-dt*(Alpha[i]+k));
        poly1[i] = -4+dt*(Alpha[i]+k)+exp(-dt*(Alpha[i]+k))*(4+3*dt*(Alpha[i]+k)+dt*dt*(Alpha[i]+k)*(Alpha[i]+k));
        poly2[i] = 4-2*dt*(Alpha[i]+k)+exp(-dt*(Alpha[i]+k))*(-4-2*dt*(Alpha[i]+k));
        poly3[i] = -4+3*dt*(Alpha[i]+k)-dt*dt*(Alpha[i]+k)*(Alpha[i]+k)+exp(-dt*(Alpha[i]+k))*(4+dt*(Alpha[i]+k));
    }


    MPI_Sendrecv(&Beta[n], 1, MPI_DOUBLE, right, 0,
                 &Beta[0], 1, MPI_DOUBLE, left, 0, comm, &status);
    MPI_Sendrecv(&Alpha[0], 1, MPI_DOUBLE, left, 1,
                 &Alpha[n], 1, MPI_DOUBLE, right, 1, comm, &status);

    if (myid == 0) {
        L3[0] = -pow((double)(1+k),(double)-3);
        phi1[0] = exp(-dt*(1+k)/2);
        phi2[0] = (1-exp(-dt*(1+k)/2))/(1+k);
        phi3[0] = exp(-dt*(1+k));
        poly1[0] = -4+dt*(1+k)+exp(-dt*(1+k))*(4+3*dt*(1+k)+dt*dt*(1+k)*(1+k));
        poly2[0] = 4-2*dt*(1+k)+exp(-dt*(1+k))*(-4-2*dt*(1+k));
        poly3[0] = -4+3*dt*(1+k)-dt*dt*(1+k)*(1+k)+exp(-dt*(1+k))*(4+dt*(1+k));
    }

    if (myid == 0) {
        C[1] = Cinit;
    }
    double C0 = Cinit;
    //C[0] = Cinit*at2cb;                 //atom-1×ª»»³Écm

    size_t simd_blocks = (n/4)*4;  // simd blocks for length n
    doublev4 vc, va, vb, vd, vphi1, vphi2, vphi3, vfc, vfa, vfb, vfd, vl3, vpoly1, vpoly2, vpoly3;
    doublev4 vdt2 = pow(dt, -2);

    if (!myid) {
        printf("[Rank %d] Computing...\n", myid);
    }

    while (t < TIME) {

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[0], &AlphaBeta[4],
            Alpha, Beta, C, C0, n, myid);
        MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);

        mpiSendRecvHalos(C, n, left, right, comm, &status);

        // Compute F_C(t)
        computeF(Fc, n, C0, k, Alpha, Beta, C);

        if (myid == 0) {
            Fc[0] = Alpha[1] * C[2] - C0 * (2*Beta[1]*C0 - 1 - k) + tot_AB[0] - tot_AB[4];
        }

        // Compute A (simd)
        //for (i = 0; i < n; i++) {
        //    A[i+1] = phi1[i]*C[i+1] + phi2[i]*Fc[i];
        //}
        for (i = 0; i < simd_blocks; i+=4) {
            simd_load(vphi1, &(phi1[i]));
            simd_load(vphi2, &(phi2[i]));
            //simd_load(vfc, &(Fc[i]));  // in small-scale tests, there is some penalty with simd_load
            simd_loadu(vfc, &(Fc[i]));
            simd_loadu(vc, &(C[i+1]));

            va = simd_vmad(vphi1, vc, vphi2 * vfc);

            simd_storeu(va, &(A[i+1]));
        }
        for (i = simd_blocks; i < n; i++) {
            A[i+1] = phi1[i]*C[i+1] + phi2[i]*Fc[i];
        }

        double A0 = A[1];
        MPI_Bcast(&A0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[1], &AlphaBeta[5],
            Alpha, Beta, A, A0, n, myid);
        MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);

        mpiSendRecvHalos(A, n, left, right, comm, &status);

        // Compute F_A(t)
        computeF(Fa, n, A0, k, Alpha, Beta, A);

        if (myid == 0) {
            Fa[0] = Alpha[1] * A[2] - A0 * (2*Beta[1]*A0 - 1 - k) + tot_AB[1] - tot_AB[5];
        }

        // Compute B (simd)
        //for (i = 0; i < n; i++) {
        //    B[i+1] = phi1[i]*C[i+1] + phi2[i]*Fa[i];
        //}
        for (i = 0; i < simd_blocks; i+=4) {
            simd_load(vphi1, &(phi1[i]));
            simd_load(vphi2, &(phi2[i]));
            //simd_load(vfa, &(Fa[i]));  // in small-scale tests, there is some penalty with simd_load
            simd_loadu(vfa, &(Fa[i]));
            simd_loadu(vc, &(C[i+1]));

            vb = simd_vmad(vphi1, vc, vphi2 * vfa);

            simd_storeu(vb, &(B[i+1]));
        }
        for (i = simd_blocks; i < n; i++) {
            B[i+1] = phi1[i]*C[i+1] + phi2[i]*Fa[i];
        }

        double B0 = B[1];
        MPI_Bcast(&B0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[2], &AlphaBeta[6],
            Alpha, Beta, B, B0, n, myid);

        MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);
        mpiSendRecvHalos(B, n, left, right, comm, &status);

        // Compute F_B(t)
        computeF(Fb, n, B0, k, Alpha, Beta, B);

        if (myid == 0) {
            Fb[0] = Alpha[1] * B[2] - B0 * (2*Beta[1]*B0 - 1 - k) + tot_AB[2] - tot_AB[6];
        }

        // Compute D
        //for (i = 0; i < n; i++) {
        //    D[i+1] = phi1[i]*A[i+1] + phi2[i]*(2*Fb[i]-Fc[i]);
        //}
        for (i = 0; i < simd_blocks; i+=4) {
            doublev4 vtwos = 2.0;
            simd_load(vphi1, &(phi1[i]));
            simd_load(vphi2, &(phi2[i]));
            //simd_load(vfb, &(Fb[i]));  // in small-scale tests, there is some penalty with simd_load
            //simd_load(vfc, &(Fc[i]));  // but why?
            simd_loadu(vfb, &(Fb[i]));
            simd_loadu(vfc, &(Fc[i]));
            simd_loadu(va, &(A[i+1]));

            vd = simd_vmad(vphi2, simd_vmss(vtwos, vfb, vfc), vphi1 * va);

            simd_storeu(vd, &(D[i+1]));
        }
        for (i = simd_blocks; i < n; i++) {
            D[i+1] = phi1[i]*A[i+1] + phi2[i]*(2*Fb[i]-Fc[i]);
        }


        double D0 = D[1];
        MPI_Bcast(&D0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[3], &AlphaBeta[7],
            Alpha, Beta, D, D0, n, myid);

        MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);
        mpiSendRecvHalos(D, n, left, right, comm, &status);

        // Compute F_D(t)
        computeF(Fd, n, D0, k, Alpha, Beta, D);

        if (myid == 0) {
            Fd[0] = Alpha[1] * D[2] - D0 * (2*Beta[1]*D0 - 1 - k) + tot_AB[3] - tot_AB[7];
        }

        // Compute C (simd, most of the time reduced here)
        //for (i = 0; i < n; i++) {
        //    C[i+1] = phi3[i] * C[i+1]
        //        + pow(dt, -2) * L3[i]
        //        * (poly1[i] * Fc[i]
        //            + poly2[i] * (Fa[i] + Fb[i])
        //            + poly3[i] *Fd[i]
        //            );
        //}

        for (i = 0; i < simd_blocks; i+=4) {
            simd_load(vphi3, &(phi3[i]));
            simd_load(vpoly1, &(poly1[i]));
            simd_load(vpoly2, &(poly2[i]));
            simd_load(vpoly3, &(poly3[i]));
            simd_load(vl3, &(L3[i]));
            simd_load(vfc, &(Fc[i]));
            simd_load(vfa, &(Fa[i]));
            simd_load(vfb, &(Fb[i]));
            simd_load(vfd, &(Fd[i]));
            simd_loadu(vc, &(C[i+1]));

            vc = simd_vmad(vphi3, vc,
                    vdt2 * vl3 * simd_vmad(vpoly1, vfc,
                        simd_vmad(vpoly2, vfa + vfb, vpoly3 * vfd)));

            simd_storeu(vc, &(C[i+1]));
        }
        for (i = simd_blocks; i < n; i++) {
            C[i+1] = phi3[i] * C[i+1]
                + pow(dt, -2) * L3[i]
                * (poly1[i] * Fc[i]
                    + poly2[i] * (Fa[i] + Fb[i])
                    + poly3[i] *Fd[i]
                    );
        }

        if (myid == 0) {
            C0 = C[1];
        }
        MPI_Bcast(&C0, 1, MPI_DOUBLE, 0, comm);

        for (i=0; i<8; i++) {
            AlphaBeta[i] = 0.0;
        }

        t += dt;
    }

    // Stop the timer
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    double tot_time = elapsed_time;
    MPI_Reduce(&elapsed_time, &tot_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (!myid) {
        printf("[Rank %d] Processes = %d, N = %d, Total time = %.2f s\n", myid, numprocs, N, tot_time);

        //for (i = 0; i < n; i++) {
        //    printf("[Rank %d] C[%d] = %e\n", myid, myid*n+i+1, C[i+1]);
        //}
    }

    free(C);
    free(A);
    free(B);
    free(D);
    free(Fa);
    free(Fb);
    free(Fc);
    free(Fd);
    free(Alpha);
    free(Beta);
    free(L3);
    free(phi1);
    free(phi2);
    free(phi3);
    free(poly1);
    free(poly2);
    free(poly3);

    MPI_Finalize();

    return 0;
}

