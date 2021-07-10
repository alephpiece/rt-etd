#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define KB 8.625e-5
#define T 823
#define PI 3.1415926536
#define TIME 1.0e5

void computeF(double *fx,
              const size_t i_start, const size_t i_end,
              const double x0, const double k,
              const double *alpha, const double *beta,
              const double *x) {
    size_t i;
    for (i = i_start; i <= i_end; i++) {
        // Original
        //fx[i-1] = beta[i-1]*x[i-1]*x0 - beta[i]*x[i]*x0 + alpha[i]*x[i+1] + k*x[i];
        // Version 2
        fx[i-1] = x0 * (beta[i-1] * x[i-1]) - x[i] * (x0 * beta[i] - k) + alpha[i] * x[i+1];
    }
}


void computeSumAlphaBeta(double *sum_alpha, double *sum_beta,
                         const double *alpha, const double *beta,
                         const double *x, const double x0,
                         const size_t N,
                         const int rank) {
    size_t i;
    for (i = 0; i < N; i++) {
        *sum_alpha += alpha[i] * x[i+1];
        *sum_beta += beta[i+1] * x[i+1] * x0;
    }

    if (rank == 0) {
        *sum_alpha -= alpha[0] * x[1];
        *sum_beta -= beta[1] * x[1] * x0;
    }
}


int main(int argc,char *argv[]) {
    // Get inputs
    size_t N = 64;
    if (argc > 1) {
        N = atoi(argv[1]);
    }

    int myid = 0, numprocs = 1;
    size_t n = N;
    size_t eqn_start = 1;
    size_t eqn_end = n;
    // Does this necessary?
    // We could compute F(0) twice for rank 0
    //if (myid == 0)
    //    eqn_start=2;

    // Does this necessary?
    // We could ensure each array involved has length n+1 with an ending 0 value.
    //if (myid==numprocs-1)
    //    eqn_end=n-1;

    // Start the timer
    printf("[Rank %d] Start the timer...\n", myid);
    clock_t start_time = clock();

    printf("[Rank %d] Allocate memory...\n", myid);
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

    double *C = (double *)calloc((n+2),sizeof(double));
    double *A = (double *)calloc((n+2),sizeof(double));
    double *B = (double *)calloc((n+2),sizeof(double));
    double *D = (double *)calloc((n+2),sizeof(double));
    double *Fa = (double *)calloc(n,sizeof(double));
    double *Fb = (double *)calloc(n,sizeof(double));
    double *Fc = (double *)calloc(n,sizeof(double));
    double *Fd = (double *)calloc(n,sizeof(double));
    double *Alpha = (double *)calloc((n+1),sizeof(double));
    double *Beta = (double *)calloc((n+1),sizeof(double));
    double *L3 = (double *)calloc(n,sizeof(double));
    double *phi1 = (double *)calloc(n,sizeof(double));
    double *phi2 = (double *)calloc(n,sizeof(double));
    double *phi3 = (double *)calloc(n,sizeof(double));
    double *poly1 = (double *)calloc(n,sizeof(double));
    double *poly2 = (double *)calloc(n,sizeof(double));
    double *poly3 = (double *)calloc(n,sizeof(double));

    //double acons = 6.02e23;
    //double at2cb = 8.385e22;
    //double cb2at = 1./at2cb;

    printf("[Rank %d] Initialize variables...\n", myid);
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

    while (t < TIME) {

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[0], &AlphaBeta[4],
            Alpha, Beta, C, C0, n, myid);
        //MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_C(t)
        computeF(Fc, eqn_start, eqn_end, C0, k, Alpha, Beta, C);

        if (myid == 0) {
            Fc[0] = Alpha[1] * C[2] - C0 * (2*Beta[1]*C0 - 1 - k) + tot_AB[0] - tot_AB[4];
        }

        // Compute A
        for (i = 1; i < n+1; i++) {
            A[i] = phi1[i-1]*C[i] + phi2[i-1]*Fc[i-1];
        }

        double A0 = A[1];
        //MPI_Bcast(&A0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[1], &AlphaBeta[5],
            Alpha, Beta, A, A0, n, myid);
        //MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_A(t)
        computeF(Fa, eqn_start, eqn_end, A0, k, Alpha, Beta, A);

        if (myid == 0) {
            Fa[0] = Alpha[1] * A[2] - A0 * (2*Beta[1]*A0 - 1 - k) + tot_AB[1] - tot_AB[5];
        }

        // Compute B
        for (i = 1; i < n+1; i++) {
            B[i] = phi1[i-1]*C[i] + phi2[i-1]*Fa[i-1];
        }

        double B0 = B[1];
        //MPI_Bcast(&B0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[2], &AlphaBeta[6],
            Alpha, Beta, B, B0, n, myid);
        //MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_B(t)
        computeF(Fb, eqn_start, eqn_end, B0, k, Alpha, Beta, B);

        if (myid == 0) {
            Fb[0] = Alpha[1] * B[2] - B0 * (2*Beta[1]*B0 - 1 - k) + tot_AB[2] - tot_AB[6];
        }

        // Compute D
        for (i = 1; i < n+1; i++) {
            D[i] = phi1[i-1]*A[i] + phi2[i-1]*(2*Fb[i-1]-Fc[i-1]);
        }

        double D0 = D[1];
        //MPI_Bcast(&D0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(&AlphaBeta[3], &AlphaBeta[7],
            Alpha, Beta, D, D0, n, myid);
        //MPI_Reduce(AlphaBeta, tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_D(t)
        computeF(Fd, eqn_start, eqn_end, D0, k, Alpha, Beta, D);

        if (myid == 0) {
            Fd[0] = Alpha[1] * D[2] - D0 * (2*Beta[1]*D0 - 1 - k) + tot_AB[3] - tot_AB[7];
        }

        // Compute C
        for (i = 1; i < n+1; i++) {
        /*
			C[i] =
                phi3[i-1] * C[i]
                + pow((double)dt,(double)-2) * L3[i-1]
                * (poly1[i-1] * Fc[i-1]
                    + poly2[i-1] * (Fa[i-1] + Fb[i-1])
                    + poly3[i-1] * Fd[i-1]
                    );
        */
            C[i] = phi3[i-1] * C[i]
                + pow(dt, -2) * L3[i-1]
                * (poly1[i-1] * Fc[i-1]
                    + poly2[i-1] * (Fa[i-1] + Fb[i-1])
                    + poly3[i-1] *Fd[i-1]
                    );
        }

        if (myid == 0) {
            C0 = C[1];
        }

        for (i=0; i<8; i++) {
            AlphaBeta[i] = 0.0;
        }

        t += dt;
    }

    // Stop the timer
    printf("[Rank %d] Stop the timer...\n", myid);
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("[Rank %d] N = %d, Total time = %.2f s\n", myid, N, elapsed_time);

    //for (i = 0; i < n; i++) {
    //    printf("[Rank %d] C[%d] = %e\n", myid, myid*n+i+1, C[i+1]);
    //}

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

    return 0;
}

