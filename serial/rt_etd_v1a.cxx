#include <omp.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

constexpr double KB = 8.625e-5;
constexpr double T = 823;
constexpr double PI = 2 * std::acos(0.0);
constexpr double TIME = 1.0e5;
constexpr double dt = 1.0e-1;
constexpr double dt_to_neg_2 = std::pow(dt, -2);

using namespace std::chrono;

void computeF(double *fx,
              const int i_start, const int i_end,
              const double x0, const double k,
              const double *alpha, const double *beta,
              const double *x) {
    for (int i = i_start; i <= i_end; i++) {
        // Original
        //fx[i-1] = beta[i-1]*x[i-1]*x0 - beta[i]*x[i]*x0 + alpha[i]*x[i+1] + k*x[i];
        // Version 2
        fx[i-1] = x0 * (beta[i-1] * x[i-1]) - x[i] * (x0 * beta[i] - k) + alpha[i] * x[i+1];
    }
}


void computeSumAlphaBeta(double sum_ab,
                         const double *alpha, const double *beta,
                         const double *x, const double x0,
                         const int N,
                         const int rank = 0) {
    double sum_alpha = 0.0;
    double sum_beta = 0.0;

    for (int i = 0; i < N; i++) {
        sum_alpha += x[i+1] * alpha[i];
        sum_beta += x[i+1] * beta[i+1];
    }
    if (rank == 0) {
        sum_alpha -= x[1] * alpha[0];
        sum_beta -= x[1] * beta[1];
    }

    sum_ab = sum_alpha - x0 * sum_beta;
}


int main(int argc,char *argv[]) {
    // Get inputs
    int N = 50;
    if (argc > 1) {
        N = std::atoi(argv[1]);
    }

    // Start the timer
    auto start_time = high_resolution_clock::now();

    //double sumalpha=0.0,sumalpha_a=0.0,sumalpha_b=0.0,sumalpha_d=0.0;
    //double sumbeta=0.0,sumbeta_a=0.0,sumbeta_b=0.0,sumbeta_d=0.0;
    double vec_sum_ab[4]= {0.0}; // sum_alpha - sum_beta for A, B, C, D
    double tot_sum_ab[4]= {0.0};
    double t = 0.0;
    constexpr double k = 0.0;
    double Alpha0,Beta0;

    double Vat = 1.205e-29;								//m^3¡£1.205e-29m^3=1.205e-23cm^3
    double EFvac = 1.77;								//¿ÕÎ»µÄÐÎ³ÉÄÜeV
    double EMvac = 1.1;                                //¿ÕÎ»µÄÇ¨ÒÆÄÜeV
    double Dvac= 1.0e-6*std::exp(-EMvac/(KB*T));    //m^2/s¡£1.0e-6*exp(-EMvac/(KB*T))m^2/s ¿ÕÎ»µÄÀ©É¢ÏµÊý,m^2/s=1.0e+4nm^2/s
    double Gama = 6.25e18;						//eV/m^2¡£6.25e+18eV/m^2,eV/m^2=1.0e-4eV/nm^2
    double Cinit = 1.0e-7;                              //atom-1 ³õÊ¼Å¨¶È
    Beta0 = Alpha0 = std::pow((double)(48*PI*PI/Vat/Vat),(double)1/3)*Dvac;

    int myid = 0, numprocs = 1;
    int left,right,tag1=1,tag2=2,tag3=3,tag4=4,tag5=5,tag6=6,tag7=7,tag8=8,tag_a=9,tag_b=10;

    int n = N;
    int eqn_start = 1;
    int eqn_end = n;
    // Does this necessary?
    // We could compute F(0) twice for rank 0
    //if (myid == 0)
    //    eqn_start=2;

    // Does this necessary?
    // We could ensure each array involved has length n+1 with an ending 0 value.
    //if (myid==numprocs-1)
    //    eqn_end=n-1;

    auto vec_C = (double *)calloc((n+2),sizeof(double));
    auto vec_A = (double *)calloc((n+2),sizeof(double));
    auto vec_B = (double *)calloc((n+2),sizeof(double));
    auto vec_D = (double *)calloc((n+2),sizeof(double));
    auto vec_Fa = (double *)calloc(n,sizeof(double));
    auto vec_Fb = (double *)calloc(n,sizeof(double));
    auto vec_Fc = (double *)calloc(n,sizeof(double));
    auto vec_Fd = (double *)calloc(n,sizeof(double));
    auto vec_r = (double *)calloc(n,sizeof(double));
    auto vec_EBvac = (double *)calloc(n,sizeof(double));
    auto vec_Alpha = (double *)calloc((n+1),sizeof(double));
    auto vec_Beta = (double *)calloc((n+1),sizeof(double));
    auto vec_L3 = (double *)calloc(n,sizeof(double));
    auto vec_phi1 = (double *)calloc(n,sizeof(double));
    auto vec_phi2 = (double *)calloc(n,sizeof(double));
    auto vec_phi3 = (double *)calloc(n,sizeof(double));
    auto vec_poly1 = (double *)calloc(n,sizeof(double));
    auto vec_poly2 = (double *)calloc(n,sizeof(double));
    auto vec_poly3 = (double *)calloc(n,sizeof(double));

    //double acons = 6.02e23;
    //double at2cb = 8.385e22;
    //double cb2at = 1./at2cb;

    for (int i=0; i<n; i++) {
        vec_r[i]=std::pow((double)(3*(myid*n+i+1)*Vat/(4*PI)),(double)1/3);			//nm
        vec_EBvac[i]=EFvac-2*Gama*Vat/vec_r[i];								//eV¿ÕÎ»È±ÏÝ½áºÏÄÜ
        vec_Alpha[i]=Alpha0*std::pow((double)(myid*n+i+1),(double)1/3)*std::exp(-vec_EBvac[i]/(KB*T));
        vec_Beta[i+1]=Beta0*std::pow((double)(myid*n+i+1),(double)1/3);				//¹«Ê½2

        vec_L3[i] = -std::pow((double)(vec_Alpha[i]+k),(double)-3);
        vec_phi1[i] = std::exp(-dt*(vec_Alpha[i]+k)/2);
        vec_phi2[i] = (1-std::exp(-dt*(vec_Alpha[i]+k)/2))/(vec_Alpha[i]+k);
        vec_phi3[i] = std::exp(-dt*(vec_Alpha[i]+k));
        vec_poly1[i] = -4+dt*(vec_Alpha[i]+k)+std::exp(-dt*(vec_Alpha[i]+k))*(4+3*dt*(vec_Alpha[i]+k)+dt*dt*(vec_Alpha[i]+k)*(vec_Alpha[i]+k));
        vec_poly2[i] = 4-2*dt*(vec_Alpha[i]+k)+std::exp(-dt*(vec_Alpha[i]+k))*(-4-2*dt*(vec_Alpha[i]+k));
        vec_poly3[i] = -4+3*dt*(vec_Alpha[i]+k)-dt*dt*(vec_Alpha[i]+k)*(vec_Alpha[i]+k)+std::exp(-dt*(vec_Alpha[i]+k))*(4+dt*(vec_Alpha[i]+k));
    }

    if (myid == 0) {
        vec_L3[0] = -std::pow((double)(1+k),(double)-3);
        vec_phi1[0] = std::exp(-dt*(1+k)/2);
        vec_phi2[0] = (1-std::exp(-dt*(1+k)/2))/(1+k);
        vec_phi3[0] = std::exp(-dt*(1+k));
        vec_poly1[0] = -4+dt*(1+k)+std::exp(-dt*(1+k))*(4+3*dt*(1+k)+dt*dt*(1+k)*(1+k));
        vec_poly2[0] = 4-2*dt*(1+k)+std::exp(-dt*(1+k))*(-4-2*dt*(1+k));
        vec_poly3[0] = -4+3*dt*(1+k)-dt*dt*(1+k)*(1+k)+std::exp(-dt*(1+k))*(4+dt*(1+k));
    }

    if (myid == 0) {
        vec_C[1] = Cinit;
    }
    double C0 = Cinit;
    //C[0] = Cinit*at2cb;                 //atom-1×ª»»³Écm

    auto computeF0 = [&](const double x0, const double *x, const double sum_ab)
    {
        return vec_Alpha[1] * x[2] - x0 * (2*vec_Beta[1]*x0 - 1 - k) + sum_ab;
    };


    while (t < TIME) {

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(vec_sum_ab[0],
            vec_Alpha, vec_Beta, vec_C, C0, n, myid);
        //MPI_Reduce(vec_sum_ab, tot_sum_ab, 4, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_C(t)
        computeF(vec_Fc, eqn_start, eqn_end, C0, k, vec_Alpha, vec_Beta, vec_C);

        if (myid == 0) {
            vec_Fc[0] = computeF0(C0, vec_C, tot_sum_ab[0]);
        }

        // Compute A
        for (int i = 1; i < n+1; i++) {
            vec_A[i] = vec_phi1[i-1]*vec_C[i] + vec_phi2[i-1]*vec_Fc[i-1];
        }

        double A0 = vec_A[1];
        //MPI_Bcast(&A0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(vec_sum_ab[1],
            vec_Alpha, vec_Beta, vec_A, A0, n, myid);
        //MPI_Reduce(vec_sum_ab, tot_sum_ab, 4, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_A(t)
        computeF(vec_Fa, eqn_start, eqn_end, A0, k, vec_Alpha, vec_Beta, vec_A);

        if (myid == 0) {
            vec_Fa[0] = computeF0(A0, vec_A, tot_sum_ab[1]);
        }

        // Compute B
        for (int i = 1; i < n+1; i++) {
            vec_B[i] = vec_phi1[i-1]*vec_C[i] + vec_phi2[i-1]*vec_Fa[i-1];
        }

        double B0 = vec_B[1];
        //MPI_Bcast(&B0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(vec_sum_ab[2],
            vec_Alpha, vec_Beta, vec_B, B0, n, myid);
        //MPI_Reduce(vec_sum_ab, tot_sum_ab, 4, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_B(t)
        computeF(vec_Fb, eqn_start, eqn_end, B0, k, vec_Alpha, vec_Beta, vec_B);

        if (myid == 0) {
            vec_Fb[0] = computeF0(B0, vec_B, tot_sum_ab[2]);
        }

        // Compute D
        for (int i = 1; i < n+1; i++) {
            vec_D[i] = vec_phi1[i-1]*vec_A[i] + vec_phi2[i-1]*(2*vec_Fb[i-1]-vec_Fc[i-1]);
        }

        double D0 = vec_D[1];
        //MPI_Bcast(&D0, 1, MPI_DOUBLE, 0, comm);

        // Compute the sum of alpha and beta
        computeSumAlphaBeta(vec_sum_ab[3],
            vec_Alpha, vec_Beta, vec_D, D0, n, myid);
        //MPI_Reduce(vec_sum_ab, tot_sum_ab, 4, MPI_DOUBLE, MPI_SUM, 0, comm);

        // Compute F_D(t)
        computeF(vec_Fd, eqn_start, eqn_end, D0, k, vec_Alpha, vec_Beta, vec_D);

        if (myid == 0) {
            vec_Fd[0] = computeF0(D0, vec_D, tot_sum_ab[3]);
        }

        // Compute C
        for (int i = 1; i < n+1; i++) {
        /*
			vec_C[i] =
                vec_phi3[i-1] * vec_C[i]
                + std::pow((double)dt,(double)-2) * vec_L3[i-1]
                * (vec_poly1[i-1] * vec_Fc[i-1]
                    + vec_poly2[i-1] * (vec_Fa[i-1] + vec_Fb[i-1])
                    + vec_poly3[i-1] * vec_Fd[i-1]
                    );
        */
            vec_C[i] = vec_phi3[i-1] * vec_C[i]
                + dt_to_neg_2 * vec_L3[i-1]
                * (vec_poly1[i-1] * vec_Fc[i-1]
                    + vec_poly2[i-1] * (vec_Fa[i-1] + vec_Fb[i-1])
                    + vec_poly3[i-1] *vec_Fd[i-1]
                    );
        }

        if (myid == 0) {
            C0 = vec_C[1];
        }

        t += dt;
    }

    // Stop the timer
    auto end_time = high_resolution_clock::now();
    auto elapsed_time = duration_cast<milliseconds>(end_time - start_time).count() / 1000.0;
    printf("[Rank %d] N = %d, Total time = %.2f s\n", myid, N, elapsed_time);

    //for (int i = 0; i < n; i++) {
    //    printf("[Rank %d] C[%d] = %e\n", myid, myid*n+i+1, vec_C[i+1]);
    //}

    free(vec_C);
    free(vec_A);
    free(vec_B);
    free(vec_D);
    free(vec_Fa);
    free(vec_Fb);
    free(vec_Fc);
    free(vec_Fd);
    free(vec_r);
    free(vec_EBvac);
    free(vec_Alpha);
    free(vec_Beta);
    free(vec_L3);
    free(vec_phi1);
    free(vec_phi2);
    free(vec_phi3);
    free(vec_poly1);
    free(vec_poly2);
    free(vec_poly3);

    return 0;
}

