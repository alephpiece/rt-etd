#include <omp.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <malloc.h>

constexpr double KB = 8.625e-5;
constexpr double T = 823;
constexpr double PI = 2 * std::acos(0.0);
constexpr double TIME = 1.0e5;
constexpr double dt = 1.0e-1;
constexpr double dt_to_neg_2 = std::pow(dt, -2);

using namespace std::chrono;


#define VEC_LENGTH 64
#define VEC_ALIGN 32
//void* alignedAlloc(size_t count, size_t size) {
//    return std::calloc(count, size);
//    //auto chunks = count / VEC_LENGTH + (count % VEC_LENGTH != 0);
//    //auto tot_size = chunks * VEC_LENGTH * size;
//    //
//    //auto p = memalign(VEC_ALIGN, tot_size);
//    //memset(p, 0, tot_size);
//
//    //return p;
//}


void computeF(double *fx,
              const int n_eqn,
              const double x0, const double k,
              const double *alpha, const double *beta,
              const double *x) {
    //asm volatile ("# computeF begin");
    for (int i = 1; i < n_eqn+1; i++) {
        // Original
        //fx[i-1] = beta[i-1]*x[i-1]*x0 - beta[i]*x[i]*x0 + alpha[i]*x[i+1] + k*x[i];
        // Version 2
        fx[i-1] = x0 * (beta[i-1] * x[i-1]) - x[i] * (x0 * beta[i] - k) + alpha[i] * x[i+1];

        // Version 3
        //fx[i-1] = x0 * (beta[i-1] * x[i-1]) + x[i] * (k - x0 * beta[i] * x[i]) + alpha[i] * x[i+1];
    }
    //asm volatile ("# computeF end");
}


void computeSumAlphaBeta(double &sum_alpha, double &sum_beta,
                         const double *x, const double x0,
                         const double *alpha, const double *beta,
                         const int N,
                         const int rank = 0) {
    for (int i = 0; i < N; i++) {
        sum_alpha += alpha[i] * x[i+1];
        sum_beta += beta[i+1] * x[i+1] * x0;
    }
    if (rank == 0) {
        sum_alpha -= alpha[0] * x[1];
        sum_beta -= beta[1] * x[1] * x0;
    }
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
    double vec_AB[8]= {0.0}; //sumalpha,sumalpha_a,sumalpha_b,sumalpha_d,sumbeta,sumbeta_a,sumbeta_b,sumbeta_d
    double vec_tot_AB[8]= {0.0};
    double t = 0.0,k=0.0;
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

    const int is_root = myid ? 0 : 1;
    const int not_root = !is_root;
    //printf("[Rank %d] is_root = %d\n", myid, is_root);

    // Does this necessary?
    // We can compute F[0] twice so that this could be eliminated.
    //int eqn_start = 1;
    //int eqn_end = n;
    //if (myid == 0)
    //    eqn_start=2;

    // Does this necessary?
    // We could ensure each array involved has length n+1 with an ending 0 value.
    //if (myid==numprocs-1)
    //    eqn_end=n-1;
    auto compute_FX1_X2 = [&](double *FX, double *X2, const double x0, const double *X1) {
        for (int i = 1; i < n+1; i++) {
            double fx = x0 * (vec_Beta[i-1] * X1[i-1]) - X1[i] * (x0 * vec_Beta[i] - k) + vec_Alpha[i] * X1[i+1];
            FX[i-1] = fx;
            X2[i] = vec_phi1[i-1]*X1[i] + vec_phi2[i-1]*fx;
        }
    };

    while (t < TIME) {

        // Compute F_C(t), A
        computeSumAlphaBeta(vec_AB[0], vec_AB[4], vec_C, C0,
            vec_Alpha, vec_Beta, n, myid);

        //computeF(vec_Fc, n, C0, k, vec_Alpha, vec_Beta, vec_C);
        compute_FX1_X2(vec_Fc, vec_A, C0, vec_C);

        //MPI_Reduce(vec_AB, vec_tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);
        double fc0 = C0 * (-2*vec_Beta[1]*C0 + 1 + k) - vec_tot_AB[4]+vec_tot_AB[0]+vec_Alpha[1]*vec_C[2];
        vec_Fc[0] = is_root * fc0 + not_root * vec_Fc[0];
        vec_A[1] = vec_phi1[0]*vec_C[1] + vec_phi2[0]*vec_Fc[0];

        // Compute F_A(t), B
        double A0 = vec_A[1];
        //MPI_Bcast(&A0, 1, MPI_DOUBLE, 0, comm);

        computeSumAlphaBeta(vec_AB[1], vec_AB[5], vec_A, A0,
            vec_Alpha, vec_Beta, n, myid);

        //computeF(vec_Fa, n, A0, k, vec_Alpha, vec_Beta, vec_A);
        compute_FX1_X2(vec_Fa, vec_B, A0, vec_A);

        //MPI_Reduce(vec_AB, vec_tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);
        double fa0 = A0 * (-2*vec_Beta[1]*A0 + 1 + k) - vec_tot_AB[5]+vec_tot_AB[1]+vec_Alpha[1]*vec_A[2];
        vec_Fa[0] = is_root * fa0 + not_root * vec_Fa[0];
        vec_B[1] = vec_phi1[0]*vec_C[1] + vec_phi2[0]*vec_Fa[0];

        // Compute F_B(t), D
        double B0 = vec_B[1];
        //MPI_Bcast(&B0, 1, MPI_DOUBLE, 0, comm);

        computeSumAlphaBeta(vec_AB[2], vec_AB[6], vec_B, B0,
            vec_Alpha, vec_Beta, n, myid);

        //computeF(vec_Fb, n, B0, k, vec_Alpha, vec_Beta, vec_B);
        // D is different
        for (int i = 1; i < n+1; i++) {
            double fx = B0 * (vec_Beta[i-1] * vec_B[i-1]) - vec_B[i] * (B0 * vec_Beta[i] - k) + vec_Alpha[i] * vec_B[i+1];
            vec_Fb[i-1] = fx;
            vec_D[i] = vec_phi1[i-1]*vec_A[i] + vec_phi2[i-1]*(2*fx - vec_Fc[i-1]);
        }

        //MPI_Reduce(vec_AB, vec_tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);
        double fb0 = B0 * (-2*vec_Beta[1]*B0 + 1 + k) - vec_tot_AB[6]+vec_tot_AB[2]+vec_Alpha[1]*vec_B[2];
        vec_Fb[0] = is_root * fb0 + not_root * vec_Fb[0];
        vec_D[1] = vec_phi1[0]*vec_A[1] + vec_phi2[0]*(2*vec_Fb[0] - vec_Fc[0]);

        // Compute F_D(t), C
        double D0 = vec_D[1];
        //MPI_Bcast(&D0, 1, MPI_DOUBLE, 0, comm);

        computeSumAlphaBeta(vec_AB[3], vec_AB[7], vec_D, D0,
            vec_Alpha, vec_Beta, n, myid);

        //computeF(vec_Fd, n, D0, k, vec_Alpha, vec_Beta, vec_D);
        // C is also different, note that the first i is 2
        for (int i = 2; i < n+1; i++) {
            double fx = D0 * (vec_Beta[i-1] * vec_D[i-1]) - vec_D[i] * (D0 * vec_Beta[i] - k) + vec_Alpha[i] * vec_D[i+1];
            vec_Fd[i-1] = fx;
            vec_C[i] = vec_phi3[i-1] * vec_C[i]
                + dt_to_neg_2 * vec_L3[i-1]
                * (vec_poly1[i-1] * vec_Fc[i-1]
                    + vec_poly2[i-1] * (vec_Fa[i-1] + vec_Fb[i-1])
                    + vec_poly3[i-1] * fx
                    );
        }

        //MPI_Reduce(vec_AB, vec_tot_AB, 8, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (is_root)
            vec_Fd[0] = D0 * (-2*vec_Beta[1]*D0 + 1 + k) - vec_tot_AB[7]+vec_tot_AB[3] + vec_Alpha[1]*vec_D[2];
        else
            vec_Fd[0] = D0 * (vec_Beta[0] * vec_D[0]) - vec_D[1] * (D0 * vec_Beta[1] - k) + vec_Alpha[1] * vec_D[2];
        vec_C[1] = vec_phi3[0] * vec_C[1]
            + dt_to_neg_2 * vec_L3[0]
            * (vec_poly1[0] * vec_Fc[0]
                + vec_poly2[0] * (vec_Fa[0] + vec_Fb[0])
                + vec_poly3[0] * vec_Fd[0]
                );

        C0 = vec_C[1];
        //MPI_Bcast(&C0, 1, MPI_DOUBLE, 0, comm);

        for (int i = 0; i < 8; i++) {
            vec_AB[i] = 0.0;
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

