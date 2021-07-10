#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define KB 8.625e-5
#define T 823
#define PI 3.141592
#define TIME 1.0e5

int main(int argc,char *argv[]){
    // Get inputs
    int N = 64;
    if (argc > 1) {
        N = atoi(argv[1]);
    }

    // Start the timer
    clock_t start_time = clock();

    int i;
	//double sumalpha=0.0,sumalpha_a=0.0,sumalpha_b=0.0,sumalpha_d=0.0;
	//double sumbeta=0.0,sumbeta_a=0.0,sumbeta_b=0.0,sumbeta_d=0.0;
	double bufAB[8]={0.0};//sumalpha,sumalpha_a,sumalpha_b,sumalpha_d,sumbeta,sumbeta_a,sumbeta_b,sumbeta_d
	double sumAB[8]={0.0};
	double dt = 1.0e-1;
	double t = 0.0,k=0.0;
	double Alpha0,Beta0;
	
	double Vat = 1.205e-29;								//m^3。1.205e-29m^3=1.205e-23cm^3
	double EFvac = 1.77;								//空位的形成能eV
	double EMvac = 1.1;                                //空位的迁移能eV
	double Dvac= 1.0e-6*exp(-EMvac/(KB*T));    //m^2/s。1.0e-6*exp(-EMvac/(KB*T))m^2/s 空位的扩散系数,m^2/s=1.0e+4nm^2/s
	double Gama = 6.25e18;						//eV/m^2。6.25e+18eV/m^2,eV/m^2=1.0e-4eV/nm^2
	double Cinit = 1.0e-7;                              //atom-1 初始浓度
	Beta0 = Alpha0 = pow((double)(48*PI*PI/Vat/Vat),(double)1/3)*Dvac;
	
	int myid = 0, numprocs = 1;
	int left,right,begin_line,end_line,tag1=1,tag2=2,tag3=3,tag4=4,tag5=5,tag6=6,tag7=7,tag8=8,tag_a=9,tag_b=10;
	int n = N;


	double *bufC = (double *)calloc((n+2),sizeof(double));
	double *bufa = (double *)calloc((n+2),sizeof(double));
	double *bufb = (double *)calloc((n+2),sizeof(double));
	double *bufd = (double *)calloc((n+2),sizeof(double));
	double *bufFa = (double *)calloc(n,sizeof(double));
	double *bufFb = (double *)calloc(n,sizeof(double));
	double *bufFc = (double *)calloc(n,sizeof(double));
	double *bufFd = (double *)calloc(n,sizeof(double));
	double *bufr = (double *)calloc(n,sizeof(double));
	double *bufEBvac = (double *)calloc(n,sizeof(double));
	double *bufAlpha = (double *)calloc((n+1),sizeof(double));
	double *bufalpha = (double *)calloc(n,sizeof(double));
	double *bufalpha_a = (double *)calloc(n,sizeof(double));
	double *bufalpha_b = (double *)calloc(n,sizeof(double));
	double *bufalpha_d = (double *)calloc(n,sizeof(double));
	double *bufBeta = (double *)calloc((n+1),sizeof(double));
	double *bufbeta = (double *)calloc(n,sizeof(double));
	double *bufbeta_a = (double *)calloc(n,sizeof(double));
	double *bufbeta_b = (double *)calloc(n,sizeof(double));
	double *bufbeta_d = (double *)calloc(n,sizeof(double));
	double *bufL3 = (double *)calloc(n,sizeof(double));
	double *bufphi1 = (double *)calloc(n,sizeof(double));
	double *bufphi2 = (double *)calloc(n,sizeof(double));
	double *bufphi3 = (double *)calloc(n,sizeof(double));
	double *bufpoly1 = (double *)calloc(n,sizeof(double));
	double *bufpoly2 = (double *)calloc(n,sizeof(double));
	double *bufpoly3 = (double *)calloc(n,sizeof(double));

	//double acons = 6.02e23;
    //double at2cb = 8.385e22;
    //double cb2at = 1./at2cb;
	if(myid==0){
		bufC[1] = Cinit;
	}
	double C0=Cinit;
	double a0,b0,d0;
	//C[0] = Cinit*at2cb;                 //atom-1转换成cm

	for(i=0;i<n;i++){
		bufr[i]=pow((double)(3*(myid*n+i+1)*Vat/(4*PI)),(double)1/3);			//nm
		bufEBvac[i]=EFvac-2*Gama*Vat/bufr[i];								//eV空位缺陷结合能
		bufAlpha[i]=Alpha0*pow((double)(myid*n+i+1),(double)1/3)*exp(-bufEBvac[i]/(KB*T));
		bufBeta[i+1]=Beta0*pow((double)(myid*n+i+1),(double)1/3);				//公式2
		
		bufL3[i] = -pow((double)(bufAlpha[i]+k),(double)-3);
		bufphi1[i] = exp(-dt*(bufAlpha[i]+k)/2);
		bufphi2[i] = (1-exp(-dt*(bufAlpha[i]+k)/2))/(bufAlpha[i]+k);
		bufphi3[i] = exp(-dt*(bufAlpha[i]+k));
		bufpoly1[i] = -4+dt*(bufAlpha[i]+k)+exp(-dt*(bufAlpha[i]+k))*(4+3*dt*(bufAlpha[i]+k)+dt*dt*(bufAlpha[i]+k)*(bufAlpha[i]+k));
		bufpoly2[i] = 4-2*dt*(bufAlpha[i]+k)+exp(-dt*(bufAlpha[i]+k))*(-4-2*dt*(bufAlpha[i]+k));
		bufpoly3[i] = -4+3*dt*(bufAlpha[i]+k)-dt*dt*(bufAlpha[i]+k)*(bufAlpha[i]+k)+exp(-dt*(bufAlpha[i]+k))*(4+dt*(bufAlpha[i]+k));
	}
	if(myid==0){
		bufL3[0] = -pow((double)(1+k),(double)-3);
		bufphi1[0] = exp(-dt*(1+k)/2);
		bufphi2[0] = (1-exp(-dt*(1+k)/2))/(1+k);
		bufphi3[0] = exp(-dt*(1+k));
		bufpoly1[0] = -4+dt*(1+k)+exp(-dt*(1+k))*(4+3*dt*(1+k)+dt*dt*(1+k)*(1+k));
		bufpoly2[0] = 4-2*dt*(1+k)+exp(-dt*(1+k))*(-4-2*dt*(1+k));
		bufpoly3[0] = -4+3*dt*(1+k)-dt*dt*(1+k)*(1+k)+exp(-dt*(1+k))*(4+dt*(1+k));
	}

	begin_line=1,end_line=n;
	if(myid==0)
		begin_line=2;
	if(myid==numprocs-1)
		end_line=n-1;
	
	while (t < TIME){
		
		//C
		for(i=0;i<n;i++){
			bufalpha[i]=bufAlpha[i]*bufC[i+1];
			bufbeta[i]=bufBeta[i+1]*bufC[i+1]*C0;
			bufAB[0] += bufalpha[i];	//bufAB[0]=sumalpha
			bufAB[4] += bufbeta[i];		//bufAB[4]=sumbeta
		}
		if(myid==0){
			bufAB[0] = bufAB[0]-bufalpha[0];
			bufAB[4] = bufAB[4]-bufbeta[0];
		}
		

		if(myid==0){
			bufFc[0] = -2*bufBeta[1]*bufC[1]*bufC[1]-sumAB[4]+sumAB[0]+bufAlpha[1]*bufC[2]+(1+k)*bufC[1];
			for(i=begin_line;i<=end_line;i++){
				bufFc[i-1] = bufBeta[i-1]*bufC[i-1]*bufC[1]-bufBeta[i]*bufC[i]*bufC[1]+bufAlpha[i]*bufC[i+1]+k*bufC[i];
			}
		}
		else if(myid==numprocs-1){
			for(i=begin_line;i<=end_line;i++){
				bufFc[i-1] = bufBeta[i-1]*bufC[i-1]*C0-bufBeta[i]*bufC[i]*C0+bufAlpha[i]*bufC[i+1]+k*bufC[i];
			}
			bufFc[n-1]=bufBeta[n-1]*bufC[n-1]*C0-bufBeta[n]*bufC[n]*C0+k*bufC[n-1];
		}
		else{
			for(i=begin_line;i<=end_line;i++){
				bufFc[i-1] = bufBeta[i-1]*bufC[i-1]*C0-bufBeta[i]*bufC[i]*C0+bufAlpha[i]*bufC[i+1]+k*bufC[i];
			}
		}

		//a
		for(i=1;i<n+1;i++){
			bufa[i] = bufphi1[i-1]*bufC[i] + bufphi2[i-1]*bufFc[i-1];
		}
		if(myid==0){
			a0=bufa[1];
		}
		for(i=0;i<n;i++){
			bufalpha_a[i]=bufAlpha[i]*bufa[i+1];
			bufbeta_a[i]=bufBeta[i+1]*bufa[i+1]*a0;
			bufAB[1] += bufalpha_a[i];	//bufAB[1]=sumalpha_a
			bufAB[5] += bufbeta_a[i];		//bufAB[5]=sumbeta_a
		}
		if(myid==0){
			bufAB[1] = bufAB[1]-bufalpha_a[0];
			bufAB[5] = bufAB[5]-bufbeta_a[0];
		}
		if(myid==0){
			bufFa[0] = -2*bufBeta[1]*bufa[1]*bufa[1]-sumAB[5]+sumAB[1]+bufAlpha[1]*bufa[2]+(1+k)*bufa[1];
			for(i=begin_line;i<=end_line;i++){
				bufFa[i-1] = bufBeta[i-1]*bufa[i-1]*bufa[1]-bufBeta[i]*bufa[i]*bufa[1]+bufAlpha[i]*bufa[i+1]+k*bufa[i];
			}
		}
		else if(myid==numprocs-1){
			for(i=begin_line;i<=end_line;i++){
				bufFa[i-1] = bufBeta[i-1]*bufa[i-1]*a0-bufBeta[i]*bufa[i]*a0+bufAlpha[i]*bufa[i+1]+k*bufa[i];
			}
			bufFa[n-1]=bufBeta[n-1]*bufa[n-1]*a0-bufBeta[n]*bufa[n]*a0+k*bufa[n-1];
		}
		else{
			for(i=begin_line;i<=end_line;i++){
				bufFa[i-1] = bufBeta[i-1]*bufa[i-1]*a0-bufBeta[i]*bufa[i]*a0+bufAlpha[i]*bufa[i+1]+k*bufa[i];
			}
		}
		
		//b
		for(i=1;i<n+1;i++){
			bufb[i] = bufphi1[i-1]*bufC[i] + bufphi2[i-1]*bufFa[i-1];
		}
		if(myid==0){
			b0=bufb[1];
		}
		for(i=0;i<n;i++){
			bufalpha_b[i]=bufAlpha[i]*bufb[i+1];
			bufbeta_b[i]=bufBeta[i+1]*bufb[i+1]*b0;
			bufAB[2] += bufalpha_b[i];	//bufAB[1]=sumalpha_a
			bufAB[6] += bufbeta_b[i];		//bufAB[5]=sumbeta_a
		}
		if(myid==0){
			bufAB[2] = bufAB[2]-bufalpha_b[0];
			bufAB[6] = bufAB[6]-bufbeta_b[0];
		}
		
		if(myid==0){
			bufFb[0] = -2*bufBeta[1]*bufb[1]*bufb[1]-sumAB[6]+sumAB[2]+bufAlpha[1]*bufb[2]+(1+k)*bufb[1];
			for(i=begin_line;i<=end_line;i++){
				bufFb[i-1] = bufBeta[i-1]*bufb[i-1]*bufb[1]-bufBeta[i]*bufb[i]*bufb[1]+bufAlpha[i]*bufb[i+1]+k*bufb[i];
			}
		}
		else if(myid==numprocs-1){
			for(i=begin_line;i<=end_line;i++){
				bufFb[i-1] = bufBeta[i-1]*bufb[i-1]*b0-bufBeta[i]*bufb[i]*b0+bufAlpha[i]*bufb[i+1]+k*bufb[i];
			}
			bufFb[n-1]=bufBeta[n-1]*bufb[n-1]*b0-bufBeta[n]*bufb[n]*b0+k*bufb[n-1];
		}
		else{
			for(i=begin_line;i<=end_line;i++){
				bufFb[i-1] = bufBeta[i-1]*bufb[i-1]*b0-bufBeta[i]*bufb[i]*b0+bufAlpha[i]*bufb[i+1]+k*bufb[i];
			}
		}
		
		//d
		for(i=1;i<n+1;i++){
			bufd[i] = bufphi1[i-1]*bufa[i] + bufphi2[i-1]*(2*bufFb[i-1]-bufFc[i-1]);
		}
		if(myid==0){
			d0=bufd[1];
		}
		for(i=0;i<n;i++){
			bufalpha_d[i]=bufAlpha[i]*bufd[i+1];
			bufbeta_d[i]=bufBeta[i+1]*bufd[i+1]*d0;
			bufAB[3] += bufalpha_d[i];	//bufAB[1]=sumalpha_a
			bufAB[7] += bufbeta_d[i];		//bufAB[5]=sumbeta_a
		}
		if(myid==0){
			bufAB[3] = bufAB[3]-bufalpha_d[0];
			bufAB[7] = bufAB[7]-bufbeta_d[0];
		}
		if(myid==0){
			bufFd[0] = -2*bufBeta[1]*bufd[1]*bufd[1]-sumAB[7]+sumAB[3]+bufAlpha[1]*bufd[2]+(1+k)*bufd[1];
			for(i=begin_line;i<=end_line;i++){
				bufFd[i-1] = bufBeta[i-1]*bufd[i-1]*bufd[1]-bufBeta[i]*bufd[i]*bufd[1]+bufAlpha[i]*bufd[i+1]+k*bufd[i];
			}
		}
		else if(myid==numprocs-1){
			for(i=begin_line;i<=end_line;i++){
				bufFd[i-1] = bufBeta[i-1]*bufd[i-1]*d0-bufBeta[i]*bufd[i]*d0+bufAlpha[i]*bufd[i+1]+k*bufd[i];
			}
			bufFd[n-1]=bufBeta[n-1]*bufd[n-1]*d0-bufBeta[n]*bufd[n]*d0+k*bufd[n-1];
		}
		else{
			for(i=begin_line;i<=end_line;i++){
				bufFd[i-1] = bufBeta[i-1]*bufd[i-1]*d0-bufBeta[i]*bufd[i]*d0+bufAlpha[i]*bufd[i+1]+k*bufd[i];
			}
		}
		
		for(i=1;i<n+1;i++){
			bufC[i] = bufphi3[i-1]*bufC[i]+pow((double)dt,(double)-2)*bufL3[i-1]*(bufpoly1[i-1]*bufFc[i-1]+bufpoly2[i-1]*(bufFa[i-1]+bufFb[i-1])+bufpoly3[i-1]*bufFd[i-1]);
		}
		if(myid==0){
			C0=bufC[1];
		}
		t = t + dt;
		for(i=0;i<8;i++){
			bufAB[i]=0.0;
		}
	}
	
    // Stop the timer
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("[Rank %d] N = %d, Total time = %.2f s\n", myid, N, elapsed_time);

    //for (i = 0; i < n; i++) {
    //    printf("[Rank %d] C[%d] = %e\n", myid, myid*n+i+1, bufC[i+1]);
    //}

	free(bufC);
	free(bufa);
	free(bufb);
	free(bufd);
	free(bufFa);
	free(bufFb);
	free(bufFc);
	free(bufFd);
	free(bufr);
	free(bufEBvac);
	free(bufAlpha);
	free(bufalpha);
	free(bufalpha_a);
	free(bufalpha_b);
	free(bufalpha_d);
	free(bufBeta);
	free(bufbeta);
	free(bufbeta_a);
	free(bufbeta_b);
	free(bufbeta_d);
	free(bufL3);
	free(bufphi1);
	free(bufphi2);
	free(bufphi3);
	free(bufpoly1);
	free(bufpoly2);
	free(bufpoly3);

	return 0;
}

