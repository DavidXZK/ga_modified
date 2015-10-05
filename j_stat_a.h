#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

const double PI = 3.1415926;
const int N = 3189;     //num of sample
const int RNUM = 250;   //random num
const int XNUM = 0;    //theta para num
const int CNUM = 6;
//para: mur mub cl ch k rou 6 theta 21
const int PNUM = CNUM + XNUM;    //para num
const int IN = 10;       //industry num
const int NN = 24;     //Z num
const int HH = 100;
const int P = 800;    //population quantity
const int M = 20;  //code size;
const int FLAG = 0;
int color[N];    //black or white 0 1 2
int indu[N];     //belong to which industry
//double X[N][XNUM];  //X变量
double Y[N][2],Z[N][NN];
double m[N][2*NN];
double t[N];     // parameter t 抽检比例
double rand1[2][RNUM],rand2[N][RNUM];
double obj[2*P];
double mp,vp,best;
double qij_mean[N],qijh_mean[N],qijl_mean[N];
double paras[PNUM];
double q0[N],W[2*NN][2*NN],rev[N][RNUM];
double F[N][RNUM],G[N][RNUM];
double q_bar[N][RNUM],q_ijr[N][RNUM],q_ijb[N][RNUM];

ofstream fout("result_ga_modified.txt");

//************************ normal distribution function*******************
double normal(double z){
	double temp = exp((-1)*z*z/2)/sqrt(2*PI);
	return temp;
}
double normcdf(double z){
	if(z>6) return 1;
	if(z<-6) return 0;
	static const double gamma = 0.231641900,
	a1 = 0.319381530,
	a2 = -0.356563782,
	a3 =  1.781477973,
	a4  = -1.821255978,
	a5  =  1.330274429;
	double k = 1.0 /(1 + fabs(z)*gamma);
	double n = k * (a1 + k * (a2 + k * (a3 + k *(a4 + k * a5))));
	n = 1 - normal(z) * n;
	if(z < 0)
		return 1.0 - n;
	return n;
}
//***************************readdata**************************
void readData(){
	cout<<"GA() in"<<endl;
	ifstream fw("WW");
	for(int i = 0;i < 2*NN;i ++){
		for(int j = 0;j < 2*NN;j ++){
			fw>>W[i][j];
		}
	}//for
	fw.close();
	//mur mub cl ch k rho east middle dindus1~9 lnasset lnyear
	//regsoe regforeign regprivate regcollective reglegal lna2 lny2 1
	//-3 10
	ifstream fin_para("paras.txt");
	for(int i=0;i<PNUM;i++){
		fin_para>>paras[i];
	}
/*	ifstream fin_range("para_range.txt");
	for(int i=0;i<PNUM;i++){
		fin_range>>para_range[i][0]>>para_range[i][1];
	}*/
	//fin_range.close();
	// range above
	mp = 1;vp = 0.8;best = 0.0;
	//ifstream fin_x("para_x.txt");    //para x
	ifstream fin_y("para_y.txt");   //para lny0 y1 color industry q0 t
	ifstream fin_z("para_z.txt");   //para x + wave1 + wave2 + q0 + 1
	for(int i=0;i<N;i++){
		/*for(int j = 0;j < XNUM;j ++){
			fin_x>>X[i][j];
		}*/
		for(int j = 0;j < NN;j ++){
			fin_z>>Z[i][j];
		}
		fin_y>>Y[i][0]>>Y[i][1]>>color[i]>>indu[i]>>q0[i]>>t[i];
		Y[i][0] /=100;
	}//for
	//fin_x.close();
	fin_y.close();
	fin_z.close();
	//读取随机正态，每次一致
	ifstream frand("randnums");
	ofstream rr("rrr");
	for(int i=0;i<2;i++){
		for(int j=0;j<RNUM;j++){
			frand>>rand1[i][j];
		}
	}
	for(int i=0;i<N;i++){
		for(int j=0;j<RNUM;j++){
			frand>>rand2[i][j];
			rr<<rand2[i][j]<<" ";
			if(abs(rand2[i][j])>=10){
				cout<<i<<" "<<j<<endl;
			}
		}
		rr<<endl;
	}
	rr.close();
	frand.close();
	cout<<"GA() out"<<endl;
}
//***************************gets******************************
void getS(double *para){
	for(int i=0;i<N;i++){
		for(int j=0;j<RNUM;j++){
			F[i][j] = normcdf(rand1[0][j] + para[0]);
			G[i][j] = normcdf(rand1[1][j] + para[1]);
			q_bar[i][j] = (F[i][j]*q0[i])/(F[i][j]*q0[i] + G[i][j]*(1-q0[i]));
			q_ijr[i][j] = F[i][j] + q_bar[i][j]*(1-F[i][j]);
			q_ijb[i][j] = q_bar[i][j]*(1-G[i][j]);
		}
	}
	double rho = para[5];
	double rhox = rho/(1.0 - rho);
	double sum,sum1,sum2;
	for(int i=0;i<N;i++){
		sum = 0,sum1 = 0,sum2 = 0;
		for(int j=0;j<RNUM;j++){
			sum += pow(q_bar[N][j],rhox);
			sum1 += pow(q_ijr[N][j],rhox);
			sum2 += pow(q_ijb[N][j],rhox);
		}
		qij_mean[i] = sum/RNUM;
		qijh_mean[i] = sum1/RNUM;
		qijl_mean[i] = sum2/RNUM;
	}//for	
	// we got qij_mean, ci.
	double rhoy = 1.0/(rho - 1.0);
	double revmean[N][2];
	for(int i=0;i<N;i++){
		sum = 0;
		double temp[RNUM];
		for(int j=0;j<RNUM;j++){
			temp[j] = exp(-(rand2[i][j]));
		}
		if(color[i]==0){
			double ch[RNUM],cl[RNUM],dh[RNUM],dl[RNUM],paih,pail;
			double sum_paih = 0,sum_pail = 0;
			for(int j=0;j<RNUM;j++){
				ch[j] = para[3]*temp[j];
				cl[j] = para[2]*temp[j];
				dh[j] = para[4]*pow(ch[j]/rho,rhoy);
				dl[j] = para[4]*pow(cl[j]/rho,rhoy);
				paih = t[i]*(ch[j]/rho - ch[j])*dh[j]*qijh_mean[i] + (1-t[i])*(ch[j]/rho - ch[j])*dh[j]*qij_mean[i];
				pail = t[i]*(cl[j]/rho - cl[j])*dl[j]*qijl_mean[i] + (1-t[i])*(cl[j]/rho - cl[j])*dl[j]*qij_mean[i];
				sum_paih += paih;
				sum_pail += pail;
			}
			if(sum_paih>sum_pail){
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dh[j]*ch[j]/rho*qijh_mean[i];
					sum += rev[i][j];
				}
			}
			else{
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dl[j]*cl[j]/rho*qijl_mean[i];
					sum += rev[i][j];
				}
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==1){        //HIGH
			for(int j=0;j<RNUM;j++){
				double ch = para[3]*temp[j];
				double h = ch/rho;
				double dh = para[4]*pow(h,rhoy);
				rev[i][j] = dh*h*qijh_mean[i];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==2){//low
			for(int j=0;j<RNUM;j++){
				double cl = para[2]*temp[j];
				double l = cl/rho;
				double dl = para[4]*pow(l,rhoy);
				rev[i][j] = dl*l*qijl_mean[i];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
		}//if
	}//for
	double TotalR[N];
	memset(TotalR,0,sizeof(double)*N);
	for(int i=0;i<RNUM;i++){
		double sums[IN];
		memset(sums,0,sizeof(sums));
		for(int j=0;j<N;j++){
			sums[indu[j]] += rev[j][i];
		}
		for(int j=0;j<N;j++){
			TotalR[j] += rev[j][i]/sums[indu[j]]; 
		}
	}//for RNUM
	for(int i=0;i<N;i++){
		if(TotalR[i]<0){
			cout<<i<<"TotalR = "<<TotalR[i]<<endl;
		}
		revmean[i][1] = TotalR[i]/RNUM;
	}
	// now we get phat == revmean
	for(int i=0;i<N;i++){
		for(int j=0;j<2;j++){
			revmean[i][j] = Y[i][j] - revmean[i][j];
		}
	}//for
	double m_mean[2*NN];
	memset(m_mean,0,sizeof(m_mean));
	for(int i=0;i<N;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<NN;k++){
				m[i][NN*j+k] = revmean[i][j]*Z[i][k];
			}
		}	
	}
	double S[2*NN][2*NN];
	for(int i=0;i<2*NN;i++){
		for(int j=0;j<2*NN;j++){
			double ss = 0;
			for(int k=0;k<N;k++){
				ss += m[k][i]*m[k][j];
			}
			S[i][j] = ss/(double)N;
		}
	}
	cout<<" @"<<S[0][0]<<"@"<<endl;
	ofstream fout("SS_a.txt");
	for(int i=0;i<2*NN;i++){
		for(int j=0;j<2*NN;j++){
			fout<<S[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
}
//***************************obj*******************************
double jstat_a(double *para){
	for(int i=0;i<N;i++){
		for(int j=0;j<RNUM;j++){
			F[i][j] = normcdf(rand1[0][j] + para[0]);
			G[i][j] = normcdf(rand1[1][j] + para[1]);
			q_bar[i][j] = (F[i][j]*q0[i])/(F[i][j]*q0[i] + G[i][j]*(1-q0[i]));
			q_ijr[i][j] = F[i][j] + q_bar[i][j]*(1-F[i][j]);
			q_ijb[i][j] = q_bar[i][j]*(1-G[i][j]);
		}
	}
	double rho = para[5];
	double rhox = rho/(1.0 - rho);
	double sum,sum1,sum2;
	for(int i=0;i<N;i++){
		sum = 0,sum1 = 0,sum2 = 0;
		for(int j=0;j<RNUM;j++){
			sum += pow(q_bar[i][j],rhox);
			sum1 += pow(q_ijr[i][j],rhox);
			sum2 += pow(q_ijb[i][j],rhox);
		}
		qij_mean[i] = sum/RNUM;
		qijh_mean[i] = sum1/RNUM;
		qijl_mean[i] = sum2/RNUM;
	}//for	
	// we got qij_mean, ci.
	double rhoy = 1.0/(rho - 1.0);
	double revmean[N][2];
	for(int i=0;i<N;i++){
		sum = 0;
		double temp[RNUM];
		for(int j=0;j<RNUM;j++){
			temp[j] = exp(-(rand2[i][j]));
		}
		if(color[i]==0){
			double ch[RNUM],cl[RNUM],dh[RNUM],dl[RNUM],paih,pail;
			double sum_paih = 0,sum_pail = 0;
			for(int j=0;j<RNUM;j++){
				ch[j] = para[3]*temp[j];
				cl[j] = para[2]*temp[j];
				dh[j] = para[4]*pow(ch[j]/rho,rhoy);
				dl[j] = para[4]*pow(cl[j]/rho,rhoy);
				paih = t[i]*(ch[j]/rho - ch[j])*dh[j]*qijh_mean[i] + (1-t[i])*(ch[j]/rho - ch[j])*dh[j]*qij_mean[i];
				pail = t[i]*(cl[j]/rho - cl[j])*dl[j]*qijl_mean[i] + (1-t[i])*(cl[j]/rho - cl[j])*dl[j]*qij_mean[i];
				sum_paih += paih;
				sum_pail += pail;
			}
			if(sum_paih>sum_pail){
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dh[j]*ch[j]/rho*qijh_mean[i];
					sum += rev[i][j];
				}
			}
			else{
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dl[j]*cl[j]/rho*qijl_mean[i];
					sum += rev[i][j];
				}
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==1){        //HIGH
			for(int j=0;j<RNUM;j++){
				double ch = para[3]*temp[j];
				double h = ch/rho;
				double dh = para[4]*pow(h,rhoy);
				rev[i][j] = dh*h*qijh_mean[i];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==2){//low
			for(int j=0;j<RNUM;j++){
				double cl = para[2]*temp[j];
				double l = cl/rho;
				double dl = para[4]*pow(l,rhoy);
				rev[i][j] = dl*l*qijl_mean[i];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
		}//if
	}//for
	double TotalR[N];
	memset(TotalR,0,sizeof(double)*N);
	for(int i=0;i<RNUM;i++){
		double sums[IN];
		memset(sums,0,sizeof(sums));
		for(int j=0;j<N;j++){
			sums[indu[j]] += rev[j][i];
		}
		for(int j=0;j<N;j++){
			TotalR[j] += rev[j][i]/sums[indu[j]]; 
		}
	}//for RNUM
	for(int i=0;i<N;i++){
		if(TotalR[i]<0){
			cout<<i<<"TotalR = "<<TotalR[i]<<endl;
		}
		revmean[i][1] = TotalR[i]/RNUM;
	}
	// now we get phat == revmean
	for(int i=0;i<N;i++){
		for(int j=0;j<2;j++){
			revmean[i][j] = Y[i][j] - revmean[i][j];
		}
	}//for
	double m_mean[2*NN];
	memset(m_mean,0,sizeof(m_mean));
	for(int i=0;i<N;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<NN;k++){
				m[i][NN*j+k] = revmean[i][j]*Z[i][k];
			}
		}	
	}
	for(int j=0;j<2*NN;j++){
        for(int i=0;i<N;i++){
	         m_mean[j] += m[i][j];
        }
	}
	for(int i=0;i<2*NN;i++){
		m_mean[i] = m_mean[i]/N;
	}
	double jstat_a = 0;
	double obj1[2*NN];
	memset(obj1,0,sizeof(obj1));
	for(int i=0;i<2*NN;i++){
		for(int j=0;j<2*NN;j++){
			obj1[i] += m_mean[j]*W[j][i]; 
		}
	}
	for(int i=0;i<2*NN;i++){
		jstat_a += obj1[i]*m_mean[i]; 
	}
	jstat_a *= N;
	return jstat_a;
}//jstat
