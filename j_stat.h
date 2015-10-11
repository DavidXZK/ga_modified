#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

using namespace std;

const double PI = 3.1415926;
const int N = 3189;     //num of sample
const int RNUM = 250;   //random num
const int XNUM = 21;    //theta para num
//para: mur mub cl ch k rou 6 theta 21
const int CNUM = 6;
const int PNUM = CNUM+XNUM;    //para num
const int IN = 10;       //industry num
const int NN = 24;     //Z num
const int FLAG = 0;
int color[N];    //black or white 0 1 2
int indu[N];     //belong to which industry
double X[N][XNUM];  //X变量
double Y[N][2],Z[N][NN];
double m[N][2*NN];
double t[N];     // parameter t 抽检比例
double rand1[2][RNUM],rand2[N][RNUM];
double qij_mean[IN],qijh_mean[IN],qijl_mean[IN];
double q0[N],W[2*NN][2*NN],rev[N][RNUM];
double F[IN][RNUM],G[IN][RNUM];
double q_bar[IN][RNUM],q_ijr[IN][RNUM],q_ijb[IN][RNUM],qr0[IN];
double para[PNUM];
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
	cout<<"jstat() in"<<endl;
	ifstream fw("WW");
	for(int i = 0;i < 2*NN;i ++){
		for(int j = 0;j < 2*NN;j ++){
			fw>>W[i][j];
		}
	}//for
	fw.close();
	ifstream fin_q0("para_q0.txt");
	for(int i=0;i<IN;i++){
		fin_q0>>qr0[i];
	}
	fin_q0.close();
	ifstream fin_x("para_x.txt");    //para x
	ifstream fin_y("para_y.txt");   //para lny0 y1 color industry q0 t
	ifstream fin_z("para_z.txt");   //para x + wave1 + wave2 + q0 + 1
	for(int i=0;i<N;i++){
		for(int j = 0;j < XNUM;j ++){
			fin_x>>X[i][j];
		}
		for(int j = 0;j < NN;j ++){
			fin_z>>Z[i][j];
		}
		fin_y>>Y[i][0]>>Y[i][1]>>color[i]>>indu[i]>>q0[i]>>t[i];
	}//for
	fin_x.close();
	fin_y.close();
	fin_z.close();
	//读取随机正态，每次一致
	ifstream frand("randnums");
	for(int i=0;i<2;i++){
		for(int j=0;j<RNUM;j++){
			frand>>rand1[i][j];
		}
	}
	for(int i=0;i<N;i++){
		for(int j=0;j<RNUM;j++){
			frand>>rand2[i][j];
		}
	}
	frand.close();
	ifstream fin_para("para.txt");
	for(int i=0;i<PNUM;i++){
		fin_para>>para[i];
	}
	cout<<"jstat() out"<<endl;
}
//***************************getS*******************************
double getS(double *para){
	for(int i=0;i<IN;i++){
		for(int j=0;j<RNUM;j++){
			F[i][j] = normcdf(rand1[0][j] + para[0]);
			G[i][j] = normcdf(rand1[1][j] + para[1]);
			q_bar[i][j] = (F[i][j]*qr0[i])/(F[i][j]*qr0[i] + G[i][j]*(1-qr0[i]));
			q_ijr[i][j] = F[i][j] + q_bar[i][j]*(1-F[i][j]);
			q_ijb[i][j] = q_bar[i][j]*(1-G[i][j]);
		}
	}
	double rho = para[5];
	double rhox = rho/(1.0 - rho);
	double sum,sum1,sum2;
	for(int i=0;i<IN;i++){
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
		double xx = 0;
		for(int j=0;j<XNUM;j++){
			xx += para[CNUM+j]*X[i][j];
		}
		sum = 0;
		if(color[i]==0){
			double ch[RNUM],cl[RNUM],dh[RNUM],dl[RNUM],paih,pail;
			double sum_paih = 0,sum_pail = 0;
			for(int j=0;j<RNUM;j++){
				double temp = exp(-(xx+rand2[i][j]));
				ch[j] = para[3]*temp;
				cl[j] = para[2]*temp;
				dh[j] = para[4]*pow(ch[j]/rho,rhoy);
				dl[j] = para[4]*pow(cl[j]/rho,rhoy);
				paih = t[i]*(ch[j]/rho - ch[j])*dh[j]*qijh_mean[indu[i]] + (1-t[i])*(ch[j]/rho - ch[j])*dh[j]*qij_mean[indu[i]];
				pail = t[i]*(cl[j]/rho - cl[j])*dl[j]*qijl_mean[indu[i]] + (1-t[i])*(cl[j]/rho - cl[j])*dl[j]*qij_mean[indu[i]];
				sum_paih += paih;
				sum_pail += pail;
			}
			if(sum_paih>sum_pail){
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dh[j]*ch[j]/rho*qijh_mean[indu[i]];
					sum += rev[i][j];
				}
			}
			else{
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dl[j]*cl[j]/rho*qijl_mean[indu[i]];
					sum += rev[i][j];
				}
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==1){        //HIGH
			for(int j=0;j<RNUM;j++){
				double temp = exp(-(xx+rand2[i][j]));
				double ch = para[3]*temp;
				double h = ch/rho;
				double dh = para[4]*pow(h,rhoy);
				rev[i][j] = dh*h*qijh_mean[indu[i]];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==2){//low
			for(int j=0;j<RNUM;j++){
				double temp = exp(-(xx+rand2[i][j]));
				double cl = para[2]*temp;
				double l = cl/rho;
				double dl = para[4]*pow(l,rhoy);
				rev[i][j] = dl*l*qijl_mean[indu[i]];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
		}//if
	}//for
	double sums[IN];
	memset(sums,0,sizeof(sums));
	for(int i = 0;i < N;i ++){
		sums[indu[i]] += revmean[i][0];
	}
	for(int i = 0;i < N;i ++){
		revmean[i][1] = revmean[i][0] / sums[indu[i]];
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
	ofstream fout("SS");
	for(int i=0;i<2*NN;i++){
		for(int j=0;j<2*NN;j++){
			fout<<S[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
}//getS
//************************************
double jstat(double *para){
	for(int i=0;i<IN;i++){
		for(int j=0;j<RNUM;j++){
			F[i][j] = normcdf(rand1[0][j] + para[0]);
			G[i][j] = normcdf(rand1[1][j] + para[1]);
			q_bar[i][j] = (F[i][j]*qr0[i])/(F[i][j]*qr0[i] + G[i][j]*(1-qr0[i]));
			q_ijr[i][j] = F[i][j] + q_bar[i][j]*(1-F[i][j]);
			q_ijb[i][j] = q_bar[i][j]*(1-G[i][j]);
		}
	}
	double rho = para[5];
	double rhox = rho/(1.0 - rho);
	double sum,sum1,sum2;
	for(int i=0;i<IN;i++){
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
		double xx = 0;
		for(int j=0;j<XNUM;j++){
			xx += para[CNUM+j]*X[i][j];
		}
		sum = 0;
		if(color[i]==0){
			double ch[RNUM],cl[RNUM],dh[RNUM],dl[RNUM],paih,pail;
			double sum_paih = 0,sum_pail = 0;
			for(int j=0;j<RNUM;j++){
				double temp = exp(-(xx+rand2[i][j]));
				ch[j] = para[3]*temp;
				cl[j] = para[2]*temp;
				dh[j] = para[4]*pow(ch[j]/rho,rhoy);
				dl[j] = para[4]*pow(cl[j]/rho,rhoy);
				paih = t[i]*(ch[j]/rho - ch[j])*dh[j]*qijh_mean[indu[i]] + (1-t[i])*(ch[j]/rho - ch[j])*dh[j]*qij_mean[indu[i]];
				pail = t[i]*(cl[j]/rho - cl[j])*dl[j]*qijl_mean[indu[i]] + (1-t[i])*(cl[j]/rho - cl[j])*dl[j]*qij_mean[indu[i]];
				sum_paih += paih;
				sum_pail += pail;
			}
			if(sum_paih>sum_pail){
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dh[j]*ch[j]/rho*qijh_mean[indu[i]];
					sum += rev[i][j];
				}
			}
			else{
				for(int j=0;j<RNUM;j++){
					rev[i][j] = dl[j]*cl[j]/rho*qijl_mean[indu[i]];
					sum += rev[i][j];
				}
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==1){        //HIGH
			for(int j=0;j<RNUM;j++){
				double temp = exp(-(xx+rand2[i][j]));
				double ch = para[3]*temp;
				double h = ch/rho;
				double dh = para[4]*pow(h,rhoy);
				rev[i][j] = dh*h*qijh_mean[indu[i]];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
			continue;
		}
		if(color[i]==2){//low
			for(int j=0;j<RNUM;j++){
				double temp = exp(-(xx+rand2[i][j]));
				double cl = para[2]*temp;
				double l = cl/rho;
				double dl = para[4]*pow(l,rhoy);
				rev[i][j] = dl*l*qijl_mean[indu[i]];
				sum += rev[i][j];
			}
			revmean[i][0] = sum/RNUM;
		}//if
	}//for
	double sums[IN];
	memset(sums,0,sizeof(sums));
	for(int i=0;i<N;i++){
		sums[indu[i]] += revmean[i][0];
	}
	for(int i=0;i<N;i++){
		revmean[i][1] = revmean[i][0] / sums[indu[i]];
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
	double jstat = 0;
	double obj1[2*NN];
	memset(obj1,0,sizeof(obj1));
	for(int i=0;i<2*NN;i++){
		for(int j=0;j<2*NN;j++){
			obj1[i] += m_mean[j]*W[j][i]; 
		}
	}
	for(int i=0;i<2*NN;i++){
		jstat += obj1[i]*m_mean[i]; 
	}
	jstat *= N;
	return jstat;
}//jstat