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
const int CNUM = 7;
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
double para_range[PNUM][2];
double paras[P][PNUM];
double paras2[2*P][PNUM];
int codeX[P][PNUM][M];
int codeX1[P][PNUM][M];
int codeX2[P][PNUM][M];

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
//********************************quick sort*****************************
void sort(double a[],int b[],int zz,int yy){    //快排,b记录的是顺序
	int z,y,zs;
	double k;
	if(zz<yy){
		z = zz;
		zs = b[zz];
		y = yy;
		k = a[z];
		do{
			while((z<y)&&(a[y]>=k)) y--;
			if(z<y){
				a[z] = a[y];
				b[z] = b[y];
				z ++;
			}
			while((z<y)&&(a[z]<=k)) z++;
			if(z < y){
				a[y] = a[z];
				b[y] = b[z];
			}
		}while(z!=y);
		a[z] = k;
		b[z] = zs;
		sort(a,b,zz,z-1);
		sort(a,b,z+1,yy);
	}
}
//***************************readdata**************************
void readData(){
	cout<<"GA() in"<<endl;
	if(FLAG == 0){
		for(int i = 0;i < 2*NN;i ++){
			for(int j = 0;j < 2*NN;j ++){
				W[i][j] = (i==j)?1:0;
			}
		}
	}else{
		ifstream fw("WW");
		for(int i = 0;i < 2*NN;i ++){
			for(int j = 0;j < 2*NN;j ++){
				fw>>W[i][j];
			}
		}//for
		fw.close();
	}
	//mur mub cl ch k rho east middle dindus1~9 lnasset lnyear
	//regsoe regforeign regprivate regcollective reglegal lna2 lny2 1
	//-3 10
	ifstream fin_range("para_range.txt");
	for(int i=0;i<PNUM;i++){
		fin_range>>para_range[i][0]>>para_range[i][1];
	}
	fin_range.close();
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
	}//for
	//fin_x.close();
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
	cout<<"GA() out"<<endl;
}
//***************************obj*******************************
double objective(double *para){
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
				//double temp = exp(-(rand2[i][j]));
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
			}else{
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
				//double temp = exp(-(rand2[i][j]));
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
				//double temp = exp(-(rand2[i][j]));
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
	double objs = 0;
	double obj1[2*NN];
	memset(obj1,0,sizeof(obj1));
	for(int i=0;i<2*NN;i++){
		for(int j=0;j<2*NN;j++){
			obj1[i] += m_mean[j]*W[j][i]; 
		}
	}
	for(int i=0;i<2*NN;i++){
		objs += obj1[i]*m_mean[i]; 
	}
	objs *= N;
	return objs;
}//obj
//************************************
void init(){
	memset(codeX,0,sizeof(int)*P*PNUM*M);
	memset(codeX1,0,sizeof(int)*P*PNUM*M);
	memset(codeX2,0,sizeof(int)*P*PNUM*M);
}//array initial

void popinit(){
	//first 6 paras we don't have the initials,so random
	for(int i=0;i<P;i++){
		for(int j=0;j<PNUM;j++){
			if(j==3){
				paras[i][j] = paras[i][2] + (para_range[j][1] - paras[i][2])*((double)rand()/RAND_MAX);
			}else{
				paras[i][j] = para_range[j][0] + (para_range[j][1] - para_range[j][0])*((double)rand()/RAND_MAX);
			}
		}
	}
	//the left 21(X)
	/*ifstream fin("para_initial.txt");
	for(int j=6;j<PNUM;j++){
		fin>>paras[0][j];
	}
	for(int i=1;i<HH;i++){
		for(int j=6;j<PNUM;j++){
			paras[i][j] = paras[0][j];
		}
	}
	fin.close();*/
}//popinit

void encoding(){//编码
	cout<<"encoding() in"<<endl;
	for(int i=0;i<P;i++){
		for(int j=0;j<PNUM;j++){
			int y = (int)((paras[i][j] - para_range[j][0])/(para_range[j][1] - para_range[j][0])*pow(2.0,M));
			int count = 0;
			while(true){
				int aa = y%2;
				y = y/2;
				codeX[i][j][count] = aa;
				count ++;
				if(y==0)
				break;
			}
		}
	}
	cout<<"encoding() out"<<endl;
}
void mutation(){
	cout<<"mu() in"<<endl;
	int flag[P];
	memset(flag,0,sizeof(flag));
	int cc = P-1;
	for(int i=0;i<P;i++){
		if(flag[i]==1) continue;
		int which = rand()%cc;//交配哪一个
		int count = -1;
		for(int j=i+1;j<P;j++){
			if(flag[j]==0){
				count ++;
			}
			else{
				continue;
			}
			if(count==which){
			    flag[j] = 1;
			    cc--;
			    for(int l =0;l<PNUM;l++){
			    	double r = (double)rand()/RAND_MAX;
			    	if(r<mp){//发生交配
						int kk = rand()%M;
						if(l==2){
							double cll1 = 0,cll2 = 0,chh1 = 0,chh2 = 0;
							for(int k=0;k<M;k++){
								if(k<=kk){
									if(codeX[j][l][k]==1){
										cll1 += pow(2.0,k);
									}
									if(codeX[i][l][k]==1){
										cll2 += pow(2.0,k);
									}
								}
								else{
									if(codeX[i][l][k]==1){
										cll1 += pow(2.0,k);
									}
									if(codeX[j][l][k]==1){
										cll2 += pow(2.0,k);
									}
								}
								if(codeX[i][3][k]==1){
									chh1 += pow(2.0,k);
								}
								if(codeX[j][3][k]==1){
									chh2 += pow(2.0,k);
								}
							}
							cll1 = cll1*(para_range[2][1] - para_range[2][0])/pow(2.0,M)+para_range[2][0];
							cll2 = cll2*(para_range[2][1] - para_range[2][0])/pow(2.0,M)+para_range[2][0];
							chh1 = chh1*(para_range[3][1] - para_range[3][0])/pow(2.0,M)+para_range[3][0];
							chh2 = chh2*(para_range[3][1] - para_range[3][0])/pow(2.0,M)+para_range[3][0];
							if((chh1>cll1)&&(chh2>cll2)){
								for(int k=0;k<M;k++){
									if(k<=kk){
										//前段交换
										codeX1[i][l][k] = codeX[j][l][k];
										codeX1[j][l][k] = codeX[i][l][k];
									}
									if(k>kk){
										codeX1[i][l][k] = codeX[i][l][k];
										codeX1[j][l][k] = codeX[j][l][k];
									}	
								}
							}
						}//if l==3
						else if(l==3){
							double cll1 = 0,cll2 = 0,chh1 = 0,chh2 = 0;
							for(int k=0;k<M;k++){
								if(k<=kk){
									if(codeX[j][l][k]==1){
										chh1 += pow(2.0,k);
									}
									if(codeX[i][l][k]==1){
										chh2 += pow(2.0,k);
									}
								}
								else{
									if(codeX[i][l][k]==1){
										chh1 += pow(2.0,k);
									}
									if(codeX[j][l][k]==1){
										chh2 += pow(2.0,k);
									}
								}
								if(codeX[i][2][k]==1){
									cll1 += pow(2.0,k);
								}
								if(codeX[j][2][k]==1){
									cll2 += pow(2.0,k);
								}
							}
							cll1 = cll1*(para_range[2][1] - para_range[2][0])/pow(2.0,M)+para_range[2][0];
							cll2 = cll2*(para_range[2][1] - para_range[2][0])/pow(2.0,M)+para_range[2][0];
							chh1 = chh1*(para_range[3][1] - para_range[3][0])/pow(2.0,M)+para_range[3][0];
							chh2 = chh2*(para_range[3][1] - para_range[3][0])/pow(2.0,M)+para_range[3][0];
							if((chh1>cll1)&&(chh2>cll2)){
								for(int k=0;k<M;k++){
									if(k<=kk){
										//前段交换
										codeX1[i][l][k] = codeX[j][l][k];
										codeX1[j][l][k] = codeX[i][l][k];
									}
									if(k>kk){
										codeX1[i][l][k] = codeX[i][l][k];
										codeX1[j][l][k] = codeX[j][l][k];
									}	
								}
							}
						}//if l==3
						else{
							for(int k=0;k<M;k++){
								if(k<=kk){
									//前段交换
									codeX1[i][l][k] = codeX[j][l][k];
									codeX1[j][l][k] = codeX[i][l][k];
								}
								if(k>kk){
									codeX1[i][l][k] = codeX[i][l][k];
									codeX1[j][l][k] = codeX[j][l][k];
								}	
							}
						}
					} 
					else{
						for(int k=0;k<M;k++){
							codeX1[i][l][k] = codeX[i][l][k];
							codeX1[j][l][k] = codeX[j][l][k];
						}
					}//if else mur
				}//for
			}//if
		}//for
	}//for
	cout<<"mu() out"<<endl;
}

void variation(){
	cout<<"va() in"<<endl;
	for(int i=0;i<P;i++){
		for(int j=0;j<PNUM;j++){
			double r = (double)rand()/RAND_MAX;
			if(r<vp){
				int k = rand()%M;
				if(codeX1[i][j][k]==0){
					codeX1[i][j][k] = 1;
				}else{
					codeX1[i][j][k] = 0;
				}
				if(j==2){
					double cll = 0,chh = 0;
					for(int l=0;l<M;l++){
						if(codeX1[i][2][l]==1){
							cll += pow(2.0,l);
						}
						if(codeX1[i][3][l]==1){
							chh += pow(2.0,l);
						}
					}
					cll = cll*(para_range[2][1] - para_range[2][0])/pow(2.0,M)+para_range[2][0];
					chh = chh*(para_range[3][1] - para_range[3][0])/pow(2.0,M)+para_range[3][0];
					if(cll>=chh){
						if(codeX1[i][j][k]==0){
							codeX1[i][j][k] = 1;
						}else{
							codeX1[i][j][k] = 0;
						}
					}
				}
				if(j==3){
					double cll = 0,chh = 0;
					for(int l=0;l<M;l++){
						if(codeX1[i][2][l]==1){
							cll += pow(2.0,l);
						}
						if(codeX1[i][3][l]==1){
							chh += pow(2.0,l);
						}
					}
					cll = cll*(para_range[2][1] - para_range[2][0])/pow(2.0,M)+para_range[2][0];
					chh = chh*(para_range[3][1] - para_range[3][0])/pow(2.0,M)+para_range[3][0];
					if(cll>=chh){
						if(codeX1[i][j][k]==0){
							codeX1[i][j][k] = 1;
						}else{
							codeX1[i][j][k] = 0;
						}
					}
				}//if
			}//if vp
		}
	}
	cout<<"va() out"<<endl;
}

void decoding(){
	cout<<"de() in"<<endl;
	for(int i=0;i<P;i++){
		for(int j=0;j<PNUM;j++){
			paras2[i][j] = paras[i][j];
			double mm = 0;
			for(int k=0;k<M;k++){
				if(codeX1[i][j][k]==1){
					mm += pow(2.0,k);
				}
			}
			mm = mm*(para_range[j][1] - para_range[j][0])/pow(2.0,M)+para_range[j][0];
			paras2[P+i][j] = mm;
		}
	}
	cout<<"de() out"<<endl;
}

void choose(){
	//选择最好的P个
	int order[2*P];
	for(int i=0;i<2*P;i++){
		order[i] = i;
	}
	sort(obj,order,0,2*P-1);
	best = obj[0];
	cout<<"best = "<<best<<endl;
	for(int i=0;i<P;i++){
		for(int j=0;j<PNUM;j++){
			paras[i][j] = paras2[order[i]][j];
			if(order[i]<P){
				for(int l=0;l<M;l++){
					codeX2[i][j][l] = codeX[order[i]][j][l];
				}
			}
			else{
				for(int l=0;l<M;l++){
					codeX2[i][j][l] = codeX1[order[i] - P][j][l];
				}
			}
		}
	}//for
	for(int i=0;i<P;i++){
		for(int j=0;j<PNUM;j++){
			for(int l=0;l<M;l++){
				codeX[i][j][l] = codeX2[i][j][l];
			}
		}
	}//for	 
}// choose
