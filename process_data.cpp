#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
const int N = 3189;
const int XNUM = 21;
double X[N][XNUM];
double Y[N][2];
int color[N];
int industry[N];
int wave[N][2];
double t[10],q0[10];
int main(int argc,char**argv){
	memset(color,0,sizeof(color));
	ifstream fin("131-14.csv");
	int temp_int,red,black;
	double temp_double;
	string line;
	getline(fin,line);
	for(int i = 0;i < N;i ++){
		int index = 0;
		fin>>temp_int;
		fin>>X[i][index++]>>X[i][index++]>>temp_int;
		fin>>Y[i][0]>>Y[i][1];
		fin>>temp_int>>temp_double;   //lnlabor09
		fin>>temp_double;
		X[i][index++] = temp_double;
		fin>>temp_double;
		X[i][index++] = log(temp_double);
		fin>>temp_int;
		fin>>wave[i][0]>>wave[i][1];
		fin>>red>>black;
		if(red ==1){
			color[i] = 1;
		}else if(black == 1){
			color[i] = 2;
		}
		for(int j=0;j<14;j++){
			fin>>X[i][index ++];
		}//for
		for(int j=0;j<30;j++){
			fin>>temp_int;
		}
		industry[i] = 0;
		for(int j=0;j<9;j++){
			if(X[i][9+j] == 1){
				industry[i] = j+1;
				break;
			}
		}//for
	}//for i
	int rednum[10],blacknum[10],insnum[10];
	memset(rednum,0,sizeof(rednum));
	memset(blacknum,0,sizeof(blacknum));
	memset(insnum,0,sizeof(insnum));
	for(int i=0;i<N;i++){
		insnum[industry[i]] ++;
		if(color[i] == 1){
			rednum[industry[i]] ++;
		}else if(color[i] == 2){
			blacknum[industry[i]] ++;
		}
	}
	for(int i=0;i<10;i++){
		t[i] = (double)(rednum[i]+blacknum[i])/(double)insnum[i];
		q0[i] = (double)rednum[i]/(rednum[i]+blacknum[i]);
	}
	ofstream fout_q0("para_q0.txt");
	for(int i=0;i<10;i++){
		fout_q0<<q0[i]<<endl;
	}
	fout_q0.close();
	//para X:east middle lnasset lnyear regsoe regforeign regprivate regcollective reglegal dindus1~dindus9
	//go to fileX:east middle dis1~9 lnasset lnyear regsoe~reglegal lnasset2 lnyear2 1
	ofstream fout_X("para_x.txt");
	for(int i = 0;i < N;i ++){
		fout_X<<X[i][0]<<" "<<X[i][1]<<" ";
		for(int j = 9;j<18;j++){
			fout_X<<X[i][j]<<" ";
		}
		fout_X<<X[i][2]<<" "<<X[i][3]<<" ";
		for(int j = 4;j < 9;j ++){
			fout_X<<X[i][j]<<" ";
		}
		fout_X<<pow(X[i][2],2)<<" "<<pow(X[i][3],2)<<" "<<1<<endl;
	}//for X
	fout_X.close();
	cout<<"X done"<<endl;
	//para Z:fileX(21) wave1 wave2 q0 1
	ofstream fout_Z("para_z.txt");
	for(int i = 0;i < N;i ++){
		fout_Z<<X[i][0]<<" "<<X[i][1]<<" ";
		for(int j = 9;j<18;j++){
			fout_Z<<X[i][j]<<" ";
		}
		fout_Z<<X[i][2]<<" "<<X[i][3]<<" ";
		for(int j = 4;j < 9;j ++){
			fout_Z<<X[i][j]<<" ";
		}
		fout_Z<<pow(X[i][2],2)<<" "<<pow(X[i][3],2)<<" ";
		fout_Z<<wave[i][0]<<" "<<wave[i][1]<<" "<<q0[industry[i]]<<" "<<1<<endl;
	}//for Z
	fout_Z.close();
	cout<<"Y done"<<endl;
	//para Y:lnasset share color industry q0 t
	ofstream fout_Y("para_y.txt");
	for(int i=0;i<N;i++){
		fout_Y<<Y[i][0]<<" "<<Y[i][1]<<" "<<color[i]<<" "<<industry[i]<<" "<<q0[industry[i]]<<" "<<t[industry[i]]<<endl;
	}//for
	fout_Y.close();
	cout<<"Z done"<<endl;
	return 0;
}
