#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

const int Q = 24*2;
double WW[Q][Q];
double WW_hat[Q][Q];
double sum[Q][Q];

void readData(){
	ifstream fin("WW");
	for(int i=0;i<Q;i++){
		for(int j=0;j<Q;j++){
			fin>>WW[i][j];
		}
	}//for
	fin.close();
}

void check(){
	for(int i=0;i<Q;i++){
		for(int j=0;j<Q;j++){
			sum[i][j] = 0;
			for(int k=0;k<Q;k++){
				sum[i][j] += WW[i][k]*WW_hat[j][k];
			}
		}
	}
	ofstream fout("II");
	for(int i=0;i<Q;i++){
		for(int j=0;j<Q;j++){
			fout<<sum[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
}

bool Gauss(double a[][Q],double b[][Q])     //求矩阵的逆
{
	double max,temp;
	double t[Q][Q];
	for(int i=0;i<Q;i++)
	{
		for(int j=0;j<Q;j++)
		{
			t[i][j] = a[i][j];
		}
	}//将矩阵a存储在t
	for(int i=0;i<Q;i++)
	{
		for(int j=0;j<Q;j++)
		{
			b[i][j] = (i==j)?1.0:0;
		}
	}// set matrix B I
	for(int i=0;i<Q;i++)
	{
		max = t[i][i];
		int k = i;
		for(int j=i+1;j<Q;j++)
		{
			if(fabs(t[j][i])>fabs(max))
			{
				max = t[j][i];
				k = j;
			}
		}
		if(k!=i)
		{
			for(int j=0;j<Q;j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;
				temp = b[i][j];
				b[i][j] = b[k][j];
				b[k][j] = temp;
			}
		}//if
		if(t[i][i] == 0)
		{
			cout<<"no exist"<<endl;
			return false;
		}
		temp = t[i][i];
		for(int j=0;j<Q;j++)
		{
			t[i][j] = t[i][j]/temp;
			b[i][j] = b[i][j]/temp;
		}
		for(int j=0;j<Q;j++)
		{
			if(j!=i)
			{
				temp = t[j][i];
				for(k=0;k<Q;k++)
				{
					t[j][k] = t[j][k] - t[i][k]*temp;
					b[j][k] = b[j][k] - b[i][k]*temp;
				}
			}
		}
	}// big for
	return true;
}
