#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;
const int N = 73379;

int indu[] = {13,14,15,18,19,21,22,24,26,27,29,30,31,34,37,39,40,42};
double q0[18];
int howmany[18],redmany[18],lmany[18];

int main(int argc,char** argv)
{
	ifstream fin("data_consumer_end");
	ofstream fout("XX");
	ofstream fout1("YY");
	ofstream fout2("ZZ");
	ofstream fout3("q0");
	ofstream fout4("bound");
	ofstream fout5("t");
	string s;
	getline(fin,s);
	cout<<s<<endl;
	double a[N][10];double b[2];int red,black;
	double temp;
	double max = 0;
	for(int i=0;i<N;i++)
	{
		fin>>a[i][0]>>a[i][1];

		fin>>b[0]>>b[1];
		fin>>a[i][2];
		fin>>temp;
		fin>>a[i][3]>>a[i][4];
		fin>>temp;
		a[i][3] = log(a[i][3]);
		if(a[i][3]>max)
			max = a[i][3];
		a[i][4] = log(a[i][4]);
		//cout<<a[i][0]<<" "<<a[i][1]<< " "<<a[i][3]<<" "<<a[i][4]<<endl;
		for(int j=5;j<10;j++)
		{
			fin>>a[i][j];
		}
		fin>>red>>black;
		int which = 0;
		fout<<a[i][0]<<" "<<a[i][1]<<" ";
		for(int j=0;j<18;j++)
		{
			if(a[i][2]==indu[j])
			{
				//cout<<a[i][2];
			//	cin>>temp;
				which = j+1;
				fout<<1<<" ";
				lmany[j] ++;
				if(red+black==1)
				{
					howmany[j] ++;
				}
				if(red==1)
				{
					redmany[j] ++;
				}
			}
			else
			{
				fout<<0<<" ";
			}
		}//for
		fout<<a[i][3]<<" "<<a[i][4]<<" ";
		for(int j=5;j<10;j++)
		{
			fout<<a[i][j]<<" ";
		}
		fout<<pow(a[i][3],2)<<" "<<pow(a[i][4],2)<<" ";
		int d  = 0;
		if(red==1)
		{
			d=1;
		}
		if(black==1)
		{
			d=2;
		}
		fout<<d<<" "<<which<<endl;
		fout1<<b[0]/10000<<" "<<b[1]<<endl;
	}
	fin.close();
	fout1.close();
	cout<<"max = "<<max<<endl;
	for(int i=0;i<18;i++)
	{

		q0[i] = (double)redmany[i]/(double)howmany[i];
	}
	for(int i=0;i<N;i++)
	{
		fout2<<a[i][4]<<" "<<pow(a[i][4],2)<<" "<<a[i][3]<<" "<<pow(a[i][3],2)<<" "<<a[i][0]<<" "<<a[i][1]<<" ";
		for(int j=5;j<10;j++)
		{
			fout2<<a[i][j]<<" ";
		}
		int f = 17;
		for(int j=0;j<17;j++)
		{
			if(a[i][2]==indu[j])
			{
				fout2<<1<<" ";
				f = j;
			}
			else
			{
				fout2<<0<<" ";
			}
		}
		fout2<<q0[f]<<" ";
		fout2<<1<<endl;
	}
	fout2.close();
	for(int i=0;i<18;i++)
	{
		fout3<<q0[i]<<endl;
	}
	fout3.close();
	int sum = 0;
	fout4<<1<<endl;
	for(int i=0;i<18;i++)
	{
		sum += lmany[i];
		fout4<<sum<<endl;
	}
	fout4.close();
	for(int i=0;i<18;i++)
	{
		double aa = (double)howmany[i]/(double)lmany[i];
		for(int j=0;j<lmany[i];j++)
		{
			fout5<< aa<<endl;
		}
	}
	fout5.close();
	return 0;

}