#include "ga_steps.h"
#include "mpi.h"
using namespace std;

int main(int argc,char**argv){
	int myid,numprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);  //current process id
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); //num of process
	MPI_Status status;
	srand((unsigned int)time(NULL));
	if(myid==0){
		cout<<numprocs<<endl;
		double r1 = (double)rand()/RAND_MAX;
		double r2 = (double)rand()/RAND_MAX;
		if(r1 == r2){
			cout<<"random error"<<endl;
			return -1;
		} 
	}//assert
	double **sendbuff = NULL;
	double *obj_mpi = NULL;
	double *slaverecv = NULL;
	double *obj_slave = NULL;
	double **sends = NULL;
	double *recvs = NULL;
	double *obj_slaves = NULL;
	double *obj_mpis = NULL;
	int n=5;
	int times = 2000;
	double beta = 0;
	int ll = P/(numprocs -1);
	readData();
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0){
		sendbuff = new double*[numprocs];
		sends = new double*[numprocs];
		obj_mpi = new double[ll];
		obj_mpis = new double[ll];
		for(int j=0;j<numprocs;j++){
			sendbuff[j] = new double[ll*PNUM];
			sends[j] = new double[ll*PNUM];
		}
		cout<<"before init"<<endl;
		init();
		popinit();
		encoding();
		for(int x = 1;x<numprocs;x++){
			for(int y=0;y<ll;y++){
				for(int k=0;k<PNUM;k++){
					sends[x][y*PNUM+k] = paras[(x-1)*ll+y][k];
				}
			}
		}//for
		for(int x=1;x<numprocs;x++){
			MPI_Send(sends[x],ll*PNUM,MPI_DOUBLE,x,99,MPI_COMM_WORLD);
		}
	}else{
		recvs = new double[ll*PNUM];
		obj_slaves = new double[ll];
		MPI_Recv(recvs,ll*PNUM,MPI_DOUBLE,0,99,MPI_COMM_WORLD,&status);
		for(int i=0;i<ll;i++){
			double temps[PNUM];
			for(int j=0;j<PNUM;j++){
				temps[j] = recvs[i*PNUM+j];
			}
			obj_slaves[i] = objective(temps);
		}
		MPI_Send(obj_slaves,ll,MPI_DOUBLE,0,98,MPI_COMM_WORLD);
	}
	if(myid==0){
		for(int x=1;x<numprocs;x++){
			MPI_Recv(obj_mpis,ll,MPI_DOUBLE,x,98,MPI_COMM_WORLD,&status);
			for(int y=0;y<ll;y++){
				obj[ll*(x-1)+y] = obj_mpis[y];
				//cout<<"obj "<<P+ll*(x-1)+y<<" "<<obj[P+ll*(x-1)+y]<<endl;
			}
		}
	}//if
	for(int j=0;j<times;j++){
		if(myid==0){
			mutation();
			variation();
			decoding();
			for(int x = 1;x<numprocs;x++){
				for(int y=0;y<ll;y++){
					for(int k=0;k<PNUM;k++){
						sendbuff[x][y*PNUM+k] = paras2[P+(x-1)*ll+y][k];
					}
				}
			}//for
			for(int x=1;x<numprocs;x++){
				MPI_Send(sendbuff[x],ll*PNUM,MPI_DOUBLE,x,97,MPI_COMM_WORLD);
			}
		}else{
			slaverecv = new double[ll*PNUM];
			obj_slave = new double[ll];
			MPI_Recv(slaverecv,ll*PNUM,MPI_DOUBLE,0,97,MPI_COMM_WORLD,&status);
			for(int i=0;i<ll;i++){
				double temps[PNUM];
				for(int j=0;j<PNUM;j++){
					temps[j] = slaverecv[i*PNUM+j];
				}
				obj_slave[i] = objective(temps);
			}
			MPI_Send(obj_slave,ll,MPI_DOUBLE,0,96,MPI_COMM_WORLD);
		}
		if(myid==0){
			for(int x=1;x<numprocs;x++){
				MPI_Recv(obj_mpi,ll,MPI_DOUBLE,x,96,MPI_COMM_WORLD,&status);
				for(int y=0;y<ll;y++){
					obj[P+ll*(x-1)+y] = obj_mpi[y];
				}
			}
			choose();
			cout<<j<<" best = "<<best<<endl;
			fout<<j<<" "<<best<<" ";
			for(int i=0;i<PNUM;i++){
				fout<<paras[0][i]<<" ";
			}
			fout<<endl;
			if(best<=beta){
				return 0;
			}
		}//if
	}//times
	fout.close();
//	cout<<"best "<<best<<endl;
	MPI_Finalize();
	return 0;
}
