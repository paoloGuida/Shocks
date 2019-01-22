/* 
 * File:   advection_F.cpp
 * Author: Paolo Guida
 */
# include <cstdlib>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <sstream>
using namespace std;
void write ( string file, int Nx, int tstep, double **u,double dx,double dy );
int main()
{
	int i;
	int j;
	int k;
	int Nx;

int tstep;
	cout<<"set the number of grid points ";
	cin>>Nx;
	cout<<"set the number of tsteps points ";
	cin>>tstep;
	double Lx=2;
        double time=1;
        
	double dx=Lx/(Nx);
	double dt=time/(tstep);
        double alpha=1;
	double **u;
        double *x;
	double *u0;
	double *u_pred;
	double *u_corr;
	int number;
	double s=alpha*dt/dx;
	cout<<"type 1 for the step function as u0 or 2 for the gaussian ";
	cin>>number;
        x=new double[Nx+1];
	u=new double*[Nx+1];
	u0=new double[Nx+1];
	u_pred=new double[Nx+1];
	u_corr=new double[Nx+1];
	cout<<s<<endl;
	for (i=0;i<=Nx;i++)
	{
		u[i]=new double[tstep];
	}
        for (i=0;i<=Nx;i++)
	{
		x[i]=i*dx;
	}

	if (number==1)
	{
        for (i=0;i<=Nx;i++)
        {
		
		if (x[i]<=0.25)
		{
			u0[i]=1;
		}
		else
			{u0[i]=0;}

                cout<<i<<" "<<x[i]<<" "<<u0[i]<<endl;
        }
        }
        else
	{
        for (i=0;i<=Nx;i++)
        {
		u0[i]=exp(-200*pow(x[i]-.5,2));
        }
	}
	for (i=0;i<=Nx;i++)
        {
		u[i][0]=u0[i]; 

	}
 	ofstream output;
	for (j=1;j<tstep;j++){/*
u[0][j+1]=-1/4*s*(s-1)*u[Nx-1][j]+1/4*s*(5-s)*u[Nx][j]-1/4*(s-1)*(s+4)*u[0][j]+1/4*s*(s-1)*u[1][j];
u[1][j+1]=-1/4*s*(s-1)*u[Nx][j]+1/4*s*(5-s)*u[0][j]-1/4*(s-1)*(s+4)*u[1][j]+1/4*s*(s-1)*u[2][j];
		for (i=2;i<Nx;i++)
		{
u[i][j+1]=-1/4*s*(s-1)*u[i-2][j]+1/4*s*(5-s)*u[i-1][j]-1/4*(s-1)*(s+4)*u[i][j]+1/4*s*(s-1)*u[i+1][j];

	         }

u[Nx][j+1]=-1/4*s*(s-1)*u[Nx-2][j]+1/4*s*(5-s)*u[Nx-1][j]-1/4*(s-1)*(s+4)*u[Nx][j]+1/4*s*(s-1)*u[0][j];*/
u[0][j]=u[0][j-1]+s*((u[Nx][j-1]-u[0][j-1])+1/4*(1-s)*(u[0][j-1]-u[Nx-1][j-1]-u[0+1][j-1]+u[Nx][j-1])/4);
u[1][j]=u[1][j-1]+s*((u[0][j-1]-u[1][j-1])+(1-s)*(u[1][j-1]-u[Nx][j-1]-u[2][j-1]+u[0][j-1])/4);
		for (i=2;i<Nx;i++)
		{
u[i][j]=u[i][j-1]+s*((u[i-1][j-1]-u[i][j-1])+(1-s)*(u[i][j-1]-u[i-2][j-1]-u[i+1][j-1]+u[i-1][j-1])/4);

		}

u[Nx][j]=u[Nx][j-1]+s*((u[Nx-1][j-1]-u[Nx][j-1])+1/4*(1-s)*(u[Nx][j-1]-u[Nx-2][j-1]-u[1][j-1]+u[Nx-1][j-1])/4);

			/*u_pred[0]=u_corr[0]+1/2*(1-s)*(u_corr[1]-u_corr[Nx])/2;
			for (i=1;i<Nx;i++)
			{
			u_pred[i]=u_corr[i]+1/2*(1-s)*(u_corr[i+1]-u_corr[i-1])/2;
			}
                   	u_pred[Nx]=u_corr[Nx]+1/2*(1-s)*(u_corr[0]-u_corr[Nx-1])/2;
			u_corr[0]=u_corr[0]+s*(u_pred[1]-u_pred[Nx]);
			for (i=1;i<Nx;i++)
			{
			u_corr[i]=u_corr[i]+s*(u_pred[i+1]-u_pred[i-1]);
			}
			u_corr[Nx]=u_corr[Nx]+s*(u_pred[0]-u_pred[Nx-1]);
			
			for (i=0;i<=Nx;i++)
			{
			u[i][j+1]=u_corr[i];
			}*/
			
   		} 
        cout<<"COURANT NUMBER="<<alpha*dt/dx<<endl;
	output.close();	
        string text2= "wave_F_";
        string text3=".txt";
        stringstream file2;
        file2<<text2<<Nx<<text3;
	string file=file2.str();
	write(file,Nx,tstep,u,dx,dt);
	return 0;
	}
	void write ( string file, int Nx, int tstep, double **u,double dx,double dy )
	{
	int i,j;
        ofstream output;
        output.open(file.c_str());
	for (i=0;i<=Nx;i++)
	{
	for (j=0;j<tstep;j++)
		
        {
                output<<u[i][j]<<" ";
        }
                output << "\n";
	}
	output.close();
	return;
}


