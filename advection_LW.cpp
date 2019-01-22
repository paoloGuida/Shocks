
/* 
 * File:   advection_LW.cpp
 * Author: Paolo Guida
 *
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
        double *u_sub;
        double *u0;
	double number;
        x=new double[Nx+1];
	u=new double*[Nx+1];
        u_sub=new double [Nx+1];
        u0=new double [Nx+1];
	double s=alpha*dt/dx;
	cout<<"type 1 for the step function as u0 or 2 for the gaussian ";
	cin>>number;
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
		
               //cout<<i<<" "<<x[i]<<" "<<u[i][0]<<endl;
        }
	}
	for (i=0;i<=Nx;i++)
        {
		u[i][0]=u0[i];
	}
 	ofstream output;

		for (j=0;j<tstep;j++)
		{	
			u[0][j+1]=u[0][j]-s/2*(u[1][j]-u[Nx][j])+pow(s,2)/2*(u[1][j]-2*u[0][j]+u[Nx][j]);
                   	
			for (i=1;i<Nx;i++)
			{
				u[i][j+1]=u[i][j]-s/2*(u[i+1][j]-u[i-1][j])+pow(s,2)/2*(u[i+1][j]-2*u[i][j]+u[i-1][j]);
			}
				u[Nx][j+1]=u[Nx][j]-s/2*(u[0][j]-u[Nx-1][j])+pow(s,2)/2*(u[0][j-1]-2*u[Nx][j]+u[Nx-1][j]);

		}
        cout<<"COURANT NUMBER="<<alpha*dt/dx<<endl;
	output.close();	
        string text2= "wave_LW_";
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

/*  for (k=0;k<=Nx;k++)
			{
				if (k==0)
			{
			        u_sub[k]=(u0[k]+u0[k])/2+alpha*dt/dx/2*(u0[k]-u0[k]);           
			}			
				else if (k==Nx
			{
			        u_sub[k]=(u0[k]+u0[k-1])/2+alpha*dt/dx/2*(u0[k-1]-u0[k]);           
			}		
			}
			for (i=1;i<=Nx;i++)
			{
                            	u[i][j+1]=u0[i]-alpha*dt/dx*(u_sub[i]-u_sub[i-1]);           
			}
                        for (i=0;i<=Nx;i++)
			{
                            	u0[i]=u[i][j+1];           
			}*/
