

/* 
 * File:   advection_UP.cpp
 * Author: Paolo Guida
 */

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
	cout<<"set the number of grid points ";
	cin>>Nx;
	double Lx=2;
        double time=1;
        int tstep=80;
	double dx=Lx/(Nx);
	double dt=time/(tstep);
        double alpha=1;
	double *error;
	double **u;
        double *x;
	double *u0;
	int number;
       	double s=alpha*dt/dx;
	cout<<"type 1 for the step function as u0 or 2 for the gaussian ";
	cin>>number;
        x=new double[Nx+1];
	u=new double*[Nx+1];
	u0=new double[Nx+1];
	error=new double[Nx+1];
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
/*	for (j=0;j<=tstep;j++)
	{
	for (i=0;i<=Nx;i++)
        	{
			
			exact[i][j]=exp(-200*pow(x[i]-a*t[i]-.5,2));
		}
	}
	}*/
	for (i=0;i<=Nx;i++)
        {
		u[i][0]=u0[i];
	}

 	ofstream output;

	for (j=0;j<tstep;j++)
		{	
		//u[0][j+1]=u[Nx][j];
		for (i=1;i<=Nx;i++)
		{
			
                        u[i][j+1]=u[i][j]-alpha*dt/dx*(u[i][j]-u[i-1][j]); 
			         
		}
	
		} 
	for (i=0;i<=Nx;i++)
		{
			error[i]=fabs(u[i][tstep]-u0[i]);
			cout<<error[i]<<endl;
		}
	
        cout<<"COURANT NUMBER="<<alpha*dt/dx<<endl;
	output.close();	
        string text2= "wave_UP_";
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
/*
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
	int Nx=64;
	double Lx=1;
        double time=1;
        int tstep=1000;
	double dx=Lx/(Nx-1);
	double dt=time/(tstep-1);
        double alpha=1;
	double **u;
        double *x;
        x=new double[Nx+1];
	u=new double*[Nx+1];
        
	for (i=0;i<=Nx;i++)
	{
		u[i]=new double[tstep+1];
	}
        for (i=0;i<=Nx;i++)
	{
		x[i]=i*dx;
                u[i][0]=0;
	}

        for (i=0;i<=Nx;i++)
        {
		u[i][0]=exp(-200*pow(x[i]-.5,2));
               // cout<<i<<" "<<x[i]<<" "<<u[i][0]<<endl;
        }

 	ofstream output;

		for (j=0;j<tstep;j++)
		{	
			for (i=1;i<=Nx;i++)
			{
                            u[i][j+1]=u[i][j]-alpha*dt/dx*(u[i][j]-u[i-1][j]);           
			}
		}
        cout<<"COURANT NUMBER="<<alpha*dt/dx<<endl;
	output.close();	
        string text2= "wave_upwind_";
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
	for (i=1;i<=Nx;i++)
	{
	for (j=1;j<=tstep;j++)
		
        {
                output<<u[i][j]<<" ";
        }
                output << "\n";
	}
  output.close();
  return;
}*/
