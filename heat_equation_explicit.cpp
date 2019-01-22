//heat equation explicit
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
#include <vector>
using namespace std;
void write ( string file, int Nx, int Nt, double print[],double dt,double dx );
int main()
{
	int Nt;
	int Nx;
	double a;
	double b;
	double*x;
	double*t;
	double alpha;
	double*u0;
	double x_0=0;
	double x_N=1;
	
	Nx=128;
	double dx=(x_N-x_0)/(Nx-1);
        x=new double[Nx];
	for (int i=0; i<Nx; i++)
	{
		x[i]=i*dx;
	}
	//time
	double t_0=0;
	double t_end=5;
	Nt=5000;
	double dt=(t_end-t_0)/(Nt-1);

    t=new double[Nt];
	for (int n=0; n<Nt; n++)
	{
		t[n]=n*dt;
			}	
	//field initialization
	u0=new double[Nx];
	for (int i=0; i<Nx; i++)
	{
		u0[i]=1;
	}
	//alpha
	alpha=1e-1;
	alpha=alpha/pow(dx,2)*dt; 
	a=1;
	b=1;
	//define solution vector
	double u[Nx][Nt];
	double *source;

	double print[Nx*Nt];
	source=new double[Nx];
			for (int i=0;i<Nx;i++)
				{
					source[i]=-cos(2*3.14*x[i]);
				}
			for (int i=0;i<Nx;i++)
				{
					u[i][0]=u0[i];
				}
	for (int n=0;n<Nt-1;n++)
		{	
			cout<<n<<endl;

			
			for (int i=0;i<Nx;i++)
				{	
					if (i==0){
					u[i][n+1]=u[i][n]+alpha*(u[i+1][n]-2*u[i][n]+a)+source[i]*dt;
					}
					else if(i==Nx-1){
					u[i][n+1]=u[i][n]+alpha*(b-2*u[i][n]+u[i-1][n])+source[i]*dt;
					}
					else{u[i][n+1]=u[i][n]+alpha*(u[i+1][n]-2*u[i][n]+u[i-1][n])+source[i]*dt;}
					
				}
				for (int i=0;i<Nx;i++)
				{
					print[i+Nx*n]=u[i][n];
					cout<<u[i][n]<<endl;
				}
		}
	string file= "explicit_64.txt";
	write (file,Nx,Nt,print,dt,dx);


		return 0;
}


void write ( string file, int Nx, int Nt, double print[],double dt,double dx )
{
  ofstream output;
  output.open(file.c_str());
  for ( int j = 0; j < Nt; j++ )
  {
	output<<j*dt<<" ";
	
 for ( int i = 0; i < Nx; i++ )
    { 
      if (j==0){output<<setw(9)<<i*dx<<" ";}
      else{output << setw(9) << print[i+j*Nx] << "  ";}
    }
    
    output << "\n";
  }
  output.close();
  return;
}
