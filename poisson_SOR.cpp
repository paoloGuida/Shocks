
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;
void write ( string file, int Nx, int Ny, double **u,double dx,double dy );
int main()
{
	int i;
	int j;
	int k;
	int Nx=256;
	int Ny=256;
	double Lx=1;
	double Ly=1;
	double dx=Lx/(Nx-1);
	double dy=Ly/(Ny-1);
	double tot=0;
	int maxIter=10000;
	double tol=1e-5;
	double absRes=0.5;
	//BCs
	double n=0;
	double s=1;
	double e=0;
	double w=0;
	double **u;
	double **uNew;
	double beta=1.2;
	u=new double*[Ny+1];
	uNew=new double*[Ny+1];
	for (j=0;j<=Ny;j++)
	{
		uNew[j]=new double[Nx+1];
		u[j]=new double[Nx+1];
	}
	for (i=0;i<=Nx;i++)
	{
		for (j=0;j<=Ny;j++)
		{
			u[i][j]=0;	
		}
	}
  /*for (i=0;i<=Nx;i++)
	{
	for (j=0;j<=Ny;j++)
		{
			cout<<u[i][j]<<" ";
		}
		cout<<endl;
	}*/
	for (j=0;j<=Ny;j++)
	{
		u[0][j]=s;
		u[Nx][j]=n;
	}
	for (i=1;i<Nx;i++)
	{
		u[i][0]=w;
		u[i][Ny]=e;
	}
	for (i=0;i<=Nx;i++)
	{
		for (j=0;j<=Ny;j++)
		{
			uNew[i][j]=u[i][j];		
		}
	}
	int iter=0;
	string residue="residue_SOR_256.txt";
 	ofstream output;
  	output.open(residue.c_str());
	for (k=0;k<maxIter;k++)
	{
		
		tot=0;
		iter=iter+1;
		output<<iter<<" ";
		cout<<iter<<" ";
		for (i=1;i<Nx;i++)
		{	
			for (j=1;j<Ny;j++)
			{
				uNew[i][j]=(u[i+1][j]+uNew[i-1][j]+u[i][j+1]+uNew[i][j-1])*beta/4+u[i][j]*(1-beta);

			}
		
		}

	for (i=0;i<=Nx;i++)
	{
	for (j=0;j<=Ny;j++)
		{
			
			tot+=fabs(uNew[i][j]-u[i][j]);
			u[i][j]=uNew[i][j];	
		}
	}
		output<<absRes<<"\n";
		cout<<absRes<<"\n";
		absRes=tot/(Nx-1)/(Ny-1);
		if (absRes<tol){break;}

	}
		output.close();		

	string file= "poisson_SOR_256.txt";
	write(file,Nx,Ny,u,dx,dy);
	return 0;
	
}
void write ( string file, int Nx, int Ny, double **u,double dx,double dy )
{
	int i,j;
  ofstream output;
  output.open(file.c_str());
	for (i=0;i<=Nx;i++)
	{
	for (j=0;j<=Ny;j++)
		{
			output<<u[i][j]<<" ";
		}
    output << "\n";
	}
  output.close();
  return;
}
