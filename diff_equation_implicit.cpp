# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <math.h>
# include <vector>
using namespace std;

vector<double>solver(const vector<double>& a,
                      const vector<double>& b,
                      const vector<double>& c,
                      vector<double>& d,int Nx);
void write ( string file, int m, int n, double print[],double dt,double dx );
int main()
{	
	int Nt;
	int Nx;
	double u_a;
	double u_b;
	double*x;
	double*t;
	double alpha;
	double*u;
	
	//space
	double x_0=0;
	double x_N=1;
	Nx=256;
	double dx=(x_N-x_0)/(Nx-1);
    	x=new double[Nx];
	for (int i=0; i<Nx; i++)
	{
		x[i]=i*dx;
	}
	//time
	u_a=1;
	u_b=1;
	double t_0=0;
	double t_end=10;
	Nt=1000;
	double dt=(t_end-t_0)/(Nt-1);
    	t=new double[Nt];
	for (int n=0; n<Nt; n++)
	{
		t[n]=n*dt;
	}	

	//alpha
	alpha=1e-1;
	alpha=alpha/pow(dx,2)*dt;
        cout<<alpha<<endl; 
   	vector<double> a(Nx-1, -alpha);
    	vector<double> b(Nx, 1+2*alpha);
    	vector<double> c(Nx-1, -alpha);

    	vector<double> d(Nx);
	vector<double> solution(Nx);
    	double f[Nx][Nt];
    	vector<double> source(Nx,0);
    	u=new double [Nx*Nt];
			for (int i=0;i<Nx;i++)
			{
				source[i]=-cos(2*3.14*x[i]);
				//cout<<source[i]<<endl;
				f[i][0]=1;
			}

	for (int n=0;n<Nt;n++)
		{
			cout<<n<<endl;
			d[0]=f[0][n]+alpha*u_a+source[0]*dt;
			d[Nx-1]=f[Nx-1][n]+alpha*u_b+source[Nx-1]*dt;
			for (int i=1;i<Nx-1;i++)
				{
					d[i]=f[i][n]+source[i]*dt;
					
					//cout<<f[i]<<endl;
				}
	for (int i=0;i<Nx;i++){cout<<d[i]<<endl;}
			cout<<n<<endl;
			solution=solver(a,b,c,d,Nx);
			
			
			for (int i=0;i<Nx;i++)
			{
				

				f[i][n+1]=solution[i];
				cout<<f[i][n+1]<<endl;
				u[i+n*Nx]=f[i][n];
			}
					

		}
			string file = "implicit_256.txt";
			write ( file, Nx, Nt, u,dt,dx);

return 0;
}
vector<double>solver(const vector<double>& a,
                      const vector<double>& b,
                      const vector<double>& c,
                      vector<double>& d,int Nx) 
 { 
	vector<double> solution(Nx);
	vector<double> b_(Nx);
	b_[0]=b[0];
    for(int i=1; i<Nx; i++){
 
        double m = a[i-1]/b_[i-1];
        b_[i] =b[i]-m*c[i-1];
        d[i] -= m*d[i-1];
    }
    solution[Nx-1] = d[Nx-1]/b_[Nx-1];
 	
    for(int i=Nx-2; i >= 0; i--){
        solution[i] = (d[i]- c[i]*solution[i+1])/b_[i];}
		for (int i=0;i<Nx;i++)
	{
		//cout<<solution[i]<<endl;
}
    return solution;
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
