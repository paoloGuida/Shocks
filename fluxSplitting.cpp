# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
using namespace std;

int main()
{
	//data
	int nx,nt;
	double a,b;
	double t0,tn;
	double dx,dt;
	double gamma;
	double epsilon;
	double sigma;


	a=-10;
	b=10;
	t0=0;
	tn=0.01;
	nx=50;
	nt=25;
	gamma=1.4;
	dx=(b-a)/(nx-1);
	dt=(tn-t0)/(nt-1);
	sigma=0.5*dt/dx;


	double rhol,ul,pl;
	double rhor,ur,pr;

	rhol=1;
	ul=0;
	pl=100000;

	rhor=0.01;
	ur=0;
	pr=1000;
	

	/*cout<<"Initial conditions U vector left"<<endl;	

	cin>>rhol;cout<<endl;
	cin>>ul;cout<<endl;
	cin>>pl;cout<<endl;

	cout<<"Initial conditions U vector right"<<endl;

	cin>>rhor;cout<<endl;
	cin>>ur;cout<<endl;
	cin>>pr;cout<<endl;*/	


	//initialize variables
	double *rho;
	double *u;
	double *p;
	double *rhou;
	double *rhoe;
	double *e;
	double *x;
	double *M;

	double *fplus1;
	double *fplus2;
	double *fplus3;
	double *fminus1;
	double *fminus2;
	double *fminus3;
	double *c;


	rho=new double [nx];
	u=new double [nx];
	p=new double [nx];
	rhou= new double [nx];
	e=new double [nx];
	x=new double [nx];
	rhoe=new double [nx];
	M=new double [nx];
	c=new double [nx];
	fplus1=new double [nx];
	fplus2=new double [nx];
	fplus3=new double [nx];

	fminus1=new double [nx];
	fminus2=new double [nx];
	fminus3=new double [nx];

	for (int i=0;i<nx;i++)
		{
			x[i]=a+i*dx;
			//cout<<x[i]<<endl;
		}

	for (int i=0;i<nx;i++)
		{
			if (x[i]<0)
			{
				rho[i]=rhol;
				u[i]=ul;
				p[i]=pl;
				rhou[i]=rho[i]*u[i];
				rhoe[i]=(p[i]/0.4+0.5*(rho[i]*pow(u[i],2)));
				c[i]=pow(p[i]*gamma/rho[i],0.5);

				M[i]=u[i]/c[i];
			}
			else
			{
				rho[i]=rhor;
				u[i]=ur;
				p[i]=pr;
				rhou[i]=rho[i]*u[i];
				rhoe[i]=p[i]/0.4+0.5*(rho[i]*pow(u[i],2));

				c[i]=pow(p[i]*gamma/rho[i],0.5);
				M[i]=u[i]/c[i];
			
			}
				cout<<c[i]<<endl;

		}

	ofstream output;
	string file="b.txt";
	output.open(file.c_str());

		for (int n=0;n<=nt;n++)
		{

		for (int i=0;i<nx;i++)
		{

			fplus1[i]=rho[i]*c[i]/4*pow((M[i]+1),2)*(1);
			fplus2[i]=rho[i]*c[i]/4*pow(M[i]+1,2)*2*c[i]/gamma*(1+(gamma-1)*M[i]/2);
			fplus3[i]=rho[i]*c[i]/4*pow(M[i]+1,2)*2*c[i]*c[i]/(gamma*gamma-1)*pow(1+(gamma-1)*M[i]/2,2); 
			//
			fminus1[i]=-rho[i]*c[i]/4*pow((M[i]-1),2)*(1);
			fminus2[i]=-rho[i]*c[i]/4*pow((M[i]-1),2)*(2*c[i]/gamma*(-1+(gamma-1)*M[i]/2));
			fminus3[i]=-rho[i]*c[i]/4*pow((M[i]-1),2)*2*c[i]*c[i]/(gamma*gamma-1)*pow((1-(gamma-1)*M[i]/2),2);
			//

		}	

		for (int i=1;i<nx-1;i++)
		{
			rho[i]=rho[i]-dt/dx*(fplus1[i]-fplus1[i-1]+fminus1[i+1]-fminus1[i]);
			rhou[i]=rhou[i]-dt/dx*(fplus2[i]-fplus2[i-1]+fminus2[i+1]-fminus2[i]);
			rhoe[i]=rhoe[i]-dt/dx*(fplus3[i]-fplus3[i-1]+fminus3[i+1]-fminus3[i]);
			//
		}
		for (int i=1;i<nx-1;i++)
		{
		    u[i]=rhou[i]/rho[i];
	            p[i]=(gamma-1)*(rhoe[i]-rho[i]*u[i]*u[i]/2);
		    c[i]=pow(p[i]*gamma/rho[i],0.5);
		    M[i]=u[i]/c[i];
		    
		     
		}
			
			//
		
		}
		for (int i=0;i<nx;i++)
		{
		    output<<x[i]<<" ";
		    output<<u[i]<<" ";
		    output<<rho[i]<<" ";
		    output<<p[i]<<" ";
		    output<<M[i]<<" ";
		    output<<endl;
		     
		}
  	output.close();

	return 0;
}


