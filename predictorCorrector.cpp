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
	double sigma;
	a=-10;
	b=10;
	t0=0;
	tn=0.01;
	nx=25;
	nt=25;
	gamma=1.4;
	double rhol,ul,pl;
	double rhor,ur,pr;
	
	rhol=1;
	ul=0;
	pl=100000;

	rhor=0.01;
	ur=0;
	pr=1000;

/*	cout<<"Initial conditions U vector left"<<endl;	

	cin>>rhol;cout<<endl;
	cin>>ul;cout<<endl;
	cin>>pl;cout<<endl;

	cout<<"Initial conditions U vector right"<<endl;

	cin>>rhor;cout<<endl;
	cin>>ur;cout<<endl;
	cin>>pr;cout<<endl;
*/
	dx=(b-a)/(nx-1);
	dt=(tn-t0)/(nt-1);
	sigma=0.5*dt/dx;
	//initialize variables
	double *rho;
	double *u;
	double *p;
	double *rhou;
	double *e;
	double *x;
	double *c;
	double epsilon;

	double *frho;
	double *frhou;
	double *fe;
	double *rhopred;
	double *rhoupred;
	double *epred;
	double *ppred;
	double *upred;
	double *M;

	rho=new double [nx];
	u=new double [nx];
	p=new double [nx];
	rhou= new double [nx];
	e=new double [nx];
	x=new double [nx];
	c=new double [nx];
	
	frho=new double [nx];	
	frhou=new double [nx];
	fe=new double [nx];
	rhopred=new double [nx];
	rhoupred=new double [nx];
	epred=new double [nx];
	ppred=new double [nx];
	upred=new double [nx];
	M=new double [nx];

	for (int i=0;i<nx;i++)
		{
			x[i]=a+i*dx;
		}
	cout<<"insert an artificial viscosity value"
	cin>>epsilon;cout<<endl;
	for (int i=0;i<nx;i++)
		{
			if (x[i]<0)
			{
				rho[i]=rhol;
				u[i]=ul;
				p[i]=pl;
				rhou[i]=rho[i]*u[i];
				e[i]=p[i]/(gamma-1)+1/2*(rho[i]*pow(u[i],2));
			}
			else
			{
				rho[i]=rhor;
				u[i]=ur;
				p[i]=pr;
				rhou[i]=rho[i]*u[i];
				e[i]=p[i]/(gamma-1)+1/2*(rho[i]*pow(u[i],2));

			}

		}
	ofstream output;
	string file="b.txt";
	output.open(file.c_str());
	for (int n=0;n<=nt;n++)
	{
		for (int i=0;i<nx;i++)
		{
			frho[i]=rho[i]*u[i];
			frhou[i]=rho[i]*pow(u[i],2)+p[i];
			fe[i]=(e[i]+p[i])*u[i];	
			cout<<frhou[i]<<endl;
			//
		}	
		for (int i=0;i<nx;i++)
		{
			rhopred[i]=0.5*(rho[i]+rho[i+1])-sigma*(-frho[i]+frho[i+1]);
			rhoupred[i]=0.5*(rhou[i]+rhou[i+1])-sigma*(-frhou[i]+frhou[i+1]);
			epred[i]=0.5*(e[i]+e[i+1])-sigma*(-fe[i]+fe[i+1]);
			upred[i]=rhoupred[i]/rhopred[i];
			ppred[i]=(gamma-1)*(epred[i]-0.5*(rhopred[i]*pow(upred[i],2)));	
/*			cout<<i<<" ";
			cout<<rhopred[i]<<" ";
			cout<<rhoupred[i]<<" ";
			cout<<upred[i]<<endl;
			*/
		}
		for (int i=0;i<nx;i++)
		{	
			frho[i]=rhopred[i]*upred[i];
			frhou[i]=rhopred[i]*pow(upred[i],2)+ppred[i];
			fe[i]=(epred[i]+ppred[i])*upred[i];
			//
		}	
		for (int i=1;i<nx-1;i++)
		{
			rho[i]=rho[i]-dt/dx*(frho[i]-frho[i-1])+epsilon*(rho[i+1]+rho[i-1]-2*rho[i]);
			rhou[i]=rhou[i]-dt/dx*(frhou[i]-frhou[i-1])+epsilon*(rhou[i+1]+rhou[i-1]-2*rhou[i]);
			u[i]=rhou[i]/rho[i];
			e[i]=e[i]-dt/dx*(fe[i]-fe[i-1])+epsilon*(e[i+1]+e[i-1]-2*e[i]);
			p[i]=(gamma-1)*(e[i]-0.5*rho[i]*u[i]*u[i]);
			c[i]=pow(p[i]*gamma/rho[i],0.5);
			cout<<u[i]<<endl;
			M[0]=0;
			M[nx]=0;
			M[i]=u[i]/c[i];
		}

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
