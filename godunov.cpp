# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <algorithm>
# include <bits/stdc++.h>

#define TOL 0.00001


using namespace std;
double starPressure(double ul,double pl,double rhol,double ur,double pr,double rhor);
struct riemann
{
	double pMid;
	double rhoMid;
	double uMid;
};

double preFun(double P,double rho,double p,double c)
{
	double flux;
	if(P<=p)      //rarefaction wave
	{
		flux=(pow(P/p,0.1429)-1)*5*c;

	}
	else	      //shock
	{
		flux=(P-p)*pow(0.8333/(rho*(0.1667*rho+P)),0.5);
		
	}	

	return flux;
}
double preFunDer(double P,double rho,double p,double c)
{
	double fluxDerivative;
	if(P<=p)      //rarefaction wave
	{
		
		fluxDerivative=(1/(rho*c))*pow(P/p,-0.8571);
	}
	else	      //shock
	{
		
		fluxDerivative=(1-0.5*(P-p)/(0.1667/rho+P))*pow(0.8333/rho/(0.1667*rho+p),0.5);
	}	
			
	return fluxDerivative;
}
riemann riemannSolver(double rhol,double rhor,double ul,double ur,double cl,double cr,double pl,double pr,double pStar,double time,double domain,int cells);

int main()
{
	//data
	int nx,nt,nv;
	double a,b;
	double t0,tn;
	double dx,dt;
	double gamma;
	double sigma;
	int cells=20;

	
	a=-10;
	b=10;

	t0=0;
	tn=0.01;
	nx=21;
	nt=100;
	nv=nx-1;
	gamma=1.4;
	double rhol,ul,pl,cl;
	double rhor,ur,pr,cr;
	
	rhol=1;
	ul=0;
	pl=100000;
	cl=sqrt(gamma*pl/rhol);

	rhor=0.125;
	ur=0;
	pr=10000;
	cr=sqrt(gamma*pr/rhor);
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
	



	rho=new double [nv];
	u=new double [nv];
	p=new double [nv];
	rhou= new double [nv];
	e=new double [nv];
	x=new double [nv];
	c=new double [nv];
	


	for (int i=0;i<nx;i++)
		{
			x[i]=a+i*dx;
		}

	for (int i=0;i<nx;i++)
		{
			if (x[i]<0)
			{
				rho[i]=rhol;
				u[i]=ul;
				p[i]=pl;
				rhou[i]=rho[i]*u[i];
				c[i]=pow(gamma*p[i]/rho[i],0.5);
				e[i]=p[i]/(gamma-1)+1/2*(rho[i]*pow(u[i],2));
			}
			else
			{
				rho[i]=rhor;
				u[i]=ur;
				p[i]=pr;
				c[i]=pow(gamma*p[i]/rho[i],0.5);
				rhou[i]=rho[i]*u[i];
				e[i]=p[i]/(gamma-1)+1/2*(rho[i]*pow(u[i],2));

			}
			

		}
	double *uMid,*rhoMid,*pMid,*flux1,*flux2,*flux3;
	double pStar;
	uMid=new double [nv];
	pMid=new double [nv];
	rhoMid=new double [nv];
	
	flux1=new double[nv];
	flux2=new double[nv];
	flux3=new double[nv];

	ofstream output;
	string file="b.txt";
	output.open(file.c_str());
	for (int n=0;n<nt;n++)
	{
		for (int i=0;i<nv;i++)
		{
			
			rho[i]=rho[i]-dt*(flux1[i]-flux1[i-1])/dx;
			rhou[i]=rhou[i]-dt*(flux2[i]-flux2[i-1])/dx;
			e[i]=e[i]-dt*(flux3[i]-flux3[i-1])/dx;
			u[i]=rhou[i]/rho[i];
			c[i]=pow(gamma*p[i]/rho[i],0.5);
			p[i]=(gamma-1)*(e[i]-0.5*(rho[i]*pow(u[i],2)));
			cout<<i<<" "<<rhou[i]<<" "<<u[i]<<" "<<rho[i]<<" "<<p[i]<<" "<<endl;
		}	
		
		cout<<"Riemann solution"<<endl;
		for (int i=0;i<nv;i++)
		{
			
			pStar=starPressure(u[i],p[i],rho[i],u[i+1],p[i+1],rho[i+1]);
			riemann variables=riemannSolver(rho[i],rho[i+1],u[i],u[i+1],c[i],c[i+1],p[i],p[i+1],pStar,dt,dx,cells);
			
			pMid[i]=variables.pMid;
			rhoMid[i]=variables.rhoMid;
			uMid[i]=variables.uMid;
			cout<<i<<" "<<uMid[i]<<" "<<rhoMid[i]<<" "<<pMid[i]<<" "<<endl;
		}
		cout<<"Fluxes"<<endl;	
		for (int i=0;i<nv;i++)
		{
			
			flux1[i]=rhoMid[i]*uMid[i];
			flux2[i]=rhoMid[i]*uMid[i]*uMid[i];
			flux3[i]=uMid[i]*(pMid[i]/(gamma-1)+1/2*(rhoMid[i]*pow(uMid[i],2))+pMid[i]);	
			cout<<i<<" "<<flux1[i]<<" "<<flux2[i]<<" "<<flux3[i]<<" "<<endl;
		}
		cout<<"Values update"<<endl;
		for (int i=1;i<nv;i++)
		{
			
			rho[i]=rho[i]-dt*(flux1[i]-flux1[i-1])/dx;
			rhou[i]=rhou[i]-dt*(flux2[i]-flux2[i-1])/dx;
			e[i]=e[i]-dt*(flux3[i]-flux3[i-1])/dx;
			u[i]=rhou[i]/rho[i];
			c[i]=pow(gamma*p[i]/rho[i],0.5);
			p[i]=(gamma-1)*(e[i]-0.5*(rho[i]*pow(u[i],2)));
			cout<<i<<" "<<rhou[i]<<" "<<u[i]<<" "<<rho[i]<<" "<<p[i]<<" "<<endl;
		}	
		
		 	


		
	}
		for (int i=0;i<nx;i++)
		{
		 	
		    	output<<x[i]<<" ";
		    	output<<u[i]<<" ";
		    	output<<rho[i]<<" ";
		    	output<<p[i]<<" ";
		    	
		output<<endl;
		}
  	output.close();
	return 0;
}
double starPressure(double ul,double pl,double rhol,double ur,double pr,double rhor)
{
	double P=3.1527e+04;
	double pStar;
	double check=0;
	double cl=pow(1.4*pl/rhol,0.5);
	double cr=pow(1.4*pr/rhor,0.5);
	double maxIter=400;
	double fluxL;
	double fluxDerivativeL;
	double fluxR;
	double fluxDerivativeR;
	double uStar;
	for(int n=0;n<=maxIter;n++)
		{
			fluxL=preFun(P,rhol,pl,cl);
			fluxDerivativeL=preFunDer(P,rhol,pl,cl);
			fluxR=preFun(P,rhor,pr,cr);
			fluxDerivativeR=preFunDer(P,rhor,pr,cr);

			pStar=P-(fluxL+fluxR+ur-ul)/(fluxDerivativeL+fluxDerivativeR);
			check=2*abs((pStar-P)/(pStar+P));
			
			if (check<TOL)
			{
				break;
			}	
			if(pStar<0)
			{
				pStar=TOL;
			}
			P=pStar;
		}

	return pStar;
}

			
riemann riemannSolver(double rhol,double rhor,double ul,double ur,double cl,double cr,double pl,double pr,double pStar,double time,double domain,int cells)
	{
		double dx=domain/(cells);
		//cout<<"dx "<<dx<<endl;
		double *xpos;
		double *u,*p,*rho;
		double j;
		xpos=new double [cells];
		u=new double [cells];
		p=new double [cells];
		rho=new double [cells];
		int mid=(cells-1)/2+1;
		double fluxR=preFun(pStar,rhor,pr,cr);
		double fluxL=preFun(pStar,rhol,pl,cl);
		double uStar=0.5*(ul+ur+fluxR-fluxL);
		double s,c,sl,sr;
		double shl,shr;
		double rhos,us,ps;
		double cml,cmr,stl,str,pml,pmr;
		j=-1;
		for (int i=0;i<cells;i++)
		{
			j=j+1;
			xpos[i]=(j-1/2)*dx;
			s=(xpos[i]-domain/2)/time;
			//cout<<i<<" "<<xpos[i]<<" "<<s<<endl;
			if (s<=uStar)
			{
				if (pStar<=pl) //left fan
				{
					shl=ul-cl;
					if (s<=shl)
					{
						rhos=rhol;
						us=ul;
						ps=pl;
					}
					else
					{
						cml=cl*pow(pStar/pl,0.1429);
						stl=uStar-cml;
						if(s>stl)
						{
							rhos=rhol*pow(pStar/pl,1/1.4);
							us=uStar;
							ps=pStar;
						}
						else
						{
							us=0.8333*(cl+0.2*ul+s);
							c=0.8333*(cl+0.2*(ul-s));
							rhos=rhol*pow(c/cl,5);
							ps=pl*pow(c/cl,7);
						}
						
					}
				}
				else
				{
					pml=pStar/pl;
					sl=ul-cl*pow(0.8571*pml+0.1429,0.5);
					
					if(s<=sl)
					{
						rhos=rhol;
						us=ul;
						ps=pl;
					}
					else
					{
						rhos=rhol*(pml+0.1667)/(pml*0.1667+1);
						us=uStar;
						ps=pStar;
					}
				}
			}
			else
			{
				if (pStar>pr)
				{
					pmr=pStar/pr;
					sr=ur+cr*pow(0.8571*pmr+0.1429,0.5);
					if(s>sr)
					{
						rhos=rhor;
						us=ur;
						ps=pr;
					}
					else
					{
						rhos=rhor*(pmr+0.1667)/(pmr*0.1667+1);
						us=uStar;
						ps=pStar;
						
					}
				}
				else
				{
					shr=ur+cr;
					if (s>=shr)
					{
						rhos=rhor;
						us=ur;
						ps=pr;
					}
					else
					{
						cmr=cr*pow(pStar/pr,0.1429);
						str=uStar+cmr;
						if (s<=str)
						{
						rhos=rhor*pow(pStar/pr,1/1.4);
						us=uStar;
						ps=pStar;
						}
						else
						{
							us=0.8333*(-cr+0.2*ur+s);
							c=0.8333*(cr-0.2*(ur-s));
							rhos=rhor*pow(c/cr,5);	
							ps=pr*pow(c/cr,7);
						}
						
					}
					
				}
				
				
			}
			
			u[i]=us;
			p[i]=ps;
			rho[i]=rhos;
			
			
		}
		/*double pAverage,rhoAverage,uAverage;
		rhoAverage=0;
		pAverage=0;
		uAverage=0;
		for (int i=0;i<cells;i++)
				{
					rhoAverage+=rho[i];
					pAverage+=p[i];
					uAverage+=u[i];
				}
				
		return riemann{pAverage/cells,rhoAverage/cells,uAverage/cells};*/
		
		return riemann{p[mid-1],rho[mid-1],u[mid-1]};
			
			
	}
		
