# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <algorithm>
# include <bits/stdc++.h>

#define TOL 1e-10
#define gamma 1.4

using namespace std;
struct expansion
{
	double p;
	double dp;
	double c;
};
struct shock
{
	double p;
	double dp;
	double c; 
};
struct riemann
{
	double pInt;
	double cInt;
	double uInt;
};
expansion expansionSolver(double u0,double c0,double p0,double uStar)
{
	double c,p,dp;
	c=c0-((1.4-1)/2)*abs(uStar-u0);
	p=p0*pow(c/c0,(2*1.4)/(1.4-1));
	dp=1.4*(p/c);
	return{p,dp,c};
}
shock shockSolver(double u0,double c0,double p0,double uStar)
{
	double m,p,dp,c;
	m=((gamma+1)/4)*(abs(uStar-u0)/c0)+sqrt(1+pow((gamma+1)/4*abs(uStar-u0)/c0,2));
	//potenziale errore
	
	p=p0*(1+(2*1.4/(1.4+1))*(pow(m,2)-1));
	dp=2*1.4*(p0/c0)*(pow(m,3)/(1+pow(m,2)));
	c=c0*sqrt((1.4+1+(1.4-1)*(p/p0))/(1.4+1+(1.4-1)*(p0/p)));
	return{p,dp,c};
}

riemann riemannSolver(double ul,double cl,double pl,double ur,double cr,double pr)
{
	
      	int iter=10;

	double delta,rl,rr,z,u,p1,dp1,c1,leftWave,rightWave,p2,dp2,c2,u1,u2;
	//cout<<"INPUT RIEMANN SOLVER "<<ul<<" "<<cl<<" "<<pl<<" "<<ur<<" "<<cr<<" "<<pr<<endl;
      	delta=(1.4-1)/2;
      	rl=cl+((1.4-1)/2)*ul;
      	rr=cr-((1.4-1)/2)*ur;
      	z=(cr/cl)*pow(pl/pr,delta/1.4);
      	u=(z*rl-rr)/(delta*(1+z));
      	//cout<<"questa u "<<u<<endl;
	for (int i=0;i<=iter;i++)
		{	
			//left
            		if (u<=ul || (u-ul)<1e-5)
            			{
            				shock shockl=shockSolver(ul,cl,pl,u);
					p1=shockl.p;
					dp1=shockl.dp;
					c1=shockl.c;
  					leftWave=1;
  					  					//cout<<"leftWave=1"<<endl;
  				}
			else
				{
					expansion expl=expansionSolver(ul,cl,pl,u);
					p1=expl.p;
					dp1=expl.dp;
					c1=expl.c;
                    			leftWave=2;
                    			  					//cout<<"leftWave=2"<<endl;
            			}
            			dp1=-dp1;
            		//right
			if (u>=ur && (u-ur)>1e-5)
            			{
            				shock shockr=shockSolver(ur,cr,pr,u);
					p2=shockr.p;
					dp2=shockr.dp;
					c2=shockr.c;
  					rightWave=1;
  					//cout<<"rightWave=1"<<endl;
  				}
			else
				{
					expansion expr=expansionSolver(ur,cr,pr,u);
					p2=expr.p;
					dp2=expr.dp;
					c2=expr.c;
                    			rightWave=2;
                    			  					//cout<<"rightWave=2"<<endl;
            			}  			
   			if(abs(1-(p1/p2))<TOL)
   			{ 
   			
   				break;
   			
   			}
   			u=u-((p1-p2)/(dp1-dp2));
      		}
      		u1=u;
      		u2=u;	
      		//cout<<"u value "<<u<<endl;
      		double u_i,p_i,c_i,Val1,Val2;
      		double A,B,C,s1,s2,s3,s4,Mach1L,Mach2L,Mach1R,Mach2R,jump1L,jump2L,jump1R,jump2R,S_L,S_R,Val3,Val4;
      	if (leftWave==1)
      	{

      		
      		A=2/(gamma+1);
       		B=-((3-gamma)/(gamma+1)*ul+u1);
       		C=ul*u1-(2/(gamma+1))*pow(cl,2)-((gamma-1)/(gamma+1))*pow(ul,2);
       		s1=(-B-sqrt(pow(B,2)-4*A*C))/(2*A);
       		s2=(-B+sqrt(pow(B,2)-4*A*C))/(2*A);
       		Mach1L=(ul-s1)/cl;
      		Mach2L=(ul-s2)/cl;
      		jump1L=1+((2*gamma/(gamma+1))*(pow(Mach1L,2)-1));
      		jump2L=1+((2*gamma/(gamma+1))*(pow(Mach2L,2)-1));
      		//cout<<"jump1L  "<<jump1L<<" jump2L "<<jump2L<<endl;
       		if (abs((p1/pl)-jump1L)<3e-5)
         	{S_L=s1;}
      		else if (abs((p1/pl)-jump2L)<3e-5)
          	{S_L=s2;}
      		else if (abs((p1/pl)-jump1L)<3e-1)
          	{S_L=s1;}
      		else if (abs((p1/pl)-jump2L)<3e-1)
         	{S_L=s2;}
    		
          	
          	if (S_L<0 && abs(S_L)>1e-5)
             	{	
             		if (u>0)
                	{
                	u_i=u1;
                	p_i=p1;
                	c_i=c1;
                	}
             		else
               		{
               		u_i=u2;
                	p_i=p2;
                	c_i=c2;
             		}
             	}
          	else if (S_L>=0 || S_L<1e-5)
                {
                u_i=ul;
                p_i=pl;
                c_i=cl;
                }
                
           }
          if (leftWave==2)
       		{
       		Val1=ul-cl;
       		Val2=u1-c1;
       		
       		if (Val2<=0 && u>=0)
          	{
          	u_i=u1;
          	p_i=p1;
          	c_i=c1;
          	}
       		else if (Val2<0 && u<0)
          	{	
          	u_i=u2;
          	p_i=p2;
          	c_i=c2;
       		}
       		if (Val1>=0)
          	{
          	u_i=ul;
          	p_i=pl;
          	c_i=cl;
          	}
      
     		if (Val1<0 && Val2>0)
          	{
          	u_i=(cl+((gamma-1)/2)*ul)/(1+((gamma-1)/2)); 
          	c_i=u_i;
          	p_i=pl*(pow((c_i/cl),(2*gamma)/(gamma-1)));
          	}
       		
       		S_L=0;
    		}
    		
    		if (rightWave==1)
       		{
       		A=2/(gamma+1);
       		B=-(((3-gamma)/(gamma+1))*ur+u2);
       		C=ur*u2-(2/(gamma+1))*pow(cr,2)-((gamma-1)/(gamma+1))*pow(ur,2);
       		s3=(-B-sqrt(pow(B,2)-4*A*C))/(2*A);
       		s4=(-B+sqrt(pow(B,2)-4*A*C))/(2*A);
       		Mach1R=(ur-s3)/cr;
     		Mach2R=(ur-s4)/cr;
      		jump1R=1+((2*gamma/(gamma+1))*(pow(Mach1R,2)-1));
      		jump2R=1+((2*gamma/(gamma+1))*(pow(Mach2R,2)-1));
            	//cout<<" jump1R  "<<jump1R<<" jump2R "<<jump2R<<endl;
      		if (abs((p2/pr)-jump1R)<1e-5)
          	{
          	S_R=s3;
          	}
      		else if (abs((p2/pr)-jump2R)<1e-5)
          	{
          	S_R=s4;
          	}
      		else if (abs((p2/pr)-jump1R)<1e-3)
          	{
          	S_R=s3;
          	}
      		else if (abs((p2/pr)-jump2R)<1e-3)
          	{
          	S_R=s4;
          	}
      
          	
          	if (S_R<=0)
             	{
             	u_i=ur;
             	c_i=cr;
             	p_i=pr;
             	}
   
    		}
    		if (rightWave==2)
       		{
       		Val4=ur+cr;
       		Val3=u2+c2;
         	if (Val4>0 && Val3<0)
            	{
            	u_i=-(cr-((gamma-1)/2)*ur)/(1+((gamma-1)/2));
            	c_i=-u_i;
            	p_i=pr*(pow(c_i/cr,(2*gamma)/(gamma-1)));
            	}
      
         	if (Val4<=0)
           	{
            	u_i=ur;
            	c_i=cr;
            	p_i=pr;
            	}
        	S_R=0;
    		}
    	return {p_i,c_i,u_i};
}
      	
int main()
{
	//data
	int nx,nt,nv;
	double a,b;
	double t0,tn;
	double dx,dt;
	double gamma;
	double sigammaa;
	int cells=20;

	
	a=-10;
	b=10;

	t0=0;
	tn=0.01;
	nx=400;
	nt=10000;

	nv=nx-1;
	gamma=1.4;
	double rhol,ul,pl,cl;
	double rhor,ur,pr,cr;
	
	rhol=1;
	ul=0;
	pl=100000;
	cl=sqrt(gamma*pl/rhol);

	rhor=0.01;
	ur=0;
	pr=1000;
	cr=sqrt(gamma*pr/rhor);
/*	cout<<"Initial conditions U vector left"<<endl;	

	cin>>rhol;cout<<endl;
	cin>>ul;cout<<endl;
	cin>>pl;cout<<endl;

	cout<<"Initial conditions U vector right"<<endl;

	cin>>rhor;cout<<endl;
	cin>>ur;cout<<endl;
	cin>>pr;cout<<endl;
*/int j;
	dx=(b-a)/(nx-1);
	dt=(tn-t0)/(nt-1);
	sigammaa=dt/dx;
	//cout<<sigammaa<<endl;
	//initialize variables
	double *rho;
	double *u;
	double *p;
	double *rhou;
	double *e;
	double *x;
	double *c;
	double epsilon;
	double *uMid,*rhoMid,*pMid,*cMid,*eMid,*flux1,*flux2,*flux3;

	uMid=new double [nv];
	pMid=new double [nv];
	rhoMid=new double [nv];
	cMid=new double[nv];
	eMid=new double[nv];
	flux1=new double[nv];
	flux2=new double[nv];
	flux3=new double[nv];


	rho=new double [nx];
	u=new double [nx];
	p=new double [nx];
	rhou= new double [nx];
	e=new double [nx];
	x=new double [nx];
	c=new double [nx];
	


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

j=0;
	ofstream output;
	string file="results.txt";
	output.open(file.c_str());
	for (int n=0;n<nt;n++)
	{

		for (int i=0;i<nv;i++)
		{


			
			riemann variables=riemannSolver(u[i],c[i],p[i],u[i+1],c[i+1],p[i+1]);
			pMid[i]=variables.pInt;
			cMid[i]=variables.cInt;
			uMid[i]=variables.uInt;
			rhoMid[i]=pMid[i]*gamma/(c[i]*c[i]);
			//cout<<i<<" "<<pMid[i]<<" "<<cMid[i]<<" "<<uMid[i]<<" "<<endl;
			eMid[i]=pMid[i]/(gamma-1)+0.5*rhoMid[i]*uMid[i]*uMid[i];
		}
		//cout<<"Fluxes"<<endl;	
		for (int i=0;i<nv;i++)
		{
			
			flux1[i]=rhoMid[i]*uMid[i];
			flux2[i]=rhoMid[i]*uMid[i]*uMid[i]+pMid[i];
			flux3[i]=uMid[i]*(eMid[i]+pMid[i]);	
			//cout<<i<<" "<<flux1[i]<<" "<<flux2[i]<<" "<<flux3[i]<<" "<<endl;
		}
		//cout<<"Values update"<<endl;
		for (int i=1;i<nx;i++)
		{
			
			rho[i]=rho[i]-dt*(flux1[i]-flux1[i-1])/dx;
			rhou[i]=rhou[i]-dt*(flux2[i]-flux2[i-1])/dx;
			e[i]=e[i]-dt*(flux3[i]-flux3[i-1])/dx;
			u[i]=rhou[i]/rho[i];
			c[i]=pow(gamma*p[i]/rho[i],0.5);
			p[i]=(gamma-1)*(e[i]-0.5*(rho[i]*pow(u[i],2)));
			//cout<<i<<" "<<rhou[i]<<" "<<u[i]<<" "<<rho[i]<<" "<<p[i]<<" "<<endl;
		}	
		u[0]=ul;
		u[nx-1]=ur;
		p[0]=pl;
		p[nv]=pr;
		rho[0]=rhol;
		rho[nv]=rhor;
		c[0]=pow(gamma*p[0]/rho[0],0.5);
		e[0]=p[0]/(gamma-1)+1/2*(rho[0]*pow(u[0],2));
		c[nv]=pow(gamma*p[nv]/rho[nv],0.5);
		e[nv]=p[nv]/(gamma-1)+1/2*(rho[nv]*pow(u[nv],2));

		
	}
	j=0;
		for (int i=0;i<nx;i++)
		{
		 	
		    	output<<x[i]<<" ";
		    	output<<u[i]<<" ";
		    	output<<rho[i]<<" ";
		    	output<<p[i]<<" ";
		    	output<<u[i]/c[i]<<" ";
		    	
		output<<endl;
		}
  	output.close();
	return 0;
}


