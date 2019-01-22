# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <algorithm>
# include <bits/stdc++.h>

#define TOL 1e-4
using namespace std;


      	
int main()
{
	//data
	int nx,nt,nv;
	double a,b;
	double t0,tn;
	double dx,dt;
	double gamma;
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

	rhor=0.125;
	ur=0;
	pr=10000;
	cr=sqrt(gamma*pr/rhor);
int j;
	dx=(b-a)/(nx-1);
	dt=(tn-t0)/(nt-1);

	//cout<<sigammaa<<endl;
	//initialize variables
	double *rho;
	double *u;
	double *p;
	double *rhou;
	double *e;
	double *x;
	double *c;
	double *h;
	double epsilon;
	double *uAv,*rhouAv,*rhoAv,*hAv,*cAv,*rhoDel,*rhorhoDel,*eDel,*f1,*f2,*f3,*flux1,*flux2,*flux3,*k1,*k2,*k3,*k4,*k5,*k6,*k7,*k8,*k9,*alpha1Av,*alpha2Av,*alpha3Av,*lambda1,*lambda2,*lambda3;

	double A[3][3];
	
	flux1=new double[nv];
	flux2=new double[nv];
	flux3=new double[nv];
	h=new double[nx];     


	rho=new double [nx];
	u=new double [nx];
	p=new double [nx];
	rhou= new double [nx];
	e=new double [nx];
	x=new double [nx];
	c=new double [nx];
                                                                                                                                                                             
	uAv=new double [nv];
	rhouAv=new double [nv];
	rhoAv=new double [nv];
	hAv= new double [nv];
	cAv=new double [nv];
	rhoDel=new double [nv];
	rhorhoDel=new double [nv];
	eDel=new double [nv];
	f1=new double [nv];
	f2=new double [nv];
	f3=new double [nv];
	k1=new double [nv];	
	k2=new double [nv];
	k3=new double [nv];
	k4=new double [nv];
	k5=new double [nv];
	k6=new double [nv];
	k7=new double [nv];
	k8=new double [nv];
	k9=new double [nv];
	alpha1Av=new double [nv];
	alpha2Av=new double [nv];
	alpha3Av=new double [nv];
	lambda1=new double [nv];
	lambda2=new double [nv];
	lambda3=new double [nv];	
	
	/*double R[3][3];
	double L[3][3];
	*/
	double Afinal[3][3];
	double prodotto[3][3];
	double trans[3][3];
	
	
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
			h[i]=(p[i]+e[i])/rho[i];
			
		}

j=0;
	ofstream output;
	string file="results.txt";
	output.open(file.c_str());
	for (int n=0;n<nt;n++)
	{
		//Averaged values
		for (int i=0;i<nv;i++)
		{
			uAv[i]=(sqrt(rho[i])*u[i]+sqrt(rho[i+1])*u[i+1])/ (sqrt(rho[i])+sqrt(rho[i+1]));
			rhoAv[i]=sqrt(rho[i]*rho[i+1]);
			hAv[i]=(sqrt(rho[i])*h[i]+sqrt(rho[i+1])*h[i+1])/ (sqrt(rho[i])+sqrt(rho[i+1]));
			cAv[i]=sqrt((gamma-1)*(hAv[i]-0.5*pow(uAv[i],2)));


		
		}
		//Gradients
		for (int i=0;i<nv;i++)
		{
			
			rhoDel[i]=rho[i+1]-rho[i];

			rhorhoDel[i]=rhou[i+1]-rhou[i];
			eDel[i]=e[i+1]-e[i];
			//cout<<rhorhoDel[i]<<endl;


		}
		//cell centers fluxes
		for (int i=0;i<nv;i++)
		{
			f1[i]=u[i]*rho[i];
			f2[i]=rho[i]*u[i]*u[i]+p[i];
			f3[i]=rho[i]*u[i]*h[i];

			lambda1[i]=abs(uAv[i]);
			lambda2[i]=abs(uAv[i]+cAv[i]);
			lambda3[i]=abs(uAv[i]-cAv[i]);
		if (abs(lambda1[i])<=TOL)
			{ lambda1[i]=(pow(lambda1[i],2)+pow(TOL,2))/2/TOL;}
			else if(abs(lambda2[i])<=TOL)
			{ lambda2[i]=(pow(lambda2[i],2)+pow(TOL,2))/2/TOL;}
			else if(abs(lambda3[i])<=TOL)
			{ lambda3[i]=(pow(lambda3[i],2)+pow(TOL,2))/2/TOL;}
			
		}	

		//eigenvectors
		for (int i=0;i<nv;i++)
		{
			k1[i]=1;
			k2[i]=uAv[i];
			k3[i]=uAv[i]*uAv[i]/2;
			k4[i]=rhoAv[i]/(2*cAv[i]);
			k5[i]=rhoAv[i]/(2*cAv[i])*(uAv[i]+cAv[i]);
			k6[i]=rhoAv[i]/(2*cAv[i])*(hAv[i]-uAv[i]*cAv[i]);
			k7[i]=-rhoAv[i]/(2*cAv[i]);
			k8[i]=-rhoAv[i]/(2*cAv[i])*(uAv[i]-cAv[i]);
			k9[i]=-rhoAv[i]/(2*cAv[i])*(hAv[i]-uAv[i]*cAv[i]);

		}
		//roe waves 
		
		for (int i=0;i<nv;i++)
		{	
			double L[3][3]={lambda1[i], 0,0,0,lambda2[i],0,0,0,lambda3[i]};
			double R[3][3]={k1[i],k4[i],k7[i],k2[i],k5[i],k8[i],k3[i],k6[i],k9[i]};
			
			    for(int m=0; m<3; m++)
			    		{
        				for(int j=0; j<3; j++)
        					{

            					 prodotto[m][j] = 0;

            				
            						for(int k=0; k<3; k++)
            						{
                						prodotto[m][j] += L[m][k] * R[k][j];
            						}
        					} 
					}
					 
        				float determinant = 0;
    
    
    					//finding determinant
   					 for(int m = 0; m < 3; m++)
   					 {
        			determinant = determinant + (R[0][m] * (R[1][(m+1)%3] * R[2][(m+2)%3] - R[1][(m+2)%3] * R[2][(m+1)%3]));
    
    					 }
    for(int m = 0; m < 3; m++)
    {
        for(int j = 0; j < 3; j++)
        {
            trans[m][j]=((R[(j+1)%3][(m+1)%3] * R[(j+2)%3][(m+2)%3]) - (R[(j+1)%3][(m+2)%3] * R[(j+2)%3][(m+1)%3]))/ determinant;

        
        }
        
    }

					for(int m=0; m<3; m++)
			    		{
        				for(int j=0; j<3; j++)
        					{
            					 Afinal[m][j] = 0;
            				
            						for(int k=0; k<3; k++)
            						{
                						Afinal[m][j] += prodotto[m][k] * trans[k][j];
            						}
        					}
					}
					/*cout<<"matrix"<<endl;
					for(int m=0; m<3; m++)
        					{
							for(int j=0; j<3; j++)
            								{
                							cout<<Afinal[m][j]<<" ";
                							
            								}
            								cout<<endl;
        					}*/
        	flux1[i]=0.5*(f1[i]+f1[i+1])-0.5*(Afinal[1][1]*rhoDel[i]+Afinal[1][2]*rhorhoDel[i]+Afinal[1][3]*eDel[i]);	
        	flux2[i]=0.5*(f2[i]+f2[i+1])-0.5*(Afinal[2][1]*rhoDel[i]+Afinal[2][2]*rhorhoDel[i]+Afinal[2][3]*eDel[i]);
        	flux3[i]=0.5*(f3[i]+f3[i+1])-0.5*(Afinal[3][1]*rhoDel[i]+Afinal[3][2]*rhorhoDel[i]+Afinal[3][3]*eDel[i]);	
        	//cout<<flux1[i]<<endl;		
		}		
		for (int i=1;i<nv;i++)
		{
			rho[i]=rho[i]-dt/dx*(flux1[i]-flux1[i-1]);
			rhou[i]=rhou[i]-dt/dx*(flux2[i]-flux2[i-1]);
			u[i]=rhou[i]/rho[i];
			e[i]=e[i]-dt/dx*(flux3[i]-flux3[i-1]);
			p[i]=(gamma-1)*(e[i]-0.5*(rho[i]*pow(u[i],2)));
			c[i]=pow(gamma*p[i]/rho[i],0.5);
			h[i]=e[i]/rho[i]+p[i]/rho[i];
		}
		


		
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


