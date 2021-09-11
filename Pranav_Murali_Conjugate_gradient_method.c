#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define H 1.0
#define B 1.0
#define nx 101
#define ny 101
#define pi 3.141592653589793
double x[nx*ny],y[nx*ny],b[nx*ny],dx,dy,beta_matrix,actual[nx][ny];
double Mp[nx*ny],Me[nx*ny],Mw[nx*ny],Mn[nx*ny],Ms[nx*ny];
double solver();
double initialization();
double output_files();
double analytical();
int main()
{
    initialization();
	solver();
	analytical();
	output_files();
    return 0;
}
double output_files()
{
	FILE *output1,*output2,*output3; 
    output1=fopen("conjugated.dat","w");
    output2=fopen("Mid_vertical_conjugated.dat","w");
    output3=fopen("Mid_horizontal_conjugated.dat","w");
    dx=H/(double)(nx-1.0);
    dy=B/(double)(ny-1.0);
    int i,j,L,n;n=nx*ny;
    fprintf( output1, "VARIABLES=\"X\",\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
    fprintf(output1,"ZONE F=POINT\n");
	fprintf( output1, "I=%d, J=%d\t\n",nx,ny );
	fprintf( output2, "VARIABLES=\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
    dx=H/(double)(nx-1.0);
    fprintf( output3, "VARIABLES=\"X\",\"coded_Temp\",\"analytical_Temp\"\n");
    dy=B/(double)(ny-1.0);
	double xp,yp;
    for(j=1;j<=nx;j++)
	{
        for(i=1;i<=ny;i++)
        {
            L=(i-1)*ny+j;
            yp=(j-1)*dy;xp=(i-1)*dx;
            fprintf(output1,"%5.8f\t%5.8f\t%5.8f\t%5.8f\n",xp,yp,x[L],actual[i][j]);
			if(xp==.5)    
			{
				fprintf(output2,"%5.8f\t%5.8f\t%5.8f\n",yp,x[L],actual[i][j]);
			} 
			if(yp==.5)
			{
				fprintf(output3,"%5.8f\t%5.8f\t%5.8f\n",xp,x[L],actual[i][j]);
			}        
        }
    }
    fclose(output1);
    fclose(output2);
    fclose(output3);
}
double solver()
{
	int i,j;int n;n=nx*ny;
	dx=H/(double)(nx-1.0);
    dy=B/(double)(ny-1.0);  
    for(i=1;i<=n;i++)
    {
        x[i]=0.0;
    }
    double p[n],r[n],rold[n],Ap[n],alpha,beta;
	double alpha_num, alpha_den, Ax[n];
    double error;
    int k;
	k=0;
    for (i=1; i<=n; i++)
    {
       	if(i==1)
       	{
       		Ax[i]=Mp[i]*x[i]+Mn[i]*x[i+1]+Me[i]*x[ny+i];
       	}
		else if((i>1)&&(i<=ny))
		{
			Ax[i]=Ms[i]*x[i-1]+ Mp[i]*x[i]
			        +Mn[i]*x[i+1]+Me[i]*x[i+ny];
		}
		else if((i>ny)&&(i<=(n-ny)))
		{
			Ax[i]=Mw[i]*x[i-ny]+ Ms[i]*x[i-1]+Mp[i]*x[i]
			        +Mn[i]*x[i+1]+Me[i]*x[i+ny];
		}
		else if((i>ny)&&(i<(n)))
		{
			Ax[i]=Mw[i]*x[i-ny]+Ms[i]*x[i-1]+Mp[i]*x[i]
			        +Mn[i]*x[i+1]+Me[i]*x[i+ny];
		}
		else if(i==n)
		{
			Ax[i]=Mw[i-1]*x[i-ny-1]+Ms[i-1]*x[i-2]
			        +Mp[i-1]*x[i-1];
		}
		p[i]=b[i]-Ax[i];
		r[i]=b[i]-Ax[i];
	}
	error=1.0;
    do
    {
        error=0;
		alpha_num=0;
		alpha_den=0;
		for (i=1; i<=n; i++)
        {
        	if(i==1)
        	{
        		Ap[i]=Mp[i]*p[i]+Mn[i]*p[i+1]+Me[i]*p[ny+i];
				alpha_num+=r[i]*r[i];     
				alpha_den+=p[i]*Ap[i];
			}
			else if((i>1)&&(i<=ny))
			{
				Ap[i]=Ms[i]*p[i-1]+Mp[i]*p[i]
				        +Mn[i]*p[i+1]+Me[i]*p[i+ny];
				alpha_num+=r[i]*r[i];  
				alpha_den+=p[i]*Ap[i];	  	
			}
			else if((i>ny)&&(i<=(n-ny)))
			{
				Ap[i]=Mw[i]*p[i-ny]+Ms[i]*p[i-1]+Mp[i]*p[i]
				        +Mn[i]*p[i+1] + Me[i]*p[i+ny];
				alpha_num+=r[i]*r[i];     
				alpha_den+=p[i]*Ap[i]; 
			}
			else if((i>(n-ny))&&(i<(n)))
			{
				Ap[i]=Mw[i]*p[i-ny]+Ms[i]*p[i-1]+Mp[i]*p[i]
				        +Mn[i]*p[i+1];
				alpha_num+=r[i]*r[i];     
				alpha_den+=p[i]*Ap[i]; 
			}
			else if(i==n)
			{
				Ap[i]=Mw[i]*p[i-ny]+Ms[n]*p[i-1]+Mp[i]*p[i];
				alpha_num+=r[i]*r[i];  
				alpha_den+=p[i]*Ap[i];
			}
		}
        alpha=alpha_num/alpha_den;
		double beta_num,beta_den;
        beta_num=0; 
		beta_den=0;
        for (i=1; i<=n; i++)
        {
            x[i]=x[i]+alpha*p[i];
			rold[i]=r[i];
			beta_den+=rold[i]*rold[i];
            r[i]=r[i]-alpha*Ap[i];
            beta_num+=r[i]*r[i];
            error+=r[i]*r[i];
        }
        beta=beta_num/beta_den;
        for (i=1; i<=n; i++)
        {
            p[i]=r[i]+beta*p[i];
        }
        k++;
        printf("%d\t\t%5.8f\n",k,error);
    }while (error>1.0e-14);   
}
double analytical()
{
	int i,j;double xp,yp;
	dx=H/(double)(nx-1.0);
    dy=B/(double)(ny-1.0);
	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			xp=(i-1)*dx;
			yp=(j-1)*dy;    
			actual[i][j]=sin(2.0*pi*xp)*sin(2.0*pi*yp);
		}
	}
}
double initialization()
{
	int i,j,L,n;
	n=nx*ny;
	double xp,yp;
	dx=H/(double)(nx-1.0);
    dy=B/(double)(ny-1.0);
	beta_matrix=dx/dy;
    for(i=1;i<=nx;i++)
    {
        for(j=1;j<=ny;j++)
        {
            L=(i-1)*ny+j; 
           if(i==1)
           {
               Mp[L]=1.0;
               Me[L]=0.0;
			   Mw[L]=0.0;
			   Mn[L]=0.0;
			   Ms[L]=0.0;
			   b[L]=0.0;

           }
		   else if(i==nx)
		   {
			   Mp[L]=1.0;
               Me[L]=0.0;
			   Mw[L]=0.0;
			   Mn[L]=0.0;
			   Ms[L]=0.0;
			   b[L]=0.0;
		   }
		   else if(j==1)
		   {
			   Mp[L]=1.0;
               Me[L]=0.0;
			   Mw[L]=0.0;
			   Mn[L]=0.0;
			   Ms[L]=0.0;
			   b[L]=0.0;
		   }
		   else if(j==ny)
		   {
			   Mp[L]=1.0;
               Me[L]=0.0;
			   Mw[L]=0.0;
			   Mn[L]=0.0;
			   Ms[L]=0.0;
			   b[L]=0.0;
		   }
           else
           {
               xp=(i-1)*dx;   
			   yp=(j-1)*dy;
               Mp[L]=-2.0*(1+beta_matrix*beta_matrix);
			   Mw[L]=beta_matrix*beta_matrix;
			   Mn[L]=1.0;
			   Ms[L]=1.0;
			   Me[L]=beta_matrix*beta_matrix;
               b[L]=-8.0*pi*pi*(sin(2.0*pi*xp)*sin(2.0*pi*yp))*dx*dx;
           }
        }
    }
}
