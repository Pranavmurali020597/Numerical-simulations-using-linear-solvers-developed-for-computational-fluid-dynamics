#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define nx 100
#define ny 100
#define l 1
#define W 1
#define pi 3.14159265
double Tnew[nx*ny],Told[nx*ny],beta,b[nx*ny],Mp[nx*ny],Mnw[nx*ny],Mw[nx*ny],Mn[nx*ny],Me[nx*ny],Ms[nx*ny],Mse[nx*ny];
double Lw[nx*ny],Lp[nx*ny],Ls[nx*ny],Un[nx*ny],Ue[nx*ny],N,dx,dy;
double ILU();
double output_files();
double analytic();
double initializtion();
double solver();
int main()
{
	initializtion();
	ILU();
	solver();
	analytic();
	output_files();
	return 0;
}
double initializtion()
{
	int i,j,L,N;N=nx*ny;
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			L=(i-1)*ny+j;
			Tnew[L]=0.0;
			Told[L]=0.0;
			b[L]=0.0;
		}
		
	}
}
double ILU()
{
	int i,j,L;N=nx*ny;
	int Nj;
	Nj=ny;dx=l/(double)(nx-1);
	dy=W/(double)(ny-1);
	beta=dx/dy;
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			L=(i-1)*ny+j;
			if(i==1)
			{
				Mp[L]=1.0;
				Mw[L]=0.0;
				Mn[L]=0.0;
				Ms[L]=0.0;
				Me[L]=0.0;
				b[L]=0.0;
			}
			else if((i==nx)||(j==1))
			{
				Mp[L]=1.0;
				Mw[L]=0.0;
				Mn[L]=0.0;
				Ms[L]=0.0;
				Me[L]=0.0;
				b[L]=0.0;
			}
			else if(j==ny)
			{
				Mp[L]=1.0;
				Mw[L]=0.0;
				Mn[L]=0.0;
				Ms[L]=0.0;
				Me[L]=0.0;
				b[L]=1.0;
			}
			else
			{
				Mp[L]=-2.0*(1+beta*beta);
				Mw[L]=beta*beta;
				Mn[L]=1.0;
				Ms[L]=1.0;
				Me[L]=beta*beta;
				b[L]=0.0;
			}
		}
	
	}

	for(L=1;L<=N;L++)
	{
			
		Lw[L]=Mw[L];
		
		Ls[L]=Ms[L];
		Lp[L]=Mp[L];
		if((L-1)>=1)
		{
			Lp[L]=Lp[L]-Ls[L]*Un[L-1];
		}
		
		if((L-Nj)>=1)
		{
			Lp[L]=Lp[L]-Lw[L]*Ue[L-Nj];
		}
		Un[L]=Mn[L]/Lp[L];
		
		Ue[L]=Me[L]/Lp[L];
	}
	for(L=1;L<=(N);L++)
	{
		Mnw[L]=0.0;
		if((L-Nj)>=1)
		{
			Mnw[L]=Lw[L]*Un[L-Nj];
		}
		
	}
	for(L=1;L<=(N);L++)
	{
		Mse[L]=0.0;
		if((L-1)>=1)
		{
		Mse[L]=Ls[L]*Ue[L-1];
		}
	}
}
double output_files()
{
	FILE *oup,*output2,*output3;double dx,dy,x,y;dx=l/(double)(nx-1);
	dy=l/(double)(ny-1);
	oup=fopen("ILU.dat","w");
	output2=fopen("ILU_vertical.dat","w");
	output3=fopen("ILU_horizontal.dat","w");
	int i,j,L;
	fprintf( oup, "VARIABLES=\"X\",\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	fprintf(oup,"ZONE F=POINT\n");
	fprintf( oup, "I=%d, J=%d\t\n", nx, ny );
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			L=(i-1)*ny+j;
			fprintf(oup,"%lf\t%lf\t%lf\t%lf\n",x,y,Tnew[L],Told[L]);
		}
			
	}
	int x1, x2,y1, y2,L1,L2;double avg1,avg2,avg,ana;
    x1=nx/2;
	x2=(nx/2)-1;
	fprintf( output2, "VARIABLES=\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	for(j=1;j<=ny;j++)
	{
		L1=(x1-1)*ny+j;
		avg1=Tnew[L1];
		L2=(x2-1)*ny+j;
		avg2=Tnew[L2];
		ana=(Told[L1]+Told[L2])/2;
		avg=(avg1+avg2)/2;
		y=(j-1)*dy;
		fprintf(output2,"%5.8f\t%5.8f\t%5.8f\n",y,avg,ana);
	}
	y1=ny/2;
	y2=(ny/2)-1;
	fprintf( output3, "VARIABLES=\"X\",\"coded_Temp\",\"analytical_Temp\"\n");
	for(i=1;i<=nx;i++)
	{
		L1=(i-1)*ny+y1;
		avg1=Tnew[L1];
		L2=(i-1)*ny+y2;
		avg2=Tnew[L2];
		ana=(Told[L1]+Told[L2])/2;
		avg=(avg1+avg2)/2;
		x=(i-1)*dx;
		fprintf(output3,"%5.8f\t%5.8f\t%5.8f\n",x,avg,ana);
	}
	fclose(oup);

}
double solver()
{
	double Y[nx*ny],e,en,ed;
	int N,i,j,L;
	N=nx*ny;
	e=1.0;
	while(e>1.0e-14)
	{
		
		Y[1]=(b[1]+Mse[1]*Told[1+(ny-1)])/Lp[1];
		for(L=2;L<=N;L++)
		{
			Y[L]=b[L]-Ls[L]*Y[L-1];
			if((L-(ny-1))>=1)
			{
				Y[L]=Y[L]+Mnw[L]*Told[L-(ny-1)];
			}
			if((L+(ny-1))<=N)
			{
				Y[L]=Y[L]+Mse[L]*Told[L+(ny-1)];
			}
			if((L-ny)>=1)
			{
				Y[L]=Y[L]-Lw[L]*Y[L-ny];
			}
			Y[L]=Y[L]/Lp[L];	
		}
		Tnew[N]=Y[N];
		for(L=N-1;L>=1;L--)
		{
			Tnew[L]=Y[L]-Un[L]*Tnew[L+1];
			if((L+ny)<=N)
			{
				Tnew[L]=Tnew[L]-Ue[L]*Tnew[L+ny];	
			}
			
		}
		en=0.0;
		ed=0.0;
		for(i=1;i<=N;i++)
		{
			en+=fabs(Tnew[i]-Told[i]);
			ed+=fabs(Tnew[i]);
		}
		e=en/ed;
		for(i=1;i<=N;i++)
		{
			Told[i]=Tnew[i];
			Y[i]=0.0;
		}printf("%lf\n",e);	
	}

}
double analytic()
{
int i,j,k,L;double dx,dy,x,y,a[nx*ny],sum[nx*ny],e,term;
	dx=l/(double)(nx-1);
	dy=l/(double)(ny-1);
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			L=(i-1)*ny+j;
			x=dx*(i-1);
			y=dy*(j-1);
			e=1.0;
			sum[L]=0.0;k=1;
			for(k=1;k<=110;k++)
			{
				sum[L]+=((pow(-1,k+1)+1)/k)*(sin(k*pi*x))*sinh(k*pi*y)/(sinh(k*pi))*2/pi;
			}
			Told[L]=sum[L];
		}
	}
}
		
		
	
