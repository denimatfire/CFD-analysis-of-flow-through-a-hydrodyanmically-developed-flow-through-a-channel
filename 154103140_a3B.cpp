#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
using namespace std;

int main()
{int n,m,i,j,Re,l,h;
   float dx,dy,b,sum,serror=10,c1,c2,c3,c4,werror=10,iterr=0,x,y,w=0.01,U=-1;

    clock_t begin,end;
    double cpu_time_used;

    ofstream f1,f2,f3,f4,f5,f6,f7;
    f1.open("plot.plt");
   
   cout<<"Assignment 3B\n";
   cout<<"Enter number of grid point along x-axis-";
   cin>>n;
   cout<<"\nEnter number of grid points along y-axis-";
   cin>>m;
   cout<<"\nEnter the Reynolds number";
   cin>>Re;
   cout<<"\nEnter the length L-";
   cin>>l;
   cout<<"\nEnter the Height H-";
   cin>>h;

  float** u=new float*[n];
   for(i=0;i<n;i++)
        u[i]= new float[m];
   float** v=new float*[n];
   for(i=0;i<n;i++)
        v[i]= new float[m];
   float** w0=new float*[n];
   for(i=0;i<n;i++)
        w0[i]=new float[m];
   float** w1=new float*[n];
   for(i=0;i<n;i++)
        w1[i]=new float[m];
    float** s1=new float*[n];
   for(i=0;i<n;i++)
        s1[i]=new float[m];
    float** s0=new float*[n];
   for(i=0;i<n;i++)
        s0[i]=new float[m];


   dx=((float)l/(n-1));
   dy=((float)h/(m-1));
   b=dx/dy;

   //u,v boundary conditions//
    for(i=0;i<n;i++){
        u[i][0]=0;
        v[i][0]=0;//bottom wall
        u[i][m-1]=0;
        v[i][m-1]=0;//top wall
    }

    for(j=0;j<m;j++){
        u[n-1][j]=U;//right wall
        v[n-1][j]=0;
    }
    for(j=0;j<(m/2);j++){
        u[0][j]=0;//left wall
        v[0][j]=0;
    }
    //initial values of stream and vorticity
    for(i=1;i<n-1;i++)
		for(j=1;j<m-1;j++)
			{
				s1[i][j]=0.0;
				s0[i][j]=0.0;
			}

	for(i=1;i<n-1;i++)
		for(j=1;j<m-1;j++)
			{
				w1[i][j]=0.0;
				w0[i][j]=0.0;
			}
			float cc=0.000001;
begin =clock();
do{
   //stream function boundary conditions
   for(i=0;i<n;i++){
    s0[i][0]=0;//bottom wall
    s0[i][m-1]=U*h;//top wall
   }
   for(j=0;j<m-1;j++){
    s0[n-1][j+1]=s0[n-1][j]+(U*dy);//right wall
   }
   for(j=0;j<(m/2);j++){
        s0[0][j]=0;//left wall
   }

   for(j=(m/2);j<m;j++){
    s0[0][j]=s0[0][j-1]+(u[0][j]*dy);//2*s0[1][j]-s0[2][j];//left  wall
   }
 
   //vorticity boundary conditions
   for(i=0;i<n;i++){
    w0[i][0]=(2*(s0[i][0]-s0[i][1]))/pow(dy,2);//bottom wall
    w0[i][m-1]=(2*(s0[i][m-1]-s0[i][m-2]))/pow(dy,2);//top wall
   }
   for(j=0;j<m;j++){
    w0[n-1][j]=(2*(s0[m-1][j]-s0[m-2][j]))/pow(dx,2);//right wall
   }
   for(j=0;j<(m/2);j++){
    w0[0][j]=(-2*(s0[1][j]-s0[0][j]))/pow(dx,2);//left wall
   }
   for(j=(m/2);j<m;j++){
    w0[0][j]=w0[1][j];//left
   }
 
   //update stream function
   sum=0;
    for(i=1;i<n-1;i++){
        for(j=1;j<m-1;j++){
            s1[i][j]=(s0[i+1][j]+s0[i-1][j]+(b*b)*(s0[i][j+1]+s0[i][j-1])+w0[i][j]*dx*dx)/(2*(1+b*b));
            sum=sum+(s1[i][j]-s0[i][j]);
            s0[i][j]=s1[i][j];
        }
    }
 
    serror=sqrt((sum*sum)/((m-1)*(n-1)));
    //update u v
    for(i=1;i<n-1;i++){
        for(j=1;j<m-1;j++){
            u[i][j]=(s0[i][j+1]-s0[i][j-1])/(2*dy);
            v[i][j]=(s0[i-1][j]-s0[i+1][j])/(2*dx);
        }
    }
    for(j=(m/2);j<m;j++){
        u[0][j]=u[1][j];
        v[0][j]=v[1][j];
    }
  
    //solve for vorticity and error
    sum=0;
    for(i=1;i<n-1;i++){
        for(j=1;j<m-1;j++){
                c1=(1-(u[i][j]*dx*Re*0.5));
                c2=(1+(u[i][j]*dx*Re*0.5));
                c3=(b*b)*(1-(v[i][j]*dy*Re*0.5));
                c4=(b*b)*(1+(v[i][j]*dy*Re*0.5));
            w1[i][j]=w0[i][j]*(1-w)+w*((c1*w0[i+1][j]+c2*w0[i-1][j]+c3*w0[i][j+1]+c4*w0[i][j-1])/(2*(1+(b*b))));
            sum=sum+(w1[i][j]-w0[i][j]);
            w0[i][j]=w1[i][j];
        }
    }
 
    werror=sqrt((sum*sum)/((m-1)*(n-1)));

    iterr++;
    cout<<"iterration no- "<<iterr<<"\tserror- "<<serror<<"\twerror- "<<werror<<"\n";
    //f4<<"\n"<<iterr<<"\t"<<serror<<"\t"<<werror;
    }while(serror>cc || werror>cc);
end = clock();
cpu_time_used =((double)(end-begin))/CLOCKS_PER_SEC;
cout<<" \ntime used "<<cpu_time_used;

    f1<<"ZONE"<<"\tI="<<n<<"\tJ="<<m<<"\n";
        for(y=0,j=0;j<m;j++,y=y+dy){
            for(x=0,i=0;i<n;i++,x=x+dx){
                f1<<x<<"\t"<<y<<"\t"<<u[i][j]<<"\t"<<v[i][j]<<"\t"<<s0[i][j]<<"\t"<<w0[i][j]<<"\n";
            }
        }

	
   }
