#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
/* TASEP(L sites ,N particles) with PBC */
/* Calculate current J vs N/L  */
/*  a.out T  */
void main(int argc, char *argv[]);
void Initial_Config(int L, int N, int *X);
int Dynamics(int L, int N, int *X);
double J(int L, int N, int T);
double Jx0(int L, int N, int T, double *v0_mean, int *x0);
double Jx0x1(int L, int N, int T, double *v0_mean, int *x0, double *v1_mean, int *x1);
void QSort(int x[ ], int left, int right);


void main(int argc, char *argv[])
{
  int  T,L,N;
  int i,t;
  int seed;
  double rho,j,v0_mean,v1_mean;
  static int x0[20][1000001],x1[20][1000001];

  char filename[100];
  FILE *fp,*fp1,*fp2;

     /*   Initialization  */
     T=atoi(argv[1]);
     L=100;
     printf("T=%d L=%d  \n",T,L);
     seed=10;
     init_genrand(seed); /* seed の設定 */

     sprintf(filename,"rho_vs_J_v0_mean_v1_mean_L%d.csv",L);
     fp=fopen(filename,"w");
     fprintf(fp,"rho,J,v0_mean,v1_mean \n");
     sprintf(filename,"t_N_vs_x0_L%d.csv",L);
     fp1=fopen(filename,"w");
     fprintf(fp1,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");
     sprintf(filename,"t_N_vs_x1_L%d.csv",L);
     fp2=fopen(filename,"w");
     fprintf(fp2,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");
     /* MC Simulation */

     for(i=0;i<20;i++){
     if(i==0){
     N=2;
     }
     else{
     N=5*i;
     }
     rho=(double)(N*1.0)/L;
     j=Jx0x1(L,N,T,&v0_mean,&x0[i][0],&v1_mean,&x1[i][0])/L;
     printf("N=%d rho=%lf J=%lf v0_mean=%lf v1_mean=%lf  \n",N,rho,j,v0_mean,v1_mean);
     fprintf(fp,"%lf,%lf,%lf,%lf \n",rho,j,v0_mean,v1_mean);
     }
     fclose(fp);

     for(t=0;t<=T;t++){
     fprintf(fp1,"%d",t);
     for(i=0;i<20;i++){
     fprintf(fp1,",%d",x0[i][t]);
     }
     fprintf(fp1,"\n");
     }
     fclose(fp1);

     for(t=0;t<=T;t++){
     fprintf(fp2,"%d",t);
     for(i=0;i<20;i++){
     fprintf(fp2,",%d",x1[i][t]);
     }
     fprintf(fp2,"\n");
     }
     fclose(fp2);

  return;
}

double Jx0x1(int L, int N, int T, double *v0_mean, int *x0, double *v1_mean, int *x1)
{
   int t,n;
   int move,x0_ini,x1_ini;
   int X[1000];
   int move_sum;
   double j;

    Initial_Config(L,N,X);
    move_sum=0;
     *x0=0;
     *x1=X[1]-X[0];

     for(t=1;t<=T;t++){
       x0_ini=X[0];
       x1_ini=X[1];
       move=Dynamics(L,N,X);
       move_sum+=move;

       if(X[0]==x0_ini){
         *(x0+t)=*(x0+t-1);
       }
       else{
         *(x0+t)=*(x0+t-1)+(X[0]-x0_ini+L)%L;
       }

       if(X[1]==x1_ini){
         *(x1+t)=*(x1+t-1);
       }
       else{
         *(x1+t)=*(x1+t-1)+(X[1]-x1_ini+L)%L;
       }

     }
     j=(double)move_sum*1.0/T;
     *v0_mean=(double)*(x0+T)/T;
     *v1_mean=(double)(*(x1+T)-*x1)/T;
return(j);
}

double Jx0(int L, int N, int T, double *v0_mean, int *x0)
{
   int t,n;
   int move,x0_ini;
   int X[1000];
   int move_sum;
   double j;

    Initial_Config(L,N,X);
    move_sum=0;
     *x0=0;
     for(t=1;t<=T;t++){
       x0_ini=X[0];
       move=Dynamics(L,N,X);
       move_sum+=move;
       if(X[0]==x0_ini){
         *(x0+t)=*(x0+t-1);
       }
       else{
         *(x0+t)=*(x0+t-1)+(X[0]-x0_ini+L)%L;
       }
     }
     j=(double)move_sum*1.0/T;
     *v0_mean=(double)*(x0+T)/T;

return(j);
}

double J(int L, int N, int T)
{
   int t,n;
   int move;
   int X[1000];
   int move_sum;
   double j;

    Initial_Config(L,N,X);
    move_sum=0;
     for(t=0;t<T;t++){
       move=Dynamics(L,N,X);
       move_sum+=move;
     }
     j=(double)move_sum*1.0/T;

return(j);
}



int Dynamics(int L, int N, int *X)
{
  int n,i,j,occupied,x_next;
  int move;
      move=0;

      for(i=0;i<N;i++){
      n=(int)(Uniform()*N);
      x_next=(*(X+n)+1)%L;
      occupied=0;
        for(j=0;j<N;j++){
          if(*(X+j)==x_next){
            occupied=1;
          }
        }
        if(occupied==0){
          *(X+n)=x_next;
          move+=1;
        }
      }

  return(move);
}


void Initial_Config(int L, int N, int *X)
{
  int i,j,x;
  int X_ini[1000];
      x=Uniform()*L;
      X_ini[0]=x;
      for(i=1;i<N;i++){
        x=Uniform()*(L-i);
      for(j=0;j<i;j++){
        if(X_ini[j]<=x){
          x+=1;
        }
      }
      X_ini[i]=x;
      }

      QSort(X_ini,0,N-1);

      for(i=0;i<N;i++){
        *(X+i)=X_ini[i];
      }
  return;
}


/*  小さい順にデータxを並べ替える */
void QSort(int x[ ], int left, int right)
{
    int i, j;
    int pivot,temp;

    i = left;
    j = right;

    pivot = x[(left + right) / 2];

    while (1) {

        while (x[i] < pivot)
            i++;

        while (pivot < x[j])
            j--;
        if (i >= j)
            break;

        temp=x[i];
        x[i]=x[j];
        x[j]=temp;
        i++;
        j--;
    }

    if (left < i - 1)
        QSort(x, left, i - 1);
    if (j + 1 <  right)
        QSort(x, j + 1, right);
}
