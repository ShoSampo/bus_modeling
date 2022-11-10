#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
/* TASEP(L sites ,N particles) with PBC */
/* Calculate current J vs N/L  */
/*  a.out T  */
void main(int argc, char *argv[]);
void Initial_Config_X(int L, int N, int *X);
void Initial_Config_Y(int L, int M, int *Y);
int Dynamics_X(int L, int N, int *X);
int Dynamics_Y(int L, int N, int M, int *X, int *Y);
double Jy0y1(int L, int N, int M, int T, double *v_x0_mean, int *x0, double *v_x1_mean, int *x1, double *v_y0_mean, int *y0, double *v_y1_mean, int *y1);
void QSort(int d[ ], int left, int right);


void main(int argc, char *argv[])
{
  int  T,L,N,M;
  int i,t;
  int seed;
  double rho_x,rho_y,j_x,j_y;
  double v_x0_mean,v_x1_mean,v_y0_mean,v_y1_mean;
  static int x0[20][1000001],x1[20][1000001],y0[20][1000001],y1[20][1000001];

  char filename[100];
  FILE *fp,*fp1,*fp2;

    /*   Initialization  */
    T=atoi(argv[1]);
    L=100;
    M=2;
    printf("T=%d L=%d  \n",T,L);
    seed=10;
    init_genrand(seed); /* seed の設定 */

    sprintf(filename,"rho_vs_J_v_y0_mean_v_y1_mean_L%d.csv",L);
    fp=fopen(filename,"w");
    fprintf(fp,"rho,J,v_y0_mean,v_y1_mean \n");
    sprintf(filename,"t_N_vs_y0_L%d.csv",L);
    fp1=fopen(filename,"w");
    fprintf(fp1,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");
    sprintf(filename,"t_N_vs_y1_L%d.csv",L);
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
      rho_x=(double)(N*1.0)/L;
      rho_y=(double)(M*1.0)/L;
      j_y=Jy0y1(L,N,M,T,&v_x0_mean,&x0[i][0],&v_x1_mean,&x1[i][0],&v_y0_mean,&y0[i][0],&v_y1_mean,&y1[i][0])/L;
      printf("N=%d rho_x=%lf J=%lf v_y0_mean=%lf v_y1_mean=%lf  \n",N,rho_x,j_y,v_y0_mean,v_y1_mean);
      fprintf(fp,"%lf,%lf,%lf,%lf \n",rho_x,j_y,v_y0_mean,v_y1_mean);
    }
    fclose(fp);

    for(t=0;t<=T;t++){
      fprintf(fp1,"%d",t);
      for(i=0;i<20;i++){
        fprintf(fp1,",%d",y0[i][t]);
      }
      fprintf(fp1,"\n");
    }
    fclose(fp1);

    for(t=0;t<=T;t++){
      fprintf(fp2,"%d",t);
      for(i=0;i<20;i++){
        fprintf(fp2,",%d",y1[i][t]);
      }
      fprintf(fp2,"\n");
    }
    fclose(fp2);

  return;
}


double Jy0y1(int L, int N, int M, int T, double *v_x0_mean, int *x0, double *v_x1_mean, int *x1, double *v_y0_mean, int *y0, double *v_y1_mean, int *y1)
{
  int t,n;
  int move_x,x0_ini,x1_ini;
  int move_y,y0_ini,y1_ini;
  int X[1000],Y[1000];
  int move_y_sum,move_x_sum;
  double j_y;

    Initial_Config_X(L,N,X);
    Initial_Config_Y(L,M,Y);
    move_x_sum=0;
     *x0=0;
     *x1=X[1]-X[0];
    move_y_sum=0;
     *y0=0;
     *y1=Y[1]-Y[0];

    for(t=1;t<=T;t++){
      x0_ini=X[0];
      x1_ini=X[1];
      move_x=Dynamics_X(L,N,X);
      move_x_sum+=move_x;

      y0_ini=Y[0];
      y1_ini=Y[1];
      move_y=Dynamics_Y(L,N,M,X,Y);
      move_y_sum+=move_y;

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

      if(Y[0]==y0_ini){
        *(y0+t)=*(y0+t-1);
      }
      else{
        *(y0+t)=*(y0+t-1)+(Y[0]-y0_ini+L)%L;
      }

      if(Y[1]==y1_ini){
        *(y1+t)=*(y1+t-1);
      }
      else{
        *(y1+t)=*(y1+t-1)+(Y[1]-y1_ini+L)%L;
      }

    }
    j_y=(double)move_y_sum*1.0/T;
    *v_y0_mean=(double)*(y0+T)/T;
    *v_y1_mean=(double)(*(y1+T)-*y1)/T;
return(j_y);
}


int Dynamics_X(int L, int N, int *X)
{
  int n, i, j, occupied, x_next;
  int move_x;
      move_x=0;

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
          move_x+=1;
        }
      }

  return(move_x);
}


int Dynamics_Y(int L,int N, int M, int *X, int *Y)
{
  int n, i, j, k, occupied, y_next;
  int move_y;
      move_y=0;

      for(i=0;i<M;i++){
        n=(int)(Uniform()*M);
        y_next=(*(Y+n)+1)%L;
        occupied=0;
        for(j=0;j<M;j++){
          if(*(Y+j)==y_next){
            occupied=1;
          }
        }
        for(k=0;k<N;k++){
          if(*(X+k)==y_next){
            occupied=1;
          }
        }
        if(occupied==0){
          *(Y+n)=y_next;
          move_y+=1;
        }
      }

  return(move_y);
}


void Initial_Config_X(int L, int N, int *X)
{
  int i,j,d;
  int X_ini[1000];
      d=Uniform()*L;
      X_ini[0]=d;
      for(i=1;i<N;i++){
        d=Uniform()*(L-i);
      for(j=0;j<i;j++){
        if(X_ini[j]<=d){
          d+=1;
        }
      }
      X_ini[i]=d;
      }

      QSort(X_ini,0,N-1);

      for(i=0;i<N;i++){
        *(X+i)=X_ini[i];
      }
  return;
}


void Initial_Config_Y(int L, int M, int *Y)
{
  int i,j,d;
  int Y_ini[1000];
      d=Uniform()*L;
      Y_ini[0]=d;
      for(i=1;i<M;i++){
        d=Uniform()*(L-i);
      for(j=0;j<i;j++){
        if(Y_ini[j]<=d){
          d+=1;
        }
      }
      Y_ini[i]=d;
      }

      QSort(Y_ini,0,M-1);

      for(i=0;i<M;i++){
        *(Y+i)=Y_ini[i];
      }
  return;
}


/*  小さい順にデータdを並べ替える */
void QSort(int d[ ], int left, int right)
{
    int i, j;
    int pivot,temp;

    i = left;
    j = right;

    pivot = d[(left + right) / 2];

    while (1) {

        while (d[i] < pivot)
            i++;

        while (pivot < d[j])
            j--;
        if (i >= j)
            break;

        temp=d[i];
        d[i]=d[j];
        d[j]=temp;
        i++;
        j--;
    }

    if (left < i - 1)
        QSort(d, left, i - 1);
    if (j + 1 <  right)
        QSort(d, j + 1, right);
}
