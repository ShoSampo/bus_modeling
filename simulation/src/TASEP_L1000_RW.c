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
void Dynamics_Coupled(int L, int N, int M, int *X, int *Y, int *move_x, int *move_y);
void J_v(int L, int N, int M, int T, double *v_x0_mean, int *x0, double *v_x1_mean, int *x1, double *v_y0_mean, int *y0, double *v_y1_mean, int *y1, double *j_x, double *j_y);
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
  FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5;

    /*   Initialization  */
    T=atoi(argv[1]);
    L=1000;
    M=2;
    printf("T=%d L=%d  \n",T,L);
    seed=10;
    init_genrand(seed); /* seed の設定 */

    sprintf(filename,"rho_vs_J_v_x0_mean_v_x1_mean_L%d.csv",L);
    fp=fopen(filename,"w");
    fprintf(fp,"rho,J,v_x0_mean,v_x1_mean \n");
    sprintf(filename,"t_N_vs_x0_L%d.csv",L);
    fp1=fopen(filename,"w");
    fprintf(fp1,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");
    sprintf(filename,"t_N_vs_x1_L%d.csv",L);
    fp2=fopen(filename,"w");
    fprintf(fp2,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");

    sprintf(filename,"rho_vs_J_v_y0_mean_v_y1_mean_L%d.csv",L);
    fp3=fopen(filename,"w");
    fprintf(fp3,"rho,J,v_y0_mean,v_y1_mean \n");
    sprintf(filename,"t_N_vs_y0_L%d.csv",L);
    fp4=fopen(filename,"w");
    fprintf(fp4,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");
    sprintf(filename,"t_N_vs_y1_L%d.csv",L);
    fp5=fopen(filename,"w");
    fprintf(fp5,"t,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95 \n");

    /* MC Simulation */

    for(i=0;i<20;i++){
      if(i==0){
        N=20;
      }
      else{
        N=50*i;
      }
      rho_x=(double)(N*1.0)/L;
      rho_y=(double)(M*1.0)/L;
      J_v(L,N,M,T,&v_x0_mean,&x0[i][0],&v_x1_mean,&x1[i][0],&v_y0_mean,&y0[i][0],&v_y1_mean,&y1[i][0],&j_x,&j_y);
      printf("N=%d rho_x=%lf J_x=%lf J_y=%lf v_x0_mean=%lf v_x1_mean=%lf v_y0_mean=%lf v_y1_mean=%lf  \n",N,rho_x,j_x,j_y,v_x0_mean,v_x1_mean,v_y0_mean,v_y1_mean);
      fprintf(fp,"%lf,%lf,%lf,%lf \n",rho_x,j_x,v_x0_mean,v_x1_mean);
      fprintf(fp3,"%lf,%lf,%lf,%lf \n",rho_x,j_y,v_y0_mean,v_y1_mean);
    }
    fclose(fp);
    fclose(fp3);

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

    for(t=0;t<=T;t++){
      fprintf(fp4,"%d",t);
      for(i=0;i<20;i++){
        fprintf(fp4,",%d",y0[i][t]);
      }
      fprintf(fp4,"\n");
    }
    fclose(fp4);

    for(t=0;t<=T;t++){
      fprintf(fp5,"%d",t);
      for(i=0;i<20;i++){
        fprintf(fp5,",%d",y1[i][t]);
      }
      fprintf(fp5,"\n");
    }
    fclose(fp5);

  return;
}


void J_v(int L, int N, int M, int T, double *v_x0_mean, int *x0, double *v_x1_mean, int *x1, double *v_y0_mean, int *y0, double *v_y1_mean, int *y1, double *j_x, double *j_y)
{
  int t,n;
  int x0_ini,x1_ini;
  int y0_ini,y1_ini;
  int X[1000],Y[1000];
  int move_x, move_y;
  int move_y_sum,move_x_sum;

    Initial_Config_X(L,N,X);
    Initial_Config_Y(L,M,Y);
    move_x=0;
    move_x_sum=0;
     *x0=0;
     *x1=X[1]-X[0];
    move_y=0;
    move_y_sum=0;
     *y0=0;
     *y1=Y[1]-Y[0];

    for(t=1;t<=T;t++){
      x0_ini=X[0];
      x1_ini=X[1];
      y0_ini=Y[0];
      y1_ini=Y[1];
      Dynamics_Coupled(L, N, M, X, Y, &move_x, &move_y);
      move_x_sum+=move_x;
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
    *j_x=(double)move_x_sum*1.0/(T*L);
    *v_x0_mean=(double)*(x0+T)/T;
    *v_x1_mean=(double)(*(x1+T)-*x1)/T;

    *j_y=(double)move_y_sum*1.0/(T*L);
    *v_y0_mean=(double)*(y0+T)/T;
    *v_y1_mean=(double)(*(y1+T)-*y1)/T;
  return;
}


void Dynamics_Coupled(int L, int N, int M, int *X, int *Y, int *move_x, int *move_y)
{
  int n, i, j, k, occupied, x_next, y_next;
     *move_x=0;
     *move_y=0;
      for(i=0;i<N+M;i++)
      {
        n=(int)(Uniform()*(N+M));
        occupied=0;
        /*粒子Yの処理*/
        if(n==0||n==1)
        {
          y_next=(*(Y+n)+1)%L;
          for(j=0;j<N;j++){
            if(*(X+j)==y_next){
              occupied=1;
            }
          }
          if(occupied==0){
            *(Y+n)=y_next;
            *move_y+=1;
          }
        }
        /*粒子Xの処理*/
        else
        {
          n-=2;
          x_next=(*(X+n)+1)%L;
          for(j=0;j<N;j++){
            if(*(X+j)==x_next){
              occupied=1;
            }
          }
          if(occupied==0){
            *(X+n)=x_next;
            *move_x+=1;
          }
        }
      }

  return;
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
