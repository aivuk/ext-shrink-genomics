#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define MAXTR 100
#define EPS 10e-28

void read_eigen (double **v, double *d, int n)
{
   int i, j, k;
   FILE *InPut;

   InPut = fopen("./EigenOut.txt", "r");

   for ( i = 0; i < n; i++){
      fscanf(InPut, "%lf", &d[i]);
   }
   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++){
         fscanf(InPut, "%lf", &v[i][j]);
      }
}

void extend (double **mat, double **v, double *d, int n)
{
   int i, j, k;
   double **aux;

   aux = (double **) malloc (n*sizeof(double *));
   for (j = 0; j < n; j++)
      aux[j] = (double *) malloc (n*sizeof(double));
   
   for ( i = 0; i < n; i++)
      if(d[i] < 0.000003) d[i] = d[i-1];

   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         mat[i][j] = 0.0;
      
   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         aux[i][j] = d[i]*v[j][i]; 

   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         for (k = 0; k < n; k++)
            mat[i][j] += v[i][k]*aux[k][j];
//   free(aux);
}

void print_mat(double **Mat, char *outfile, int Ne, int tr)
{
   int i, j;
   FILE *OutPut;
   OutPut = fopen (outfile, "w");

   for (i = 0; i < Ne; i++){
      for ( j = 0; j < tr-1; j++){
         fprintf(OutPut, "%lf\t", Mat[i][j]);
      }
      fprintf(OutPut, "%lf\n", Mat[i][j]);
   }
}

void bend (double **mat_ext, int tr)
{
   int j;
   double **v, *d;
   
   v = (double **) malloc (tr*sizeof(double *));
   for (j = 0; j < tr; j++)
      v[j] = (double *) malloc (tr*sizeof(double));

   d = (double *) malloc (tr*sizeof(double));

   read_eigen (v, d, tr);
   extend (mat_ext, v, d, tr);
}

main (){
   int i, j, k, l, Ne, tr;
   double **Covar;

   scanf("%d", &tr);
   
   Covar = (double **) malloc (tr*sizeof(double));
   for ( j = 0; j < tr; j++)
      Covar[j] = (double *) malloc (tr*sizeof(double));

   bend(Covar, tr);
   print_mat(Covar, "Cov_Ext.T.csv", tr, tr);
}
