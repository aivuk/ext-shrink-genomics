#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define MAXTR 100
#define EPS 10e-28

double normal(double mu, double sig)
{
   double x1, x2, w, y1;
   
   do {
      x1 = 2.0 * ((double)(rand()) / (RAND_MAX + 1.)) - 1.0;
      x2 = 2.0 * ((double)(rand()) / (RAND_MAX + 1.)) - 1.0;
      w = x1 * x1 + x2 * x2;
   } while ( w >= 1.0 );

   w = sqrt( (-2.0 * log( w ) ) / w );
   y1 = x1 * w;

   return( mu + y1 * sig );
}

void print_eigen (double **v, double *d, int n)
{
   int i, j, k;
   FILE *OutPut;

   OutPut = fopen("./EigenOut.txt", "w");

   for ( i = 0; i < n; i++){
      if(i<n-1)fprintf(OutPut, "%lf ", d[i]);
      else fprintf(OutPut, "%lf\n\n", d[i]);
   }
      
   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++){
         if(j<n-1) fprintf(OutPut, "%lf ", v[i][j]);
         else fprintf(OutPut, "%lf\n", v[i][j]);
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
      if(d[i] < 0.00001) d[i] = d[i-1];

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
}

void eigsrt (double *d, double ***v, int n)
{
   int i, k, j;
   double p;
   for (i = 0; i<n-1; i++){
      p = d[k=i];
      for (j = i; j<n; j++)
         if (d[j] >= p) p=d[k=j];
      if( k!=i ){
         d[k] = d[i];
         d[i] = p;
         for ( j = 0; j<n; j++){
            p = (*v)[j][i];
            (*v)[j][i] =(*v)[j][k];
            (*v)[j][k] = p;
         }
      }
   }
}

void rot (double **a, double s, double tau, int i, int j, int k, int l)
{
   double g = a[i][j];
   double h = a[k][l];
   a[i][j] = g - s*(h+g*tau);
   a[k][l] = h + s*(g-h*tau);
}

void jacobi (double **a, double **v, double *d, int n)
{
   int nrot = 0, i, j, ip, iq;
   double tresh, theta, tau, t, sm, s, h, g, c;
   double *b, *z;

   z = (double *) malloc (n*sizeof(double));
   b = (double *) malloc (n*sizeof(double));

   for (ip = 0; ip < n; ip++){
      for (iq = 0; iq < n; iq++)
         v[ip][iq] = 0.;
   }
   for (ip = 0; ip < n; ip++){
      v[ip][ip] = 1.;
      d[ip] = a[ip][ip];
      b[ip] = d[ip];
      z[ip] = 0.;
   }
   for ( i = 1; i <= 50; i++) {
      sm = 0.0;
      for (ip = 0; ip < n-1; ip++){
         for (iq = ip+1; iq<n; iq++)
            sm += fabs(a[ip][iq]);
      }
      if (sm < EPS){
         eigsrt (d, &v, n);
         return;
      }
      if (i < 4)
         tresh = 0.2*sm/(n*n);
      else
         tresh = 0.0;
      for (ip = 0; ip < n-1; ip++){
         for (iq = ip+1; iq<n; iq++){
            g = 100.0*fabs(a[ip][iq]);
            if ( i > 4 && g <= EPS*fabs(d[ip]) && g <= EPS*fabs(d[iq]))
               a[ip][iq] = 0.0;
            else if (fabs(a[ip][iq]) > tresh){
               h = d[iq]-d[ip];
               if (g <= EPS*fabs(h))
                  t = (a[ip][iq])/h;
               else {
                  theta = 0.5*h/(a[ip][iq]);
                  t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta));
                  if (theta < 0.0) t = -t;
               }
               c = 1.0/sqrt(1+t*t);
               s = t*c;
               tau = s/(1.0+c);
               h = t*a[ip][iq];
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               a[ip][iq] = 0.0;
               for ( j = 0; j < ip; j++)
                  rot(a, s, tau, j, ip, j, iq);
               for ( j = ip+1; j < iq; j++)
                  rot(a, s, tau, ip, j, j, iq);
               for ( j = iq+1; j < n; j++)
                  rot(a, s, tau, ip, j, iq, j);
               for ( j = 0; j < n; j++)
                  rot(v, s, tau, j, ip, j, iq);
               ++nrot;
            }
         }
      }
      for (ip=0; ip<n; ip++){
         b[ip] += z[ip];
         d[ip] = b[ip];
         z[ip] = 0.0;
      }
   }
   printf("\nNÃO COPNVERIGUEWD QARARRO\n");
}

void bend (double **mat, double **mat_ext, int tr)
{
   int j;
   double **v, *d;
   
   v = (double **) malloc (tr*sizeof(double *));
   for (j = 0; j < tr; j++)
      v[j] = (double *) malloc (tr*sizeof(double));

   d = (double *) malloc (tr*sizeof(double));

   jacobi(mat, v, d, tr);
   print_eigen (v, d, tr);
}

void cholesky (double **orig, double **chol, int tr, int dept)
{
   int i, j, k;

   FILE *error;
start:
   for (i = 0; i < tr; i++) {
      chol[i][i] = orig[i][i];
      for (k = 0; k < i; k++)
         chol[i][i] -= chol[k][i]*chol[k][i];
      if (chol[i][i] <= 0) {
         bend(orig, orig, tr);
         error = fopen("./std.error.txt", "a");
         fprintf(error, "\nERROR: non-positive definite matrix!\n");
         dept++;
         fprintf(error, "\nproblem from %d %f %d\n", i, chol[i][i], dept);
         fclose(error);
         goto start;
      }
      chol[i][i] = sqrt(chol[i][i]);

      for (j = i+1; j < tr; j++) {
         chol[i][j] = orig[i][j];
         for (k = 0; k < i; k++)
            chol[i][j] -= chol[k][i]*chol[k][j];
         chol[i][j] /= chol[i][i];
         chol[j][i] = chol[i][j];
      }
   }
}

void mat_cov (double **data, int Ne, int tr, double **Covar, double **Corr, double **Star)
{
   int i, j, k, l;
   double *x, *y, *wxy, Var_s=0.0, sum_var = 0, sum_Cov2 = 0 ;
   double correlation, covariance, lambda;
   double yt, xt;
   double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0, epx = 0.0, epy = 0.0;

   x = (double *) malloc (Ne*sizeof(double));
   y = (double *) malloc (Ne*sizeof(double));
   wxy = (double *) malloc (tr*sizeof(double));

   for (i = 0; i < tr; i++){
      for (j = 0; j <= i; j++){
         sxy = 0.0;
         sxx = 0.0;
         syy = 0.0; 
         epx = 0.0;
         epy = 0.0;
         ax = 0.0;
         ay = 0.0; 
         Var_s = 0.0;
         for(k = 0; k < Ne; k++){
            x[k] = data[k][i];
            y[k] = data[k][j];
            ax += x[k];
            ay += y[k];
         }
         ax /= Ne;
         ay /= Ne;
         for(k = 0; k < Ne; k++){
            xt = x[k] - ax;
            yt = y[k] - ay;
            sxx += xt*xt;
            syy += yt*yt;
            wxy[k] = xt*yt;
            sxy += wxy[k];
            epx += xt;
            epy += yt;
         }
         covariance = (sxy - epx*epy/Ne)/(Ne-1);
         sxy /= Ne;
         for(k = 0; k < Ne; k++)
            Var_s += (wxy[k] - sxy)*(wxy[k] - sxy);

         sxx = (sxx - epx*epx/Ne)/(Ne-1);
         syy = (syy - epy*epy/Ne)/(Ne-1);
         correlation = covariance/sqrt(sxx*syy); 
         Corr[i][j] = correlation;
         Corr[j][i] = correlation;
         Covar[i][j] = covariance;   
         Covar[j][i] = covariance;
         if(i!=j){
            sum_var += Var_s*Ne/((Ne-1)*(Ne-1)*(Ne-1));
            sum_Cov2 += covariance*covariance;
         }
      }
   }
   lambda = sum_var/sum_Cov2;
   if (lambda > 1) lambda = 1;
   if (lambda < 0) lambda = 0;
   for (i = 0; i < tr; i++){
      for (j = 0; j < tr; j++){
         if(i!=j) Star[i][j] = Covar[i][j] * (1-lambda);
         else Star[i][j] = Covar[i][j];
      }
   }
}

double solve(double *left, double **chol, double *right, int tr)
{
   int i, k;
   double result = 0.;

   for (i = 0; i < tr; i++) {
      for (k = 0; k < i; k++) {
         right[i] -= right[k]*chol[k][i];
      }
      right[i] /= chol[i][i];
   }
   
   for (i = tr-1; i >= 0; i--){
      for (k = i+1; k < tr; k++) {
         right[i] -= right[k]*chol[i][k];
      }
      right[i] /= chol[i][i];
   }
   
   for (i = 0; i < tr; i++){ 
      result += right[i]*left[i];
   }

   return result/tr;
}

void print_mat(double **Mat, char *outfile, int Ne, int tr)
{
   int i, j;
   FILE *OutPut;
   OutPut = fopen (outfile, "w");

   for (i = 0; i < Ne; i++){
      for ( j = 0; j < tr-1; j++){ fprintf(OutPut, "%lf\t", Mat[i][j]);
      }
      fprintf(OutPut, "%lf\n", Mat[i][j]);
   }
}

main (){
   int i, j, k, l, Ne, tr;
   double **data, **Covar, **Corr, **Star, **Cov_Ext, **Corr_Ext;

   FILE *Dados;

   Dados = fopen("dados.txt", "r");

   fscanf(Dados, "%d", &Ne); // numero de observaçoes
   fscanf(Dados, "%d", &tr); // numero de dimensoes dos dados

   data = (double **) malloc (Ne*sizeof(double));
   for ( j = 0; j < Ne; j++)
      data[j] = (double *) malloc (tr*sizeof(double));
   
   Covar = (double **) malloc (tr*sizeof(double));
   for ( j = 0; j < tr; j++)
      Covar[j] = (double *) malloc (tr*sizeof(double));
   
   Corr = (double **) malloc (tr*sizeof(double));
   for ( j = 0; j < tr; j++)
      Corr[j] = (double *) malloc (tr*sizeof(double));
   
   Star = (double **) malloc (tr*sizeof(double));
   for ( j = 0; j < tr; j++)
      Star[j] = (double *) malloc (tr*sizeof(double));
   
   
   for ( k = 0; k < Ne; k++)
      for (j = 0; j < tr; j++)
         fscanf(Dados, "%lf", &data[k][j]);
   
   mat_cov (data, Ne, tr, Covar, Corr, Star);

   print_mat(Covar, "Covar.csv", tr, tr);
   print_mat(Corr, "Corr.csv", tr, tr);
   print_mat(Star, "Star.csv", tr, tr);

}
