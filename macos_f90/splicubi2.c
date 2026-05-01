#include <stdio.h>
#include <stdlib.h>

#define SREAL double

/* --- Fortran to C Linkage Mapping (Windows and Linux with Intel ifx -names uppercase) --- */
#if defined(_WIN32) || defined(_WIN64)
#  define splie2f2c_ SPLIE2F2C_
#  define splin2f2c_ SPLIN2F2C_
#endif

/* --- Function Prototypes --- */
void splie2(SREAL[],SREAL[],SREAL**,int,int,SREAL**);
void splin2(SREAL[],SREAL[],SREAL**,SREAL**,
            int,int,SREAL,SREAL,int,int,SREAL*,SREAL*,SREAL*);
void spline(SREAL[],SREAL[],int,SREAL,SREAL,SREAL[]);
void splint(SREAL[],SREAL[],SREAL[],int,SREAL,SREAL*);
SREAL *vector(int,int);
double **dmatrix(int,int,int,int);
void free_vector(SREAL*,int,int);
void errorf(char*);

/* --- Global arrays for grid surfaces --- */
double*** ya_arr = NULL;  
double*** y2a_arr = NULL; 

/* --- 1. Fortran Interface: Initialization --- */
void splie2f2c_(SREAL x1a[], SREAL x2a[], SREAL* ya_f, int* m, int* n, int* srf_id)
{
  SREAL **ya, **y2a;
  int ir, ic;
  const short nGridSrf = 256;
  static char first_entry = 1;

  if (first_entry) {
    ya_arr = (double***) calloc(nGridSrf, sizeof(double**));
    y2a_arr = (double***) calloc(nGridSrf, sizeof(double**));
    first_entry = 0;
  }

  ya = dmatrix(1, *m, 1, *n); 
  y2a = dmatrix(1, *m, 1, *n);

  for (ir = 1; ir <= *m; ++ir) {
    for (ic = 1; ic <= *n; ++ic) {
      ya[ir][ic] = ya_f[(ir - 1) * (*n) + ic];
    }
  }

  splie2(x1a, x2a, ya, *m, *n, y2a);
  ya_arr[*srf_id] = ya;   
  y2a_arr[*srf_id] = y2a; 
}

/* --- 2. Fortran Interface: Interpolation --- */
void splin2f2c_(SREAL x1a[], SREAL x2a[], int* m, int* n,
                SREAL* x1, SREAL* x2, SREAL* fn, SREAL* dfdx1, SREAL* dfdx2,
                int* i0, int* j0, int* srf_id)
{
  splin2(x1a, x2a, ya_arr[*srf_id], y2a_arr[*srf_id], *m, *n, *x1, *x2, *i0, *j0,
         fn, dfdx1, dfdx2);
}

/* --- Internal Spline Logic --- */

void splie2(SREAL x1a[], SREAL x2a[], SREAL** ya, int m, int n, SREAL** y2a)
{
  int j;
  for (j = 1; j <= m; j++) {
    spline(x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j]);
  }
}

void splin2(SREAL x1a[], SREAL x2a[], SREAL** ya, SREAL** y2a,
            int m, int n, SREAL x1, SREAL x2, int i0, int j0,
            SREAL *y, SREAL* dfdx1, SREAL* dfdx2)
{
  int i, j;
  static SREAL *ytmp1, *yytmp1, *ytmp2, *yytmp2, *ya_tmp, *y2a_tmp;
  double A, B;
  static char first_entry = 1;

  if (first_entry) {
    ytmp1 = vector(1, n);  yytmp1 = vector(1, n);
    ytmp2 = vector(1, m);  yytmp2 = vector(1, m);
    ya_tmp = vector(1, m); y2a_tmp = vector(1, m);
    first_entry = 0;
  }

  for (j = 1; j <= m; j++) {
    splint(x2a, ya[j], y2a[j], n, x2, &yytmp1[j]);
  }

  spline(x1a, yytmp1, m, 1.0e30, 1.0e30, ytmp1); 
  splint(x1a, yytmp1, ytmp1, m, x1, y);

  A = (x1a[i0 + 1] - x1) / (x1a[i0 + 1] - x1a[i0]); 
  B = 1.0 - A;
  *dfdx1 = (yytmp1[i0 + 1] - yytmp1[i0]) / (x1a[i0 + 1] - x1a[i0])
           - (3.0 * A * A - 1.0) / 6.0 * (x1a[i0 + 1] - x1a[i0]) * ytmp1[i0]
           + (3.0 * B * B - 1.0) / 6.0 * (x1a[i0 + 1] - x1a[i0]) * ytmp1[i0 + 1];

  for (j = 1; j <= n; j++) {
    for (i = 1; i <= m; i++) { 
        ya_tmp[i] = ya[i][j]; 
        y2a_tmp[i] = y2a[i][j]; 
    }
    splint(x1a, ya_tmp, y2a_tmp, m, x1, &yytmp2[j]);
  }
  spline(x2a, yytmp2, n, 1.0e30, 1.0e30, ytmp2);

  A = (x2a[j0 + 1] - x2) / (x2a[j0 + 1] - x2a[j0]); 
  B = 1.0 - A;
  *dfdx2 = (yytmp2[j0 + 1] - yytmp2[j0]) / (x2a[j0 + 1] - x2a[j0])
           - (3.0 * A * A - 1.0) / 6.0 * (x2a[j0 + 1] - x2a[j0]) * ytmp2[j0]
           + (3.0 * B * B - 1.0) / 6.0 * (x2a[j0 + 1] - x2a[j0]) * ytmp2[j0 + 1];
}

void spline(SREAL x[], SREAL y[], int n, SREAL yp1, SREAL ypn, SREAL y2[])
{
    int i, k;
    SREAL p, qn, sig, un, *u;
    u = vector(1, n - 1);
    if (yp1 > 0.99e30) y2[1] = u[1] = 0.0;
    else {
      y2[1] = -0.5;
      u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
    }
    for (i = 2; i <= n - 1; i++) {
      sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
      u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30) qn = un = 0.0;
    else {
      qn = 0.5;
      un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
    }
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);
    for (k = n - 1; k >= 1; k--) y2[k] = y2[k] * y2[k + 1] + u[k];
    free_vector(u, 1, n - 1);
}

void splint(SREAL xa[], SREAL ya[], SREAL y2a[], int n, SREAL x, SREAL* y)
{
    int klo, khi, k;
    SREAL h, b, a;
    klo = 1; khi = n;
    while (khi - klo > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) errorf("Bad XA input to routine SPLINT");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}

SREAL *vector(int nl, int nh)
{
    SREAL *v;
    v = (SREAL *)malloc((unsigned) (nh - nl + 1) * sizeof(SREAL));
    if (!v) errorf("allocation failure in vector()");
    return v - nl;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;
    m = (double **) malloc((unsigned) (nrh - nrl + 1) * sizeof(double*));
    if (!m) errorf("allocation failure 1 in dmatrix()");
    m -= nrl;
    for(i = nrl; i <= nrh; i++) {
      m[i] = (double *) malloc((unsigned) (nch - ncl + 1) * sizeof(double));
      if (!m[i]) errorf("allocation failure 2 in dmatrix()");
      m[i] -= ncl;
    }
    return m;
}

void free_vector(SREAL *v, int nl, int nh)
{
   free((char*) (v + nl));
}

void errorf(char *msg)
{
   fprintf(stdout, "%s", msg); exit(0);
}