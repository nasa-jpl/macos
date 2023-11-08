//
// splicubi2.c
// For 2D grid cubic spline interpolations, adapted from Numerical Recipes
// C routines.
// John Z. Lou, Jet Propulsion Laboratory
// 10/2010
//

// *** Warning: maximum number of grid surfaces currently supported: 256;
//     need to change the constant below to increase the limit.           

#include <stdio.h>
#include <stdlib.h>

#define SREAL double

void splie2(SREAL[],SREAL[],SREAL**,int,int,SREAL**);
void splin2(SREAL[],SREAL[],SREAL**,SREAL**,
            int,int,SREAL,SREAL,int,int,SREAL*,SREAL*,SREAL*);
void spline(SREAL[],SREAL[],int,SREAL,SREAL,SREAL[]);
void splint(SREAL[],SREAL[],SREAL[],int,SREAL,SREAL*);
SREAL *vector(int,int);
double **dmatrix(int,int,int,int);
void free_vector(SREAL*,int,int);
void errorf(char*);

void splie2f2c_(SREAL[],SREAL[],SREAL*,int*,int*,int*);
void splin2f2c_(SREAL[],SREAL[],int*,int*,
                SREAL*,SREAL*, SREAL*,SREAL*,SREAL*,int*,int*,int*);

#if 0
int main()
{
  int m,n,i,j;
  SREAL *x1arr,*x2arr,*y1m,fn,dfdx1,dfdx2;
  SREAL x1,x2,dx,x,y;
  int i0,j0,srf_id,iray;

  m=n=99; dx=1.0/98.0;
  x1arr=vector(1,m); x2arr=vector(1,n);
  for (i=1;i<=m; i++) x1arr[i]=1.0+(i-1)*dx;
  for (i=1;i<=n; i++) x2arr[i]=1.0+(i-1)*dx;
  y1m=vector(1,m*n);
  
  // Define a conic y1m: z = x^2+y^2
  for (i=1;i<=m;i++) {
    y=1.0+(i-1)*dx;    
    for (j=1;j<=n;j++) {
      x=1.0+(j-1)*dx;
      y1m[(i-1)*n+j]=x*x+y*y;
    }
  }

  srf_id=1;
  i0=48, x1=1.0+(i0-1+0.3)*dx;
  j0=21, x2=1.0+(j0-1+0.5)*dx;
  //
  splie2f2c_(x1arr,x2arr,y1m,&m,&n,&srf_id);
  splin2f2c_(x1arr,x2arr,&m,&n,&x1,&x2,&fn,&dfdx1,&dfdx2,&i0,&j0,&srf_id);
  //
  printf("**** x1,x2 = %lf, %lf\n",x1,x2);
  //printf("**** CS computed fn = %lf, exact fn = %lf\n",fn,y1m[(i0-1)*n+j0]);
  printf("**** CS computed fn = %lf, exact fn = %lf\n",fn,x1*x1+x2*x2);
  printf("**** CS computed dfdx1 = %lf, exact dfdx1 = %lf\n",dfdx1,2*x1);
  printf("**** CS computed dfdx2 = %lf, exact dfdx2 = %lf\n",dfdx2,2*x2);
  return 0;
} // main - testing only
#endif


double*** ya_arr;  // global array of matrix, storing grid y values for cubic splines
double*** y2a_arr; // global array of matrix, storing second derivatives of cubic splines


// *** This function is called from Fortran code in MACOS ***
// This routine is called ONLY ONCE for each grid surface during MACOS grid surface init.
// ya_f, of 1-D array, is allocated in Fortran code. It contains 2D grid surface data.
void splie2f2c_(SREAL x1a[],SREAL x2a[],SREAL* ya_f,int* m,int* n,int* srf_id)
{
  SREAL **ya,**y2a;
  int ir,ic;
  const short nGridSrf=256;
  static char first_entry=1;

  if (first_entry) {
    ya_arr=(double***) calloc(nGridSrf,sizeof(double**));
    y2a_arr=(double***) calloc(nGridSrf,sizeof(double**));
    first_entry=0;
  }

  ya=dmatrix(1,*m,1,*n); y2a=dmatrix(1,*m,1,*n);

  for (ir=1; ir<=*m; ++ir) {
    for (ic=1; ic<=*n; ++ic) {
      ya[ir][ic]=ya_f[(ir-1)*(*n)+ic];
    }
  }

  splie2(x1a,x2a,ya,*m,*n,y2a);
  ya_arr[*srf_id]=ya;   // save ya for this grid surface
  y2a_arr[*srf_id]=y2a; // save y2a for this grid surface
} // splie2_f2c, wrapper routine providing interface from Fortran to C


void splie2(SREAL x1a[],SREAL x2a[],SREAL** ya,int m,int n,SREAL** y2a)
{
  int j;

  for (j=1;j<=m;j++) {
    spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
  }
} // splie2


// *** This function is called from Fortran code in MACOS ***
// Given p={x1,x2}, it computes cubic spline fn at p and partial derivatives
// dfdx1, dfdx2 at p. 
void splin2f2c_(SREAL x1a[],SREAL x2a[],int* m,int* n,
                SREAL* x1,SREAL* x2, SREAL* fn,SREAL* dfdx1,SREAL* dfdx2,
                int* i0,int* j0,int* srf_id)
{
  splin2(x1a,x2a,ya_arr[*srf_id],y2a_arr[*srf_id],*m,*n,*x1,*x2,*i0,*j0,
         fn,dfdx1,dfdx2);
} // splin2f2c_


//
// Remember the convention used in these routines:
// i0,x1 are in vertical/'y' direction; j0,x2 are in horizontal/'x' direction
//
void splin2(SREAL x1a[],SREAL x2a[],SREAL** ya,SREAL** y2a,
            int m,int n,SREAL x1,SREAL x2,int i0,int j0,
            SREAL *y,SREAL* dfdx1,SREAL* dfdx2)
{
  int i,j;
  static SREAL *ytmp1,*yytmp1,*ytmp2,*yytmp2,*ya_tmp,*y2a_tmp;
  double A,B;
  static char first_entry=1;

  if (first_entry) {
    ytmp1=vector(1,n); yytmp1=vector(1,n);
    ytmp2=vector(1,m); yytmp2=vector(1,m);
    ya_tmp=vector(1,m); y2a_tmp=vector(1,m);
    first_entry=0;
  }

  if (0 & !first_entry) {
    printf("***** splin2: ya = %x, y2a = %x\n",ya,y2a);
  }

  for (j=1;j<=m;j++) {
    splint(x2a,ya[j],y2a[j],n,x2,&yytmp1[j]); // yytmp1 is row spline at x2
  }

  spline(x1a,yytmp1,m,1.0e30,1.0e30,ytmp1); // ytmp1 is 2nd derivative
  splint(x1a,yytmp1,ytmp1,m,x1,y);


  // Now comput dfdx1 at (x1,x2)
  A=(x1a[i0+1]-x1)/(x1a[i0+1]-x1a[i0]); B=1.0-A;
  *dfdx1=(yytmp1[i0+1]-yytmp1[i0])/(x1a[i0+1]-x1a[i0])
          -(3.0*A*A-1.0)/6.0*(x1a[i0+1]-x1a[i0])*ytmp1[i0]
          +(3.0*B*B-1.0)/6.0*(x1a[i0+1]-x1a[i0])*ytmp1[i0+1];

  // To compute dfdx2, need first compute interpolated values and 2nd derivatives along  
  // 'n' direction with fixed 'x1'.
  for (j=1;j<=n;j++) {
    for (i=1;i<=m;i++) { ya_tmp[i]=ya[i][j]; y2a_tmp[i]=y2a[i][j]; }
    splint(x1a,ya_tmp,y2a_tmp,m,x1,&yytmp2[j]); // yytmp2 is column spline at x1
  }
  spline(x2a,yytmp2,n,1.0e30,1.0e30,ytmp2); // ytmp2 is 2nd derivative
  
  // With yytmp2 and ytmp2 computed, now compute dfdx2
  A=(x2a[j0+1]-x2)/(x2a[j0+1]-x2a[j0]); B=1.0-A;
  *dfdx2=(yytmp2[j0+1]-yytmp2[j0])/(x2a[j0+1]-x2a[j0])
          -(3.0*A*A-1.0)/6.0*(x2a[j0+1]-x2a[j0])*ytmp2[j0]
          +(3.0*B*B-1.0)/6.0*(x2a[j0+1]-x2a[j0])*ytmp2[j0+1];
} // splin2


void spline(SREAL x[],SREAL y[],int n,SREAL yp1,SREAL ypn,SREAL y2[])
{
    int i,k;
    SREAL p,qn,sig,un,*u;

    u=vector(1,n-1);
    if (yp1 > 0.99e30)
      y2[1]=u[1]=0.0;
    else {
      y2[1] = -0.5;
      u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
      qn=un=0.0;
    else {
      qn=0.5;
      un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)
      y2[k]=y2[k]*y2[k+1]+u[k];
    free_vector(u,1,n-1);
} // spline


void splint(SREAL xa[],SREAL ya[],SREAL y2a[],int n,SREAL x,SREAL* y)
{
    int klo,khi,k;
    SREAL h,b,a;

    klo=1;
    khi=n;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k;
      else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) errorf("Bad XA input to routine SPLINT");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
} // splint


SREAL *vector(int nl,int nh)
{
    SREAL *v;

    v=(SREAL *)malloc((unsigned) (nh-nl+1)*sizeof(SREAL));
    if (!v) errorf("allocation failure in vector()");
    return v-nl;
} // vector


double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
    int i;
    double **m;

    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) errorf("allocation failure 1 in dmatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
      if (!m[i]) errorf("allocation failure 2 in dmatrix()");
      m[i] -= ncl;
    }
    return m;
} // dmatrix


void free_vector(SREAL *v,int nl,int nh)
{
   free((char*) (v+nl));
} // free_vector


void errorf(char *msg)
{
   fprintf(stdout,"%s",msg); exit(0);
}  // errorf

