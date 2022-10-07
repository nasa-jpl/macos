//
// fk.c
//
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>  // for memcpy()

#define parentProc  1
#define childProc   0
#define max_np  31  // up to 31 child procs  
int fd[max_np][2],nproc,nimg;

void fkprocs_(int *np, int *myip, int *mypart, int *nimg_in)
{
  pid_t  pid;
  int ip,i,ii,part;

  if (*np <= 1 || *np > max_np) {
    nimg = (*nimg_in)/2;
    *np = 1, *myip = 1, *mypart = parentProc;
    nproc = *np;
    return;
  }

  nimg = (*nimg_in)/2;

  if (*np>nimg) *np = nimg;  // load balancing
  nproc = *np;
  //printf(" --fkprocs_2: *np = %i\n", *np);

  pid = 1;
  part = 1; // default to parent
  for (ip=1; ip<=(*np)-1; ++ip) {
    if (pid != 0) {
      pipe(fd[ip]);
      pid = fork();
    }
    if (pid == 0) {
       close(fd[ip][0]);  // close read fd
       part=childProc;
       break;
    }
    else {
      close(fd[ip][1]);  // close write fd
      part=parentProc;
    }
  }

  *myip = ip;  // parant has ip = np
  *mypart = part;
}  // fkprocs_


void pipeio_(int *part, int *ip, double *pix, int *ppsize)
{
  static char first_entry=1;
  int i, p, pi, pj, psize, psize2;
  static double *data;
  FILE *fp;

  if (nproc==1) {
    // Summary
    printf("\n ** Processes used: %i\n",nproc);
    printf(" ** Images generated and composed: %i\n",nimg);
    return;
  }

  psize = *ppsize, psize2=psize*psize;

  if (*part==parentProc) {
    first_entry = 0;
    data = (double*) calloc(psize2,sizeof(double));

    //printf(" **pipeio: party = parent\n");
    for (i=1; i<=nproc-1; ++i) {
      read(fd[i][0],data,psize2*sizeof(double));

      // Combine images
      for (p=0; p<psize2; ++p) {
        pix[p] += data[p];
        //pix[p] = (double)(i%10);  % test only
      }
    }
    free(data);
#if 0
    // Dump final combined image
    fp = fopen("DFSimg.txt","w");
    for (pi=0; pi<psize; ++pi) {
      for (pj=0; pj<psize; ++pj)
        fprintf(fp,"%e ",pix[(pi-1)*psize+pj]);
      fprintf(fp,"\n");
    }
    fclose(fp);
#endif
    // Summary
    printf("\n ** Processes used: %i\n",nproc);
    printf(" ** Images generated and composed: %i\n",nimg);
   }
  else {
    // Child
    write(fd[*ip][1],pix,psize2*sizeof(double));
    exit(0);  // no reason to be alive any longer
  } 
} // pipeio_


void recvintarr_(int *ip, int *ibuf, int *len)
{
   read(fd[*ip][0], ibuf, (*len)*sizeof(int)); 
}  // recvintarr_

void sendintarr_(int *ip, int *ibuf, int *len)
{
   write(fd[*ip][1], ibuf, (*len)*sizeof(int));      
}  // sendintarr_


void recvdwdxloc_(int *ip, double *dbuf, int *buflen)
{ 
  read(fd[*ip][0],dbuf,(*buflen)*sizeof(double));
} // recvdwdxloc_


void senddwdxloc_(int *ip, double *dbuf, int *buflen)
{
  write(fd[*ip][1],dbuf,(*buflen)*sizeof(double));
} // senddwdxloc_


void cdblcopy_(double *dest, double *src, int *len)
{
  memcpy((void*)dest, (void*)src, (*len)*sizeof(double));
}  // cdblcopy_


void writedwdx_(double *buff, int *len)
{
  FILE *fp;
  fp =fopen("dwdx.bin","wb");
  printf("*** writedwdx: len = %i\n", *len);
  fwrite((void*) buff, sizeof(double), *len, fp);
  fclose(fp);
}


void writeidx_(int *buf, int *len, int *nrow)
{
  FILE *fp;
  fp =fopen("colidx.bin","wb");
  fwrite((void*) buf, sizeof(int), *len, fp);
  fwrite((void*) nrow, sizeof(int), 1, fp);
  fclose(fp); 
}


void quitchild_()
{
  exit(0);
}
