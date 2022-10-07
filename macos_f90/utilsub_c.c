/***********************************************************************
        Begin file utilsub_c.c 
     +----------------------------------------------------------------+
     |  Copyright (C) 1995-2007, California Institute of Technology.  |
     |  U.S. Government Sponsorship Is Acknowledged.                  |
     +----------------------------------------------------------------+
***********************************************************************/

//#ifdef MACOS_CMD

//
// C utility routines for MACOS
//

#include <stdio.h>
#include <unistd.h> 
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#ifdef SUNOS
/* sws for solaris */
#include <dirent.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#define RxFile  0
#define MaxFileCnt  1000
#define MaxFnLen    132

//getcwd(dirname,256); 
void files_sort(char[][MaxFnLen],int);
int match_ext(char*, char*);

void shellcmd_(char *cmd, int *clen)
{
  cmd[*clen]='\0';
  system(&cmd[1]);  // exclude '!' character
}

void
pwd_()
{
  char dirname[256];
  getcwd(dirname,256); // get current working directory 
  printf("  Current Folder: \n  %s\n\n",dirname);
}

void 
chgdir_(char *dirname, int* clen)
{
  char curr_dirname[256], cmd[256];
  dirname[*clen]='\0'; 
  sprintf(cmd,"cd %s",dirname);
  //system(cmd);
  chdir(dirname);
  getcwd(curr_dirname,256); 
  printf("  Current Folder: \n  %s\n\n",curr_dirname);
}


void
starteditor_(char *cmd,int* cmd_len,char* filename,int* clen)
{ 
  char edit_cmd[128];
  cmd[*cmd_len]='\0', filename[*clen]='\0';
  sprintf(edit_cmd,"%s %s",cmd,filename);
  //printf("  filename = %s\n",filename);
  system(edit_cmd);
}


void
startmatlab_(char *cmd,int* cmd_len,char* arg,int* clen)
{
  char ml_cmd[128];
  cmd[*cmd_len]='\0', arg[*clen]='\0';
  sprintf(ml_cmd,"%s %s",cmd,arg);
  system(ml_cmd);
}


#include <sys/types.h>
#include <dirent.h>
#include <string.h>

char files[MaxFileCnt][MaxFnLen];
int file_cnt=0;

void
lsfiles_(int* fileType, char* cmdArg, int* clen)
{
  char dirname[256], cbuf[3][256], cmd[130];
  DIR *dirp;
  struct dirent *dp;
  int slen,max_slen[3],tot_files,mm;

  getcwd(dirname,256);
  if (*fileType==0) 
    printf("  MACOS Prescription File(s) with File Id in Current Folder\n  %s:\n",dirname);
  else
    printf("  Files in Current Folder: \n  %s\n",dirname);
  
  if (*clen>0) cmdArg[*clen]='\0';  

  if (*fileType==2) {
    if (*clen>0) sprintf(cmd,"ls -l %s",cmdArg);
    else sprintf(cmd,"ls -l");
    system(cmd);
    return;
  }

  dirp = opendir("."); file_cnt=0;
  while (dirp) {
    if ((dp = readdir(dirp)) != NULL) {
      if (*clen==0 && dp->d_name[0]=='.') continue;

      if (*clen>0) {
        if (strcmp(cmdArg,"*")==0 || strcmp(cmdArg,"*.*")==0) ;
        else if (strcmp(cmdArg,".*") ==0) {
	  if (dp->d_name[0] != '.') continue; 
         }
        else if (!match_ext(dp->d_name,cmdArg)) continue;
      }

      slen=strlen(dp->d_name);
      if (*fileType==RxFile) {
        if (slen>=3 &&
            (strncmp((char*) &(dp->d_name[slen-3]),".in",3)==0)) {
	  strcpy(files[file_cnt],dp->d_name);
	  file_cnt++;
        }
       }
      else {
	strcpy(files[file_cnt],dp->d_name);
#ifdef SUNOS
        /* solaris doesn't have DT_DIR */
        {
          struct stat *s;
          stat(dp->d_name, s);
          if (s->st_mode & S_IFDIR) { // directory
            files[file_cnt][strlen(dp->d_name)]='/';
            files[file_cnt][strlen(dp->d_name)+1]='\0';
          }
        }
#else
        if (dp->d_type==DT_DIR) { // directory
          files[file_cnt][strlen(dp->d_name)]='/';
          files[file_cnt][strlen(dp->d_name)+1]='\0';
        }
#endif
        file_cnt++;
      }
     }
    else
      break;
  } 

  closedir(dirp);
  if (!file_cnt) {
    if (*clen==0) 
      printf("  No prescription file found in the folder!\n");
    else 
      printf("  No file found in the folder!\n");
    return;
  }
  tot_files=file_cnt;
 
  // Sort filenames 
  files_sort(files,tot_files);


#if 0
  printf("pfiles = %s\n", pfiles);
  pfiles += MaxFnLen;
  printf("pfiles = %s\n", pfiles);
#endif

  // Align columns when print filenames
  max_slen[0]=0; max_slen[1]=0; max_slen[2]=0;
  file_cnt=0;
  while (file_cnt<tot_files) {
    slen=strlen(files[file_cnt]);
    if (slen>max_slen[file_cnt%3]) max_slen[file_cnt%3]=slen;
    ++file_cnt;
  } 
  //
  for (mm=0; mm<3; ++mm) {
    max_slen[mm]+=1, cbuf[mm][max_slen[mm]]='\0';
  }

  // Display files 
  file_cnt=0;
  while (file_cnt<tot_files) {
    slen=strlen(files[file_cnt]);  
    if (*fileType==RxFile) {
      if (slen>=3 && 
          (strncmp((char*) &(files[file_cnt][slen-3]),".in",3)==0)) {
        mm=file_cnt%3;
        memset(cbuf[mm],' ',max_slen[mm]);
        memcpy(cbuf[mm],files[file_cnt],strlen(files[file_cnt]));
        printf("  %3i  %s  ", ++file_cnt,cbuf[mm]);
        if ((file_cnt%3)==0) printf("\n");
      }
     }
    else {
      mm=file_cnt%3;
      memset(cbuf[mm],' ',max_slen[mm]);
      memcpy(cbuf[mm],files[file_cnt],strlen(files[file_cnt]));
      printf("  %3i  %s  ", ++file_cnt,cbuf[mm]);
      if ((file_cnt%3)==0) printf("\n");
   }
  } // while 
  printf("\n");
} // lsfiles_


void getrxfn_(int *id, char* rxfn, int *status)
{
  if ((*id)>0 && (*id)<=file_cnt) {
     memcpy(rxfn,(char*)(files[*id-1]),strlen(files[*id-1])-3);
     rxfn[strlen(files[*id-1])-2]=' ';
     *status=1;
   }
  else {
    *status=0;  // failure
  }
} // getrxfn_


void files_sort(char fils[][MaxFnLen], int nfils)
{
  int i,j;
  char tmp[MaxFnLen];

  for (i=0; i<nfils-1; ++i)
    for (j=i+1; j<nfils; ++j) {
      if (strcmp(fils[i],fils[j])>0) {
	strcpy(tmp,fils[i]); strcpy(fils[i],fils[j]);
        strcpy(fils[j],tmp);
      }
    }
} // files_sort


int match_ext(char* src, char* tgt)
{
  int i;

  i=0;
  while (src[i] != '\0') {
    if (src[i] != '.') ++i;
    else break;
  }
  if (src[i]=='.') {
    return (!strcmp(&src[i],&tgt[1]));
  }
  return 0; 
} // match_ext


void cstdout_(char msg[],int *len)
{
  if (*len) {
    msg[*len]='\0'; printf("%s",msg);
   }
  else printf("\n"); 
}  // cstdout_

//#endif
