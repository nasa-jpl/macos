/***********************************************************************
        Begin file utilsub_c.c 
     +----------------------------------------------------------------+
     |  Copyright (C) 1995-2026, California Institute of Technology.  |
     |  U.S. Government Sponsorship Is Acknowledged.                  |
     +----------------------------------------------------------------+
***********************************************************************/

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h> 

/* --- 1. Cross-platform Header & Mapping Logic --- */
#if defined(MSWIN) || defined(_WIN32) || defined(_WIN64)
#   include <io.h>
#   include <process.h>
#   include <direct.h>
#   define getcwd _getcwd
#   define chdir  _chdir
#   define F_OK 0
    /* Windows: ifx /assume:underscore produces UPPERCASE with trailing underscore */
#   define FTN_NAME(lower, upper) upper##_
#elif defined(LNXOS) || defined(__linux__)
#   include <unistd.h>
#   include <dirent.h>
#   include <sys/types.h>
/* Linux Intel ifx: default lowercase with trailing underscore */
#   define FTN_NAME(lower, upper) lower##_
#else
#   include <unistd.h>
#   include <dirent.h>
#   include <sys/types.h>
    /* Other (gfortran): lowercase with trailing underscore */
#   define FTN_NAME(lower, upper) lower##_
#endif

#define MaxFileCnt  1000
#define MaxFnLen    132

/* Global state for file listing */
char files[MaxFileCnt][MaxFnLen];
int file_cnt = 0;

/* Internal Utility Prototypes */
void files_sort(char[][MaxFnLen], int);
int match_ext(char*, char*);

/* --- 2. OS Bridge Functions --- */

/* Executes a system command */
void FTN_NAME(shellcmd, SHELLCMD)(char *cmd, int *clen)
{
    cmd[*clen] = '\0';
    /* Offset by 1 if the Fortran string carries a control character/space */
    system(&cmd[0]);  
}

/* Prints the Current Working Directory */
void FTN_NAME(pwd, PWD)()
{
    char dirname[256];
    if (getcwd(dirname, 256)) 
        printf("  Current Folder: \n  %s\n\n", dirname);
}

/* Changes the Current Working Directory */
void FTN_NAME(chgdir, CHGDIR)(char *dirname, int* clen)
{
    char curr_dirname[256];
    dirname[*clen] = '\0'; 
    if (chdir(dirname) == 0) {
        if (getcwd(curr_dirname, 256)) 
            printf("  Current Folder: \n  %s\n\n", curr_dirname);
    } else {
        perror("  Error changing directory");
    }
}

/* Launches a text editor */
void FTN_NAME(starteditor, STARTEDITOR)(char *cmd, int* cmd_len, char* filename, int* clen)
{ 
    char edit_cmd[512];
    cmd[*cmd_len] = '\0';
    filename[*clen] = '\0';
    sprintf(edit_cmd, "%s %s", cmd, filename);
    system(edit_cmd);
}

/* --- 3. Directory Listing (LSFILES) --- */

void FTN_NAME(lsfiles, LSFILES)(int* fileType, char* cmdArg, int* clen)
{
#if defined(MSWIN) || defined(_WIN32)
    /* Windows uses a different directory API (FindNextFile). 
       For now, we provide a basic 'dir' bridge for compatibility. */
    printf("  Executing directory listing via Windows Shell...\n");
    if (*fileType == 0) system("dir /B *.in");
    else system("dir /B");
    return;
#else
    char dirname[256];
    DIR *dirp;
    struct dirent *dp;
    int slen;

    getcwd(dirname, 256);
    printf("  Listing Files in: %s\n", dirname);
    
    if (*clen > 0) cmdArg[*clen] = '\0';  

    dirp = opendir("."); 
    file_cnt = 0;
    while (dirp && (dp = readdir(dirp)) != NULL) {
        if (file_cnt >= MaxFileCnt) break;
        if (*clen == 0 && dp->d_name[0] == '.') continue;
        
        slen = strlen(dp->d_name);
        /* Filter for .in files (RxFile) */
        if (*fileType == 0) {
            if (slen >= 3 && (strcmp(&dp->d_name[slen-3], ".in") == 0)) {
                strncpy(files[file_cnt++], dp->d_name, MaxFnLen);
            }
        } else {
            strncpy(files[file_cnt++], dp->d_name, MaxFnLen);
        }
    } 
    if (dirp) closedir(dirp);

    files_sort(files, file_cnt);
    for (int i = 0; i < file_cnt; i++) {
        printf("  %3i  %s  ", i + 1, files[i]);
        if ((i + 1) % 3 == 0) printf("\n");
    }
    printf("\n");
#endif
}

/* Returns a filename from the internal list based on index ID */
void FTN_NAME(getrxfn, GETRXFN)(int *id, char* rxfn, int *status)
{
    if ((*id) > 0 && (*id) <= file_cnt) {
        int len = strlen(files[*id-1]);
        memcpy(rxfn, files[*id-1], len);
        rxfn[len] = '\0';
        *status = 1;
    } else {
        *status = 0; 
    }
}

/* Standard output bridge */
void FTN_NAME(cstdout, CSTDOUT)(char msg[], int *len)
{
    if (*len > 0) {
        char* buf = malloc(*len + 1);
        if (buf) {
            memcpy(buf, msg, *len);
            buf[*len] = '\0';
            printf("%s", buf);
            free(buf);
        }
    } else {
        printf("\n"); 
    }
}

/* --- 4. Internal Helpers --- */

void files_sort(char fils[][MaxFnLen], int nfils)
{
    char tmp[MaxFnLen];
    for (int i = 0; i < nfils - 1; ++i) {
        for (int j = i + 1; j < nfils; ++j) {
            if (strcmp(fils[i], fils[j]) > 0) {
                strcpy(tmp, fils[i]);
                strcpy(fils[i], fils[j]);
                strcpy(fils[j], tmp);
            }
        }
    }
}

/*added by Luis Marchen*/
#if defined(MSWIN) || defined(_WIN32) || defined(_WIN64) || defined(LNXOS) || defined(__linux__)
#include <time.h>

void FDATE_(char *date, int len)
{
    time_t t = time(NULL);
    char *s = ctime(&t);
    memcpy(date, s, 24);
}

float DTIME_(float *time)
{
    float user = (float)clock() / CLOCKS_PER_SEC;
    time[0] = user;
    time[1] = 0.0f;
    return user;
}

#endif

/* GETENV*/
void GETENV_(char *name, char *value, int name_len, int value_len)
{
    char tmp[256];
    int len = name_len < 255 ? name_len : 255;
    memcpy(tmp, name, len);
    tmp[len] = '\0';
    /* trim trailing spaces (Fortran pads strings with spaces) */
    while (len > 0 && tmp[len-1] == ' ') tmp[--len] = '\0';
    char *result = getenv(tmp);
    if (result) {
        int rlen = strlen(result);
        if (rlen > value_len) rlen = value_len;
        memcpy(value, result, rlen);
        memset(value + rlen, ' ', value_len - rlen);
    } else {
        memset(value, ' ', value_len);
    }
}