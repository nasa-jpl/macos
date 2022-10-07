/***********************************************************************
        Begin file mhist.c
     +----------------------------------------------------------------+
     |  Copyright (C) 1995-2008, California Institute of Technology.  |
     |  U.S. Government Sponsorship Is Acknowledged.                  |
     +----------------------------------------------------------------+
***********************************************************************/

//
// Utility functions for MACOS command history and recall 
// John Z. Lou, Jet Propulsion Laboratory
// Last updated: 08/2008 
//

#ifdef READLINE_LIBRARY

#if defined (HAVE_CONFIG_H)
#include <config.h>
#endif

#include <stdio.h>
#include <sys/types.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#else 
extern void exit();
#endif

#ifdef READLINE_LIBRARY
#  include "readline.h"
#  include "history.h"
#else
#  include <readline/readline.h>
#  include <readline/history.h>
#endif

extern HIST_ENTRY **history_list();
extern void fntprint();
void show_history(HIST_ENTRY**);

int tot_hist=0;

void
mhist_(char* mp, char *cmd)
{
  char *temp, *prompt, cbuf[16], cbf[4];
  int done, non_empty, cmd_id_rerun, cmd_name_rerun, in_hist; 
  HIST_ENTRY **list;

  // Test Calling a Fortran routine
  //fntprint_();  // works!

  if (mp[0]==' ') 
    prompt = (char*) NULL;
  else {
    cbuf[0] = ' ';
    strncpy(&cbuf[1],mp,5); cbuf[6]='\0';
    prompt = strcat(cbuf,"> ");
  }
  temp = (char *) NULL;
  done = 0;

  while (!done) {
      register int i, j, k;
      char *temp2;
      done = 1;
      temp = readline(prompt);

      // Test for EOF. 
      if (!temp) {
	fprintf (stderr,"  **mhist: Invalid command input, quit!\n");
	exit (1);
      }

      non_empty = 0;
      for (i=0; i<strlen(temp); ++i) {
        if (temp[i] == '%' && non_empty == 0) break;
        if (temp[i] != ' ') { 
          temp = &temp[i], non_empty = 1;
          break; 
        } 
      }
      //printf(" **mhist: temp =%s\n",temp);

      cmd_id_rerun = (strncmp(temp,"!",1)==0) &&
                     (strlen(temp)>1) && isdigit(temp[1]);

      cmd_name_rerun = 0;
      if ((strncmp(temp,"!",1)==0) && (strlen(temp)>1)) {
 	list=history_list();
	if (list=history_list()) {
	  for (i=0; list[i]; i++) {
	    if (strncmp(&temp[1],list[i]->line,strlen(&temp[1])) == 0)
	      temp2=list[i]->line, cmd_name_rerun = 1; 
          } // for
	}	
      }

      // If there is anything on the line, print and remember it. 

      if (non_empty && !cmd_id_rerun && !cmd_name_rerun) {
	  //fprintf (stderr, "%s\r\n", temp);
	  in_hist = 0; 
	  if (list=history_list()) {
            for (i=0; list[i]; i++) {
	      if (strcmp(temp,list[i]->line) == 0) { 
                in_hist = 1; break;
	      } 
	    }
	  }
	  if (!in_hist) { add_history(temp); ++tot_hist; }
      }

      // Check for `command' that we handle.
      if (strlen(temp)>=3) {
	for (i=0; i<3; ++i) cbf[i]=toupper(temp[i]);
       }
      else cbf[0] = ' ';
      

      if (strncmp(cbf,"HIS",3)==0) {
          list = history_list();
	  show_history(list);
#if 0
	  if (list) {
            for (i=0; list[i]; i++)
              fprintf (stderr, "    %d: %s\r\n", i, list[i]->line);
          } 
#endif
          done = 0;
       }
      else if (strncmp(cbf,"CLE",3)==0) {
          clear_history(); tot_hist=0;
	  fprintf (stderr, "    Command history cleared!\n");
          done = 0;
       }	
      else if (cmd_id_rerun) {
	//rl_stuff_char('h');
	j = atoi((char*)&temp[1]);
	if (j>=1 && j<=tot_hist && (list=history_list())) {
	  j--; // starting from 0 to match cmd index in history
          for (i=0; list[i]; i++)
	    if (j==i) {
	      if (strlen(list[i]->line)>=3) {
                for (k=0; k<3; ++k) 
	          cbf[k]=toupper(list[i]->line[k]);
	       }
	      else cbf[0] = ' ';
	      fprintf(stderr, "   ** Run command: %s\n",list[i]->line);
	      if (strncmp(cbf,"HIS",3) == 0) {
	         show_history(list);
	         done = 0;
	       }
	      else {
	        strncpy(cmd,list[i]->line,strlen(list[i]->line));
	      }
	      break;
	    }  // j==i 
	  }
	 else fprintf(stderr, "   ** Command not found in history!\n");
       }  // cmd_id_return
      else if (cmd_name_rerun) {
	  if (strlen(temp2)>=3) {
	    for (k=0; k<3; ++k) cbf[k]=toupper(temp2[k]);
	  }
	  if (strncmp(cbf,"HIS",3) == 0) {
	    show_history(list);
            done = 0;
	   }
	  else {
	    strncpy(cmd,temp2,strlen(temp2));
	  }
       }
      else {
	  strncpy(cmd,temp,strlen(temp));
      }
      free(temp);
    } // while
} // mhist_


void
show_history(HIST_ENTRY **lst)
{
  lst = history_list();
  if (lst) {
    register int i;
    for (i=0; lst[i]; i++)
      fprintf (stderr, "    %d: %s\r\n", i+1, lst[i]->line);
  }
} // show_history


void
maddh_(char *cmd)
{
  int in_his, non_empty, i;
  HIST_ENTRY **lst;

  non_empty = 0;
  for (i=0; i<strlen(cmd); ++i) {
    if (cmd[i] == '%' && non_empty == 0) {
      // first non-blank char is '%', return
      return;
    }  
    if (cmd[i] != ' ') { 
      cmd = &cmd[i], non_empty = 1; 
      break; 
    }
  }

  if (non_empty) {
    in_his = 0;
    if (lst=history_list()) {
      for (i=0; lst[i]; i++) {
        if (strcmp(cmd,lst[i]->line) == 0) {
          in_his = 1; break;
        }
      }
    }
    if (!in_his) { add_history(cmd); ++tot_hist; }
  }
}  // maddh_
 
#endif
 
