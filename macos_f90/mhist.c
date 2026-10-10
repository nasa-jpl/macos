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

/* ------------------------------------------------------------------
 * Sub-prompt cache (always built, regardless of READLINE_LIBRARY).
 *
 * Lets Fortran ACCEPT routines pass the actual prompt text down to
 * readline (when present) so readline renders it and manages cursor
 * state authoritatively.  Replaces the old fragile pattern of
 * Fortran-side WRITE for the prompt + readline-with-empty-prompt for
 * input, which was a recurring source of overwrite/blank-line bugs.
 *
 *   Fortran: CALL set_sub_prompt(promptStr)   [sets cache]
 *            CALL READ_LOH(...)                [normal sub-prompt read]
 *            -- mhist_ sees mp[0]==' ' and uses the cached string --
 *
 * Defined OUTSIDE the READLINE_LIBRARY guard so non-readline builds
 * (smacos / smacos_dvr) still have the symbol.  In those builds the
 * cache is set by Fortran but never consulted; harmless.
 * ------------------------------------------------------------------ */

#include <string.h>      /* memcpy, for the setter below */
#include <stdio.h>       /* fputs/fflush, for print_sub_prompt_ below */

char sub_prompt_buf[256] = "";

void
set_sub_prompt_ (const char *s, int slen)
{
  int n = (slen > 0 && slen < (int)sizeof(sub_prompt_buf) - 2)
              ? slen
              : (int)sizeof(sub_prompt_buf) - 2;
  /* Trim trailing blanks (Fortran strings are blank-padded, plus our
   * caller may have used ICLEN which strips trailing spaces). */
  while (n > 0 && s[n-1] == ' ') n--;
  if (n < 0) n = 0;
  memcpy(sub_prompt_buf, s, n);
  /* Always append a single trailing space so readline's cursor sits
   * one column right of the prompt, matching the legacy
   * ' ',A,'[',A,']: ' format that the old WRITE-the-prompt code
   * produced. */
  sub_prompt_buf[n] = ' ';
  sub_prompt_buf[n+1] = '\0';
}

void
clear_sub_prompt_ (void)
{
  sub_prompt_buf[0] = '\0';
}

/* Non-readline fallback: render the cached sub-prompt before a bare
 * Fortran READ(*).  Readline builds never need this (readline itself
 * renders sub_prompt_buf inside mhist_); builds WITHOUT readline
 * (Windows, fresh clones lacking the pre-built libreadline.a) fall
 * through to a plain READ in PARSE_LOH, which prints nothing -- the
 * "pert freezes" class: the program blocks on input with the prompt
 * invisible.  Called from PARSE_LOH's #else branch. */
void
print_sub_prompt_ (void)
{
  if (sub_prompt_buf[0] != '\0') {
    fputs(sub_prompt_buf, stdout);
    fflush(stdout);
  }
}

#ifdef READLINE_LIBRARY

#if defined (HAVE_CONFIG_H)
#include <config.h>
#endif

#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <unistd.h>      /* isatty(), STDIN_FILENO */

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

/* Weak reference to giza_process_events — resolved at link time only when
 * the Giza library is linked (pgplot target with USE_GIZA=ON).  When linked
 * against classic PGPLOT, the symbol is NULL and the hook becomes a no-op. */
extern void giza_process_events(void) __attribute__((weak));

static int
_macos_rl_event_hook (void)
{
  if (giza_process_events) giza_process_events();
  return 0;
}

int tot_hist=0;

void
mhist_(char* mp, char *cmd)
{
  char *temp, *prompt, cbuf[16], cbf[4];
  int done, non_empty, cmd_id_rerun, cmd_name_rerun, in_hist; 
  HIST_ENTRY **list;

  // Test Calling a Fortran routine
  //fntprint_();  // works!

  if (mp[0]==' ') {
    /* Sub-prompt.  Use the cached prompt string set by the Fortran
     * caller via set_sub_prompt_().  Empty string fallback if the
     * caller forgot to set it (should not normally happen). */
    prompt = (sub_prompt_buf[0] != '\0') ? sub_prompt_buf : "";
  }
  else {
    cbuf[0] = ' ';
    strncpy(&cbuf[1],mp,5); cbuf[6]='\0';
    prompt = strcat(cbuf,"> ");
  }
  temp = (char *) NULL;
  done = 0;

  /* Let Giza repaint plot windows while we wait at the prompt */
  rl_event_hook = _macos_rl_event_hook;

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

      /* Non-TTY transcript echo.
       *
       * On piped/redirected stdin the kernel doesn't echo, and
       * readline's own non-TTY echo path is uneven.  When this was a
       * sub-prompt (cache was set), echo "<resp>\n" so the transcript
       * records what was typed and the next prompt starts on a fresh
       * line.  No-op on TTY -- readline owns rendering and the
       * cached prompt lets it manage cursor authoritatively, so no
       * manual newline is needed (was the source of cosmetic blank
       * lines in earlier iterations).
       *
       * After this, clear the sub-prompt cache so a stale entry
       * doesn't leak into the next call.
       */
      if (mp[0] == ' ') {
        if (!isatty(STDIN_FILENO)) {
          printf("%s\n", temp);
          fflush(stdout);
        }
        sub_prompt_buf[0] = '\0';
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

#endif  /* READLINE_LIBRARY */
