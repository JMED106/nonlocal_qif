/********************************/
/* Library for processing files */
/********************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "file.h"
#include "utils.h"

/* Usage:  */
/*************************************************/
/*   FILE *fio[6];				 */
/* 						 */
/*   t_file File;				 */
/*   printf("\nHola soy %s\n","JM");		 */
/*   File.open = OpenFile;			 */
/*   File.open(&((fio)[0]),"archibo_con_b","s"); */
/*************************************************/

t_file LoadFileLibrary(char **filenames, char **modes) {
  t_file file;
  int i;
  if(*filenames != NULL) {
    for(i = 0; filenames[i] != NULL; i++);
    file.n = i;
    file.filenames = (char**) malloc(file.n*sizeof(char*));
    file.mode = (char**) malloc(file.n*sizeof(char*));
    for(i = 0; i < file.n; i++) {
      (file.filenames)[i] = filenames[i];
      (file.mode)[i] = modes[i];
    }
  }
  file.open = OpenFile;
  file.multiopen = OpenMultipleFiles;
  file.closeall = CloseFiles;
  return file;
}

/* === FUNCTION  OpenFile ====================
 * Description:  Opens a file
 *   Variables:  FILE pointer, name and mode (read, ...)
 * ======================================= */

void OpenFile(FILE **files,char filename[], char mode[]) {
  char cmd[200];

  sprintf(cmd,"%s",filename);
  *files = fopen(cmd,mode);
}

/* === FUNCTION  OpenMultipleFiles ====================
 * Description:  Opens multiple files (as OpenFile)
 *   Variables:  Array length needed (?)
 * ======================================= */

void OpenMultipleFiles(FILE *(*files)[], t_file *file) {
  char cmd[200];
  int i;

  for(i=0;i<file->n ;i++ ) {
    sprintf(cmd,"%s",(file->filenames)[i]);
    (*files)[i] = fopen(cmd,(file->mode)[i]);
  }
}

/* === FUNCTION  CloseFile ====================
 * Description:  Closes a file
 *   Variables:  FILE
 * ======================================= */

void CloseFile(FILE **file) {
  fclose(*file);
}

/* === FUNCTION  CloseFiles ====================
 * Description:  Closes multiple files
 *   Variables:  FILE, t_file
 * ======================================= */

void CloseFiles(FILE *(*files)[],t_file *Tfile) {
  do {
    fclose((*files)[Tfile->n -1]);
    Tfile->n--;
  } while(Tfile->n > 0);    
}
