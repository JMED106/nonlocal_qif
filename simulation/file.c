/********************************/
/* Library for processing files */
/********************************/

#include<stdio.h>
#include<string.h>

#include "file.h"

/* Usage:  */
/*************************************************/
/*   FILE *fio[6];				 */
/* 						 */
/*   t_file File;				 */
/*   printf("\nHola soy %s\n","JM");		 */
/*   File.open = OpenFile;			 */
/*   File.open(&((fio)[0]),"archibo_con_b","s"); */
/*************************************************/

t_file LoadFileLibrary(void) {
t_file file;
 file.open = OpenFile;
 file.multiopen = OpenMultipleFiles;
 return file;
}

/* === FUNCTION  OpenFile ====================
 * Description:  Opens a file
 *   Variables:  FILE pointer, name and mode (read, ...)
 * ======================================= */

void OpenFile(FILE **files, char *filename, char *mode) {
  char cmd[200];

  sprintf(cmd,"%s",filename);
  *files = fopen(cmd,mode);
}

/* === FUNCTION  OpenMultipleFiles ====================
 * Description:  Opens multiple files (as OpenFile)
 *   Variables:  Array length needed (?)
 * ======================================= */

void OpenMultipleFiles(FILE *(*files)[], char **filenames, char **mode, int n) {
  char cmd[200];
  int i;
  for(i=0;i<n ;i++ ) {
    sprintf(cmd,"%s",filenames[i]);
    (*files)[i] = fopen(cmd,mode[i]);
  }
}
