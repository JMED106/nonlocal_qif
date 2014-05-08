#ifndef __FILE__H__
#define __FILE__H__



struct file {
  void (*open)(FILE **, char *filename, char *mode);
  void (*multiopen)(FILE *(*)[], char **filenames, char **mode,int n);
};
typedef struct file t_file;

void OpenFile(FILE **files, char *filename, char *mode);
void OpenMultipleFiles(FILE *(*files)[], char **filenames, char **mode, int n);
t_file LoadFileLibrary(void);


#endif
