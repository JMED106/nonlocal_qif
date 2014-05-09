#ifndef __FILE__H__
#define __FILE__H__



struct file {
  char **filenames;
  char **mode;
  int n;
  void (*open)(FILE **, char filename[], char mode[]);
  void (*multiopen)(FILE *(*)[], struct file *);
  void (*closeall)(FILE *(*)[], struct file *);
};
typedef struct file t_file;


void OpenFile(FILE **files, char filename[], char mode[]);
void OpenMultipleFiles(FILE *(*files)[], t_file *file);
void CloseFile(FILE **file);
void CloseFiles(FILE *(*files)[],t_file *Tfile);
t_file LoadFileLibrary(char **filenames, char **modes);


#endif
