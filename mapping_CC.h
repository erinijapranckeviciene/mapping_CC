#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_help();
char * read_line(FILE * file_name);
char * word_next(char * line);
int hash(char * c, int length, int type);
int char2int(char * line);
double char2double(char * line);
