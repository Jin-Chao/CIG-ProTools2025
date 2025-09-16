#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define getName(var)  #var
#define getNameWithPostfix(var, num) getName(var) #num

#define CLEAR_DEFAULT 1024

int varinfo_file_open(char *fname);
int varinfo_file_close();
int varinfo_file_print(char *vname, void* begin, void* end);

/*return the duration in seconds */
double CalElapsedTime(struct timeval* tv_begin, struct timeval* tv_end);

