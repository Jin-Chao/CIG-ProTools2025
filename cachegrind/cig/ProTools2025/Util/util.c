#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *varinfo_file = NULL;

int
varinfo_file_open(char *fname)
{
	varinfo_file = fopen(fname, "w");
    if(varinfo_file != NULL)
		return 0;

	return -1;
}

int
varinfo_file_print(char *vname, void* begin, void* end)
{
	if(varinfo_file != NULL && vname != NULL && begin != NULL && end != NULL)
	{
		fprintf(varinfo_file, "%s %p %p\n", vname, begin, end);
		return 0;
	}
	return -1;
}


int
varinfo_file_close()
{
	if(varinfo_file != NULL)
    {
		fflush(varinfo_file);
		fsync(fileno(varinfo_file));
		fclose(varinfo_file);
    }

	return 0;
}

/*return the duration in seconds */
double
CalElapsedTime(struct timeval* tv_begin, struct timeval* tv_end)
{
    double ela_time;
    long ela_secs, ela_usecs;

    if(tv_end->tv_usec >= tv_begin->tv_usec)
    {
        ela_secs  = tv_end->tv_sec  - tv_begin->tv_sec;
        ela_usecs = tv_end->tv_usec - tv_begin->tv_usec;
    } else {
        ela_secs  = tv_end->tv_sec  - tv_begin->tv_sec - 1;
        ela_usecs = tv_end->tv_usec - tv_begin->tv_usec + 1000000;
    }

    ela_usecs += ela_secs * 1000000;
    ela_time = (double)ela_usecs / 1000000;

    return ela_time;
}
