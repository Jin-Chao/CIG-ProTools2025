#define _GNU_SOURCE
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

#define SIZE_DEFAULT 1048576

#define getName(var)  #var

int64_t *A = NULL;
int64_t *B = NULL;
int64_t *C = NULL;
int64_t *D = NULL;
int64_t *matrix_fill= NULL;

//initialization
/*
void initialization(int size)
{
    register int i;

    for (i = 0; i < size; i++)
    {
        A[i] = i*i;
        B[i] = i*i;
        C[i] = i*i;
        D[i] = i*i;
    }

}*/

void initialization(int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        A[i] = i*i;
    }

    for (i = 0; i < size; i++)
    {
        B[i] = i*i;
    }

    for (i = 0; i < size; i++)
    {
        C[i] = i*i;
    }

    for (i = 0; i < size; i++)
    {
        D[i] = i*i;
    }
}

void clear_cache()
{
    register int i, j;

    for (i = 0; i < CLEAR_DEFAULT; i++)
    {
        for (j = 0; j < CLEAR_DEFAULT; j++)
        {
            matrix_fill[i*CLEAR_DEFAULT + j] = i*j;
        }
    }
}

void verify_results(int size)
{
   // printf("Res: %ld\n", res);
}

void conflict(int size)
{
    register int i;

    for (i = 0; i < size; i++)
    {
       A[i] = B[i] + C[i];
    }
}

void print_usage(char **argv)
{
    printf("Usage:  %s -s size -v var_info \n", *argv);
}

int main(int argc, char *argv[])
{
    int size = SIZE_DEFAULT;

    struct timeval begin,end;
    double duration;

    char *var_info = NULL;
    int cha;
    while ((cha = getopt (argc, argv, "s:d:v:h")) != -1)
	{
        switch (cha)
		{
            case 's':
				size  = atol(optarg);
				break;
            case 'v':
				var_info = optarg;
				break;
            case 'h':
            default:
				print_usage(argv);
				exit(0);
		}
	}

    if (size < 1) {
        printf("ERROR: Invalid matrix size: %d\n", size);
        exit(EXIT_FAILURE);
    }

    printf("size: %d\n", size);

    A = (int64_t*)aligned_alloc(64, sizeof(int64_t)*size);
    B = (int64_t*)aligned_alloc(64, sizeof(int64_t)*size);
    C = (int64_t*)aligned_alloc(64, sizeof(int64_t)*size);
    D = (int64_t*)aligned_alloc(64, sizeof(int64_t)*size);

    matrix_fill = (int64_t*)aligned_alloc(64, sizeof(int64_t)*CLEAR_DEFAULT*CLEAR_DEFAULT);

    if( A == NULL || B == NULL || C == NULL || D == NULL || matrix_fill == NULL)
    {
        printf("ERROR: %s cannot allocate memory!n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if(var_info != NULL)
    {
        FILE *var_file = NULL;

        var_file = fopen(var_info, "w");
        if(var_file != NULL)
        {
            fprintf(var_file, "%s %p %p\n", getName(A), A, A + size);
            fprintf(var_file, "%s %p %p\n", getName(B), B, B + size);
            fprintf(var_file, "%s %p %p\n", getName(C), C, C + size);
            fprintf(var_file, "%s %p %p\n", getName(D), D, D + size);
            fclose(var_file);
        }
    }

    //initialization
    initialization(size);

    //clear cache
    clear_cache();

    gettimeofday(&begin, NULL);
    conflict(size);
    gettimeofday(&end, NULL);
    verify_results(size);

    duration = CalElapsedTime(&begin, &end);
    printf("Time: %f with separated loops \n", duration);

    //clear cache again
    clear_cache();
	//
	//release memory
    if(A != NULL)
    {
		free(A);
		A = NULL;
    }

    if(B != NULL)
    {
		free(B);
		B = NULL;
    }

    if(C != NULL)
    {
		free(C);
		C = NULL;
    }

    if(D != NULL)
    {
		free(D);
		D = NULL;
    }

    if(matrix_fill != NULL)
    {
		free(matrix_fill);
		matrix_fill = NULL;
    }

    return 0;
}
