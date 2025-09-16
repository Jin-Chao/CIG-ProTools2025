#define _GNU_SOURCE
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

#define SIZE_DEFAULT 1024
#define BLOCK_SIZE_DEFAULT 8

int32_t *A = NULL;
int32_t *B = NULL;
int32_t *C = NULL;
int32_t *C_d = NULL;
int32_t *matrix_fill= NULL;

void initialization(int size)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            A[i*size + j] = i*j;
            B[i*size + j] = i*j;
        }
    }

	memset(C, 0, size*size*sizeof(int32_t));
	memset(C_d, 0, size*size*sizeof(int32_t));
}

void clear_cache()
{
    int i, j;
    for (i = 0; i < CLEAR_DEFAULT; i++)
    {
        for (j = 0; j < CLEAR_DEFAULT; j++)
        {
            matrix_fill[i*CLEAR_DEFAULT + j] = i*j;
        }
    }
}

void multiply_simple(int size)
{
    int i, j, k;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
        	for (k = 0; k < size; k++)
        	{
            	C[i*size+j] += A[i*size+k] * B[k*size+j];
			}
        }
    }
}

void multiply_blocked(int size, int block)
{
    int i, j, k;
    int bi, bj, bk;
    for (bi = 0; bi < size; bi += block)
    {
        for (bj = 0; bj < size; bj += block)
        {
        	for (bk = 0; bk < size; bk += block)
        	{
				for(i = bi; i < bi+block; i++)
				{
					for(j = bj; j < bj+block; j++)
					{
						register int32_t sum = 0; 
						for(k = bk; k < bk+block; k++)
						{
            				sum += A[i*size+k] * B[k*size+j];
						}
						C_d[i*size + j] += sum;
					}
				}
			}
        }
    }
}

int verify_results(int size)
{
    int i, j;
    int non_equal = 0;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
           	if(C[i*size+j] != C_d[i*size+j])
            {
				printf("Error: i[%d], j[%d]\n", i, j);
				non_equal = 1;
            }
        }
    }

	return non_equal;
}

void print_usage(char **argv)
{
    printf("Usage:  %s -s size -b block -v var_info \n", *argv);
}

int main(int argc, char *argv[])
{
    int size = SIZE_DEFAULT;
	int block = BLOCK_SIZE_DEFAULT;

    struct timeval begin,end;
    double duration;

    char *var_info = NULL;
    int cha;
    while ((cha = getopt (argc, argv, "s:b:v:h")) != -1)
	{
        switch (cha)
		{
            case 's':
				size  = atol(optarg);
				break;
            case 'b':
				block  = atol(optarg);
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

    if (size <= 0 || block <= 0 || size < block) {
        printf("ERROR: Invalid matrix size: %d or block size: %d\n", size, block);
        exit(EXIT_FAILURE);
    }

    printf("size: %d, block: %d\n", size, block);

    A = (int32_t*)aligned_alloc(64, sizeof(int32_t)*size*size);
    B = (int32_t*)aligned_alloc(64, sizeof(int32_t)*size*size);
    C = (int32_t*)aligned_alloc(64, sizeof(int32_t)*size*size);
    C_d = (int32_t*)aligned_alloc(64, sizeof(int32_t)*size*size);

    matrix_fill = (int32_t*)aligned_alloc(64, sizeof(int32_t)*CLEAR_DEFAULT*CLEAR_DEFAULT);

    if( A == NULL || B == NULL || C == NULL || C_d == NULL || matrix_fill == NULL)
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
            fprintf(var_file, "%s %p %p\n", getName(A), A, A + size*size);
            fprintf(var_file, "%s %p %p\n", getName(B), B, B + size*size);
            fprintf(var_file, "%s %p %p\n", getName(C), C, C + size*size);
            fprintf(var_file, "%s %p %p\n", getName(C_d), C_d, C_d + size*size);

            fclose(var_file);
        }
    }

    //initialization
    initialization(size);

    //clear cache
    clear_cache();

    gettimeofday(&begin, NULL);
    multiply_simple(size);
    gettimeofday(&end, NULL);

    duration = CalElapsedTime(&begin, &end);
    printf("Time: %f with the naive execution\n", duration);

    //clear cache
    clear_cache();

    gettimeofday(&begin, NULL);
    multiply_blocked(size, block);
    gettimeofday(&end, NULL);

    duration = CalElapsedTime(&begin, &end);
    printf("Time: %f with the blocked execution\n", duration);

	//verfication
    verify_results(size);

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

    if(C_d != NULL)
    {
		free(C_d);
		C_d = NULL;
    }

    if(matrix_fill != NULL)
    {
		free(matrix_fill);
		matrix_fill = NULL;
    }

    return 0;
}
