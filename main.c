#include <stdio.h>
#include <stdlib.h>

#define N 4

void build_matrix(size_t n, int index[n], int col_id[n], int values[n])	{
	int count = 0;
    index[0] = 0;
    
	for (size_t i = 0; i < n; ++i)	{
        for (size_t j = 0; j < n; ++j)	{
        	if (i > 0) { 
            	values[count] = -1;
				col_id[count] = (i-1) * n + j; 
                count++;
            }
			
			if (i < N - 1)	{
				values[count] = -1;
				col_id[count] = (i+1) * n + j;
                count++;
            }
			
			if (j > 0)	{
				values[count] = -1;
                col_id[count] = i * n + (j-1);
				count++;
			}
			
			if (j < N - 1)	{
				values[count] = -1;
				col_id[count] = i * n + (j+1);
                count++;
			}

            values[count] = 4;
			col_id[count] = i * n + j;
            count++;

			index[i * n + j + 1] = count;
        }
    }
}

void jacobi()	{

}

void gs()	{

}

int main(void)	{
	size_t size = N * N;
    int *index = malloc(sizeof(int) * size);
    int *col_id = malloc(sizeof(int) * size * N);
    int *values = malloc(sizeof(int) * size * N);
	
	if (!index || !col_id || !values)	{
		perror("Error allocation CSR");
		free(index);
		free(col_id);
		free(values);
		exit(EXIT_FAILURE);
	}

	build_matrix(N, index, col_id, values);
	
	for (size_t i = 0; i < size; ++i)	{
		printf("%d ", index[i]);
	}

	printf("\n");

	for (size_t i = 0; i < size * N; ++i)	{
		printf("(%d %d), ", values[i], col_id[i]);
	}

	printf("\n");
	
	free(values);
	free(col_id);
	free(index);

    return 0;
}

