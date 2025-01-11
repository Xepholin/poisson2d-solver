#include <stdio.h>
#include <stdlib.h>

void build_matrix(size_t n, int index[n*n], int col_id[5 * (n*n) - 4*n], int values[5 * (n*n) - 4*n])	{
	int count = 0;
    index[0] = 0;
    
	for (size_t i = 0; i < n; ++i)	{
        for (size_t j = 0; j < n; ++j)	{
        	if (i > 0) { 
            	values[count] = -1;
				col_id[count] = (i-1) * n + j; 
                count++;
            }
			
			if (i < n - 1)	{
				values[count] = -1;
				col_id[count] = (i+1) * n + j;
                count++;
            }
			
			if (j > 0)	{
				values[count] = -1;
                col_id[count] = i * n + (j-1);
				count++;
			}
			
			if (j < n - 1)	{
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

void jacobi(size_t n, int a[n], int threshold, size_t max_iter)	{
	//
}

void gs(size_t n, int a[n], int threshold, size_t max_iter)	{
	//
}

int main(int argc, char *argv[])	{
	if (argc != 2)	{
		perror("expected 1 argument (dim)");
		exit(EXIT_FAILURE);
	}

	size_t n = strtol(argv[1], NULL, 10);
	
	if (n < 2)	{
		perror("first argument dim > 1");
		exit(EXIT_FAILURE);
	}

	size_t size = n * n;
	size_t total = 5 * size - 4 * n;
	int *index = aligned_alloc(32, sizeof(int) * size);
    int *col_id = aligned_alloc(32, sizeof(int) * total);
    int *values = aligned_alloc(32, sizeof(int) * total);
	
	if (!index || !col_id || !values)	{
		perror("Error allocation CSR");
		free(index);
		free(col_id);
		free(values);
		exit(EXIT_FAILURE);
	}

	build_matrix(n, index, col_id, values);
	
	for (size_t i = 0; i < size; ++i)	{
		printf("%d ", index[i]);
	}

	printf("\n");

	for (size_t i = 0; i < total ; ++i)	{
		printf("(%d %d), ", values[i], col_id[i]);
	}

	printf("\n");
	
	free(values);
	free(col_id);
	free(index);

    return 0;
}

