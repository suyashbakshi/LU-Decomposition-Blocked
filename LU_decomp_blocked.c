//@Author Suyash Bakshi
//Code structure for LU_Decomposition() credit: https://www.nersc.gov/users/software/programming-models/openmp/openmp-tasking/openmp-tasking-example-lu/

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

#define size 512
#define b_size 8
#define num_blocks size/b_size

float A[size*size];
void display();

void row_func(int i, int j)
{
	for (int ii = i*b_size; ii < (i*b_size)+(b_size-1); ii++)
	{
		for(int jj = ii+1; jj < b_size; jj++)
		{
			for(int kk = j*b_size; kk < (j*b_size)+b_size; kk++)
			{
				A[(jj*size) + kk] = A[(jj*size) + kk] - (A[(jj*size) + ii] * A[(ii*size) + kk]);
			}
		}
	}
}

void col_func(int i, int j)
{
	for(int ii = i*b_size; ii < (i*b_size)+b_size; ii++)
	{
		for(int jj = j*b_size; jj < (j*b_size)+b_size; jj++)
		{
			A[(jj*size) + ii] /= A[(ii*size) + ii];

			for(int kk = ii + 1; kk < (i*b_size)+b_size; kk++)
			{
				A[(jj*size) + kk] = A[(jj*size) + kk] - (A[(jj*size) + ii] * A[(ii*size) + kk]);
			}
		}

	}
}

void inner_func(int i, int j, int k)
{
	for(int ii = i*b_size; ii < (i*b_size)+b_size; ii++)
	{
		for(int jj = j*b_size; jj < (j*b_size)+b_size; jj++)
		{
			for(int kk = k*b_size; kk < (k*b_size)+b_size; kk++)
			{
				A[(jj*size) + kk] = A[(jj*size) + kk] - (A[(jj*size) + ii] * A[(ii*size) + kk]);
			}
		}
	}
}


void diag_func(int i){
	for(int ii = i*b_size; ii < (i*b_size)+b_size-1; ii++)
	{
		for(int jj = ii+1; jj < (i*b_size)+b_size; jj++)
		{
			A[(jj*size) + ii] /= A[(ii*size) + ii];

			for(int kk = ii+1; kk < (i*b_size) + b_size; kk++)
			{
				A[(jj*size) + kk] = A[(jj*size) + kk] - (A[(jj*size) + ii] * A[(ii*size) + kk]);
			}
		}
	}
}

void init(){
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			A[(i*size) + j] = rand() % 200 + 2;
		}
	}

}

void display(){
	printf("---------------------\n");
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			printf("%f ", A[(i*size + j)]);
		}
		printf("\n");
	}
}

void LU_Decomposition(){
	for(int i=0; i<num_blocks; i++) {

    	diag_func(i);
	    
	    for(int j=i+1; j<num_blocks; j++) {
	        row_func(i, j);
	    }
 
	    for(int j=i+1; j<num_blocks; j++) 
	    {
	        col_func(i, j);
	     
          	for(int k=i+1; k<num_blocks; k++) 
	        {
	            inner_func(i, j, k);
	        }
	    }
	}
}

void main(){

	init();
	// printf("Before\n");
	// display();

	LU_Decomposition();

	// printf("After\n");
	// display();
}
