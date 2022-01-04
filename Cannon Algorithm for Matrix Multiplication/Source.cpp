#include <iostream>
#include <mpi.h>
#include<math.h>
#include<stdio.h>
using namespace std;

int ** left_Shift(int **mat, int BlockSize, int size)
{
	int temp;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < BlockSize; j++)
		{
			temp = mat[i][0];

			for (int k = 0; k < size; k++)
			{
				if (k == size - 1)
				{
					continue;
				}
				mat[i][k] = mat[i][k + 1];
			}
			mat[i][size - 1] = temp;
		}
	}
	return mat;
}
int ** Up_Shift(int **mat, int BlockSize, int size)
{
	int temp;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < BlockSize; j++)
		{
			temp = mat[0][i];
			for (int k = 0; k < size; k++)
			{
				if (k == size - 1)
				{
					continue;
				}
				mat[k][i] = mat[k + 1][i];
			}
			mat[size - 1][i] = temp;
		}
	}
	return mat;
}
int **MatrixMultiplication(int **sub_c,int **sub_a,int **sub_b,int BlockSize){
	for (int i = 0; i < BlockSize; ++i)
		for (int j = 0; j < BlockSize; ++j)
			for (int k = 0; k < BlockSize; ++k)
			{
				sub_c[i][j] += sub_a[i][k] * sub_b[k][j];
			}
	return sub_c;
}
void spliting(int size,int BlockSize,int** MAT) {
	MPI_Request request;
	int counter = 0, row = 0, col = 0;
	while (counter < size)
	{
		for (int i = row; i < row + BlockSize; i++)
		{
			for (int j = col; j < col + BlockSize; j++)
			{
				MPI_Isend(&MAT[i][j], 1, MPI_INT, counter, counter, MPI_COMM_WORLD, &request);
			}
		}
		if (col%BlockSize == 0)
		{
			col += BlockSize;
		}
		counter++;
		if (counter%BlockSize == 0)
		{
			row += BlockSize;
			col = 0;
		}
	}
}
void AdditionalOfDoubleMat(int **mat_c,int **mat_Ctemp,int blocksize){
	for (int i = 0; i < blocksize; i++)
	{
		for (int j = 0; j < blocksize; j++)
		{
			mat_c[i][j] += mat_Ctemp[i][j];
		}
	}
}

int main(int argc, char** argv)
{
	/*cout << "Enter the N you want to create a N*N Matrix: ";
	cin >> num;
	cout << endl;
	cout << "Enter Number of the Processes that you want to make Calculate this multiplication:";
	cin >> p;
	cout << endl;*/
	MPI_Init(NULL, NULL);
	int num = 4, p = 4;
	int size;
	int myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;
	MPI_Request request;
	int intRoot = sqrt(p);
	double check = num % intRoot;
	if ((sqrt(p) - floor(sqrt(p))) != 0) {
		printf("[ERROR] Number of processes must be a perfect square!\n");
		MPI_Finalize();
		return 0;
	}
	if (check != 0) {
		printf("[ERROR] Number of rows/columns not divisible by %d!\n", intRoot);
		MPI_Finalize();
		return 0;
	}
	int temp;
	int **MAT_A = (int **)malloc(num * sizeof(int));
	int **MAT_B = (int **)malloc(num * sizeof(int));
	int **MAT_C = (int **)malloc(num * sizeof(int));
	int **MAT_Result = (int **)malloc(num * sizeof(int));
	int BlockSize = num / sqrt(p);
	int **sub_a = (int **)malloc(BlockSize * sizeof(int));
	int **sub_b = (int **)malloc(BlockSize * sizeof(int));
	int **sub_c = (int **)malloc(BlockSize * sizeof(int));
	int **sub_Ctemp = (int **)malloc(BlockSize * sizeof(int));
	for (int i = 0; i < num; i++)
	{
		MAT_A[i] = (int *)malloc(num * sizeof(int*));
		MAT_B[i] = (int *)malloc(num * sizeof(int*));
		MAT_C[i] = (int *)malloc(num * sizeof(int*));
		MAT_Result[i] = (int *)malloc(num * sizeof(int*));
	}
	for (int i = 0; i < BlockSize; i++)
	{
		sub_Ctemp[i] = (int *)malloc(BlockSize * sizeof(int*));
		sub_a[i] = (int *)malloc(BlockSize * sizeof(int*));
		sub_c[i] = (int *)malloc(BlockSize * sizeof(int*));
		sub_b[i] = (int *)malloc(BlockSize * sizeof(int*));
	}
	//Enter 
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++) {
				MAT_B[i][j] = rand() % 30;
				MAT_A[i][j] =  rand() % 30;
				MAT_C[i][j] = 0;
				MAT_Result[i][j] = 0;
		}
	}
	//Display the Original Matrix A & B
	if (myrank == 0) {
		cout << "Matrix A: " << endl;
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				cout << MAT_A[i][j] << "	";
			}
			cout << endl;
		}
		cout << "Matrix B: " << endl;
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				cout << MAT_B[i][j] << "	";
			}
			cout << endl;
		}
	}
	for (int i = 0; i < BlockSize; i++)
	{
		for (int j = 0; j < BlockSize; j++)
		{
			sub_c[i][j] = 0;
			sub_Ctemp[i][j] = 0;
		}
	}
	int shiftCount = 0;
	//left shifting Matrix A in (Ali)gnment
	for (int i = 0; i < num; i++)
	{
		shiftCount = i / BlockSize;
		//if shiftcount==0 we won't shift
		if (shiftCount == 0)
		{
			continue;
		}
		for (int j = 0; j < shiftCount*BlockSize; j++)
		{
			temp = MAT_A[i][0];
			for (int k = 0; k < num; k++)
			{
				MAT_A[i][k] = MAT_A[i][k + 1];
			}
			MAT_A[i][num - 1] = temp;
		}
	}
	//UP Shift Matrix B in (Ali)gnment
	for (int i = 0; i < num; i++)
	{
		shiftCount = i / BlockSize;
		//if shiftcount==0 we won't shift
		if (shiftCount == 0)
		{
			continue;
		}
		for (int j = 0; j < shiftCount*BlockSize; j++)
		{
			temp = MAT_B[0][i];
			for (int k = 0; k < num; k++)
			{
				if (k == num - 1)
				{
					continue;
				}
				MAT_B[k][i] = MAT_B[k + 1][i];
			}
			MAT_B[num - 1][i] = temp;
		}
	}
	//Spliting Matrix A and distribute it to ranks
	if (myrank == 0) {
		spliting(num, BlockSize, MAT_A);
	}
	for (int i = 0; i < BlockSize; i++)
	{
		for (int j = 0; j < BlockSize; j++)
		{
			MPI_Recv(&sub_a[i][j], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		}
	}
	//Spliting Matrix B and distribute it to ranks
	if (myrank == 0) {
		spliting(num, BlockSize, MAT_B);
	}
	for (int i = 0; i < BlockSize; i++)
	{
		for (int j = 0; j < BlockSize; j++)
		{
			MPI_Recv(&sub_b[i][j], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		}
	}
	sub_c = MatrixMultiplication(sub_c, sub_a, sub_b, BlockSize);
	for (int k = 0; k < BlockSize-1; k++)
	{
		for (int i = 0; i < BlockSize; i++)
		{
			for (int j = 0; j < BlockSize; j++)
			{
				sub_Ctemp[i][j] = 0;
			}
		}
		MAT_A = left_Shift(MAT_A, BlockSize, num);
		MAT_B = Up_Shift(MAT_B, BlockSize, num);
		if (myrank == 0)
			spliting(num, BlockSize, MAT_A);
		for (int i = 0; i < BlockSize; i++)
		{
			for (int j = 0; j < BlockSize; j++)
			{
				MPI_Recv(&sub_a[i][j], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
			}
		}
		if (myrank == 0)
			spliting(num, BlockSize, MAT_B);
		for (int i = 0; i < BlockSize; i++)
		{
			for (int j = 0; j < BlockSize; j++)
			{
				MPI_Recv(&sub_b[i][j], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
			}
		}
		sub_Ctemp = MatrixMultiplication(sub_Ctemp, sub_a, sub_b, BlockSize);
		AdditionalOfDoubleMat(sub_c, sub_Ctemp, BlockSize);
	}
	//Collect the Result Matrix C in Rank 0 and Display it	
	for (int i = 0; i < p; i++)
	{
		if (myrank == i) 
		{
				for (int i = 0; i < BlockSize; i++)
				{
					for (int j = 0; j < BlockSize; j++)
					{
						MPI_Isend(&sub_c[i][j], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
					}
				}		
		}
	}
	if (myrank==0)
	{
		int row=0, col=0,count=0;
		for (int k = 0; k < p; k++)
		{
			for (int i = 0; i < BlockSize; i++)
			{
				for (int j = 0; j < BlockSize; j++)
				{
					MPI_Recv(&MAT_C[i][j], 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
					MAT_Result[i + row][j + col] = MAT_C[i][j];
				}
			}
			col += BlockSize;
			count++;
			if (count%BlockSize == 0)
			{
				row += BlockSize;
			    col = 0;
		    }
		}
		cout << "Final Matrix C (The result): " << endl;
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				cout<<MAT_Result[i][j] <<"	";
			}
			cout << endl;
		}
	}
	MPI_Finalize();
	return 0;
}

