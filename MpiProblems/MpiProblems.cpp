// MpiProblems.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "mpi.h"
#include "vector"
#include <time.h>
#include <stdlib.h>
#include "Complex.h"
using namespace std;

int process_rank;
//bool isPrime(int n)
//{
//	for (int i = 2; i * i <= n; i++)
//	{
//		if (n % i == 0)
//		{
//			return false;
//		}
//	}
//	return true;
//}
//void PointToPoint(int argc, char** argv)
//{
//	int rank;
//	int size;
//	MPI_Status status;
//	int result;
//	MPI_Request requests[100];
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	if (rank == 0)
//	{
//		int numbers[17] = { 2,3,3,3,3,4,5,19,22,25,24,50,44,33,21,7,13 };
//		int n = 17;
//		int mode = n % (size - 1);
//		int countOfSend = n / (size - 1);
//
//		for (int i = 1; i < size; i++)
//		{
//			MPI_Send(&countOfSend, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
//			MPI_Isend(&numbers[(i - 1) * countOfSend], countOfSend, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[i]);
//		}
//
//		for (int i = 1; i < size; i++)
//		{
//			int recevCount;
//
//
//			MPI_Recv(&recevCount, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
//
//			int* subPrimes = new int[recevCount];
//
//			MPI_Recv(subPrimes, recevCount, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
//			for (int j = 0; j < recevCount; j++)
//			{
//				printf("%d\n", subPrimes[j]);
//			}
//		}
//		for (int i = n - mode; i < n; i++)
//		{
//			if (isPrime(numbers[i]))
//			{
//				printf("%d\n", numbers[i]);
//			}
//		}
//	}
//	else
//	{
//		int recevCount;
//		MPI_Irecv(&recevCount, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[rank]);
//		MPI_Wait(&requests[rank], MPI_STATUS_IGNORE);
//		int* slave = new int[recevCount];
//		MPI_Recv(slave, recevCount, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//		int* sendingArray = new int[recevCount];
//		int sendCount = 0;
//
//
//		for (int i = 0; i < recevCount; i++)
//		{
//			if (isPrime(slave[i]))
//			{
//				sendingArray[sendCount] = slave[i];
//				sendCount++;
//				//	printf("%d\n", sendingArray[sendCount++]);
//
//			}
//		}
//		MPI_Send(&sendCount, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
//		MPI_Send(sendingArray, sendCount, MPI_INT, 0, 2, MPI_COMM_WORLD);
//	}
//
//	MPI_Finalize();
//
//}
//void CollectivSending(int argc, char** argv)
//{
//	int rank;
//	int size;
//	MPI_Status status;
//	int result;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	int* recBuffer;
//
//	int numbers[100][150];
//	recBuffer = new int[101];
//	int count;
//
//	if (rank == 0)
//	{
//		srand(time(NULL) + rank);
//
//		//int numbers[100][150];
//		for (int i = 0; i < 100; i++)
//		{
//			for (int j = 0; j < 150; j++)
//			{
//				numbers[i][j] = rand() % 10 + 1;
//			}
//		}
//
//	}
//
//	MPI_Scatter(&numbers[rank][100], 100, MPI_INT, recBuffer, 100, MPI_INT, 0, MPI_COMM_WORLD);
//	for (int i = 0; i < 100; i++)
//	{
//		printf("%d\n", recBuffer[i]);
//		printf("Rankk   ");
//		printf("%d\n", rank);
//	}
//
//	MPI_Finalize();
//}

void AddComplex(Complex* invec, Complex* incountvc, int* len, MPI_Datatype* dtype)
{
	for (int i = 0; i < *len; i++)
	{
		printf("AddComplex Before rank = %d in %d %d out %d %d\n", process_rank, invec[i].real, invec[i].imag,
			incountvc[i].real, incountvc[i].imag);
		incountvc[i] = incountvc[i] + invec[i];
		printf("AddComplex rank = %d in %d %d out %d %d\n", process_rank, invec[i].real, invec[i].imag,
			incountvc[i].real, incountvc[i].imag);
	}

}

void ComplexNumbers(int argc, char** argv)
{
	int size;
	MPI_Status status;
	//Complex result[6];
	Complex correct_result(0, 0);
	//Complex data[6] = { Complex(2, 3),Complex(4, 5), Complex(8, 3), Complex(10, 2), Complex(1, 4), Complex(5, 7) };
	double data[6] = { 2,3, 4, 5, 6, 7 };
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

	int n = 4;
	Complex result[4];
	Complex numbers[4] = { Complex(2,3), Complex(4,5), Complex(1,7), Complex(3, 4) };
	MPI_Datatype complex_type, complex_array_type;
	int blocklens[2];
	MPI_Datatype types[2];
	MPI_Aint indices[2];
	blocklens[0] = blocklens[1] = 1;
	types[0] = MPI_DOUBLE;
	types[1] = MPI_DOUBLE;
	MPI_Get_address(&numbers[0].real, &indices[0]);
	MPI_Get_address(&numbers[0].imag, &indices[1]);
	indices[1] = indices[1] - indices[0];
	indices[0] = 0;
	MPI_Type_create_struct(2, blocklens, indices, types, &complex_type);
	MPI_Type_commit(&complex_type);
	//MPI_Type_contiguous(2, MPI_DOUBLE, &complex_array_type);
	//MPI_Type_commit(&complex_array_type);
	MPI_Op op;
	MPI_Op_create((MPI_User_function*)AddComplex, 1, &op);
	MPI_Reduce(&numbers, &result, 4, complex_type, op, 0, MPI_COMM_WORLD);
	if (process_rank == 0)
	{
		for (int i = 0; i < 4; i++)
		{
			printf("%f + %fi\n", result[i].real, result[i].imag);
		}
	}
	MPI_Op_free(&op);
	//MPI_Type_free(&complex_array_type);
	MPI_Type_free(&complex_type);
}

void CollectiveGaterV(int argc, char** argv)
{
	int size;
	int* recBuffer = nullptr;
	int* displ = nullptr;
	int* rcounts = nullptr;
	int	n = 3;

	int mat[9] = { 3,4,14,
				   1,8,9,
				   2,12,5 };



	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	if (process_rank == 0)
	{
		recBuffer = (int*)malloc(6 * sizeof(int));

		displ = (int*)malloc(n * sizeof(int));
		rcounts = (int*)malloc(n * sizeof(int));

		displ[0] = 0;
		//rcounts[0] = 1;                                           //Global Lower Dioganal Matrix Elements;
		//rcounts[0] = 3;                                             //Global Upper Dioganal Matrix Elements;
		rcounts[0] = 1;                                             //Secondary Lower Dioganal Matrix Elements;
		for (size_t i = 1; i < size; i++)
		{
			//rcounts[i] = i + 1;  
			//displ[i] = rcounts[i - 1] + displ[i - 1];             //Global Lower Dioganal Matrix Elements;

			//rcounts[i] = n - i;
			//displ[i] = rcounts[i - 1] + displ[i - 1];               //Global Upper Dioganal Matrix Elements;	

			rcounts[i] = i + 1;
			displ[i] = rcounts[i - 1] + displ[i - 1];               //Global Upper Dioganal Matrix Elements;

			printf("%d\n", displ[i]);
			printf("%d %d\n", rcounts[i], process_rank);
			printf("\n", "----------------");
		}
	}

	//MPI_Gatherv(&mat[process_rank * n], process_rank + 1, MPI_INT, recBuffer, rcounts, displ, MPI_INT, 0, MPI_COMM_WORLD);                  //Global Lower Dioganal Matrix Elements;
	//MPI_Gatherv(&mat[process_rank * n + process_rank], n - process_rank, MPI_INT, recBuffer, rcounts, displ, MPI_INT, 0, MPI_COMM_WORLD);   //Global Upper Dioganal Matrix Elements;
	MPI_Gatherv(&mat[process_rank * n + n - (process_rank + 1)], process_rank + 1, MPI_INT, recBuffer, rcounts, displ, MPI_INT, 0, MPI_COMM_WORLD);     //Secondary Lower Dioganal Matrix Elements;

	if (process_rank == 0)
	{
		for (int j = 0; j < 6; j++)
			printf("%d\n", recBuffer[j]);
	}

	MPI_Finalize();
}

void CollectiveScatterV(int argc, char** argv)
{
	int size;
	int* displ = nullptr;
	int* scounts = nullptr;
	int	n = 3;

	int* mat = nullptr;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	int* recBuffer = (int*)malloc(3 * sizeof(int));

	if (process_rank == 0)
	{
		mat = new int[9]
		{
		  3,4,14,
		  1,8,9,
		  2,12,5
		};

		displ = (int*)malloc(n * sizeof(int));
		scounts = (int*)malloc(n * sizeof(int));

		displ[0] = 0;
		scounts[0] = 1;
		for (size_t i = 1; i < size; i++)
		{
			scounts[i] = i + 1;
			displ[i] = i * n;
		}

	}

	MPI_Scatterv(mat, scounts, displ, MPI_INT, recBuffer, process_rank + 1, MPI_INT, 0, MPI_COMM_WORLD);

	for (int j = 0; j < 3; j++)
		printf("%d %d\n", recBuffer[j], process_rank);

	MPI_Finalize();
}

void MpiPackAndUnpackExample(int argc, char** argv)
{
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

	int position = 0;
	int i;
	char* c;
	char* buffer[10];
	MPI_Status status;

	if (process_rank == 0)
	{
		i = 15;
		c = new char[2]{ 'd','c' };

		MPI_Pack(&i, 1, MPI_INT, buffer, 10, &position, MPI_COMM_WORLD);
		printf("%d %d\n", position, process_rank);

		MPI_Pack(c, 2, MPI_CHAR, buffer, 10, &position, MPI_COMM_WORLD);
		printf("%d %d\n", position, process_rank);
		MPI_Send(buffer, 10, MPI_PACKED, 1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(buffer, 10, MPI_PACKED, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Unpack(buffer, 10, &position, &i, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Unpack(buffer, 10, &position, &c, 2, MPI_CHAR, MPI_COMM_WORLD);

		printf("Variable %d Rank : %d\n", i, process_rank);
		printf("Array %.2s Rank : %d\n", &c, process_rank);
	}

	MPI_Finalize();
}

void MPI_CONTIGUOUS(int argc, char** argv)
{
	int size;
	int* displ = nullptr;
	int* scounts = nullptr;
	int	n = 3;

	int* mat = nullptr;
	int* recBuffer = (int*)malloc(sizeof(int) * n);

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	MPI_Datatype newType;
	MPI_Type_contiguous(3, MPI_INT, &newType);
	MPI_Type_commit(&newType);

	if (process_rank == 0)
	{
		mat = new int[9]
		{
		  3,4,14,
		  1,8,9,
		  2,12,5
		};

	}
	MPI_Scatter(&mat[0], 1, newType, recBuffer, 3, MPI_INT, 0, MPI_COMM_WORLD);

	for (int j = 0; j < 3; j++)
		printf("%d %d\n", recBuffer[j], process_rank);

	MPI_Type_free(&newType);
	MPI_Finalize();

}
void MPI_Vector(int size, char** argv)
{
	int	n = 3;
	int* recBuffer = (int*)malloc(sizeof(int) * n);
	int* mat = nullptr;
	MPI_Datatype newType;
	MPI_Type_vector(n, 1, n, MPI_INT, &newType);
	MPI_Type_commit(&newType);
	if (rank == 0)
	{
		mat = new int[n * n]
		{
		  3,4,14,
		  1,8,9,
		  2,12,5
		};
		for (int i = 0; i < size - 1; i++)
			MPI_Send(&mat[i], 1, newType, i + 1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(recBuffer, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int j = 0; j < 3; j++)
			printf("%d %d\n", recBuffer[j], rank);
	}
	MPI_Type_free(&newType);
}
void MPI_Indexed(int size, char** argv)
{
	int	n = 3;
	int* recBuffer = (int*)malloc(sizeof(int) * n);
	int* mat = nullptr;
	if (process_rank == 0)
	{
		MPI_Datatype newType;
		int blockLengths[3], disp[3];
		for (int i = 0; i < n; i++)
		{
			blockLengths[i] = n - i;
			disp[i] = n * i + i;
		}
		MPI_Type_indexed(n, blockLengths, disp, MPI_INT, &newType);
		MPI_Type_commit(&newType);
		mat = new int[n * n]
		{
		  3,4,14,
		  1,8,9,
		  2,12,5
		};
		MPI_Send(&mat[0], 1, newType, 1, 0, MPI_COMM_WORLD);
		MPI_Type_free(&newType);
	}
	else
	{
		int count = n * (n + 1) / 2;// upper // lower -> n*(n-1)
		MPI_Recv(recBuffer, count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int j = 0; j < count; j++)
			printf("%d %d\n", recBuffer[j], rank);
	}
}
void Transposition(double*& B, int size)
{
	double temp = 0.0;
	for (int i = 0; i < size; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			temp = B[i * size + j];
			B[i * size + j] = B[j * size + i];
			B[j * size + i] = temp;
		}
	}
}
void RandInit(double* pMatrix, int size)
{
	srand(10);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			pMatrix[i * size + j] = rand() % 10;
	}
}
double* A;
double* B;
double* C;

void PrintMatrix(double* pMatrix, int size)
{
	printf("\n////////////////////////////");
	for (int i = 0; i < size; i++)
	{
		printf("\n");
		for (int j = 0; j < size; j++)
			printf("%7.4f ", pMatrix[i * size + j]);
	}
}
void MatrixMultiplication(int size, int rank)
{
	int matrix_size = 4;
	double temp;
	MPI_Status Status;
	int ProcPartSize = matrix_size / size;
	int ProcPartElem = ProcPartSize * matrix_size;
	double* bufA = new double[ProcPartElem];
	double* bufB = new double[ProcPartElem];
	double* bufC = new double[ProcPartElem];
	if (rank == 0)
	{
		A = new double[matrix_size * matrix_size];
		B = new double[matrix_size * matrix_size];
		C = new double[matrix_size * matrix_size];
		RandInit(A, matrix_size);
		RandInit(B, matrix_size);
		PrintMatrix(A, matrix_size);
		PrintMatrix(B, matrix_size);
		Transposition(B, matrix_size);
	}
	MPI_Scatter(A, ProcPartElem, MPI_DOUBLE, bufA, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(B, ProcPartElem, MPI_DOUBLE, bufB, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	temp = 0.0;
	for (int i = 0; i < ProcPartSize; i++)
	{
		for (int j = 0; j < ProcPartSize; j++)
		{
			for (int k = 0; k < matrix_size; k++)
				temp += bufA[i * matrix_size + k] * bufB[j * matrix_size + k];
			bufC[i * matrix_size + j + ProcPartSize * rank] = temp;
			temp = 0.0;
		}
	}
	int NextProc; int PrevProc;
	for (int p = 1; p < size; p++)
	{
		NextProc = rank + 1;
		if (rank == size - 1) NextProc = 0;
		PrevProc = rank - 1;
		if (rank == 0) PrevProc = size - 1;
		MPI_Sendrecv_replace(bufB, ProcPartElem, MPI_DOUBLE, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);
		temp = 0.0;
		for (int i = 0; i < ProcPartSize; i++)
		{
			for (int j = 0; j < ProcPartSize; j++)
			{
				for (int k = 0; k < matrix_size; k++)
					temp += bufA[i * matrix_size + k] * bufB[j * matrix_size + k];
				int index = 0;
				if (rank - p >= 0)
					index = rank - p;
				else
					index = (size - p + rank);
				bufC[i * matrix_size + j + index * ProcPartSize] = temp;
				temp = 0.0;
			}
		}
	}
	MPI_Gather(bufC, ProcPartElem, MPI_DOUBLE, C, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0) PrintMatrix(C, matrix_size);
	delete[]bufA;
	delete[]bufB;
	delete[]bufC;
}

int main(int argc, char** argv)
{
	//PointToPoint(argc, argv);
	//CollectivSending(argc, argv);

	//Complex c1(1, 2);
	//Complex c2(3, 5);
	//Complex c = c1 + c2;
	//printf("%d %d", c.RealPart, c.InRealPart);
	//ComplexNumbers(argc, argv);
	//CollectiveGaterV(argc, argv);
	//CollectiveScatterV(argc, argv);

	//MpiPackAndUnpackExample(argc, argv);
	//MPI_CONTIGUOUS(argc, argv);
	MPI_Vector(argc, argv);
	MPI_Indexed(argc, argv);

}
