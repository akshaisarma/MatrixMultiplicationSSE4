/*
MATRIX MATRIX Multiplication using SSE4
Restrictions : Only for Square Matrices of any positve sizes which are integer multiples of 4
Instructions : Change N defined below: Every variable is calculated off N
*/


#include <iostream>
#include <tchar.h>
#include <windows.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

const int N = 1020;
float A[N][N];
float B[N][N];
float C[N][N];
float SC[N][N];


void main()
{
	srand ((unsigned int)time(NULL));
	//intialize matrices
	int i, j, k;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			A[i][j] = (float) rand()/ (float) 13101;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			B[i][j] = (float) rand()/ (float) 27010;

	//Query Performance
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	double t;
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);

	if (QueryPerformanceCounter((LARGE_INTEGER *)&ctr1)!= 0)
	{
		// Matrix Matrix Multiplication Standard
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
			{
				C[i][j] = 0;
				for (k =0; k < N; k++)
					C[i][j] += A[i][k] * B [k][j];
			}
		// Finish timing the code.
		QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
		t = double ((ctr2 - ctr1)* 1000000)/ freq;
		cout << "Matrix Matrix (Square) Standard C Multiplication of a "<<N<<" by "<<N<<
				" matrix with \n" <<N*N<< " entries each took : "<<t<< " us.\n"<<endl;
	}
	ctr1 = 0; ctr2 = 0;
	if (QueryPerformanceCounter((LARGE_INTEGER *)&ctr1)!= 0)
	{
		//Convert B to Column Major Format or take its Transpose making it easier to multiply. Adds to total time so in Performance.
		float temp;
		for (i = 0; i < N; i++)
			for (j = i+1; j < N; j++)
			{
				temp = B[i][j];
				B[i][j] = B[j][i];
				B[j][i] = temp;
			}
		// Matrix Matrix Multiplication SSE
		float * pA = *A;
		float * pB = *B;
		float * pC = *SC;
		int NChunk = N/4;
		int NRow   = N*4;
		__asm
		{
			mov esi, dword ptr [pA]
			mov edi, dword ptr [pC]

			mov ecx, N
			L3:
				push ecx
				mov ecx, N
				mov ebx, dword ptr [pB]			// We will be  going through every column of B (NChunk * 16 = 1 Column)
				L2:
					subss xmm1,xmm1				// clear xmm2 - Accumulator of Dot Product chunks
					push esi					// save esi
					push ecx
					mov ecx, NChunk
					L1:
						movups xmm0, [esi]
						dpps xmm0, [ebx], 0xFF	// Dot Product Instruction SSE4
						addss xmm1, xmm0
						add ebx, 16
						add esi, 16
					Loop L1
					movss dword ptr[edi], xmm1  // We will have one entry in result
					add edi, 4
					pop ecx
					pop esi						// restore esi because still in the same row in A
				Loop L2							// We will have one row in result
				add esi, NRow
				pop ecx
			Loop L3								// We will have N rows in result (== done)
		}
		// Finish timing the code.
		QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
		t = double ((ctr2 - ctr1)* 1000000)/ freq;
		cout << "Matrix Matrix (Square) SSE Multiplication of a "<<N<<" by "<<N<<
				" matrix with \n" <<N*N<< " entries each took : "<<t<< " us.\n"<<endl;

		bool check = TRUE;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				if (fabs(C[i][j] - SC[i][j]) > 0.01) //Equality won't work well with floats
				{
					check = FALSE;
					cout<< C[i][j]<< " "<<SC[i][j]<<endl;
					break;
				}
		if (check)
			cout << "\nThe two computed matrices are equal." <<endl<<endl;
		else
			cout << "\nThe two matrices are not the same."<<endl;
	}
		/*
		//Remove to see results for reasonable sizes
		cout << "Standard C" <<"            "<< "SSE4 "<<endl<<endl;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				cout << C[i][j] <<"            "<< SC[i][j]<<endl;
		*/
}
