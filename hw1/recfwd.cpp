#include <iostream>
#include <omp.h>
using namespace std;

#define V 4
#define M 2
#define NA 100000

#define X11 xi, xj
#define X12 xi, xj + newlen
#define X21 xi + newlen, xj
#define X22 xi + newlen, xj + newlen

#define U11 ui, uj
#define U12 ui, uj + newlen
#define U21 ui + newlen, uj
#define U22 ui + newlen, uj + newlen

#define V11 vi, vj
#define V12 vi, vj + newlen
#define V21 vi + newlen, vj
#define V22 vi + newlen, vj + newlen

#define PARALLEL 1

void iterFW (int X[][V], int length, int xi, int xj, int ui, int uj, int vi, int vj)
{
	for (int k = uj; k < uj + length; k++)
	{
		for (int i = xi; i < xi + length; i++)
		{
			for (int j = xj; j < xj + length; j++)
			{
				if (X[i][k] + X[k][j] < X[i][j])
					X[i][j] = X[i][k] + X[k][j];
			}
		}
	}
}

void DFW (int X[][V], int length, int xi, int xj, int ui, int uj, int vi, int vj) {
	if (length == M)
		iterFW(X, length, X11, U11, V11);
	else {
		int newlen = length / 2;
#ifdef PARALLEL
		#pragma omp parallel num_threads(4)
		{
			if (omp_get_thread_num()==0)	
				DFW(X, newlen, X11, U11, V11);
			else if (omp_get_thread_num()==1)	
				DFW(X, newlen, X12, U11, V12);
			else if (omp_get_thread_num()==2)	
				DFW(X, newlen, X21, U21, V11);
			else 	
				DFW(X, newlen, X22, U21, V12);

		}
		#pragma omp barrier

		#pragma omp parallel num_threads(4)
		{
			if (omp_get_thread_num()==0)	
				DFW(X, newlen, X11, U12, V21);
			else if (omp_get_thread_num()==1)	
				DFW(X, newlen, X12, U12, V22);
			else if (omp_get_thread_num()==2)	
				DFW(X, newlen, X21, U22, V21);
			else 	
				DFW(X, newlen, X22, U22, V22);
		}
		#pragma omp barrier
#else
		DFW (X, newlen, X11, U11, V11);
		DFW (X, newlen, X12, U11, V12);
		DFW (X, newlen, X21, U21, V11);
		DFW (X, newlen, X22, U21, V12);
		
		DFW (X, newlen, X11, U12, V21);
		DFW (X, newlen, X12, U12, V22);
		DFW (X, newlen, X21, U22, V21);
		DFW (X, newlen, X22, U22, V22);

#endif
	}
}

void CFW (int X[][V], int length, int xi, int xj, int ui, int uj, int vi, int vj) {
	if (length == M)
		iterFW(X, length, X11, U11, V11);
	else {
		int newlen = length / 2;
#ifdef PARALLEL
		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				CFW(X, newlen, X11, U11, V11);
			else
				CFW(X, newlen, X21, U21, V11);
		}
		
		#pragma omp barrier

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				DFW(X, newlen, X12, U11, V12);
			else
				DFW(X, newlen, X22, U21, V12);
		}		
		#pragma omp barrier

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				CFW(X, newlen, X12, U12, V22);
			else
				CFW(X, newlen, X22, U22, V22);
		}
		
		#pragma omp barrier

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				DFW(X, newlen, X11, U12, V21);
			else
				DFW(X, newlen, X21, U22, V12);
		}
		
		#pragma omp barrier
#else
		CFW (X, newlen, X11, U11, V11);
		CFW (X, newlen, X21, U21, V11);
		
		DFW (X, newlen, X12, U11, V12);
		DFW (X, newlen, X22, U21, V12);
		
		CFW (X, newlen, X12, U12, V22);
		CFW (X, newlen, X22, U22, V22);
		
		DFW (X, newlen, X11, U12, V21);
		DFW (X, newlen, X21, U22, V12);

#endif
	}
}

void BFW (int X[][V], int length, int xi, int xj, int ui, int uj, int vi, int vj) {
	if (length == M)
		iterFW(X, length, X11, U11, V11);
	else {
		int newlen = length / 2;
#if PARALLEL

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				BFW(X, newlen, X11, U11, V11);
			else
				BFW(X, newlen, X12, U11, V12);
		}
		
		#pragma omp barrier

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				DFW(X, newlen, X21, U21, V11);
			else
				DFW(X, newlen, X22, U21, V12);
		}
		
		#pragma omp barrier

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				BFW(X, newlen, X21, U22, V21);
			else
				BFW(X, newlen, X22, U22, V22);
		}
		
		#pragma omp barrier

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				DFW(X, newlen, X11, U12, V21);
			else
				DFW(X, newlen, X12, U12, V22);
		}
		
		#pragma omp barrier
#else
		BFW (X, newlen, X11, U11, V11);
		BFW (X, newlen, X12, U11, V12);
		
		DFW (X, newlen, X21, U21, V11);
		DFW (X, newlen, X22, U21, V12);
		
		BFW (X, newlen, X21, U22, V21);
		BFW (X, newlen, X22, U22, V22);
		
		DFW (X, newlen, X11, U12, V21);
		DFW (X, newlen, X12, U12, V22);
#endif
	}
}

void AFW (int X[][V], int length, int xi, int xj, int ui, int uj, int vi, int vj) {
	if (length == M)
		iterFW(X, length, X11, U11, V11);
	else {
		int newlen = length / 2;
#if PARALLEL
		AFW (X, newlen, X11, U11, V11);

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				BFW(X, newlen, X12, U11, V12);
			else
				CFW(X, newlen, X21, U21, V11);
		}
		
		#pragma omp barrier

		DFW (X, newlen, X22, U21, V12);
		AFW (X, newlen, X22, U22, V22);

		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num()!=0)	
				BFW(X, newlen, X21, U22, V21);
			else
				CFW(X, newlen, X12, U12, V22);
		}
		
		#pragma omp barrier

		DFW (X, newlen, X11, U12, V21);
#else
		AFW (X, newlen, X11, U11, V11);
		BFW (X, newlen, X12, U11, V12);
		CFW (X, newlen, X21, U21, V11);
		
		DFW (X, newlen, X22, U21, V12);
		AFW (X, newlen, X22, U22, V22);
		BFW (X, newlen, X21, U22, V21);
		CFW (X, newlen, X12, U12, V22);
		
		DFW (X, newlen, X11, U12, V21);

#endif
	}
}

void printSolution(int X[][V])
{
	for (int i = 0; i < V; i++)
	{
		for (int j = 0; j < V; j++)
		{
			if (X[i][j] == NA)
				cout << "       " << "-";
			else
				cout << "       "   << X[i][j];
		}
		cout << "\n";
	}
}

int main(int argc, char* argv[])
{
	int X[V][V] = { {0,   5,  NA, 10},
		{NA, 0,   3, NA},
		{NA, NA, 0,   1},
		{NA, NA, NA, 0}
	};

	// Print the solution
	AFW(X, V, 0, 0, 0, 0, 0, 0);
	printSolution(X);
	return 0;
}
