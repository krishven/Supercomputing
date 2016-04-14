#include <iostream>
using namespace std;
#define D 7

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


void aLoop (int **X, int length, int xi, int xj, int ui, int uj, int vi, int vj)
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

void DFW (int **X,int length, int *tileSize, int d, int xi, int xj, int ui, int uj, int vi, int vj) {
	int r = tileSize[d];
	if (r >= length ) {
		aLoop(X, length, X11, U11, V11);
	} else {
		int newlen = length / r;

		for (int k=0;k<r;k++) {

			for (int i=0;i<r;i++) 
				for (int j=0;j<r;j++)
					DFW(X, newlen, tileSize, d+1, xi + i*newlen, xj + j*newlen, 
				ui + i*newlen, uj + k*newlen, vi + k*newlen, vj + j*newlen);		
		}
	}
}

void CFW (int **X,int length, int *tileSize, int d, int xi, int xj, int ui, int uj, int vi, int vj) {
	int r = tileSize[d];
	if (r >= length ) {
		aLoop(X, length, X11, U11, V11);
	} else {
		int newlen = length / r;

		for (int k=0;k<r;k++) {

			for (int i=0;i<r;i++) {
				CFW(X, newlen, tileSize, d+1, xi + i*newlen, xj + k*newlen, 
				ui + i*newlen, uj + k*newlen, vi + k*newlen, vj + k*newlen);
			}

			for (int i=0;i<r;i++) 
				for (int j=0;j<r;j++)
					DFW(X, newlen, tileSize, d+1, xi + i*newlen, xj + j*newlen, 
				ui + i*newlen, uj + k*newlen, vi + k*newlen, vj + j*newlen);		
		}
	}
}

void BFW (int **X,int length, int *tileSize, int d, int xi, int xj, int ui, int uj, int vi, int vj) {
	int r = tileSize[d];
	if (r >= length ) {
		aLoop(X, length, X11, U11, V11);
	} else {
		int newlen = length / r;

		for (int k=0;k<r;k++) {
			for (int j=0;j<r;j++) {
				BFW(X, newlen, tileSize, d+1, xi + k*newlen, xj + j*newlen, 
				ui + k*newlen, uj + k*newlen, vi + k*newlen, vj + j*newlen);
			}

			for (int i=0;i<r;i++) 
				for (int j=0;j<r;j++)
					DFW(X, newlen, tileSize, d+1, xi + i*newlen, xj + j*newlen, 
				ui + i*newlen, uj + k*newlen, vi + k*newlen, vj + j*newlen);		
		}
	}
}


void AFW (int **X,int length, int *tileSize, int d, int xi, int xj, int ui, int uj, int vi, int vj) {
	int r = tileSize[d];
	if (r >= length ) {
		aLoop(X, length, X11, U11, V11);
	} else {
		int newlen = length / r;

		for (int k=0;k<r;k++) {
			AFW(X,newlen, tileSize, d+1, xi + k*newlen, xj + k*newlen, 
				ui + k*newlen, uj + k*newlen, vi + k*newlen, vj + k*newlen);

			for (int j=0;j<r;j++) {
				BFW(X, newlen, tileSize, d+1, xi + k*newlen, xj + j*newlen, 
				ui + k*newlen, uj + k*newlen, vi + k*newlen, vj + j*newlen);
				int i = j;
				CFW(X, newlen, tileSize, d+1, xi + i*newlen, xj + k*newlen, 
				ui + i*newlen, uj + k*newlen, vi + k*newlen, vj + k*newlen);
			}

			for (int i=0;i<r;i++) 
				for (int j=0;j<r;j++)
					DFW(X, newlen, tileSize, d+1, xi + i*newlen, xj + j*newlen, 
				ui + i*newlen, uj + k*newlen, vi + k*newlen, vj + j*newlen);		
		}
	}
}

void printSolution(int **X,int length)
{
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            cout.width(8);
            if (X[i][j] == 100000)
                cout << "-";
            else
                cout << X[i][j];
        }
        cout << "\n";
    }
}

int main(int argc, char* argv[]) {
	int *tileSize = new int[D];
	
	for (int i=0;i<D;i++) {
		tileSize[i] = 1;
	}

	int n = 4;
	
	int **X = new int*[n];
    for (int i = 0; i < n; ++i)
        X[i] = new int[n];


    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                X[i][j] = 0;
            else
                X[i][j] = i + j;

    AFW(X, n, tileSize, 0, 0, 0, 0, 0, 0, 0);
    printSolution(X, n);
    for ( int i = 0 ; i < n; i++ )
    {
        delete X[i];
    }
    delete[] X;  
    delete[] tileSize;      
}