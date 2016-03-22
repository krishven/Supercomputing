#include <iostream>
#include <cilk.h>
#include <cilkview.h>
using namespace std;

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

void iterFW (int **X, int length, int xi, int xj, int ui, int uj, int vi, int vj)
{
    for (int k = uj; k < uj + length; k++)
    {
        cilk_for (int i = xi; i < xi + length; i++)
        {
            cilk_for (int j = xj; j < xj + length; j++)
            {
                if (X[i][k] + X[k][j] < X[i][j])
                    X[i][j] = X[i][k] + X[k][j];
            }
        }
    }
}

void DFW (int **X, int length, int m, int xi, int xj, int ui, int uj, int vi, int vj) {
    if (length == m)
        iterFW(X, length, X11, U11, V11);
    else {
        int newlen = length / 2;
        cilk_spawn DFW (X, newlen, m, X11, U11, V11);
        cilk_spawn DFW (X, newlen, m, X12, U11, V12);
        cilk_spawn DFW (X, newlen, m, X21, U21, V11);
        DFW (X, newlen, m, X22, U21, V12);
        cilk_sync;
        cilk_spawn DFW (X, newlen, m, X11, U12, V21);
        cilk_spawn DFW (X, newlen, m, X12, U12, V22);
        cilk_spawn DFW (X, newlen, m, X21, U22, V21);
        DFW (X, newlen, m, X22, U22, V22);
        cilk_sync;
    }
}

void CFW (int **X, int length, int m, int xi, int xj, int ui, int uj, int vi, int vj) {
    if (length == m)
        iterFW(X, length, X11, U11, V11);
    else {
        int newlen = length / 2;
        cilk_spawn CFW (X, newlen, m, X11, U11, V11);
        CFW (X, newlen, m, X21, U21, V11);
        cilk_sync;
        cilk_spawn DFW (X, newlen, m, X12, U11, V12);
        DFW (X, newlen, m, X22, U21, V12);
        cilk_sync;
        cilk_spawn CFW (X, newlen, m, X12, U12, V22);
        CFW (X, newlen, m, X22, U22, V22);
        cilk_sync;
        cilk_spawn DFW (X, newlen, m, X11, U12, V21);
        DFW (X, newlen, m, X21, U22, V12);
        cilk_sync;

    }
}

void BFW (int **X, int length, int m, int xi, int xj, int ui, int uj, int vi, int vj) {
    if (length == m)
        iterFW(X, length, X11, U11, V11);
    else {
        int newlen = length / 2;
        cilk_spawn BFW (X, newlen, m, X11, U11, V11);
        BFW (X, newlen, m, X12, U11, V12);
        cilk_sync;
        cilk_spawn DFW (X, newlen, m, X21, U21, V11);
        DFW (X, newlen, m, X22, U21, V12);
        cilk_sync;
        cilk_spawn BFW (X, newlen, m, X21, U22, V21);
        BFW (X, newlen, m, X22, U22, V22);
        cilk_sync;
        cilk_spawn DFW (X, newlen, m, X11, U12, V21);
        DFW (X, newlen, m, X12, U12, V22);
        cilk_sync;
    }
}

void AFW (int **X, int length, int m, int xi, int xj, int ui, int uj, int vi, int vj) {
    if (length == m)
        iterFW(X, length, X11, U11, V11);
    else {
        int newlen = length / 2;
        AFW (X, newlen, m, X11, U11, V11);
        cilk_spawn BFW (X, newlen, m, X12, U11, V12);
        CFW (X, newlen, m, X21, U21, V11);
        cilk_sync;
        DFW (X, newlen, m, X22, U21, V12);
        AFW (X, newlen, m, X22, U22, V22);
        cilk_spawn BFW (X, newlen, m, X21, U22, V21);
        CFW (X, newlen, m, X12, U12, V22);
        cilk_sync;
        DFW (X, newlen, m, X11, U12, V21);
    }
}

// void printSolution(int **X)
// {
//     for (int i = 0; i < V; i++)
//     {
//         for (int j = 0; j < V; j++)
//         {
//             if (X[i][j] == NA)
//                 cout << "       " << "-";
//             else
//                 cout << "       "   << X[i][j];
//         }
//         cout << "\n";
//     }
// }

int cilk_main(int argc, char* argv[])
{
    int n = 4096;
    int m = n/2;

    while (m > 1) {
        int **X = new int*[n];
        for (int i = 0; i < n; ++i)
            X[i] = new int[n];

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i == j)
                    X[i][j] = 0;
                else
                    X[i][j] = i + j;

        cilk::cilkview cv;
        cv.start();
        AFW(X, n, m, 0, 0, 0, 0, 0, 0);
        cv.stop();
        cv.dump("results", false);
        cout << "M: " << m << " " << cv.accumulated_milliseconds() / 1000.f << " seconds" << endl;

        for ( int i = 0 ; i < n; i++ )
        {
            delete X[i];
        }
        delete[] X;
        m = m / 2;
    }

    return 0;
}
