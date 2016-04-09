#include <iostream>
#include <cilk/cilk.h>
#include <math.h>
#include <limits.h>

using namespace std;

void merge(double arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 =  r - m;

    int *L = new int[n1];
    int *R = new int[n2];

    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
    delete[] L;
    delete[] R;
}
void sort(double arr[], int l, int r) {
    int i, j;
    int key;
    for (i = l + 1; i <= r; i++)
    {
        key = arr[i];
        j = i - 1;
        while (j >= l && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void mergeSort(double arr[], int l, int r, int n, int m)
{
    if (n <= m) {
        sort(arr, l, r);
    }
    else
    {
        int mid = l + (r - l) / 2;

        mergeSort(arr, l, mid, mid - l + 1, m);
        mergeSort(arr, mid + 1, r, r - mid, m);
        //cilk_sync;

        merge(arr, l, mid, r);
    }
}

void cilksort(double *a, int size) {
    int m = 128;
    //printf("Before Sorting:\n");
    //printArray(a, size);
    mergeSort(a, 0, size - 1, size, m);
    //cout << "Serail sort done" <<endl;
    //printf("After Sorting:\n");
    //printArray(a, size);
    return;
}