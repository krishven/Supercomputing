#include <iostream>
// #include <cilk.h>
// #include <cilkview.h>
using namespace std;

void merge(uint64_t arr[], uint64_t l, uint64_t m, uint64_t r)
{
    uint64_t i, j, k;
    uint64_t n1 = m - l + 1;
    uint64_t n2 =  r - m;
 
    uint64_t *L = new uint64_t[n1];
    uint64_t *R = new uint64_t[n2];
 
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
 
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
 
void mergeSort(uint64_t arr[], uint64_t l, uint64_t r)
{
    if (l < r)
    {
        uint64_t m = l+(r-l)/2;
 
        mergeSort(arr, l, m);
        mergeSort(arr, m+1, r);
 
        merge(arr, l, m, r);
    }
}

void printArray(uint64_t A[], uint64_t size)
{
    uint64_t i;
    for (i=0; i < size; i++)
        cout << A[i] << "\t";
    cout <<"\n";
}
 
int main()
{
    uint64_t n=8192*128;
    uint64_t limit=8192*8192;

    while (n <= limit) {
        
        uint64_t *arr = new uint64_t[n];
        
        for (uint64_t i = 0; i < n; i++)
            arr[i] =   1;

        //cout << "Given array is \n";
        //printArray(arr, n);
        // cilk::cilkview cv;
        // cv.start();
        mergeSort(arr, 0, n-1);
        // cv.stop();
        // cv.dump("results", false);
        // cout << "V: " << n << " " << cv.accumulated_milliseconds() / 1000.f << " seconds" << endl;
        cout << n << "\nSorted array is \n";
        // printArray(arr, n);
        delete[] arr;
        n = n * 2;
    }
    
    return 0;
}