#include <iostream>
#include <cilk/cilk.h>
#include <stdint.h>
#include <algorithm>
#include<stdio.h>
using namespace std;

uint64_t M2 = 1024;
uint64_t M3 = 16;

void swap(uint64_t &a1, uint64_t &a2) {
    uint64_t t = a1;
    a1 = a2;
    a2 = t;
}

void serial_merge(double T[], uint64_t p1, uint64_t r1, uint64_t p2, uint64_t r2, double A[], uint64_t p3, uint64_t r3)
{
    uint64_t i = p1, j = p2, k = p3;

    while (i <= r1 && j <= r2 && k <= r3) {
        if (T[i] <= T[j]) {
            A[k] = T[i];
            i++;
        }
        else {
            A[k] = T[j];
            j++;
        }
        k++;
    }
    while (i <= r1 && k <= r3) {
        A[k] = T[i];
        k++;
        i++;
    }
    while (j <= r2 && k <= r3) {
        A[k] = T[j];
        k++;
        j++;
    }
}

uint64_t binary_search(double k, double T[], uint64_t l, uint64_t r) {
    uint64_t m;
    r = max(l, r + 1);

    while (l < r)
    {
        m = l + (r - l) / 2;
        if (T[m] <= k && T[m + 1] > k)
            return m + 1;
        if (T[m] < k)
            l = m + 1;
        else
            r = m;
    }
    return r;
}

void sort(double arr[], uint64_t l, uint64_t r) {
    uint64_t i, key, j;
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

void par_merge(double T[], uint64_t p1, uint64_t r1, uint64_t p2, uint64_t r2, double A[], uint64_t p3)
{
    uint64_t n1, n2, r3;


    n1 = r1 - p1 + 1, n2 = r2 - p2 + 1, r3 = p3 + n1 + n2 - 1;

    if (n1 + n2 <= M2)
        serial_merge(T, p1, r1, p2, r2, A, p3, r3);
    else {
        if (n1 < n2) {
            swap(p1, p2);
            swap(r1, r2);
            swap(n1, n2);
        }
        uint64_t q1 = (p1 + r1) / 2;
        uint64_t q2 = binary_search(T[q1], T, p2, r2);
        uint64_t q3 = p3 + (q1 - p1) + (q2 - p2);

        A[q3] = T[q1];

        cilk_spawn par_merge(T, p1, q1 - 1, p2, q2 - 1, A, p3);
//      par_merge(T, p1, q1 - 1, p2, q2 - 1, A, p3);
        par_merge(T, q1 + 1, r1, q2 , r2, A, q3 + 1);
        cilk_sync;
    }
}

void par_merge_sort_pm(double A[], double T[], uint64_t p, uint64_t r)
{
    uint64_t i, size, q;
    uint64_t n = r - p + 1;

    if (n <= M3)
        sort(A, p, r);
    else {
        q = (p + r) / 2;

        cilk_spawn par_merge_sort_pm(A, T, p, q);
//      par_merge_sort_pm(A, p, q);
        par_merge_sort_pm(A, T, q + 1, r);
        cilk_sync;

        for (i = p; i <= r; i++)
            T[i] = A[i];

        par_merge(T, p, q, q + 1, r, A, p);
    }
}

void cilksort(double *A, int n) {
    double *T = new double[n];
//  cout << "gooooo";
    par_merge_sort_pm(A, T, 0, n - 1);

    return;
}