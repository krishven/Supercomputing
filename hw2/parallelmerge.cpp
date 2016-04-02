#include <iostream>
//#include<stdio.h>
using namespace std;

#define M2 2
#define M3 2
#define N 2048

int A[N], T1[N];

void swap(int &a1, int &a2) {
	int t = a1;
	a1 = a2;
	a2 = t;
}

void serial_merge(int T[], int p1, int r1, int p2, int r2, int A[], int p3, int r3)
{
	int i = p1, j = p2, k = p3;

//	printf("%d:%d \t %d:%d \t %d:%d\n",p1,r1,p2,r2,p3,r3);
	while (i <= r1 && j <= r2) {
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
	while (i <= r1) {
		A[k] = T[i];
		k++;
		i++;
	}
	while (j <= r2) {
		A[k] = T[j];
		k++;
		j++;
	}
}

int binary_search(int k, int T[], int l, int r) {
	int m;

	while (l <= r)
	{
		m = l + (r - l) / 2;
		if (T[m] == k)
			return m;
		if (T[m] < k)
			l = m + 1;
		else
			r = m - 1;
	}
	return 0;
}

void sort(int arr[], int l, int r) {
	int i, key, j;
	for (i = l; i <= r; i++)
	{
		key = arr[i];
		j = i - 1;
		while (j >= 0 && arr[j] > key)
		{
			arr[j + 1] = arr[j];
			j = j - 1;
		}
		arr[j + 1] = key;
	}
}

void par_merge(int T[], int p1, int r1, int p2, int r2, int A[], int p3)
{
	int n1,n2,r3;
	
	if(p1<0 || r1 < 0 || p2 < 0 || r2 < 0 || p3 < 0)
		return;

	n1 = r1 - p1 + 1, n2 = r2 - p2 + 1, r3 = p3 + n1 + n2 - 1;

//	printf("%d:%d \t %d:%d \t %d:%d\n",p1,r1,p2,r2,p3,r3);
	
	if (n1 + n2 <= M2)
		serial_merge(T, p1, r1, p2, r2, A, p3, r3);
	else {
		if (n1 < n2) {
			swap(p1, p2);
			swap(r1, r2);
			swap(n1, n2);
		}
		int q1 = (p1 + r1) / 2;
		int q2 = binary_search(T[q1], T, p2, r2);
		int q3 = p3 + (q1 - p1) + (q2 - p2);
		
		if(q1<0||q2<0||q3<0)
			return;		
//		printf("q1 :%d \t q2:%d\t q3:%d\t",q1,q2,q3);

		A[q3] = T[q1];
		par_merge(T, p1, q1 - 1, p2, q2 - 1, A, p3);
		par_merge(T, q1 + 1, r1, q2 + 1, r2 - 1, A, q3);
	}
}

void par_merge_sort_pm(int A[], int p, int r)
{
	int i;
	int n = r - p + 1;

	if (n <= M3)
		sort(A, p, r);
	else {
		int q = (p + r) / 2;
		par_merge_sort_pm(A, p, q);
		par_merge_sort_pm(A, q + 1, r);
		for (i = p; i <= r; i++)
			T1[i] = A[i];
		par_merge(T1, p, q, q + 1, r, A, p);
	}
}

int main() {
	int n = N,i;

	for (i = 0; i < n; i++) {
		A[i] = n - i;
		T1[i] = n - i;
	}

	//cout << n << "\nArray is \n";
	for (i = 0; i < n; i++)
		cout << A[i] << "\t";
	cout << "\n";

	par_merge_sort_pm(A, 0, n - 1);

	
	for (i = 0; i < n; i++)
		cout << A[i] << "\t";
	cout << "\n";

	return 0;
}