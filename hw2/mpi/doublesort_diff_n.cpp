#include <mpi.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <iostream>
#include <ctime>

#define N 67108864
using namespace std;

extern void qsort(double *a, int n);

typedef struct
{
	int size;
	double array[0];
} bucket_t;

int compare( const void * n1, const void * n2)
{
	return (*(double*)n1 - * (double*)n2);
}

static int binary_search(double num, double *arr, int l, int r) {

	int mid;
	if (num > arr[r - 1])
		return r - 1;

	while (l < r) {
		mid = l + ((r - l) / 2);
		if (arr[mid] <= num && arr[mid + 1] > num)
			return mid;

		else if (arr[mid] > num) {
			r = mid ;
		}
		else {
			l = mid ;
		}

	}

	return -1;

}
static void scatterData(double *arr, uint64_t arrSize, int rank, int p) {

	int root;
	MPI_Comm comm;
	double *recvBuf, *finalBuf;
	uint64_t i;
	int q = 3;

	int *recvCounts, *disp, *recvNum;
	double globalPivot[p - 1];
	int sendCount[p];

	if (rank == 0) {

		int i;
		arr = new double[arrSize];
		finalBuf = new double[arrSize];
		recvCounts = new int[p];
		recvNum = new int[p * p];
		disp = new int[p];

		for (i = 0; i < arrSize; i++)
			arr[i] = N - i;

		for (uint64_t i = 0; i < arrSize; i++) {
            uint64_t x = rand()%arrSize;
            uint64_t y = rand()%arrSize;
            double temp = arr[x];
            arr[x] = arr[y];
            arr[y] = temp;
        }
	
	}

	int lastRankSize = (arrSize / p) + (arrSize % p);

	int tmpSize = (rank == p - 1) ? lastRankSize : (arrSize / p);

	if (arrSize % p == 0)
		recvBuf = new double[arrSize / p];
	else {
		recvBuf = new double[lastRankSize];
	}
	sendCount[rank] = tmpSize;


	MPI_Gather(&sendCount[rank], 1, MPI_INT, recvCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		disp[0] = 0;
		for (i = 1; i < p; i++) {
			disp[i] = disp[i - 1] + recvCounts[i - 1];

		}
	}
	MPI_Scatterv( arr, recvCounts, disp, MPI_DOUBLE, recvBuf, (arrSize / p) + (arrSize % p), MPI_DOUBLE,
	              0, MPI_COMM_WORLD);

	qsort(recvBuf, tmpSize, sizeof(double), compare);

	int factor = tmpSize / q;

	int index = factor;
	double indexArr[(q - 1)*p];

	i = 0;
	while (index < tmpSize) {
		indexArr[i++] = recvBuf[index];
		index += factor;
	}


	if (rank != 0) {
		MPI_Send(indexArr, q - 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	else {
		int j = q - 1;

		for (i = 1; i < p; i++)
		{
			MPI_Recv(&indexArr[j], q - 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			j = j + q - 1;
		}
		qsort(indexArr, (q - 1)*p, sizeof(double), compare);

		factor =  q - 1;
		index = factor;
		i = 0;
		while (index < (q - 1)*p) {
			globalPivot[i++] = indexArr[index];
			index += factor;
		}

		for (i = 1; i < p; i++) {
			MPI_Send(&globalPivot, p - 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}

	if (rank) {
		MPI_Recv(&globalPivot, p - 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}

	bucket_t ** buckets = (bucket_t **) malloc(sizeof(bucket_t*)*p);

	for (i = 0; i < p; i++) {
		if (i == rank) {
			buckets[i] = (bucket_t *) malloc(sizeof(bucket_t) + sizeof(double) * arrSize);

		}
		else {

			buckets[i] = (bucket_t *) malloc(sizeof(bucket_t) + sizeof(double) * tmpSize);

		}
		buckets[i]->size = 0;
	}

	int prev_index = 0;
	int j;

	for (i = 0 ; i < p - 1; i++) {


		index = binary_search(globalPivot[i], recvBuf, prev_index, tmpSize);

		for (j = prev_index; j <= index; j++) {
			buckets[i]->array[(int)buckets[i]->size] = recvBuf[j];
			buckets[i]->size++;
		}
		if (index != -1)
			prev_index = index + 1;
	}

	for (j = prev_index; j < tmpSize; j++) {
		buckets[p - 1]->array[(int)buckets[i]->size] = recvBuf[j];
		buckets[p - 1]->size++;
	}

	const int nitems = 2;
	int blocklengths[2] = {1, lastRankSize};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Datatype mpi_bucket_type;
	MPI_Aint     offsets[2];

	offsets[0] = offsetof(bucket_t, size);
	offsets[1] = offsetof(bucket_t, array);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_bucket_type);
	MPI_Type_commit(&mpi_bucket_type);

	MPI_Request * requests = (MPI_Request *) malloc(sizeof(MPI_Request) * p);
	for (j = 0; j < p; j++)
	{
		if (j != rank) {
			MPI_Isend(buckets[j], 1, mpi_bucket_type, j, 0, MPI_COMM_WORLD, &requests[j]);


		}
	}
	bucket_t *tempBucket = (bucket_t *)malloc(sizeof(bucket_t) + sizeof(double) * lastRankSize);

	for (j = 0; j < p; j++)
	{
		if (j != rank) {

			MPI_Recv(tempBucket, 1, mpi_bucket_type, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);



			memcpy(&(buckets[rank]->array[(int)buckets[rank]->size]), &tempBucket->array[0],
			       sizeof(double)*tempBucket->size);
			buckets[rank]->size += tempBucket->size;


		}
	}

	qsort(buckets[rank]->array, buckets[rank]->size, sizeof(double), compare);

#if 1

	MPI_Gather(&buckets[rank]->size, 1, MPI_INT, recvCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		disp[0] = 0;
		for (i = 1; i < p; i++) {
			disp[i] = disp[i - 1] + recvCounts[i - 1];
		}
	}

	MPI_Gatherv(buckets[rank]->array, buckets[rank]->size, MPI_DOUBLE,
	            finalBuf, recvCounts, disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		// for (i = 0; i < arrSize; i++)
		// 	cout << finalBuf[i] << "\n";

		// cout << endl;
		delete[] arr;
		delete[] finalBuf;
	}

	delete[] recvBuf;
	

	for (i = 0; i < p; i++) {
        delete buckets[i];
	}
	delete[] buckets;
	delete tempBucket;
#endif

}

int main (int argc, char* argv[])
{
	int rank, size;
	double *arr;
	uint64_t n = 536870912;
	
	MPI_Init (&argc, &argv);			
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	printf( "Hello world from process %d of %d\n", rank, size );	
	uint64_t limit = 536870912;
    uint64_t diff = 53686988;
    
  //  while(n<=limit) {
		const clock_t begin_time = clock();

		scatterData(arr, n, rank, size);

		cout << n <<"\t" << float(clock () - begin_time)/ CLOCKS_PER_SEC<<endl;
	    n = n + diff;
    //}
		
	MPI_Finalize();	
	return 0;
}