#include <mpi.h>
#define N 8192*32
extern "C++" void qsort(int a[], int n);
int main( int argc, char *argv[ ])
{
	int arr[N];
	for(int i=0;i<N;i++)
		arr[i] = N-i;
	MPI_Init( &argc, &argv );
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	qsort(arr,N);
	for(int i=0;i<N;i++)
		printf("%d ",arr[i]);
	MPI_Finalize( );
	return 0;
}
