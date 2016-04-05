
#include <mpi.h>
extern "C++" int swap_cpp( int a[], int n);
int main( int argc, char *argv[ ])
{
int arr[4]={0,1,2,3};
int n=5,r=10;
MPI_Init( &argc, &argv );
int rank;
MPI_Comm_rank( MPI_COMM_WORLD, &rank );
//printf("n=%d,r=%d\n",n,r);
swap_cpp(arr,4);
//printf("n=%d,r=%d\n",n,r);
//printf( "C( %d, %d ) = %d\n", 30, 15 + rank,
//nCr_CPP( 30, 15 + rank ) );
//nCr_CPP( arr, 15 + rank ) );
MPI_Finalize( );
return 0;
}
