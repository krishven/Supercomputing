/* C Example */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 24


typedef struct
{
    int * array;
    int index;
}bucket_t;

int compare( const void * n1, const void * n2)
{
    return (*(int*)n1 - *(int*)n2);
}

static void scatterData(int *arr, int arrSize, int rank, int p) {

	int root;
	MPI_Comm comm;
	int *recvBuf,i;
	int q =3;
	int globalPivot[p-1];
	    
	if(rank == 0) {

		int i;
	    arr = (int *) malloc(sizeof(int)*arrSize);

	    srand(time(NULL));
	    
	    for(i = 0; i < arrSize; i++)
	        arr[i] = N-i;
	}
	recvBuf = malloc(sizeof(int)*arrSize/p);
	
	MPI_Scatter(arr, arrSize/p, MPI_INT, recvBuf, arrSize/p,
                    MPI_INT, 0, MPI_COMM_WORLD);

	
	qsort(recvBuf,arrSize/p,sizeof(int),compare);

    int factor = (arrSize/p)/q;
    int index = factor;
    int indexArr[(q-1)*p];
    int qbuf[q-1];
    i=0;
    while(index < arrSize/p) {
    	indexArr[i++] = recvBuf[index];
    	index += factor;
    }

    if(rank!=0){
    	MPI_Send(indexArr,q-1,MPI_INT,0,0,MPI_COMM_WORLD);
    }

    else {
    	int j=2;
//		printf("%dth node %d:%d\n",0,indexArr[0],indexArr[1]);
    	for(i=1;i<p;i++)
    	{
        	MPI_Recv(&indexArr[j],q-1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//    		printf("%dth node %d:%d\n",i,indexArr[j],indexArr[j+1]);
    		j=j+2;
    	}
	    qsort(indexArr,(q-1)*p,sizeof(int),compare);

/*	    for(i = 0; i < (q-1)*p; i++)
    	    printf("Array[%d] = %d\n",i,indexArr[i]);
*/
    	factor =  q-1;
    	index = factor;
    	i=0;
    	while(index < (q-1)*p) {
    		globalPivot[i++] = indexArr[index]; 
    		index += factor;
    	}

	  	for(i=1;i<p;i++) {
	  		MPI_Send(globalPivot,p-1,MPI_INT,i,0,MPI_COMM_WORLD);
	  	}
	}

    if(rank) {
    	int pivot[p-1];
    	MPI_Recv(&pivot,p-1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	//printf("Rank:%d pivot[0]:%d\n",rank,pivot[0]);
    }


/*
    bucket_t ** buckets = (bucket_t **) malloc(sizeof(bucket_t*)*p);

    for(i = 0; i < p; i++) {
        buckets[int] = (bucket_t *) malloc(sizeof(bucket_t));
        buckets[i]->array = (int*) malloc(sizeof(int)*arrSize/p);
        buckets[i]->index = 0;
    }
*/
} 
int main (int argc, char* argv[])
{
	int rank, size;
	int arr[N];
	int i;
	
	for (i=0;i<N;i++)
		arr[i]=N-i;

	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
	printf( "Hello world from process %d of %d\n", rank, size );
	
	scatterData(arr,N,rank,size);
	
	MPI_Finalize();
	return 0;
}
