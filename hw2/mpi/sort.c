/* C Example */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#define N 128


typedef struct
{
    int size;
    int array[0];
}bucket_t;

typedef struct
{
    int size;
    int x[2];
}bucket_s;

int compare( const void * n1, const void * n2)
{
    return (*(int*)n1 - *(int*)n2);
}

static int binary_search(int num,int *arr,int l, int r){

	int mid;
	

	if(num > arr[r-1])
		return r-1;

	while(l < r){
		mid = l + ((r-l)/2);
		if(arr[mid]<=num && arr[mid+1]>num)
			return mid;

		else if(arr[mid] > num){
			r = mid ;
		}
		else{
			l = mid ;
		}

	}

	return -1;

}
static void scatterData(int *arr, int arrSize, int rank, int p) {

	int root;
	MPI_Comm comm;
	int *recvBuf,i,*finalBuf;
	int q =3;
	int *recvCounts,*disp;
	int globalPivot[p-1];
	
	if(rank == 0) {

		int i;
	    arr = (int *) malloc(sizeof(int)*arrSize);
		finalBuf = (int *) malloc(sizeof(int)*arrSize);
		recvCounts = (int *) malloc(p*sizeof(int));
		disp = (int *) malloc(p*sizeof(int));
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
	  		MPI_Send(&globalPivot,p-1,MPI_INT,i,0,MPI_COMM_WORLD);
	  	}
	}

    if(rank) {
    	MPI_Recv(&globalPivot,p-1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	//printf("Rank:%d pivot[0]:%d\n",rank,pivot[0]);
    }


    bucket_t ** buckets = (bucket_t **) malloc(sizeof(bucket_t*)*p);

    for(i = 0; i < p; i++) {
    	if(i == rank) {
    		buckets[i] = (bucket_t *) malloc(sizeof(bucket_t)+sizeof(int)*arrSize);
    	}
    	else {
    		buckets[i] = (bucket_t *) malloc(sizeof(bucket_t)+sizeof(int)*arrSize/p);

    	}
        buckets[i]->size = 0;
    }

    int prev_index=0;
    int j;

/*    if(rank==0){
    	for(i=0;i<arrSize/p;i++)
    		printf("sdsa%d\n",recvBuf[i]);
    }
*/
    for(i = 0 ; i < p-1; i++) {
    	//printf("b4: Rank :%d i:%d pivot:%d\n",rank,i,globalPivot[i]);
    	index = binary_search(globalPivot[i],recvBuf,prev_index, arrSize/p);
    	//printf("Rank :%d index val:%d\n",rank,index);
    	for(j = prev_index; j<=index; j++){
    		buckets[i]->array[buckets[i]->size] = recvBuf[j];
    		buckets[i]->size++;
    	}
    	if(index!=-1)
    		prev_index = index+1;
    }

	for(j = prev_index; j<arrSize/p; j++){
		buckets[p-1]->array[buckets[i]->size] = recvBuf[j];
		buckets[p-1]->size++;
	}

    const int nitems=2;
    // int blocklengths[2] = {1,arrSize/p};
    // MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype mpi_bucket_type;
    // MPI_Aint     offsets[2];

    // offsets[0] = offsetof(bucket_t, size);
    // offsets[1] = offsetof(bucket_t, array);
    MPI_Type_contiguous(1+(arrSize/p),MPI_INT,&mpi_bucket_type);
    // MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_bucket_type);
    MPI_Type_commit(&mpi_bucket_type);

	printf("sizeof bucket:%d:%d\n",sizeof(bucket_t)+sizeof(int)*arrSize/p,arrSize/p);
	#if 1
    for(j=0;j<p;j++)
	{
		if(j!=rank) {
			//MPI_Send(buckets[j]->array,arrSize/p,MPI_INT,j,0,MPI_COMM_WORLD);
			MPI_Send(buckets[j],1,mpi_bucket_type,j,0,MPI_COMM_WORLD);
			//MPI_Send(val,1,mpi_bucket_type,j,0,MPI_COMM_WORLD);
		}
	}


	bucket_t *tempBucket = malloc(sizeof(bucket_t)+sizeof(int)*arrSize/p);
	for(j=0;j<p;j++)
	{
		if(j!=rank) {
			//MPI_Recv(&tempArray[pos],arrSize/p,MPI_INT,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(tempBucket,1,mpi_bucket_type,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			memcpy(&(buckets[rank]->array[buckets[rank]->size]),&tempBucket->array[0],sizeof(int)*
				tempBucket->size);
			buckets[rank]->size += tempBucket->size;
		}
	}
	

    qsort(buckets[rank]->array,buckets[rank]->size,sizeof(int),compare);

/*
	printf("rank:%d\t",rank);
	for(j=0;j<buckets[rank]->size;j++){
		printf("%d\t",buckets[rank]->array[j]);
	}
	
	printf("\n");
*/
	MPI_Gather(&buckets[rank]->size, 1, MPI_INT, recvCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank==0){
		disp[0] = 0;
  		for(i=1;i<p;i++){
			disp[i] = disp[i-1] + recvCounts[i-1];
		}
	}

	// if(rank==0){
	// 	for(i=0;i<p;i++){
	// 		printf("%d\t",recvCounts[i]);
	// 	}
	// }
	MPI_Gatherv(buckets[rank]->array, buckets[rank]->size, MPI_INT,
            finalBuf, recvCounts, disp, MPI_INT, 0, MPI_COMM_WORLD);
	//printf("buckets rank size:%d:%d\n",buckets[rank]->size,arrSize);
	//MPI_Gather( buckets[rank]->array, buckets[rank]->size, MPI_INT, rbuf, buckets[rank]->size, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank==0) {
		for(i=0;i<arrSize;i++){
			printf("%d\t",finalBuf[i]);
		}
		printf("\n");
	}


	/*
	printf("Rank %d tmpArray\t", rank);
	for(i=0;i<arrSize;i++) {
		printf("%d\t", tempArray[i]);
	}

	printf("\n");
	*/
	
	// MPI_Type_free(&mpi_bucket_type);
 //    MPI_Finalize();
	
/*	for(i=0;i<tempBucket.size;i++){
		printf("%d\t", tempBucket.array[i]);
	}
	*/
	// printf("\n");


	
 //    for(i = 0 ; i < p; i++) {
 //    	printf("rank :%d bucket :%d\t",rank,i);
 //    	for(j=0;j<buckets[i]->size;j++){
 //    		printf("%d\t",buckets[i]->array[j]);
 //    	}
 //    	printf("\n");
	// }
	#endif
}

int main (int argc, char* argv[])
{
	int rank, size;
	int *arr;
	int i;
	/*
	for (i=0;i<N;i++)
		arr[i]=N-i;
*/
	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
	printf( "Hello world from process %d of %d\n", rank, size );
	
	scatterData(arr,N,rank,size);
	
	MPI_Finalize();
	return 0;
}
