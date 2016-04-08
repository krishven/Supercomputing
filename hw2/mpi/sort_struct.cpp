/* C Example */
#include <mpi.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <iostream>
#define N 8192
using namespace std;

typedef struct 
{
    int size;
    double array[0];
}bucket_t;

// typedef struct

// {
//     int size;
//     int x[2];
// }bucket_s;

int compare( const void * n1, const void * n2)
{

    return (*(double*)n1 - *(double*)n2);
}

static int binary_search(double num,double *arr,int l, int r){

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
static void scatterData(double *arr, uint64_t arrSize, int rank, int p) {

	int root;
	MPI_Comm comm;
	double *recvBuf,*finalBuf;
	uint64_t i;
	int q = 3;
	//double *recvCounts,*disp;
	int *recvCounts,*disp,*recvNum;
	double globalPivot[p-1];
	int sendCount[p];
	
	if(rank == 0) {

		int i;
	 //    arr = (double *) malloc(sizeof(double)*arrSize);
		// finalBuf = (double *) malloc(sizeof(double)*arrSize);
		// recvCounts = malloc(p*sizeof(int));
		// disp = malloc(p*sizeof(int));

	    //srand(time(NULL));
	    
	    arr = new double[arrSize];
	    finalBuf = new double[arrSize];
	    recvCounts = new int[p];
	    recvNum = new int[p*p];
	    disp = new int[p];

	    for(i = 0; i < arrSize; i++)
	        arr[i] = N-i;
	}

    // for(i = 0; i < arrSize/p; i++){
    // 	//printf("Array[%d]: %lf\n",i,arr[i]);
    // 	cout << "Array" << i <<":"<<arr[i]<<endl;
    // }
	//recvBuf = malloc(sizeof(double)*arrSize/p);
	int lastRankSize = (arrSize/p)+(arrSize%p);
	
	int tmpSize = (rank == p-1)?lastRankSize:(arrSize/p);
	
	if(arrSize%p==0)
		recvBuf = new double[arrSize/p];
	else{	
		recvBuf = new double[lastRankSize];
	}

	


	sendCount[rank] = tmpSize;

//	cout <<"rank" << rank <<"\t" << sendCount[rank] << "\t";

	// recvBuf = new double[arrSize/p];
	// MPI_Scatter(arr, arrSize/p, MPI_DOUBLE, recvBuf, arrSize/p,
 //                     MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Gather(&sendCount[rank], 1, MPI_INT, recvCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank==0){
		disp[0] = 0;
  		for(i=1;i<p;i++){
			disp[i] = disp[i-1] + recvCounts[i-1];
			//cout << disp[i] << "\t";
		}
	}



	MPI_Scatterv( arr, recvCounts, disp, MPI_DOUBLE, recvBuf, (arrSize/p) + (arrSize%p), MPI_DOUBLE, 
                                                               0, MPI_COMM_WORLD);
	// cout <<"rank" << rank <<"\t";
	// if(rank == p-1){
	// 	for(i=0;i<lastRankSize;i++){
	// 		cout << recvBuf[i] <<"\t";
	// 	}
	// }
	// else {
	// 	for(i=0;i<(arrSize/p);i++){
	// 		cout << recvBuf[i] <<"\t";
	// 	}		
	// }
	// cout <<endl;


	qsort(recvBuf,tmpSize,sizeof(double),compare);


	// cout << "Rank" <<rank;
 //    for(i = 0; i < tmpSize; i++)
 //    	cout << "Array" << i <<":"<<recvBuf[i]<<endl;	    
 //    cout<<endl;
//#if 1


	int factor = tmpSize/q;
	
	// if(rank == p-1)?factor
	// 	int factor = (lastRankSize)/q;
	// }
	// else
	// 	int factor = (arrSize/p)/q;

    int index = factor;
    double indexArr[(q-1)*p];
    //int qbuf[q-1];
    i=0;
    while(index < tmpSize) {
    	indexArr[i++] = recvBuf[index];
    	index += factor;
    }

 //    cout <<"Rank" << rank <<"\t";
 //    for(i = 0; i < (q-1); i++)
	//     cout << indexArr[i] << "\t";
	// cout <<endl;


    if(rank!=0){
    	MPI_Send(indexArr,q-1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }

    else {
    	int j=q-1;
//		printf("%dth node %d:%d\n",0,indexArr[0],indexArr[1]);
    	for(i=1;i<p;i++)
    	{
        	MPI_Recv(&indexArr[j],q-1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//    		printf("%dth node %d:%d\n",i,indexArr[j],indexArr[j+1]);
    		j=j+q-1;
    	}


	    qsort(indexArr,(q-1)*p,sizeof(double),compare);


		// cout <<"Rank" << rank <<"\t";
		// for(i = 0; i < (q-1)*p; i++)
		//     cout << indexArr[i] << "\t";
		// cout <<endl;

    	factor =  q-1;
    	index = factor;
    	i=0;
    	while(index < (q-1)*p) {
    		globalPivot[i++] = indexArr[index]; 
    		index += factor;
    	}

		for(i=1;i<p;i++) {
	  		MPI_Send(&globalPivot,p-1,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
	  	}
	}

    if(rank) {
    	MPI_Recv(&globalPivot,p-1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	// printf("Rank:%d pivot[0]:%d\n",rank,globalPivot[0]);
    }

	// cout <<"Rank" << rank <<"\t";
	// for(i = 0; i < p-1; i++)
	//     cout << globalPivot[i] << "\t";
	// cout <<endl;



    bucket_t ** buckets = (bucket_t **) malloc(sizeof(bucket_t*)*p);

    for(i = 0; i < p; i++) {
    	if(i == rank) {
    		buckets[i] = (bucket_t *) malloc(sizeof(bucket_t)+sizeof(double)*arrSize);
    		//memset(&buckets[i]->array,-1,sizeof(double)*arrSize);
    	}
    	else {

    		buckets[i] = (bucket_t *) malloc(sizeof(bucket_t)+sizeof(double)*tmpSize);
    		//memset(&buckets[i]->array,-1,sizeof(double)*arrSize/p);
    	}
        buckets[i]->size = 0;
    }

    int prev_index=0;
    int j;

    	 //for(i=0;i<arrSize/p;i++)
    	// 	printf("sdsa%d\n",recvBuf[i]);
		// cout << "Rank" <<rank <<"\t";
  //   	for(i = 0; i < arrSize/p; i++)
  //   		cout << "RecvBuf" << i <<":"<<recvBuf[i]<<endl;	    

  //   	cout<<endl;

    for(i = 0 ; i < p-1; i++) {
    	
    	//printf("b4: Rank :%d i:%d pivot:%d\n",rank,i,globalPivot[i]);
    	index = binary_search(globalPivot[i],recvBuf,prev_index, tmpSize);
    	//printf("Rank :%d index val:%d\n",rank,index);
    	for(j = prev_index; j<=index; j++){
    		buckets[i]->array[(int)buckets[i]->size] = recvBuf[j];
    		buckets[i]->size++;
    	}
    	if(index!=-1)
    		prev_index = index+1;
    }

	for(j = prev_index; j<tmpSize; j++){
		buckets[p-1]->array[(int)buckets[i]->size] = recvBuf[j];
		buckets[p-1]->size++;
	}


		// cout << "Rank" <<rank <<"\t";
  //   	for(i = 0; i < p; i++)
  //   		cout << "BucketSize" << i <<":"<<buckets[i]->size<<endl;	    
  //  		cout<<endl;

    const int nitems=2;
    int blocklengths[2] = {1,lastRankSize};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Datatype mpi_bucket_type;
    MPI_Aint     offsets[2];

     offsets[0] = offsetof(bucket_t, size);
     offsets[1] = offsetof(bucket_t, array);
    //MPI_Type_contiguous(1+(arrSize/p),MPI_DOUBLE,&mpi_bucket_type);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_bucket_type);
    MPI_Type_commit(&mpi_bucket_type);

	//printf("sizeof bucket:%d\n",sizeof(bucket_t)+sizeof(int)*arrSize/p);


	//cout << buckets[rank]->size << endl;


	
//#if 1
	// MPI_Gather(&bucketSize, p, MPI_INT, recvNum, p, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Request * requests = (MPI_Request *) malloc(sizeof(MPI_Request) * p);
    for(j=0;j<p;j++)
	{
		if(j!=rank) {
//			MPI_Send(buckets[j]->array,arrSize/p,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
//			printf("rankb:%d:%d\n",rank,arrSize/p);	
//			MPI_Send(buckets[j]->array,buckets[j]->size,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
			MPI_Isend(buckets[j],1,mpi_bucket_type,j,0,MPI_COMM_WORLD,&requests[j]);
//			printf("ranka:%d\n",rank);
			//MPI_Send(val,1,mpi_bucket_type,j,0,MPI_COMM_WORLD);
		}
	}



	bucket_t *tempBucket = (bucket_t *)malloc(sizeof(bucket_t)+sizeof(double)*lastRankSize);
//	cout <<"Rank"<<rank<<"\t";
	for(j=0;j<p;j++)
	{
		if(j!=rank) {
			//MPI_Recv(&tempArray[pos],arrSize/p,MPI_DOUBLE,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(tempBucket,1,mpi_bucket_type,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			// // MPI_Recv(&(buckets[rank]->array[(int)buckets[rank]->size]),arrSize/p,MPI_DOUBLE,j,0,
			// 			MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			memcpy(&(buckets[rank]->array[(int)buckets[rank]->size]),&tempBucket->array[0],
				sizeof(double)*tempBucket->size);
			buckets[rank]->size += tempBucket->size;
			//buckets[rank]->size += rnum[rank+(j*p)];
//			cout << buckets[rank]->size <<"\t";
		}
	}
	//cout <<endl;

    qsort(buckets[rank]->array,buckets[rank]->size,sizeof(double),compare);



	// cout << "rank:" << rank <<"\t";
	// for(j=0;j<buckets[rank]->size;j++){
	// 	cout << buckets[rank]->array[j] <<"\t";
	// }
	// cout << endl;	

	// printf("\n");
#if 1
//#if 1
	// MPI_Gather(&buckets[rank]->size, 1, MPI_DOUBLE, recvCounts, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&buckets[rank]->size, 1, MPI_INT, recvCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank==0){
		disp[0] = 0;
  		for(i=1;i<p;i++){
			disp[i] = disp[i-1] + recvCounts[i-1];
		}
	}

	// if(rank==0){
	//  	for(i=0;i<p;i++){
	// 		cout << recvCounts[i] << "\t";
	//  	}
	// }
	// cout <<buckets[rank]->sizeof<<endl;


//#if 0
	MPI_Gatherv(buckets[rank]->array, buckets[rank]->size, MPI_DOUBLE,
            	finalBuf, recvCounts, disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//printf("buckets rank size:%d:%d\n",buckets[rank]->size,arrSize);
	//MPI_Gather( buckets[rank]->array, buckets[rank]->size, MPI_DOUBLE, rbuf, buckets[rank]->size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	 if(rank==0) {
	    for(i = 0; i < arrSize; i++)
    		cout <<finalBuf[i]<<"\t";	    
	 	//printf("\n");
	 	cout << endl;
	 }
#endif

}

int main (int argc, char* argv[])
{
	int rank, size;
	double *arr;
	uint64_t n = N;
	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
	printf( "Hello world from process %d of %d\n", rank, size );
	
	scatterData(arr,n,rank,size);
	
	MPI_Finalize();
	return 0;
}
