#include <iostream>
#include <omp.h>
#include <cstdlib>
//#include <cilk.h>
//#include <cilkview.h>


#define NA 100000

void printSolution(int **graph,int n);

void iterFW(int **graph, int n)
{
    //cilk::cilkview cv;
    //cv.start();
    for (int k = 0; k < n; k++)
    {
        #pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            #pragma omp parallel for
            for (int j = 0; j < n; j++)
            {
                if (graph[i][k] + graph[k][j] < graph[i][j])
                    graph[i][j] = graph[i][k] + graph[k][j];
            }
        }
    }
   // cv.stop();
    //cv.dump("results", false);
    //std::cout << cv.accumulated_milliseconds() / 1000.f << " seconds" << std::endl;
    //printSolution(graph,n);
}

void printSolution(int **graph,int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (graph[i][j] == NA)
                std::cout << "       " << "-";
            else
                std::cout << "       "   << graph[i][j];
        }
        std::cout << "\n";
    }
}

int main(int argc, char* argv[])
{
/*
    int graph[V][V] = { {0,   5,  NA, 10},
        {NA, 0,   3, NA},
        {NA, NA, 0,   1},
        {NA, NA, NA, 0}
    };
*/
    int n=8192;

   omp_set_num_threads(16);

   while(n>1)
 {	
    int **graph = new int*[n];
    for (int i = 0; i < n; ++i)
        graph[i] = new int[n];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                graph[i][j] = 0;
            else
                graph[i][j] = i+j;

    // Print the solution
    double start = omp_get_wtime();
    iterFW(graph,n);
    std::cout <<n << "\t" << omp_get_wtime()-start << "\n";

    for ( int i = 0 ; i < n; i++ )
    {
        delete graph[i];
    }
    delete[] graph;    
    n=n/2;	
}
    return 0;
}
