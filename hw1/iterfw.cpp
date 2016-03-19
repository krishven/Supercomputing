#include <iostream>
#include <omp.h>
//#include <cilk.h>
//#include <cilkview.h>

#define V 4
#define NA 100000

void printSolution(int graph[][V]);

void iterFW(int graph[][V])
{
    //cilk::cilkview cv;
    //cv.start();
    double start = omp_get_wtime();
    for (int k = 0; k < V; k++)
    {
        #pragma omp parallel for
        for (int i = 0; i < V; i++)
        {
            #pragma omp parallel for
            for (int j = 0; j < V; j++)
            {
                if (graph[i][k] + graph[k][j] < graph[i][j])
                    graph[i][j] = graph[i][k] + graph[k][j];
            }
        }
    }
    std::cout << "Time diff" << omp_get_wtime()-start <<"\n";
   // cv.stop();
    //cv.dump("results", false);
    //std::cout << cv.accumulated_milliseconds() / 1000.f << " seconds" << std::endl;
    printSolution(graph);
}

void printSolution(int graph[][V])
{
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
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
    int graph[V][V] = { {0,   5,  NA, 10},
        {NA, 0,   3, NA},
        {NA, NA, 0,   1},
        {NA, NA, NA, 0}
    };

    // Print the solution
    iterFW(graph);
    return 0;
}