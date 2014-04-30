/*
 * File:   main.cpp
 * Author: artur
 *
 * Created on 28 kwiecień 2012, 15:48
 */

#include <cstdlib>
#include<cstdio>

#include<set>
#include "dijkstra.h"
#include <cmath>

using namespace std;

const int infty = 1000000000; // limit
int MAX_NUM_VERTEX=pow(2.0,16.0);
int verticesNum;
int edgesNum;
vector< vector< pair<int, double> > > comparisionMatrix;
vector<double> pathWeight;
vector<int> smallestPath;
int *prev;

void printPath(int dest);


struct cmp // the condition of the mound
{
    // if a is smaller than b

    bool operator() (const int &a, const int &b) {
        if (pathWeight[a] < pathWeight[b]) return true;
        if (pathWeight[a] > pathWeight[b]) return false;
        return a<b;
    }
};

set<int, cmp> kopiec; // ;-)

void prim_dijkstra(int s) //s starting point
{
    int v, u;
    double c;

    pathWeight.clear();
    pathWeight.resize(MAX_NUM_VERTEX, infty);
    pathWeight[s] = 0;

	prev = new int[MAX_NUM_VERTEX];
	smallestPath.clear();


    kopiec.clear();
    for (int i = 0; i < MAX_NUM_VERTEX; i++) { //throw vertices, k to the mound
        kopiec.insert(i);
        prev[i]=-1;
      //  smallestPath[i]=-1;
    }

    while (!kopiec.empty()) //go to the mound effect
    {
        u = *(kopiec.begin()); // weź wierzchołek najbliżej drzewa MST
        kopiec.erase(kopiec.begin());

        for (int i = 0; i < comparisionMatrix[u].size(); i++) //połączenia danego wierzchołka i koszty
        {
            v = comparisionMatrix[u][i].first;

            c = comparisionMatrix[u][i].second;
            if (pathWeight[u] + c < pathWeight[v]) // w alg. Prima jest tutaj c < pathWeight[v]
            {
                // uaktualniamy wagę wierzchołka v - poprawnośc przez indukcję, dla co raz większego drzewa mst wybieramy krawędzie o najmiejszej wadze
                kopiec.erase(kopiec.find(v));
                pathWeight[v] = pathWeight[u] + c; // w alg. Prima jest tutaj pathWeight[v] = c ;
                kopiec.insert(v);
                prev[v]=u;
                //	mstClosingVertices[v] = u; // domykam krawędź
            }

        }
    }

}


double initiateDijkstra(int numVertice,int numEdges,bool directed,int source,int destination,int maxVertexNum) {
    int a, b, s, g;
    double c;
	
    FILE *input = fopen("test.txt", "r+");
    verticesNum = numVertice;
    edgesNum=numEdges;
	MAX_NUM_VERTEX=maxVertexNum;
    //printf("Vertice Num %d , Edge Num %d\n",verticesNum,edgesNum);
    comparisionMatrix.resize(MAX_NUM_VERTEX); // first edge � will turn
/*	double **temp=new double*[MAX_NUM_VERTEX];
	for (int i = 0; i < MAX_NUM_VERTEX; i++)
		temp[i] = new double[MAX_NUM_VERTEX];
	
	for(int i=0;i<MAX_NUM_VERTEX;i++){
		for(int j=0;j<MAX_NUM_VERTEX;j++){
		temp[i][j]=-1;
		}
	}
	for (int i = 0; i < edgesNum; i++) {
		fscanf(input, "%d %d %lf", &a, &b, &c); // c = Cost edge from a to b
		temp[a][b]=c;
	}
	for(int i=0;i<MAX_NUM_VERTEX;i++){
		for(int j=0;j<MAX_NUM_VERTEX;j++){
			if(temp[i][j]!=-1 && temp[j][i]!=-1){
				 comparisionMatrix[i].push_back(make_pair(j, c)); // directed not addressed
			}
		}
	}*/
    for (int i = 0; i < edgesNum; i++) {
        fscanf(input, "%d %d %lf", &a, &b, &c); // c = Cost edge from a to b
        comparisionMatrix[a].push_back(make_pair(b, c)); // directed not addressed
        if(!directed)
        comparisionMatrix[b].push_back(make_pair(a, c)); //
    }

    s=source;
    g=destination;
   // fscanf(input, "%d %d", &s, &g); // sTART gOAL
    prim_dijkstra(s); // true false targets not addressed

    printf("[ALgorithm DIJKSTRA] : The cost of the minimum path from %d to %d %s is %lf\n", s, g, (directed)?"directed":"", pathWeight[g]); // index 5 czyli v nr 6
    printPath(g);
    printf("\n");
	fclose(input);

	/*for(int i = 0; i < MAX_NUM_VERTEX; ++i) {
    delete [] temp[i];
	}
	delete [] temp;*/
    return pathWeight[g];
}

/*
 * Prints the shortest path from the source to dest.
 *
 * dijkstra(int) MUST be run at least once BEFORE
 * this is called
 */
void printPath(int dest) {
	static int index=0;
	if (prev[dest] != -1)
		printPath(prev[dest]);
	printf("\n%d ", dest);
//	smallestPath[index]=dest;
	smallestPath.push_back(dest);
	index++;
}
std::vector<int> getShortestPath(){

	return smallestPath;
}

