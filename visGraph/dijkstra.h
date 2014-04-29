/*
 * Dijkstra.h
 *
 *  Created on: Sep 12, 2013
 *      Author: nut
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include <iostream>
#include<vector>

using namespace std;

double initiateDijkstra(int numVertice,int numEdges,bool directed,int source,int destination,int maxVertexNum);
void printPath(int dest);
void prim_dijkstra(int s);
vector<int> getShortestPath();

#endif /* DIJKSTRA_H_ */
