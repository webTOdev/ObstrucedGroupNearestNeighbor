/*
 * Dijkstra.h
 *
 *  Created on: Sep 12, 2013
 *      Author: nut
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include <iostream>

double initiateDijkstra(int numVertice,int numEdges,bool directed,int source,int destination,int maxVertexNum);
void printPath(int dest);
void prim_dijkstra(int s);
int * getShortestPath();

#endif /* DIJKSTRA_H_ */