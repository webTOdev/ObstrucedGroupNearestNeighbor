//Added by Nusrat

#include "ognngnn.h"

//#include <math.h>
#include <cstdio>
#include <fstream>
#include <memory>
#include <new>
#include <string>
#include <vector>
#include <algorithm>

#include "../func/gendef.h"
#include "../heap/heap.h"
#include "../visGraph/obstacleController.h"
#include "../visGraph/VisibilityGraph.h"
#include "../visGraph/VisibilityGraphController.h"
#include "odist.h"

//const int infinity = 1000000000; 

//VisibilityGraphController* vg;

void OGNN_GNN::ognnUsingEGNN(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function) {

	printf("\n----------------------------------------Searching for k-GNN---------------------------------------\n");
	
	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(queryPoints, k, kNearestNeighbor,numOfQueryPoints);
	
	printf("\nk- Group Nearest Neighbor of (%f,%f),(%f,%f),(%f,%f) is ", queryPoints[0][0],queryPoints[0][1],
			queryPoints[1][0],queryPoints[1][1],queryPoints[2][0],queryPoints[2][1]);
	for(int index=0;index<k;index++){
		printf("(%f,%f) dist %lf\t", kNearestNeighbor[index][0], kNearestNeighbor[index][1],kNearestNeighbor[index][2]);
	}

	//Add the query points in the vis graph -- **only once**
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	ObstructedDistance* obstructedDistance= new ObstructedDistance();
	obstructedDistance->writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	obstructedDistance->constructInitialVisGraph(initialVisGraph);

		//This vector will store the obstructed group NN
	std::vector < MyStruct > ognn_sorted,egnn_sorted;
	printf("\n----------------------------------------Searching for k-Obstructed GNN---------------------------------------\n");
	float *kNN_point,*q;
	
	for(int i=0;i<k;i++){
		//To sort the aggregate euclidean distance of p and Q
		kNN_point = new float[3];
		kNN_point[0] = kNearestNeighbor[i][0];
		kNN_point[1] = kNearestNeighbor[i][1];
		kNN_point[2] = kNearestNeighbor[i][2];
		printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
		egnn_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
		//Compute aggObstructedDistance
		double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
		ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
		printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);	
	}

	std::sort(egnn_sorted.begin(), egnn_sorted.end(), less_than_key());
	std::sort(ognn_sorted.begin(), egnn_sorted.end(), less_than_key());
	delete kNN_point;
	delete obstructedDistance;
	delete initialVisGraph;
}






