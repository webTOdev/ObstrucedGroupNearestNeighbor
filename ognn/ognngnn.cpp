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
void print(std::vector < MyStruct >& ognn_sorted,std::vector < MyStruct >& egnn_sorted,int k);
void centroidOfQ(Point2D queryPoints[], int numOfQueryPoints, Point2D centroid[]);

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

	std::sort(egnn_sorted.begin(), egnn_sorted.end(), more_than_key());
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
	print(egnn_sorted,ognn_sorted,k);
	double dkmax_e=egnn_sorted[0].distance;
	double dkmax_o=ognn_sorted[0].distance;

	if(dkmax_e != dkmax_o){
		double nearestNeighbor[3];
		while(1){
		//Find next Group Nearest Neighbour
			kNN_point = new float[3];
			rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor,queryPoints,numOfQueryPoints);
			kNN_point[0]=nearestNeighbor[0];
			kNN_point[1]=nearestNeighbor[1];
			kNN_point[2]=nearestNeighbor[2];
			egnn_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
			printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
			printf("\ndkmax = %lf \n",dkmax_o);	
			//Euclidean distance is greater so no need to check obs dist
			if(kNN_point[2]<=dkmax_o){
				double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
				ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
				//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);
				printf("New OGNN found\n");
				dkmax_o=aggObsDist;
			}
			else
				break;

		}

		std::sort(egnn_sorted.begin(), egnn_sorted.end(), more_than_key());
		std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
		print(egnn_sorted,ognn_sorted,k);
	}


	delete kNN_point;
	delete obstructedDistance;
	delete initialVisGraph;
	delete rt_dataPoints->kGNNHeap;
}

void OGNN_GNN::ognnSumUsingNN(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function) {

	float centroid[1][2];
	centroidOfQ(queryPoints,numOfQueryPoints,centroid);
	printf("\n----------------------------------------Searching for k-NN---------------------------------------\n");
	
	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(centroid, k, kNearestNeighbor,1);
	
	printf("\nk- Nearest Neighbor of centroid (%f,%f)is ", centroid[0][0],centroid[0][1]);
	for(int index=0;index<k;index++){
		printf("(%f,%f) dist %lf\t", kNearestNeighbor[index][0], kNearestNeighbor[index][1],kNearestNeighbor[index][2]);
	}

	//Add the query points in the vis graph -- **only once**
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	ObstructedDistance* obstructedDistance= new ObstructedDistance();
	obstructedDistance->writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	obstructedDistance->constructInitialVisGraph(initialVisGraph);

		//This vector will store the obstructed group NN
	std::vector < MyStruct > ognn_sorted,enn_centroid_sorted;
	printf("\n----------------------------------------Searching for k-Obstructed GNN---------------------------------------\n");
	float *kNN_point,*q;
	
	for(int i=0;i<k;i++){
		//To sort the aggregate euclidean distance of p and Q
		kNN_point = new float[3];
		kNN_point[0] = kNearestNeighbor[i][0];
		kNN_point[1] = kNearestNeighbor[i][1];
		kNN_point[2] = kNearestNeighbor[i][2];
		//printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
		enn_centroid_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
		//Compute aggObstructedDistance
		double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
		ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
		printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);	
	}

	std::sort(enn_centroid_sorted.begin(), enn_centroid_sorted.end(), more_than_key());
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
	print(enn_centroid_sorted,ognn_sorted,k);
	double dkmax_e=enn_centroid_sorted[0].distance;
	double dkmax_o=ognn_sorted[0].distance;
	
		double nearestNeighbor[3];
		while(dkmax_e < dkmax_o/numOfQueryPoints){
			printf("dkmax_e %lf  dkmax_o %lf",dkmax_e,dkmax_o/numOfQueryPoints);

		//Find next Group Nearest Neighbour
			kNN_point = new float[3];
			rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor,centroid,1);
			kNN_point[0]=nearestNeighbor[0];
			kNN_point[1]=nearestNeighbor[1];
			kNN_point[2]=nearestNeighbor[2];
			enn_centroid_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
			printf("\nEuclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
			dkmax_e = kNN_point[2];
			double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
			if(aggObsDist<=dkmax_o){			
				ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
				//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);
				printf("New OGNN found\n");
				dkmax_o=aggObsDist;
			}
		}

		std::sort(enn_centroid_sorted.begin(), enn_centroid_sorted.end(), more_than_key());
		std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
		print(enn_centroid_sorted,ognn_sorted,k);


	delete kNN_point;
	delete obstructedDistance;
	delete initialVisGraph;
	delete rt_dataPoints->kGNNHeap;
}

void centroidOfQ(Point2D queryPoints[], int numOfQueryPoints, Point2D centroid[]){
	float x=0.0;
	float y=0.0;
	for(int i=0;i<numOfQueryPoints;i++){
		x+=queryPoints[i][0];
		y+=queryPoints[i][1];
	}
	centroid[0][0]=x/numOfQueryPoints;
	centroid[0][1]=y/numOfQueryPoints;

}

void print(std::vector < MyStruct >& egnn_sorted,std::vector < MyStruct >& ognn_sorted,int k){
	printf("\nInitial K-Euclidean Group Nearest Neighbor\n");
	float* q;
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=egnn_sorted[index].queryPoints[0];
			q[1]=egnn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],egnn_sorted[index].distance);

			delete q;

	}
	printf("\nInitial K-Obstructed Group Nearest Neighbor\n");
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=ognn_sorted[index].queryPoints[0];
			q[1]=ognn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],ognn_sorted[index].distance);

			delete q;

	}
}






