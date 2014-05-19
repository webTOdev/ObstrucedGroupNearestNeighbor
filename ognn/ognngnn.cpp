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
#include "odistcentroid.h"

//const int infinity = 1000000000; 

//VisibilityGraphController* vg;
void print(std::vector < MyStruct >& ognn_sorted,std::vector < MyStruct >& egnn_sorted,int k);

void OGNN_GNN::ognnUsingEGNN(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function) {

	totalNumberOfPRetrieved=0.0;
	
	//printf("\n----------------------------------------Searching for k-GNN---------------------------------------\n");
	
	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(queryPoints, k, kNearestNeighbor,numOfQueryPoints,function);
	totalNumberOfPRetrieved+=k;
	
	/*printf("\nk- Group Nearest Neighbor of (%f,%f),(%f,%f),(%f,%f) is ", queryPoints[0][0],queryPoints[0][1],
			queryPoints[1][0],queryPoints[1][1],queryPoints[2][0],queryPoints[2][1]);
	for(int index=0;index<k;index++){
		printf("(%f,%f) dist %lf\t", kNearestNeighbor[index][0], kNearestNeighbor[index][1],kNearestNeighbor[index][2]);
	}*/

	//Add the query points in the vis graph -- **only once**
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	//ObstructedDistance* obstructedDistance= new ObstructedDistance();
	ObstructedDistanceCentroid* obstructedDistance= new ObstructedDistanceCentroid();
	obstructedDistance->writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	Clock sw1;
	sw1.start();
	obstructedDistance->constructInitialVisGraph(initialVisGraph);		
	sw1.stop();
	visGraphConsTime+=sw1.getDiff();


		//This vector will store the obstructed group NN
	std::vector < MyStruct > ognn_sorted,egnn_sorted;
	//printf("\n----------------------------------------Searching for k-Obstructed GNN---------------------------------------\n");
	float *kNN_point,*q;
	
	for(int i=0;i<k;i++){
		//To sort the aggregate euclidean distance of p and Q
		kNN_point = new float[3];
		kNN_point[0] = kNearestNeighbor[i][0];
		kNN_point[1] = kNearestNeighbor[i][1];
		kNN_point[2] = kNearestNeighbor[i][2];
		//printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
		egnn_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
		//Compute aggObstructedDistance
		double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
		ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
		//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);	
	}

	std::sort(egnn_sorted.begin(), egnn_sorted.end(), more_than_key());
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
	//print(egnn_sorted,ognn_sorted,k);
	double dkmax_e=egnn_sorted[0].distance;
	double dkmax_o=ognn_sorted[0].distance;

	if(dkmax_e != dkmax_o){
		double nearestNeighbor[3];
		while(1){
		//Find next Group Nearest Neighbour
			kNN_point = new float[3];
			rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor,queryPoints,numOfQueryPoints,function);
			totalNumberOfPRetrieved++;
			kNN_point[0]=nearestNeighbor[0];
			kNN_point[1]=nearestNeighbor[1];
			kNN_point[2]=nearestNeighbor[2];
			egnn_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
			//printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
			//printf("\ndkmax = %lf \n",dkmax_o);	
			//Euclidean distance is greater so no need to check obs dist
			if(kNN_point[2]<=dkmax_o){
				double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
				ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
				//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);
				//printf("New OGNN found\n");
				dkmax_o=aggObsDist;
			}
			else
				break;

		}

	}
	
	std::sort(egnn_sorted.begin(), egnn_sorted.end(), less_than_key());
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());
	if(function==0){
		printf("********************OGNN-GNN-SUM*******************************\n");
	/*	FILE * outputFile1;
		outputFile1 = fopen("Result/ognnOutput", "a+");
		fprintf(outputFile1,"\n********************OGNN-GNN-SUM*******************************\n");
		fclose(outputFile1);*/

	}
	if(function==1){
		printf("********************OGNN-GNN-MAX*******************************\n");
		/*FILE * outputFile1;
		outputFile1 = fopen("Result/ognnOutput", "a+");
		fprintf(outputFile1,"\n********************OGNN-GNN-MAX*******************************\n");
		fclose(outputFile1);*/
	}
	print(egnn_sorted,ognn_sorted,k);
	visGraphConsTime+=obstructedDistance->visGraphConsTime;
	shortestPathCalcTime+=obstructedDistance->shortestPathCalcTime;

	delete kNN_point;
	delete obstructedDistance;
	delete initialVisGraph;
	delete rt_dataPoints->kGNNHeap;
}

void OGNN_GNN::ognnSumUsingNN(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function) {

	totalNumberOfPRetrieved=0.0;
	float centroid[1][2];
	centroidOfQ(queryPoints,numOfQueryPoints,centroid);
	//printf("\n----------------------------------------Searching for k-NN---------------------------------------\n");
	
	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(centroid, k, kNearestNeighbor,1,function);
	totalNumberOfPRetrieved+=k;
	/*printf("\nk- Nearest Neighbor of centroid (%f,%f)is ", centroid[0][0],centroid[0][1]);
	for(int index=0;index<k;index++){
		printf("(%f,%f) dist %lf\t", kNearestNeighbor[index][0], kNearestNeighbor[index][1],kNearestNeighbor[index][2]);
	}*/

	//Add the query points in the vis graph -- **only once**
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	//ObstructedDistance* obstructedDistance= new ObstructedDistance();
	ObstructedDistanceCentroid* obstructedDistance= new ObstructedDistanceCentroid();
	obstructedDistance->writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	Clock sw1;
	sw1.start();
	obstructedDistance->constructInitialVisGraph(initialVisGraph);
	sw1.stop();
	visGraphConsTime+=sw1.getDiff();

		//This vector will store the obstructed group NN
	std::vector < MyStruct > ognn_sorted,enn_centroid_sorted;
	//printf("\n----------------------------------------Searching for k-Obstructed GNN---------------------------------------\n");
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
		//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);	
	}

	std::sort(enn_centroid_sorted.begin(), enn_centroid_sorted.end(), more_than_key());
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
	//print(enn_centroid_sorted,ognn_sorted,k);
	double dkmax_e=enn_centroid_sorted[0].distance;
	double dkmax_o=ognn_sorted[0].distance;
	
		double nearestNeighbor[3];
		while(dkmax_e < dkmax_o/numOfQueryPoints){
			//printf("dkmax_e %lf  dkmax_o %lf",dkmax_e,dkmax_o/numOfQueryPoints);

		//Find next Group Nearest Neighbour
			kNN_point = new float[3];
			rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor,centroid,1,function);
			totalNumberOfPRetrieved++;
			kNN_point[0]=nearestNeighbor[0];
			kNN_point[1]=nearestNeighbor[1];
			kNN_point[2]=nearestNeighbor[2];
			enn_centroid_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
			//printf("\nEuclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
			dkmax_e = kNN_point[2];
			double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
			if(aggObsDist<=dkmax_o){			
				ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
				//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);
				//printf("New OGNN found\n");
				dkmax_o=aggObsDist;
			}
		}

		std::sort(enn_centroid_sorted.begin(), enn_centroid_sorted.end(), less_than_key());
		std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());
		printf("********************OGNN-CENTROID-NN-SUM*******************************\n");
		/*FILE * outputFile1;
		outputFile1 = fopen("Result/ognnOutput", "a+");
		fprintf(outputFile1,"\n********************OGNN-CENTROID-NN-SUM*******************************\n");
		fclose(outputFile1);*/
		print(enn_centroid_sorted,ognn_sorted,k);

	visGraphConsTime+=obstructedDistance->visGraphConsTime;
	shortestPathCalcTime+=obstructedDistance->shortestPathCalcTime;
	delete kNN_point;
	delete obstructedDistance;
	delete initialVisGraph;
	delete rt_dataPoints->kGNNHeap;
}

void OGNN_GNN::ognnMaxUsingNN(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function) {

	totalNumberOfPRetrieved=0.0;
	float centroid[1][2];
	centroidOfQ(queryPoints,numOfQueryPoints,centroid);
	//printf("\n----------------------------------------Searching for k-NN---------------------------------------\n");
	
	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(centroid, k, kNearestNeighbor,1,function);
	totalNumberOfPRetrieved+=k;
	
	//printf("\nk- Nearest Neighbor of centroid (%f,%f)is ", centroid[0][0],centroid[0][1]);
	/*for(int index=0;index<k;index++){
		printf("(%f,%f) dist %lf\t", kNearestNeighbor[index][0], kNearestNeighbor[index][1],kNearestNeighbor[index][2]);
	}*/

	//Add the query points in the vis graph -- **only once**
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
//	ObstructedDistance* obstructedDistance= new ObstructedDistance();
	ObstructedDistanceCentroid* obstructedDistance= new ObstructedDistanceCentroid();
	obstructedDistance->writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	Clock sw1;
	sw1.start();
	obstructedDistance->constructInitialVisGraph(initialVisGraph);
	sw1.stop();
	visGraphConsTime+=sw1.getDiff();

		//This vector will store the obstructed group NN
	std::vector < MyStruct > ognn_sorted,enn_centroid_sorted;
	//printf("\n----------------------------------------Searching for k-Obstructed GNN---------------------------------------\n");
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
		//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);	
	}

	std::sort(enn_centroid_sorted.begin(), enn_centroid_sorted.end(), more_than_key());
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), more_than_key());
	//print(enn_centroid_sorted,ognn_sorted,k);
	double dkmax_e=enn_centroid_sorted[0].distance;
	double dkmax_o=ognn_sorted[0].distance;
	
		double nearestNeighbor[3];
		while(dkmax_e < dkmax_o){
			//printf("dkmax_e %lf  dkmax_o %lf",dkmax_e,dkmax_o);

		//Find next Group Nearest Neighbour
			kNN_point = new float[3];
			rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor,centroid,1,function);
			totalNumberOfPRetrieved++;
			kNN_point[0]=nearestNeighbor[0];
			kNN_point[1]=nearestNeighbor[1];
			kNN_point[2]=nearestNeighbor[2];
			enn_centroid_sorted.push_back(MyStruct(kNN_point[2], kNN_point));
			//printf("\nEuclidean Distance is %lf for p %f,%f\n", kNN_point[2],kNN_point[0],kNN_point[1]);	
			dkmax_e = kNN_point[2];
			if(pointInsideTheIntersectionOfCircle(kNN_point,queryPoints,numOfQueryPoints,dkmax_o)){
				double aggObsDist = obstructedDistance->computeAggObstructedDistance(initialVisGraph,kNN_point,queryPoints,numOfQueryPoints,rt_obstacle,function);
				if(aggObsDist<=dkmax_o){			
					ognn_sorted.push_back(MyStruct(aggObsDist, kNN_point));
					//printf("\nAggregate Obstructed Distance is %lf for p %f,%f\n", aggObsDist,kNN_point[0],kNN_point[1]);
					//printf("New OGNN found\n");
					dkmax_o=aggObsDist;
				}
			}
		}

		std::sort(enn_centroid_sorted.begin(), enn_centroid_sorted.end(), less_than_key());
		std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());
		printf("********************OGNN-CENTROID-NN-MAX*******************************\n");
		//FILE * outputFile1;
		//outputFile1 = fopen("Result/ognnOutput", "a+");
		//fprintf(outputFile1,"\n********************OGNN-CENTROID-NN-MAX*******************************\n");
		//fclose(outputFile1);
		print(enn_centroid_sorted,ognn_sorted,k);

	visGraphConsTime+=obstructedDistance->visGraphConsTime;
	shortestPathCalcTime+=obstructedDistance->shortestPathCalcTime;
	delete kNN_point;
	delete obstructedDistance;
	delete initialVisGraph;
	delete rt_dataPoints->kGNNHeap;
}


void print(std::vector < MyStruct >& egnn_sorted,std::vector < MyStruct >& ognn_sorted,int k){
	float* q;
	/*printf("\nInitial K-Euclidean Group Nearest Neighbor\n");
	
	for(int index=0;index<egnn_sorted.size();index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=egnn_sorted[index].queryPoints[0];
			q[1]=egnn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],egnn_sorted[index].distance);

			delete q;

	}*/

	FILE * outputFile1;
	outputFile1 = fopen("Result/ognnOutput", "a+");

	printf("\nInitial K-Obstructed Group Nearest Neighbor\n");
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=ognn_sorted[index].queryPoints[0];
			q[1]=ognn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],ognn_sorted[index].distance);
			fprintf(outputFile1,"k=%d, (%f,%f) distance %lf\n",index,q[0], q[1],ognn_sorted[index].distance);

			delete q;

	}

		fclose(outputFile1);
}

bool OGNN_GNN::pointInsideTheIntersectionOfCircle(float* p, Point2D queryPoints[],int numOfQueryPoints,double radius){
	bool inside=true;
	for(int i=0;i<numOfQueryPoints;i++){
		double dist = getDistanceBetweenTwoPoints(p,queryPoints[i]);
		if(dist>radius){
			inside=false;
		}
	}
	return inside;

}






