//Added by Nusrat

#include "ognn_gnn.h"

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

const int infinity = 1000000000; 

VisibilityGraphController* vg;




double OGNN::computeEuclideanGroupDistance(float* kNN_point,std::vector < MyStruct >& queryPoints_sorted,int numOfQueryPoints,Point2D queryPoints[] ){
	
		double aggregateEuclideanDist=0;
		//std::vector < MyStruct > queryPoints_sorted;
		for(int j=0;j<numOfQueryPoints;j++){
			float *q = new float[2];
			q[0]=queryPoints[j][0];
			q[1]=queryPoints[j][1];
			double distance = getDistanceBetweenTwoPoints(kNN_point, q);
			aggregateEuclideanDist+=distance;
			//printf("\nEuclidean Distance is %lf for q %f,%f", distance,q[0],q[1]);
			queryPoints_sorted.push_back(MyStruct(distance, q));
		}

		return aggregateEuclideanDist;
}



void OGNN_GNN::ognnMultiPointApproach(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][2], RTree* rt_obstacle,
		RTree* rt_dataPoints) {

	printf("\n----------------------------------------Searching for k-GNN---------------------------------------\n");
	double nearestNeighbor[2];
	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(queryPoints, k, kNearestNeighbor,
			numOfQueryPoints);
	
	printf("\nk- Group Nearest Neighbor of (%f,%f),(%f,%f),(%f,%f) is ", queryPoints[0][0],queryPoints[0][1],
			queryPoints[1][0],queryPoints[1][1],queryPoints[2][0],queryPoints[2][1]);
	for(int index=0;index<k;index++){
		printf("(%f,%f)\t",
			kNearestNeighbor[index][0], kNearestNeighbor[index][1]);
	}
	/*rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor);
	printf("\n Next k- Group Nearest Neighbor of (%f,%f)",nearestNeighbor[0],nearestNeighbor[1]);
	*/

	//This vector will store the obstructed group NN
	std::vector < MyStruct > ognn_sorted,egnn_sorted;
	printf("\n----------------------------------------Searching for k-Obstructed GNN---------------------------------------\n");
	float *kNN_point,*q;
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	for(int i=0;i<k;i++){
		//To sort the aggregate euclidean distance of p and Q
		std::vector < MyStruct > queryPoints_sorted;
		kNN_point = new float[2];
		kNN_point[0] = kNearestNeighbor[i][0];
		kNN_point[1] = kNearestNeighbor[i][1];
		//Compute  the aggregate euclidean distance of p and Q
		double aggregateEuclideanDist = computeEuclideanGroupDistance(kNN_point,queryPoints_sorted,numOfQueryPoints,queryPoints);		
		
		printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", aggregateEuclideanDist,kNN_point[0],kNN_point[1]);
		egnn_sorted.push_back(MyStruct(aggregateEuclideanDist, kNN_point));

		std::sort(queryPoints_sorted.begin(), queryPoints_sorted.end(), more_than_key());
		printf("Euclidean Distance %lf \n",queryPoints_sorted[0].distance);

		vector<vector<pair<double,double>>> sp;
	/*	for(int i=0;i<numOfQueryPoints;i++){
			shortestPathMatrixEuclidean(i,kNN_point[0]kNN_point[1],sp);
		}*/
				
		double obstructedDistance = computeObsDistanceForTheFarthestQueryPoint(kNN_point,queryPoints_sorted,rt_obstacle,i,queryPoints,initialVisGraph  );
		
		//Initialize the aggregate distance with the distance between p and first q 
		double aggregateObstructedDistance=obstructedDistance;	
		//first element in vector is already considered so start loop from j=1
		aggregateObstructedDistance=computeObsDistanceForOtherQueryPoints(kNN_point,queryPoints_sorted,initialVisGraph,
												   aggregateObstructedDistance,numOfQueryPoints,rt_obstacle);
		
		//printf("\nAggregate Obstructed Distance of data point p(%f,%f) is %lf\n",kNN_point[0],kNN_point[1], aggregateObstructedDistance);
		//The data point is reachable from all the query points distance !=infinity
		ognn_sorted.push_back(MyStruct(aggregateObstructedDistance, kNN_point));
		//Now remove the Data point from Vis graph as obstructed group distance calculation is finished
		removeDataPointFromVG(initialVisGraph,kNN_point);
	}
	//Sort the Euclidean GNN
	std::sort(egnn_sorted.begin(), egnn_sorted.end(), less_than_key());
	//Sort the Obstructed GNN
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());
	/*printf("\nInitial K-Euclidean Group Nearest Neighbor\n");
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=egnn_sorted[index].queryPoints[0];
			q[1]=egnn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],egnn_sorted[index].distance);

	}
	printf("\nInitial K-Obstructed Group Nearest Neighbor\n");
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=ognn_sorted[index].queryPoints[0];
			q[1]=ognn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],ognn_sorted[index].distance);

	}

	int i=k-1;
	double maxObsDistance=ognn_sorted[k-1].distance;
	while(1){
		//Find next Group Nearest Neighbour
		i++;
		kNN_point = new float[2];
		rt_dataPoints->retrieve_kth_BFN_GNNQ(nearestNeighbor,queryPoints,numOfQueryPoints);
		kNN_point[0]=nearestNeighbor[0];
		kNN_point[1]=nearestNeighbor[1];
		std::vector < MyStruct > queryPoints_sorted;
		//Compute  the aggregate euclidean distance of p and Q
		double aggregateEuclideanDist = computeEuclideanGroupDistance(kNN_point,queryPoints_sorted,numOfQueryPoints,queryPoints);
		printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", aggregateEuclideanDist,kNN_point[0],kNN_point[1]);
		if(aggregateEuclideanDist>maxObsDistance){
			printf("\n-----------------------Breaking from while loop!!-----------------------\n");
			printf("\nAggregate Euclidean Distance %lf of data point p(%f,%f) is Greater than Max obstructed distance %lf\n",aggregateEuclideanDist,kNN_point[0],kNN_point[1], maxObsDistance);
			break;
		}
				
		//printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", aggregateEuclideanDist,kNN_point[0],kNN_point[1]);
		egnn_sorted.push_back(MyStruct(aggregateEuclideanDist, kNN_point));

		std::sort(queryPoints_sorted.begin(), queryPoints_sorted.end(), more_than_key());
		double obstructedDistance = computeObsDistanceForTheFarthestQueryPoint(kNN_point,queryPoints_sorted,rt_obstacle,i,queryPoints,initialVisGraph  );
		
		//Initialize the aggregate distance with the distance between p and first q 
		double aggregateObstructedDistance=obstructedDistance;	
		//first element in vector is already considered so start loop from j=1
		aggregateObstructedDistance=computeObsDistanceForOtherQueryPoints(kNN_point,queryPoints_sorted,initialVisGraph,
												   aggregateObstructedDistance,numOfQueryPoints,rt_obstacle);
		if(aggregateObstructedDistance<maxObsDistance)
		{
			maxObsDistance=aggregateObstructedDistance;
			printf("\nAggregate Obstructed Distance %lf of data point p(%f,%f) is closer than Max %lf\n",aggregateObstructedDistance,kNN_point[0],kNN_point[1],maxObsDistance );
			ognn_sorted.push_back(MyStruct(aggregateObstructedDistance, kNN_point));
		}
		//Now remove the Data point from Vis graph as obstructed group distance calculation is finished
		removeDataPointFromVG(initialVisGraph,kNN_point);
	}

	//Sort the Euclidean GNN
	std::sort(egnn_sorted.begin(), egnn_sorted.end(), less_than_key());

	printf("\nFinal K-Euclidean Group Nearest Neighbor\n");
		for(int index=0;index<egnn_sorted.size();index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=egnn_sorted[index].queryPoints[0];
			q[1]=egnn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],egnn_sorted[index].distance);
	}

	//Sort the Obstructed GNN
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());

	printf("\nFinal K-Obstructed Group Nearest Neighbor\n");
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=ognn_sorted[index].queryPoints[0];
			q[1]=ognn_sorted[index].queryPoints[1];

			printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],ognn_sorted[index].distance);
	}

*/
}





void OGNN_GNN::constructInitialVisGraph(VisibilityGraph* initialVisGraph) {

	vector<Obstacle*> obsList;
	ifstream iFile(VISGRAPHFILE);
	string line;
	Obstacle* obs;
	//ObstacleController* obsController= new ObstacleController();

	/* While there is still a line. */
	while (getline(iFile, line)) {
		/* Printing goes here. */
		//cout << line << endl;
		//Keeping all the polygon string in a vector so that next time we get which obstacles are already added in VisGraph
		obstacleString.push_back(line);
		obs = createObstacle(line);
		obsList.push_back(obs);
	}

	initialVisGraph->setObstacle(obsList);
	vg = new VisibilityGraphController(initialVisGraph);
	//Construct Vis Graph
	vector<Line*> visEdges = vg->constructVisGraph();


	iFile.close();

}



