//Added by Nusrat

#include "ognn.h"

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

class Line;
class Obstacle;

char *VISGRAPHFILE = "Datasets/visibilityGraphPolygons.txt";
vector<string> obstacleString;
VisibilityGraphController* vg;


void generateMBR(float* mbrRange, float* queryPoints, double distance) {
	mbrRange[0] = queryPoints[0] - distance;
	mbrRange[1] = queryPoints[0] + distance;
	mbrRange[2] = queryPoints[1] - distance;
	mbrRange[3] = queryPoints[1] + distance;
}

bool visGraphContainsPoly(string poly) {
	if (std::find(obstacleString.begin(), obstacleString.end(), poly)
			!= obstacleString.end()) {
		return true;
	}

	return false;
}

int drawAndWriteFileVisEdges(vector<Line*> visEdges) {
	//Remove existing test.txt file
	if (remove("test.txt") != 0)
		perror("Error deleting file");
	else
		puts("File successfully deleted");
	/*FILE *fp;
	fp=fopen("test.txt", "w");
	fclose(fp);
*/
	//Dijkstra algorithm needs to create a vector of size n , where n is the id of the vertex
	int maxNumberOfVertex=-1;
	for (int i = 0; i < visEdges.size(); i++) {
		fileWrite(visEdges[i]->a, visEdges[i]->b);
		if(visEdges[i]->a->id > maxNumberOfVertex){
			maxNumberOfVertex=visEdges[i]->a->id;
		}else
			if(visEdges[i]->b->id > maxNumberOfVertex){
			maxNumberOfVertex=visEdges[i]->b->id;
		}
	}
	return maxNumberOfVertex+1;
}

void removeDataPointFromVG(VisibilityGraph* initialVisGraph,float* q){

	char buffer[1024];
	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f%s", "polygon((",
				q[0],q[1], q[0],q[1], "))");
	initialVisGraph = vg->removeDataPointFromVisGraph(initialVisGraph,initialVisGraph->searchObsWithString(buffer));
}

void addDataPointFromVG(VisibilityGraph* initialVisGraph,float* q){

	char buffer[1024];
	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f%s", "polygon((",
				q[0],q[1], q[0],q[1], "))");
	Obstacle* obs = createObstacle(buffer);
	initialVisGraph = vg->addNewObstacleForIncrementalVisGraph(initialVisGraph, obs);
}

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

double OGNN::computeObsDistanceForOtherQueryPoints(float* kNN_point,std::vector < MyStruct >& queryPoints_sorted,VisibilityGraph* initialVisGraph,
												   double oldAggregateObstructedDistance,int numOfQueryPoints,RTree* rt_obstacle){
		double aggregateObstructedDistance=oldAggregateObstructedDistance;
		//first element in vector is already considered so start loop from j=1
		for(int j=1;j<numOfQueryPoints;j++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			float *q = new float[2];
			q[0]=queryPoints_sorted[j].queryPoints[0];
			q[1]=queryPoints_sorted[j].queryPoints[1];
			//Add the query point from Visibility graph
			addDataPointFromVG(initialVisGraph,q);		

			double obstructedDistance = computeObstructedDistance(initialVisGraph, kNN_point, q, rt_obstacle,obstacleString);
			printf("\nObstructed Distance is %lf\n", obstructedDistance);
			//Remove this query point from VisibilityGraph as computation is finished
			removeDataPointFromVG(initialVisGraph,q);

			aggregateObstructedDistance+=obstructedDistance;
			delete q;
	}
		return aggregateObstructedDistance;
}
double OGNN::computeObsDistanceForTheFarthestQueryPoint(float* kNN_point,std::vector < MyStruct >& queryPoints_sorted, RTree* rt_obstacle,int i,Point2D queryPoints[],
												VisibilityGraph* initialVisGraph  ){
//Imagine a rectangle over knn data point with distance
		float* mbrRange = new float[4];
		generateMBR(mbrRange, kNN_point, queryPoints_sorted[0].distance);
		//Finding all the obstacle within the MBR with distance
		SortedLinList* res_list = obsInRange(rt_obstacle, mbrRange);
		float *q=new float[2];
		q[0]=queryPoints_sorted[0].queryPoints[0];
		q[1]=queryPoints_sorted[0].queryPoints[1];
		//Write then in a file so that Visibility graph can be constructed
		writeDataPointPolygonInFile(kNN_point, q, res_list);
		//Only construct the graph at the beginning
		if(i==0){			
			//Construct Initial Visibility Graph and calculate obstructed path distance
			constructInitialVisGraph(initialVisGraph);

		}//Otherwise only update the graph
		else{
			updateVisGraphFromFile(initialVisGraph);		
		}
		double obstructedDistance = computeObstructedDistance(initialVisGraph, kNN_point, q, rt_obstacle,obstacleString);
		//printf("\nObstructed Distance is %lf\n", obstructedDistance);
		//Remove this query point from VisibilityGraph as computation is finished
		removeDataPointFromVG(initialVisGraph,q);
		delete mbrRange;
		delete q;

		return  obstructedDistance;
}

void OGNN::ognnMultiPointApproach(Point2D queryPoints[], int numOfQueryPoints,
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
		double obstructedDistance = computeObsDistanceForTheFarthestQueryPoint(kNN_point,queryPoints_sorted,rt_obstacle,i,queryPoints,initialVisGraph  );
		
		//Initialize the aggregate distance with the distance between p and first q 
		double aggregateObstructedDistance=obstructedDistance;	
		//first element in vector is already considered so start loop from j=1
		aggregateObstructedDistance=computeObsDistanceForOtherQueryPoints(kNN_point,queryPoints_sorted,initialVisGraph,
												   aggregateObstructedDistance,numOfQueryPoints,rt_obstacle);
		
		//printf("\nAggregate Obstructed Distance of data point p(%f,%f) is %lf\n",kNN_point[0],kNN_point[1], aggregateObstructedDistance);
		ognn_sorted.push_back(MyStruct(aggregateObstructedDistance, kNN_point));
		//Now remove the Data point from Vis graph as obstructed group distance calculation is finished
		removeDataPointFromVG(initialVisGraph,kNN_point);
	}
	//Sort the Euclidean GNN
	std::sort(egnn_sorted.begin(), egnn_sorted.end(), less_than_key());
	//Sort the Obstructed GNN
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());
	printf("\nInitial K-Euclidean Group Nearest Neighbor\n");
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


}


double OGNN::computeObstructedDistance(VisibilityGraph* initialVisGraph,
		float* p, float* q, RTree* rt_obstacle,
		vector<string> obstacleString) {
	int maxVertexNum = drawAndWriteFileVisEdges(initialVisGraph->edges);
	double shortestPathDistance = initialVisGraph->findShortestPath(p[0], p[1],
			q[0], q[1],maxVertexNum);
	//printf("\nFirst Obstructed Distance found is %lf\n", shortestPathDistance);
	float* mbrRange = new float[4];
	bool newObsDiscovered=false;
	while (1) {
		generateMBR(mbrRange, q, shortestPathDistance);
		//Finding all the obstacle within the MBR with distance
		SortedLinList* res_list = obsInRange(rt_obstacle, mbrRange);

		//Discover new Obstacles and add in Vis Graph
		Linkable *cur = res_list->get_first();
		//If no new Obstacle discovered break the loop and this is the final shortest obstacle path Distance
		if (cur == NULL)
			break;
		vector<string> newObstcles;
		while (cur != NULL) {
			char buffer[1024];

			_snprintf(buffer, sizeof(buffer),
					"%s%f %f,%f %f,%f %f,%f %f,%f %f%s", "polygon((",
					cur->bounces[0], cur->bounces[2], cur->bounces[1],
					cur->bounces[2], cur->bounces[1], cur->bounces[3],
					cur->bounces[0], cur->bounces[3], cur->bounces[0],
					cur->bounces[2], "))");

			//If the Graph does not contain this obstacles then add it
			if (!visGraphContainsPoly(buffer)) {
				printf("\nNew Polygon discovered %s\n", buffer);
				newObstcles.push_back(buffer);
				//Update Vis Graph with this new Polygon
				Obstacle* obs = createObstacle(buffer);
				initialVisGraph = vg->addNewObstacleForIncrementalVisGraph(
						initialVisGraph, obs);
				obstacleString.insert(obstacleString.end(),newObstcles.begin(),newObstcles.end());
				newObsDiscovered=true;
			}

			cur = res_list->get_next();
		}

		double newShortestPathDistance=shortestPathDistance;
		if(newObsDiscovered){
			int maxVertexNum = drawAndWriteFileVisEdges(initialVisGraph->edges);
			double newShortestPathDistance = initialVisGraph->findShortestPath(p[0],
								p[1], q[0], q[1],maxVertexNum);
			printf("\nNew Obstructed Distance found is %lf", newShortestPathDistance);
			newObsDiscovered=false;
		}
		//If in two iteration shortest path distance does not change then this is the final shortest obstacle path Distance
		if(newShortestPathDistance==shortestPathDistance){
			printf("\nObstructed distance remains same in two consecutive iterations so breaking while loop.");
			break;
		}
		else
			shortestPathDistance=newShortestPathDistance;

	}
/*
	for(int i=0;i<obstacleString.size();i++){
		printf("\nVector has %s",obstacleString[i].c_str());
	}
*/
	return shortestPathDistance;
}

void OGNN::constructInitialVisGraph(VisibilityGraph* initialVisGraph) {

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

void OGNN::updateVisGraphFromFile(VisibilityGraph* initialVisGraph) {

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
		initialVisGraph = vg->addNewObstacleForIncrementalVisGraph(
						initialVisGraph, obs);
		
	}

	
	iFile.close();

}

SortedLinList* OGNN::obsInRange(RTree* rt_obstacle, float *mbr) {

	SortedLinList *res_list = new SortedLinList();
	rt_obstacle->rangeQuery(mbr, res_list);
	//printf("\n%f %f %f %f",mbr[0],mbr[1],mbr[2],mbr[3]);
	printf("\nRange Query returned %d entries\n", res_list->get_num());
	//res_list->print();

	return res_list;

}

void OGNN::writePointInFile(double* p){
	FILE * input = fopen(VISGRAPHFILE, "w");
	if (input == NULL) {
		printf("Cannot open file %s", VISGRAPHFILE);
	}
	fprintf(input, "polygon((");
	fprintf(input, "%f %f,%f %f", p[0], p[1], p[0], p[1]);
	fprintf(input, "))\n");

	fclose(input);
}

void OGNN::writeDataPointPolygonInFile(float* p, float* q,
		SortedLinList *res_list) {
	
	FILE * input = fopen(VISGRAPHFILE, "w");
	if (input == NULL) {
		printf("Cannot open file %s", VISGRAPHFILE);
	}

	fprintf(input, "polygon((");
	fprintf(input, "%f %f,%f %f", p[0], p[1], p[0], p[1]);
	fprintf(input, "))");
	fprintf(input, "\npolygon((");
	fprintf(input, "%f %f,%f %f", q[0], q[1], q[0], q[1]);
	fprintf(input, "))");



	Linkable *cur = res_list->get_first();

	while (cur != NULL) {
		/*printf("%f %f %f %f\n", cur->bounces[0], cur->bounces[1],
		 cur->bounces[2], cur->bounces[3]);*/
		fprintf(input, "\npolygon((");
		fprintf(input, "%f %f,", cur->bounces[0], cur->bounces[2]);
		fprintf(input, "%f %f,", cur->bounces[1], cur->bounces[2]);
		fprintf(input, "%f %f,", cur->bounces[1], cur->bounces[3]);
		fprintf(input, "%f %f,", cur->bounces[0], cur->bounces[3]);
		fprintf(input, "%f %f", cur->bounces[0], cur->bounces[2]);
		fprintf(input, "))");

		cur = res_list->get_next();
	}

	fclose(input);
}

