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


void generateMBR(float* mbrRange, Point2D queryPoints, double distance) {
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

void OGNN::ognnMultiPointApproach(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][2], RTree* rt_obstacle,
		RTree* rt_dataPoints) {

	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(queryPoints, k, kNearestNeighbor,
			numOfQueryPoints);
	printf("k- Group Nearest Neighbor of (%f,%f),(%f,%f),(%f,%f) is (%f,%f),(%f,%f),(%f,%f)\n", queryPoints[0][0],queryPoints[0][1],
			queryPoints[1][0],queryPoints[1][1],queryPoints[2][0],queryPoints[2][1],
			kNearestNeighbor[0][0], kNearestNeighbor[0][1],kNearestNeighbor[1][0], kNearestNeighbor[1][1],kNearestNeighbor[2][0], kNearestNeighbor[2][1]);
	//init a heap that stores the obstructed kGNN
	Heap *heap_kGNN_obstructed = new Heap();
	heap_kGNN_obstructed->init(rt_dataPoints->dimension);

	vector<Point2D*> dataPoints_temp, dataPoints_temp2;
	std::vector < MyStruct > ognn_sorted;
	
	float *kNN_point,*q;
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	for(int i=0;i<k;i++){
		std::vector < MyStruct > queryPoints_sorted;
		kNN_point = new float[2];
		kNN_point[0] = kNearestNeighbor[i][0];
		kNN_point[1] = kNearestNeighbor[i][1];
		double aggregateEuclideanDist=0;
		for(int j=0;j<numOfQueryPoints;j++){
			q = new float[2];
			q[0]=queryPoints[j][0];
			q[1]=queryPoints[j][1];
			double distance = getDistanceBetweenTwoPoints(kNN_point, q);
			aggregateEuclideanDist+=distance;
			//printf("\nEuclidean Distance is %lf for q %f,%f", distance,q[0],q[1]);
			queryPoints_sorted.push_back(MyStruct(distance, q));
		}
		printf("\nAggregate Euclidean Distance is %lf for p %f,%f\n", aggregateEuclideanDist,kNN_point[0],kNN_point[1]);

		std::sort(queryPoints_sorted.begin(), queryPoints_sorted.end(), more_than_key());
		//Imagine a rectangle over knn data point with distance
		float* mbrRange = new float[4];

		generateMBR(mbrRange, kNN_point, queryPoints_sorted[0].distance);
		//Finding all the obstacle within the MBR with distance
		SortedLinList* res_list = obsInRange(rt_obstacle, mbrRange);
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
		printf("\nObstructed Distance is %lf\n", obstructedDistance);
		//Remove this query point from VisibilityGraph as computation is finished
		removeDataPointFromVG(initialVisGraph,queryPoints_sorted[0].queryPoints);
		//Initialize the aggregate distance with the distance between p and first q 
		double aggregateObstructedDistance=obstructedDistance;	
		//first element in vector is already considered so start loop from j=1
		for(int j=1;j<numOfQueryPoints;j++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=queryPoints_sorted[j].queryPoints[0];
			q[1]=queryPoints_sorted[j].queryPoints[1];
			//Add the query point from Visibility graph
			addDataPointFromVG(initialVisGraph,queryPoints_sorted[j].queryPoints);		

			double obstructedDistance = computeObstructedDistance(initialVisGraph, kNN_point, q, rt_obstacle,obstacleString);
			printf("\nObstructed Distance is %lf\n", obstructedDistance);
			//Remove this query point from VisibilityGraph as computation is finished
			removeDataPointFromVG(initialVisGraph,queryPoints_sorted[j].queryPoints);

			aggregateObstructedDistance+=obstructedDistance;

		}
		printf("\nAggregate Obstructed Distance is of data point p(%f,%f) is %lf\n",kNN_point[0],kNN_point[1], aggregateObstructedDistance);
		ognn_sorted.push_back(MyStruct(aggregateObstructedDistance, kNN_point));
		removeDataPointFromVG(initialVisGraph,kNN_point);
	}
	std::sort(ognn_sorted.begin(), ognn_sorted.end(), less_than_key());
//	printf("k- Obs Group Nearest Neighbor of (%f,%f),(%f,%f),(%f,%f) is ",queryPoints[0][0],queryPoints[0][1],
//			queryPoints[1][0],queryPoints[1][1],queryPoints[2][0],queryPoints[2][1]);
		printf("\nk- Obs Group Nearest Neighbor\n");
	for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			q = new float[2];
			q[0]=ognn_sorted[index].queryPoints[0];
			q[1]=ognn_sorted[index].queryPoints[1];

			printf("(%f,%f) distance %lf \t,",q[0], q[1],ognn_sorted[index].distance);

		}
		
		
		
	


}

void OGNN::onnMultiPointApproach(Point2D queryPoints,
		double kNearestNeighbor[][2], RTree* rt_obstacle,
		RTree* rt_dataPoints) {

	double *kNN = new double[3];
	//kNearestNeighbour holds the kNN Euclidean
	rt_dataPoints->Point_BFN_NNQ(queryPoints, kNN);
	//init a heap that stores the obstructed kGNN
	Heap *heap_kNN_obstructed = new Heap();
	heap_kNN_obstructed->init(rt_dataPoints->dimension);

	vector<Point2D*> dataPoints_temp, dataPoints_temp2;
	float *firstKNN = new float[2];
	firstKNN[0] = kNN[0];
	firstKNN[1] = kNN[1];
	double distance = getDistanceBetweenTwoPoints(firstKNN, queryPoints);
	printf("\nEuclidean Distance is %lf", distance);
	//Imagine a rectangle over queryPoint with distance
	float* mbrRange = new float[4];
	generateMBR(mbrRange, queryPoints, distance);
	//Finding all the obstacle within the MBR with distance
	SortedLinList* res_list = obsInRange(rt_obstacle, mbrRange);
	//Write then in a file so that Visibility graph can be constructed
	writeDataPointPolygonInFile(firstKNN, queryPoints, res_list);
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	//Construct Initial Visibility Graph and calculate obstructed path distance
	constructInitialVisGraph(initialVisGraph);
	double obstructedDistance = computeObstructedDistance(initialVisGraph, firstKNN, queryPoints,rt_obstacle, obstacleString);
	printf("\nObstructed Distance is %lf\n", obstructedDistance);

}

double OGNN::computeObstructedDistance(VisibilityGraph* initialVisGraph,
		Point2D p, Point2D q, RTree* rt_obstacle,
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

void OGNN::writePointInFile(Point2D p){
	FILE * input = fopen(VISGRAPHFILE, "w");
	if (input == NULL) {
		printf("Cannot open file %s", VISGRAPHFILE);
	}
	fprintf(input, "polygon((");
	fprintf(input, "%f %f,%f %f", p[0], p[1], p[0], p[1]);
	fprintf(input, "))\n");

	fclose(input);
}

void OGNN::writeDataPointPolygonInFile(Point2D p, Point2D q,
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

