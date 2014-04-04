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

void OGNN::ognnMultiPointApproach(Point2D queryPoints[], int numOfQueryPoints,
		int k, double kNearestNeighbor[][2], RTree* rt_obstacle,
		RTree* rt_dataPoints) {

	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(queryPoints, k, kNearestNeighbor,
			numOfQueryPoints);
	//init a heap that stores the obstructed kGNN
	Heap *heap_kGNN_obstructed = new Heap();
	heap_kGNN_obstructed->init(rt_dataPoints->dimension);

	vector<Point2D*> dataPoints_temp, dataPoints_temp2;

}

void generateMBR(float* mbrRange, Point2D queryPoints, double distance) {
	mbrRange[0] = queryPoints[0] - distance;
	mbrRange[1] = queryPoints[0] + distance;
	mbrRange[2] = queryPoints[1] - distance;
	mbrRange[3] = queryPoints[0] + distance;
}

bool visGraphContainsPoly(string poly) {
	if (std::find(obstacleString.begin(), obstacleString.end(), poly)
			!= obstacleString.end()) {
		return true;
	}

	return false;
}

void OGNN::onnMultiPointApproach(Point2D queryPoints,
		double kNearestNeighbor[][2], RTree* rt_obstacle,
		RTree* rt_dataPoints) {

	double *kNN;
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

	double shortestPathDistance = initialVisGraph->findShortestPath(p[0], p[1],
			q[0], q[1]);
	printf("\nFirst Obstructed Distance found is %lf\n", shortestPathDistance);
	float* mbrRange = new float[4];
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

			snprintf(buffer, sizeof(buffer),
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
			}

			cur = res_list->get_next();
		}

		double newShortestPathDistance = initialVisGraph->findShortestPath(p[0],
								p[1], q[0], q[1]);
		printf("\nNew Obstructed Distance found is %lf", newShortestPathDistance);
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
		cout << line << endl;
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

SortedLinList* OGNN::obsInRange(RTree* rt_obstacle, float *mbr) {

	SortedLinList *res_list = new SortedLinList();
	rt_obstacle->rangeQuery(mbr, res_list);
	printf("\nRange Query returned %d entries\n", res_list->get_num());
	res_list->print();

	return res_list;

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

