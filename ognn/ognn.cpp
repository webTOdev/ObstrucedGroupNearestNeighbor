//Added by Nusrat

#include "../rtree/rtree.h"
#include "ognn.h"
#include <math.h>
#include "../linlist/linlist.h"
#include "../ognn_utility.h"
#include <fstream>
#include "../visGraph/obstacleController.h"
#include "../visGraph/VisibilityGraphController.h"
#include "../visGraph/VisibilityGraph.h"

char *VISGRAPHFILE = "Datasets/visibilityGraphPolygons.txt";

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
	//Imagine a rectangle over queryPoint with distance
	float* mbrRange = new float[4];
	mbrRange[0] = queryPoints[0] - distance;
	mbrRange[1] = queryPoints[0] + distance;
	mbrRange[2] = queryPoints[1] - distance;
	mbrRange[3] = queryPoints[0] + distance;
	//Finding all the obstacle within the MBR with distance
	SortedLinList* res_list = obsInRange(rt_obstacle, mbrRange);
	//Write then in a file so that Visibility graph can be constructed
	writeDataPointPolygonInFile(firstKNN, queryPoints, res_list);
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	//Construct Initial Visibility Graph and calculate obstructed path distance
	constructInitialVisGraph(initialVisGraph);
	computeObstructedDistance(initialVisGraph,firstKNN,queryPoints,rt_obstacle);

}

double OGNN::computeObstructedDistance(VisibilityGraph* initialVisGraph,
		Point2D p, Point2D q, RTree* rt_obstacle) {

	double shortestPathDistance = initialVisGraph->findShortestPath(p[0],p[1],q[0],q[1]);
	printf("\nFirst Obstructed Distance found is %lf\n",shortestPathDistance);
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
		obs = createObstacle(line);
		obsList.push_back(obs);
	}

	initialVisGraph->setObstacle(obsList);
	VisibilityGraphController* vg = new VisibilityGraphController(initialVisGraph);
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

