//Added by Nusrat

#include "../rtree/rtree.h"
#include "ognn.h"


void OGNN::ognnMultiPointApproach(Point2D queryPoints[],int numOfQueryPoints,int k,double kNearestNeighbor[][2],RTree* rt_obstacle,RTree* rt_dataPoints){

	//kNearestNeighbour holds the kGNN Euclidean
	rt_dataPoints->Point_BFN_kGNNQ(queryPoints,k,kNearestNeighbor,numOfQueryPoints);
	//init a heap that stores the obstructed kGNN
	Heap *heap_kGNN_obstructed = new Heap();
	heap_kGNN_obstructed->init(rt_dataPoints->dimension);

	vector<Point2D*> dataPoints_temp,dataPoints_temp2;




}

void OGNN::onnMultiPointApproach(Point2D queryPoints,double kNearestNeighbor[][2],RTree* rt_obstacle,RTree* rt_dataPoints) {

	double *kNN;
	//kNearestNeighbour holds the kNN Euclidean
	rt_dataPoints->Point_BFN_NNQ(queryPoints, kNN);
	//init a heap that stores the obstructed kGNN
	Heap *heap_kNN_obstructed = new Heap();
	heap_kNN_obstructed->init(rt_dataPoints->dimension);

	vector<Point2D*> dataPoints_temp, dataPoints_temp2;

	std::cout << "Length of array = " << (sizeof(kNN)/sizeof(*kNN)) << std::endl;
}

double OGNN::computeObstructedDistance(VisibilityGraph* initialVisGraph,Point2D p,Point2D q,RTree* rt_obstacle){
	double shortestPathDistance = -1;
	return shortestPathDistance;
}






