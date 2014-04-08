//Added By Nusrat
#ifndef OGNN_H_
#define OGNN_H_

#include "../func/gendef.h"
#include "../heap/heap.h"
#include "../visGraph/VisibilityGraph.h"
#include "../linlist/linlist.h"
#include "../ognn_utility.h"
#include "../rtree/rtree.h"

#include <iostream>


#include <vector>
using namespace std;

class OGNN
{
public:
	void ognnMultiPointApproach(Point2D queryPoints[],int numOfQueryPoints,int k,double kNearestNeighbor[][2],RTree* rt_obstacle,RTree* rt_dataPoints);
	

	double computeObstructedDistance(VisibilityGraph* initialVisGraph,float* p,float* q,RTree* rt_obstacle,vector<string> obstacleString);
	SortedLinList* obsInRange(RTree* rt_obstacle,float *mbr);
	void writeDataPointPolygonInFile(float* p,float* q,SortedLinList *res_list);
	void constructInitialVisGraph(VisibilityGraph* initialVisGraph);
	void writePointInFile(double* p);
	void updateVisGraphFromFile(VisibilityGraph* initialVisGraph);
	double computeEuclideanGroupDistance(float* kNN_point,std::vector < MyStruct >& queryPoints_sorted,int numOfQueryPoints,Point2D queryPoints[] );
	double computeObsDistanceForTheFarthestQueryPoint(float* kNN_point,std::vector < MyStruct >& queryPoints_sorted, RTree* rt_obstacle,int i,Point2D queryPoints[],
												VisibilityGraph* initialVisGraph  );
	double computeObsDistanceForOtherQueryPoints(float* kNN_point,std::vector < MyStruct >& queryPoints_sorted,VisibilityGraph* initialVisGraph,
												   double oldAggregateObstructedDistance,int numOfQueryPoints,RTree* rt_obstacle);
};


#endif /* OGNN_H_ */
