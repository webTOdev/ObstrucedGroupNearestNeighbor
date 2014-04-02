//Added By Nusrat
#ifndef OGNN_H_
#define OGNN_H_

#include "../func/gendef.h"
#include "../heap/heap.h"
#include "../visGraph/VisibilityGraph.h"


#include <iostream>


#include <vector>
using namespace std;

class OGNN
{
public:
	void ognnMultiPointApproach(Point2D queryPoints[],int numOfQueryPoints,int k,double kNearestNeighbor[][2],RTree* rt_obstacle,RTree* rt_dataPoints);
	void onnMultiPointApproach(Point2D queryPoints,double kNearestNeighbor[][2],RTree* rt_obstacle,RTree* rt_dataPoints);

	double computeObstructedDistance(VisibilityGraph* initialVisGraph,Point2D p,Point2D q,RTree* rt_obstacle);
};


#endif /* OGNN_H_ */
