//Added By Nusrat
#ifndef OGNN_GNN_H_
#define OGNN_GNN_H_

#include "../func/gendef.h"
#include "../heap/heap.h"
#include "../visGraph/VisibilityGraph.h"
#include "../linlist/linlist.h"
#include "../ognn_utility.h"
#include "../rtree/rtree.h"

#include <iostream>


#include <vector>
using namespace std;

class OGNN_GNN
{
public:
	void ognnMultiPointApproach(Point2D queryPoints[],int numOfQueryPoints,int k,double kNearestNeighbor[][2],RTree* rt_obstacle,RTree* rt_dataPoints);
	
	void constructInitialVisGraph(VisibilityGraph* initialVisGraph);
	
};


#endif /* OGNN_GNN_H_ */
