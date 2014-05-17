//Added By Nusrat
#ifndef OGNNGNN_H_
#define OGNNGNN_H_

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
	void ognnUsingEGNN(Point2D queryPoints[],int numOfQueryPoints,int k,double kNearestNeighbor[][3],RTree* rt_obstacle,RTree* rt_dataPoints,int function);
	void ognnSumUsingNN(Point2D queryPoints[], int numOfQueryPoints,int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function);
	bool pointInsideTheIntersectionOfCircle(float* p, Point2D queryPoints[],int numOfQueryPoints,double radius);
	void ognnMaxUsingNN(Point2D queryPoints[], int numOfQueryPoints,int k, double kNearestNeighbor[][3], RTree* rt_obstacle,
		RTree* rt_dataPoints,int function);
	int totalNumberOfPRetrieved;
	
};


#endif /* OGNNGNN_H_ */
