#ifndef ODIST_H_
#define ODIST_H_

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
#include "../linlist/linlist.h"
#include "../ognn_utility.h"
#include "../rtree/rtree.h"

const int infinity = 1000000000; 



class ObstructedDistance
{
public:
	double computeAggObstructedDistance(VisibilityGraph* initialVisGraph,
		float* p, Point2D queryPoints[],int numOfQueryPoints, RTree* rt_obstacle,
		vector<string> obstacleString);
	bool doesLineAndObstcaleIntersects(tLinestring ls,tPolygon p);

	//LineString creation 
	tLinestring createLS(double x1,double y1, double x2, double y2);
	//Polygon creation from string
	tPolygon createPolygon(double obstacle[5]);
	void writePolygonInFile(double obstacle[5]);
	void writePointAndQueryPointsInFile(float* p,Point2D queryPoints[],int numOfQueryPoints);
	void writePointInFile(float* p);
	void constructInitialVisGraph(VisibilityGraph* initialVisGraph);
	double computeObstructedDistance(VisibilityGraph* initialVisGraph,float* p, float* q,vector<int> shortestPath);
	int drawAndWriteFileVisEdges(vector<Line*> visEdges);
	void replaceObsDist(std::vector < MyStruct >& dist_O_p_qi,float* q,double obsDist);


};


#endif /* ODIST_H_ */