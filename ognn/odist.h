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
#include "../visGraph/pointHandler.h"

const int infinity = 1000000000; 



class ObstructedDistance
{
public:
	char *VISGRAPH_FILE ;
	VisibilityGraphController* vgController;
	vector<string> obstacleList;
	//function=0 sum , function=1 max
	double computeAggObstructedDistance(VisibilityGraph* initialVisGraph,
		float* p, Point2D queryPoints[],int numOfQueryPoints, RTree* rt_obstacle,int function);
	bool doesLineAndObstcaleIntersects(tLinestring ls,tPolygon p);

	//LineString creation 
	tLinestring createLS(double x1,double y1, double x2, double y2);
	//Polygon creation from string
	tPolygon createPolygon(double obstacle[5]);
	void writePolygonInFile(double obstacle[5]);
	void writePointAndQueryPointsInFile(float* p,Point2D queryPoints[],int numOfQueryPoints);
	void writePointInFile(float* p);
	void constructInitialVisGraph(VisibilityGraph* initialVisGraph);
	double computeObstructedDistance(VisibilityGraph* initialVisGraph,float* p, float* q,vector<int>& shortestPath);
	int drawAndWriteFileVisEdges(vector<Line*> visEdges);
	void replaceObsDist(std::vector < MyStruct >& dist_O_p_qi,float* q,double obsDist);
	bool checkIntersectionWithSP(float* q,vector< MyShortestPath>& shortestPath_p_qi,double obs[5],
		VisibilityGraph* initialVisGraph);
	void addOrReplaceSP(float* q,vector< MyShortestPath>& shortestPath_p_qi,vector<int>& shortestPath);
	void removeDataPointFromVG(VisibilityGraph* initialVisGraph,float* q);
	void removePointAndQueryPointsFromVisGraph(VisibilityGraph* initialVisGraph,
				float* p, Point2D queryPoints[],int numOfQueryPoints);
	bool visGraphContainsPoly(char buffer[1024]);
	void writeQueryPointsInFile(Point2D queryPoints[],int numOfQueryPoints);
	void addDataPointInVG(VisibilityGraph* initialVisGraph,float* q);

	long double visGraphConsTime;
	long double shortestPathCalcTime;

	ObstructedDistance(){
		visGraphConsTime=0.0;
		shortestPathCalcTime=0.0;
		VISGRAPH_FILE = "Datasets/visibilityGraphPolygons.txt";
	}

};


#endif /* ODIST_H_ */