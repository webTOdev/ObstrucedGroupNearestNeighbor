#ifndef ODISTCENTROID_H_
#define ODISTCENTROID_H_

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

const int infinty = 1000000000; 



class ObstructedDistanceCentroid
{
public:
	long double visGraphConsTime;
	long double shortestPathCalcTime;
	char *VISGRAPH_FILE ;
	VisibilityGraphController* vgController;
	vector<string> obstacleList;
	std::vector< vector < MyStruct > > L_R;
	std::vector< vector < Point* > > L_C;
	std::vector< vector < Point* > > L_N;

	double computeAggObstructedDistance(VisibilityGraph* initialVisGraph,
		float* p, Point2D queryPoints[],int numOfQueryPoints, RTree* rt_obstacle,int function);
	void initialize(std::vector<MyStruct>& dist_O_p_qi,Point2D queryPoints[],int numOfQueryPoints);
	bool visGraphContainsPoly(char buffer[1024]);
	void addNewObstacleInVisGraph(double* obs,VisibilityGraph* initialVisGraph);
	void writeQueryPointsInFile(Point2D queryPoints[],int numOfQueryPoints);
	void constructInitialVisGraph(VisibilityGraph* initialVisGraph) ;
	void writePointInFile(float* p);
	bool doesLineAndObstcaleIntersects(tLinestring ls,tPolygon p);
	tLinestring createLS(double x1,double y1, double x2, double y2);
	bool isVisible(float *v,float* p,VisibilityGraph* initialVisGraph);
	void relax(float *v,float *q,float* p,VisibilityGraph* initialVisGraph,std::vector<MyStruct>& dist_O_p_qi,double dist_o_v_q);
	bool isRealDistanceFor_qFound(float* q,float* p,std::vector<MyStruct>& dist_O_p_qi,double threshold,int q_index);
	void addVerticesOfObsInRangeInLc(vector<string>& obsInRange,VisibilityGraph* initialVisGraph,int q_index);
	double computeObstructedDistance(VisibilityGraph* initialVisGraph,float* p, float* q,vector<int>& shortestPath);
	int drawAndWriteFileVisEdges(vector<Line*> visEdges);


	ObstructedDistanceCentroid(){
		visGraphConsTime=0.0;
		shortestPathCalcTime=0.0;
		VISGRAPH_FILE = "Datasets/visibilityGraphPolygons.txt";
	}
};


#endif /* ODISTCENTROID_H_ */