//Added By Nusrat
#ifndef OGNN_H_
#define OGNN_H_

#include "../func/gendef.h"
#include "../heap/heap.h"

#include <vector>
using namespace std;

class OGNN
{
public:
	void ognnMultiPointApproach(Point2D o[],int k,double _rslt[][2],RTree* rt_obstacle,RTree* rt_dataPoints);
};


#endif /* OGNN_H_ */
