/*
 * PointHandler.cpp
 *
 *  Created on: Sep 12, 2013
 *      Author: nut
 */

#include "pointHandler.h"
#include <limits>

bool AreSame(double a, double b)
{
	double diff =1.0 / 100000 ;
    return  fabs(a - b) < diff ;
}
Point* searchPointByCoord(vector<Point*> pointList, double px, double py) {

	for (int i = 0; i < pointList.size(); i++) {
		if (AreSame(pointList[i]->x,px) && AreSame(pointList[i]->y,py)) {
			return pointList[i];
		}

	}
	return NULL;
}



Point * getPointById(vector<Point*> pointList, int id) {
	for (int i = 0; i < pointList.size(); i++) {
		if (pointList[i]->id == id) {
			return pointList[i];
		}

	}
	return NULL;
}
//-------------------------------------------------------------------------------
//  Distance Btw 2 Points
//-------------------------------------------------------------------------------
double distance(Point * a, Point * b) {
	return sqrt(pow(b->x - a->x, 2.0) + pow(b->y - a->y, 2.0));
}

