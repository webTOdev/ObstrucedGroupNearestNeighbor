/*
 * ognn_utility.cpp


 *
 *  Created on: Apr 2, 2014
 *      Author: nut
 */

#include "func/gendef.h"
#include <math.h>
double getDistanceBetweenTwoPoints(Point2D p,Point2D q){

	double dist = sqrt((p[0]-q[0])*(p[0]-q[0])+(p[1]-q[1])*(p[1]-q[1]));
	printf("\nDistance between (%f,%f) and (%f,%f) is %f",p[0],p[1],q[0],q[1],dist);
	return dist;
}


