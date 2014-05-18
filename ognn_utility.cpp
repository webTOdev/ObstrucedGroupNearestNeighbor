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
	//printf("\nDistance between p(%f,%f) and q(%f,%f) is %f",p[0],p[1],q[0],q[1],dist);
	return dist;
}

void centroidOfQ(Point2D queryPoints[], int numOfQueryPoints, Point2D centroid[]){
	float x=0.0;
	float y=0.0;
	for(int i=0;i<numOfQueryPoints;i++){
		x+=queryPoints[i][0];
		y+=queryPoints[i][1];
	}
	centroid[0][0]=x/numOfQueryPoints;
	centroid[0][1]=y/numOfQueryPoints;

}



