/*
 * obstacleHelper.h
 *
 *  Created on: Oct 2, 2013
 *      Author: nut
 */

#ifndef OBSTACLEHELPER_H_
#define OBSTACLEHELPER_H_

#include "boostHelper.h"
#include <string>
#include "obstacles.h"




Obstacle* createObstacle(string str);
vector<Point*> getVertices(Obstacle* o);
template<typename Vertex>
void get_coordinates(Vertex const& p);
template<typename Edge>
void get_segments(Edge const& s);






#endif /* OBSTACLEHELPER_H_ */
