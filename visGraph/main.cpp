/*
 * main.cpp
 *
 *  Created on: Oct 2, 2013
 *      Author: nut
 */
#include <iostream>
#include "obstacles.h"
#include "boostHelper.h"
#include "obstacleController.h"
#include "drawing.h"
#include "utility.h"
#include <vector>
#include "VisibilityGraphController.h"
#include "VisibilityGraph.h"
#include "graphutility.h"
#include "dijkstra.h"
#include <time.h>

clock_t startTime;
void findShortestPath(VisibilityGraph* vg, double sourceX, double sourceY,
		double destX, double destY);
void drawObs(Obstacle* obs, Point* ori);
void drawAndWriteFileVisEdges(vector<Line*> visEdges);


int mainst() {

	startTime = clock();
	vector<Obstacle*> obsList;
	//Create a list of test obstacle
	//ObstacleController* obsController= new ObstacleController();
	Obstacle* obs = createObstacle(
			"polygon((40 140,40 300,100 280,120 200,100 140,40 140))");
	obsList.push_back(obs);

	obs = createObstacle("polygon((200 40,280 40,320 200,180 190,200 40))");
	obsList.push_back(obs);

	obs =createObstacle("polygon((30 330,180 350,120 450,30 330))");
	obsList.push_back(obs);

	obs = createObstacle("polygon((400 320,520 320,500 520,400 530,400 320))");
	obsList.push_back(obs);

	obs = createObstacle("polygon((40 20,40 20))"); //For data point
	obsList.push_back(obs);

	obs =createObstacle("polygon((440 20,440 20))"); //For data point
	obsList.push_back(obs);

	obs = createObstacle("polygon((600 600,600 600))"); //For data point
	obsList.push_back(obs);

	//obs=createObstacle("polygon((250 400,320 400,320 500,250 400))");
	//obsList.push_back(obs);

	/*Obstacle* obs=createObstacle("polygon((200 20,250 20,250 200,200 200,200 20))");
	 obsList.push_back(obs);

	 obs=createObstacle("polygon((100 150,150 150,150 50,100 50,100 150))");
	 obsList.push_back(obs);
	 */

	//Create the initial Visibility graph
	VisibilityGraph* visGraph = new VisibilityGraph(obsList);
	VisibilityGraphController* vg = new VisibilityGraphController(visGraph);
	vector<Line*> visEdges = vg->constructVisGraph();
	vector<Line*> obsSide = visGraph->obsSides;

	//Incremental Visibility Graph , I am adding a new obstacle in the existing visibility  graph
	//obs=createObstacle("polygon((650 400,720 400,720 550,650 400))");//For Obs
	obs = createObstacle("polygon((250 400,320 400,320 500,250 400))");	//For Obs
	visGraph = vg->addNewObstacleForIncrementalVisGraph(visGraph,obs);
	obs=createObstacle("polygon((650 400,720 400,720 550,650 400))");//For Obs
	visGraph = vg->addNewObstacleForIncrementalVisGraph(visGraph,obs);
	visGraph = vg->removeDataPointFromVisGraph(visGraph,visGraph->searchObsWithString("polygon((600 600,600 600))"));

	//Necessary Drawing
	visGraph->print();
	visEdges = visGraph->edges;
	obsSide = visGraph->obsSides;
	vector<Obstacle*> visObs = visGraph->obstacles;
	//Draw the Obstacle Sides id in the image
		for (int i = 0; i < obsSide.size(); i++) {
			Line* l = obsSide[i];
			drawText(((l->a->x) + (l->b->x)) / 2, ((l->a->y) + (l->b->y)) / 2,
					l->id, WHITE);
			//drawEdge(l,WHITE);
		}

	//Write the edges in a file so that Dijkstra can use it
	drawAndWriteFileVisEdges(visEdges);

	//This point is not used , i created it for testign purpose
	Point* ori = new Point(440, 120);
	for (int i = 0; i < visObs.size(); i++) {
		//obsList[i]->print();
		drawObs(visObs[i], ori);

	}

	//Find the shortest path from two points using Dijkstra
	//findShortestPath(visGraph, 40, 140, 600, 600);
	findShortestPath(visGraph, 40, 20, 720, 550);
	//findShortestPath(visGraph,400,530 , 650, 400);
	/* Calculate the time */
	printf("Time elapsed: %f s\n",
			((double) clock() - startTime) / CLOCKS_PER_SEC);

	displayImage();


	return 0;
}

void drawAndWriteFileVisEdges(vector<Line*> visEdges) {
	//Remove existing test.txt file
	if (remove("test.txt") != 0)
		perror("Error deleting file");
	else
		puts("File successfully deleted");

	for (int i = 0; i < visEdges.size(); i++) {
		drawEdge(visEdges[i], BLUE);
		fileWrite(visEdges[i]->a, visEdges[i]->b);
	}
}
void drawObs(Obstacle* o, Point* ori) {

	//ObstacleController* obsController= new ObstacleController();
	vector<Point*> vertexList = getVertices(o);
	int size = vertexList.size();
	CImg<double> points(size, 2);
	double ps[size * 2];
	int i = 0;
	for (std::vector<Point*>::iterator it = vertexList.begin();
			it < vertexList.end(); ++it) {
		ps[i] = (*it)->x;
		i++;
		ps[i] = (*it)->y;
		i++;
		drawCircle((*it)->x, (*it)->y, 1.5, WHITE);
		drawText((*it)->x, (*it)->y, (*it)->id, RED);
		//Draw Lines for angle sorting
		// drawLine((*it)->x,(*it)->y,ori->x,ori->y);
	}

	double *iterator = ps;
	cimg_forX(points,i)
	{
		points(i, 0) = *(iterator++)*scale;
		points(i, 1) = *(iterator++)*scale;
	}
	drawPolygon(points);

}

void findShortestPath(VisibilityGraph* vg, double sourceX, double sourceY,
		double destX, double destY) {
	int numOfPoints = vg->nodes.size();
	int numOfEdges = vg->edges.size();
	int sourcePointId = searchPointByCoord(vg->nodes, sourceX, sourceY)->id;
	int destPointId = searchPointByCoord(vg->nodes, destX, destY)->id;
	printf("\nFinding shortest path from %d -> %d\n", sourcePointId,
			destPointId, numOfEdges, numOfPoints);
	Point* start;
	Point* goal;
	initiateDijkstra(numOfPoints, numOfEdges, false, sourcePointId,
			destPointId);
	int *shortestPath = getShortestPath();
	int i = 0;
	//Print the Shortest Path
	printf("The Shortest Path is :");
	while (shortestPath[i] != -1) {
		printf("%d ", shortestPath[i]);

		if (shortestPath[i + 1] != -1) {
			start = getPointById(vg->nodes, shortestPath[i]);
			goal = getPointById(vg->nodes, shortestPath[i + 1]);
			// Visualize:
			drawCircle(start->x, start->y, 2, GREEN);
			drawCircle(goal->x, goal->y, 2, GREEN);
			drawLine(start->x, start->y, goal->x, goal->y, GREEN);
		} else {
			goal = getPointById(vg->nodes, destPointId);
			writeText((goal->x - 5), (goal->y + 15), "Dest", GREEN);
		}

		if (i == 0) {
			writeText((start->x - 5), (start->y + 15), "Source", GREEN);
		}

		i++;
	}
}

