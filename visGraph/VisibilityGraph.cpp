/*
 * VisibilityGraph.cpp
 *
 *  Created on: Oct 4, 2013
 *      Author: nut
 */

#include "VisibilityGraph.h"
#include "pointHandler.h"
#include "dijkstra.h"
int MAX=100000;
bool lineHasPoint(tPoint p,tLinestring l);
VisibilityGraph::VisibilityGraph() {
	// TODO Auto-generated constructor stub

}

void VisibilityGraph::setObstacle(vector<Obstacle*> o){

	obstacles=o;
		for(int i=0;i<MAX;i++){
			pointVSPoints[i][0]=-1;
			pointVSPoints[i][1]=-1;
			pointVSEdges[i][0]=-1;
			pointVSEdges[i][1]=-1;

		}
		for(int i=0;i<o.size();i++){
			vector<Line*> es=o[i]->getEdges();
			obsSides.insert(obsSides.end(),es.begin(),es.end());
			vector<Point*> ps = o[i]->getVertices();
			nodes.insert(nodes.end(),ps.begin(),ps.end());

			for(int j=0;j<es.size();j++){
				Line* l=es[j];
				if(pointVSPoints[l->a->id][0]==-1){
					pointVSPoints[l->a->id][0]=l->b->id;
					pointVSEdges[l->a->id][0]=l->id;
				}
				else{
					pointVSPoints[l->a->id][1]=l->b->id;
					pointVSEdges[l->a->id][1]=l->id;
				}
				if(pointVSPoints[l->b->id][0]==-1){
					pointVSPoints[l->b->id][0]=l->a->id;
					pointVSEdges[l->b->id][0]=l->id;
				}
				else{
					pointVSPoints[l->b->id][1]=l->a->id;
					pointVSEdges[l->b->id][1]=l->id;
				}
			}


		}
}

VisibilityGraph::VisibilityGraph(vector<Obstacle*> o){
	setObstacle(o);
}


Obstacle* VisibilityGraph::addObstacle(Obstacle* o){
	obstacles.push_back(o);


		vector<Line*> es=o->getEdges();
		obsSides.insert(obsSides.end(),es.begin(),es.end());
		vector<Point*> ps = o->getVertices();
		nodes.insert(nodes.end(),ps.begin(),ps.end());

		for(int j=0;j<es.size();j++){
			Line* l=es[j];
			if(pointVSPoints[l->a->id][0]==-1){
				pointVSPoints[l->a->id][0]=l->b->id;
				pointVSEdges[l->a->id][0]=l->id;
			}
			else{
				pointVSPoints[l->a->id][1]=l->b->id;
				pointVSEdges[l->a->id][1]=l->id;
			}
			if(pointVSPoints[l->b->id][0]==-1){
				pointVSPoints[l->b->id][0]=l->a->id;
				pointVSEdges[l->b->id][0]=l->id;
			}
			else{
				pointVSPoints[l->b->id][1]=l->a->id;
				pointVSEdges[l->b->id][1]=l->id;
			}
		}

		return o;

}
void VisibilityGraph::print(){
		std::cout<<"Total # of Obstacles :"<<obstacles.size()<<std::endl;
		std::cout<<"Total # of Obstacle Edges :"<<obsSides.size()<<std::endl;
		std::cout<<"Total # of Edges :"<<edges.size()<<std::endl;
		std::cout<<"Total # of Nodes :"<<nodes.size()<<std::endl;
		std::cout<<"Obstacles :"<<std::endl;
		for(int i=0;i<obstacles.size();i++){
			obstacles[i]->print();
		}
		std::cout<<"Visible Edges :"<<std::endl;
		for(int j=0;j<edges.size();j++){
			edges[j]->print();
		}
		std::cout<<"Obstacle Sides :"<<std::endl;
		for(int j=0;j<obsSides.size();j++){
			obsSides[j]->print();
		}

		std::cout<<"Graph Nodes :"<<std::endl;
		for(int j=0;j<nodes.size();j++){
			nodes[j]->print();
		}


}

double VisibilityGraph::value()
{

	return 0.0;
}

int* VisibilityGraph::getOtherEndOfThisPoint(Point* p){
	//std::cout<<"Other Ends of p "<<p->id<<" ->"<<pointVSPoints[p->id][0]<<","<<pointVSPoints[p->id][1]<<std::endl;
	int* otherEnds=pointVSPoints[p->id];
	return otherEnds;

}
int* VisibilityGraph::getEdgesOfThisPoint(Point* p){
	int* e=pointVSEdges[p->id];
	return e;

}


bool lineHasPoint(tPoint p,tLinestring l){
	bool has=false;
	 BOOST_FOREACH(tPoint const& pt, l){
		 if(get<0>(p)==get<0>(pt) && get<1>(p)==get<1>(pt)){
		 		has=true;
		 		continue;
		 }
	 }

	return has;

}

VisibilityGraph::~VisibilityGraph() {
	// TODO Auto-generated destructor stub
}

void VisibilityGraph::removeEdgeFromVisGraph(Line* l){
	std::vector<Line*>::iterator it;
	for (it = edges.begin(); it != edges.end() /* not hoisted */; /* no increment */)
	{
	    if ((*it)->a->id == l->a->id && (*it)->b->id == l->b->id)
	    {
	        it = edges.erase(it);
	    }
	    else
	    {
	        ++it;
	    }
	}
}

void VisibilityGraph::removeObsSideFromVisGraph(Line* l){
	std::vector<Line*>::iterator it;
	for (it = obsSides.begin(); it != obsSides.end() /* not hoisted */; /* no increment */)
	{
	    if ((*it)->a->id == l->a->id && (*it)->b->id == l->b->id)
	    {
	        it = obsSides.erase(it);
	    }
	    else
	    {
	        ++it;
	    }
	}
}

void VisibilityGraph::removeNodeFromVisGraph(Point* a){
	std::vector<Point*>::iterator it;
	for (it = nodes.begin(); it != nodes.end() /* not hoisted */; /* no increment */)
	{
	    if ((*it)->id == a->id)
	    {
	        it = nodes.erase(it);
	    }
	    else
	    {
	        ++it;
	    }
	}
}

void VisibilityGraph::removeObstacleFromVisGraph(Obstacle* obs){
	std::vector<Obstacle*>::iterator it;
	for (it = obstacles.begin(); it != obstacles.end() /* not hoisted */; /* no increment */)
	{
	    if ((*it)->id == obs->id)
	    {
	        it = obstacles.erase(it);
	    }
	    else
	    {
	        ++it;
	    }
	}
}

bool VisibilityGraph::edgeExists(Line* l){
	for(int i=0;i<edges.size();i++){
		if(edges[i]->a->id == l->a->id && edges[i]->b->id == l->b->id){
			return true;
		}
	}
	return false;
}

Line* VisibilityGraph::findEdgeWithPoints(Point* a,Point* b){
	for(int i=0;i<edges.size();i++){
		if(edges[i]->a->id == a->id && edges[i]->b->id == b->id){
			return edges[i];
		}
	}
	return NULL;
}

vector<Line*> VisibilityGraph::findEdgesWithThisPoint(Point* a){
	vector<Line*> lines;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i]->a->id == a->id) {
			lines.push_back(edges[i]);
		} else if (edges[i]->b->id == a->id) {
			lines.push_back(edges[i]);
		}
	}
	return lines;
}
void  VisibilityGraph::insertEdgeInVisGraph(vector<Line*> lines){
	edges.insert(edges.end(),lines.begin(),lines.end());
}

double VisibilityGraph::findShortestPath(double sourceX, double sourceY,
		double destX, double destY,int maxVertexNum,vector<int>& shortestPath ) {
	int numOfPoints = nodes.size();
	int numOfEdges = edges.size();
	int sourcePointId = searchPointByCoord(nodes, sourceX, sourceY)->id;
	int destPointId = searchPointByCoord(nodes, destX, destY)->id;
	/*printf("\nFinding shortest path from %d -> %d\n", sourcePointId,
			destPointId, numOfEdges, numOfPoints);*/
	Point* start;
	Point* goal;
	double shortestDistance=initiateDijkstra(numOfPoints, numOfEdges, false, sourcePointId,
			destPointId,maxVertexNum);
	shortestPath = getShortestPath();
	int i = 0;
	//Print the Shortest Path
	//printf("The Shortest Path is :");
	/*while (i<shortestPath.size()) {
		printf("%d ", shortestPath[i]);
		i++;
	}*/

	return shortestDistance;

	}


