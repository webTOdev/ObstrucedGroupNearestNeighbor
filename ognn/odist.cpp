#include "odist.h"

using namespace std;

char *VISGRAPH_FILE = "Datasets/visibilityGraphPolygons.txt";
VisibilityGraphController* vgController;
vector<string> obstacleList;

//Checks whether the line and obstacle intersect
bool ObstructedDistance::doesLineAndObstcaleIntersects(tLinestring ls,tPolygon p){
	std::vector<turn_info> turns;
	bg::detail::get_turns::no_interrupt_policy policy;
	bg::get_turns<false, false, bg::detail::overlay::assign_null_policy>(ls, p, turns, policy);
	bool intersect=false;
    if(turns.size()>1){
    	intersect=true;
    }
	return intersect;
}

//LineString creation 
tLinestring ObstructedDistance::createLS(double x1,double y1, double x2, double y2){
	 tLinestring edge = boost::assign::tuple_list_of(x1, y1)(x2, y2);
	// std::cout << "hull: " << dsv(edge) << std::endl;
	 return edge;
}
//Polygon creation from string
tPolygon ObstructedDistance::createPolygon(double obstacle[5]){
	char buffer[1024];

	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f,%f %f,%f %f,%f %f%s", "polygon((",
					obstacle[0], obstacle[2], obstacle[1],
					obstacle[2], obstacle[1], obstacle[3],
					obstacle[0], obstacle[3], obstacle[0],
					obstacle[2], "))");
	//printf("%s\n",buffer);			
	tPolygon poly;

	bg::read_wkt(buffer,poly);
	bg::correct(poly);
	// std::cout << "hull: " << dsv(poly) << std::endl;
	return poly;
}

void ObstructedDistance::removeDataPointFromVG(VisibilityGraph* initialVisGraph,float* q){

	char buffer[1024];
	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f%s", "polygon((",
				q[0],q[1], q[0],q[1], "))");
	initialVisGraph = vgController->removeDataPointFromVisGraph(initialVisGraph,initialVisGraph->searchObsWithString(buffer));
}

void ObstructedDistance::removePointAndQueryPointsFromVisGraph(VisibilityGraph* initialVisGraph,
		float* p, Point2D queryPoints[],int numOfQueryPoints){

	removeDataPointFromVG(initialVisGraph,p);
	for(int i=0;i<numOfQueryPoints;i++){
		removeDataPointFromVG(initialVisGraph,queryPoints[i]);
	}


}

void ObstructedDistance::addDataPointInVG(VisibilityGraph* initialVisGraph,float* q){

	char buffer[1024];
	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f%s", "polygon((",
				q[0],q[1], q[0],q[1], "))");
	Obstacle* obs = createObstacle(buffer);
	initialVisGraph = vgController->addNewObstacleForIncrementalVisGraph(initialVisGraph, obs);
	obstacleList.push_back(buffer);
}

bool ObstructedDistance::visGraphContainsPoly(char buffer[1024]) {
	
	printf("%s - ",buffer);
	for(int i=0;i<obstacleList.size();i++){
		printf("%s \n",obstacleList[i].c_str());
		if (obstacleList[i].c_str() == buffer)
		return true;
	}

	return false;
}

double ObstructedDistance::computeAggObstructedDistance(VisibilityGraph* initialVisGraph,
		float* p, Point2D queryPoints[],int numOfQueryPoints, RTree* rt_obstacle,int function){

	double dist_OG=-1;
	vector<float*> l_Q;
	vector<double*> obsInRange;
	float *q;
	std::vector < MyStruct > dist_O_p_qi;
	std::vector < MyShortestPath > shortestPath_p_qi;
	for(int i=0;i<numOfQueryPoints;i++){
		q = new float[2];
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];
		//Compute  the aggregate euclidean distance of p and Q
		double euclideanDist = getDistanceBetweenTwoPoints(p, q);		
		
		printf("\n Euclidean Distance between p %lf,%lf is %lf and q %f,%f\n", p[0],p[1],euclideanDist,q[0],q[1]);
		dist_O_p_qi.push_back(MyStruct(euclideanDist, q));

		//delete q;
	
	}
	std::sort(dist_O_p_qi.begin(), dist_O_p_qi.end(), more_than_key());
	bool firstIteration=true;
	bool pAdded=false;
	if (remove(VISGRAPH_FILE) != 0)
		perror("Error deleting file ");
	//For once
	writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	constructInitialVisGraph(initialVisGraph);
	//Rectangle_BFN_NNQ will be called only once, so handling it from outside
	double obstacle[5];
	rt_obstacle->Rectangle_BFN_NNQ(p, obstacle);
	printf("Nearest Obstacle of (%f,%f), is (%f,%f),(%f,%f) dist %lf\n", p[0], p[1],
			obstacle[0], obstacle[2],obstacle[1], obstacle[3],obstacle[4]);
	double dmax= dist_O_p_qi[0].distance;
	if(obstacle[4]<dmax){
			obsInRange.push_back(obstacle);
	}
	//retrieve_kth_BFN_Rectangle_NNQ retrieves an obstacle whose value> dmax which has to be considered in next round
	double* extraObs = new double[5];
	do{
		dmax= dist_O_p_qi[0].distance;
		printf("dmax %lf\n",dmax);	

		while(1){
			double nObstacle[5];
			if(! firstIteration){
				if(extraObs[4]<dmax){
					obsInRange.push_back(extraObs);
				}
				//If the last obstacle gave a bigger value than current dmax than no need to check for the next
				else
					break;
			}
			rt_obstacle->retrieve_kth_BFN_Rectangle_NNQ(nObstacle,p);
			printf("Next Nearest Obstacle is (%f,%f),(%f,%f) dist %lf\n",
				nObstacle[0], nObstacle[2],nObstacle[1], nObstacle[3],nObstacle[4]);
			if(nObstacle[4]<dmax){				
				obsInRange.push_back(nObstacle);
			}
			else {
				extraObs = new double[5];
				for(int j=0;j<5;j++)
				{
					extraObs[j]=nObstacle[j];
				}
				break;
			}
			
		}


		vector<int> shortestPath;
		if (remove(VISGRAPH_FILE) != 0)
		perror("Error deleting file");
		//Clear the query points list
		l_Q.clear();

		for(int j=0;j<obsInRange.size();j++){
			double* obs=obsInRange[j];
			bool obsAdded=false;
			for(int k=0;k<numOfQueryPoints;k++){
				tPolygon tPoly=createPolygon(obs);
				tLinestring lineS;
				bool intersect=false;
				if(firstIteration){
				    lineS=createLS(p[0],p[1],queryPoints[k][0],queryPoints[k][1]);
					intersect = doesLineAndObstcaleIntersects(lineS,tPoly);
				}
				else{
					intersect = checkIntersectionWithSP(queryPoints[k],shortestPath_p_qi,obs,initialVisGraph);
				}

				//printf("Intersect %d q=%d\n",intersect,k);
				if(intersect){
					l_Q.push_back(queryPoints[k]);
					//p will be added at the first iteration and removed after the algo terminates
					if(!pAdded){
						addDataPointInVG(initialVisGraph,p);
						pAdded=true;
					}
				
					char buffer[1024];

					_snprintf(buffer, sizeof(buffer),
					"%s%f %f,%f %f,%f %f,%f %f,%f %f%s", "polygon((",
						obs[0], obs[2], obs[1],
						obs[2], obs[1], obs[3],
						obs[0], obs[3], obs[0],
						obs[2], "))");

						if(! obsAdded){
							if (!visGraphContainsPoly(buffer)) {
								Obstacle* newObs = createObstacle(buffer);
								initialVisGraph = vgController->addNewObstacleForIncrementalVisGraph(
									initialVisGraph, newObs);
								obstacleList.push_back(buffer);
							}
							obsAdded=true;

						}
						
				}
				
			}
		}

		//initialVisGraph->print();

		for(int j=0;j<l_Q.size();j++){

			vector<int> shortestPath;
			double obsDist = computeObstructedDistance(initialVisGraph,p,l_Q[j],shortestPath);
			//Update obstructed dist_O(p,q_i)
			replaceObsDist(dist_O_p_qi,l_Q[j],obsDist);
			//shortestPath_p_qi.push_back(MyShortestPath(shortestPath, l_Q[j]));
			addOrReplaceSP(l_Q[j],shortestPath_p_qi,shortestPath);
		}

		for(int j=0;j<dist_O_p_qi.size();j++){
			printf("\n(%f,%f) has distance %f\n",dist_O_p_qi[j].queryPoints[0],dist_O_p_qi[j].queryPoints[1],dist_O_p_qi[j].distance);
		}
		firstIteration=false;
		//Sort the obstructed distance as some distance may be updated by now
		std::sort(dist_O_p_qi.begin(), dist_O_p_qi.end(), more_than_key());
		//Only new obstacle retrived will be checked in next iteration
		//Other obstacles are already considered and added in the vis graph
		obsInRange.clear();

	}


	while(! l_Q.empty());
	dist_OG=0.0;
	if(function==0)
		for(int i=0;i<numOfQueryPoints;i++){
			dist_OG+=dist_O_p_qi[i].distance;
		}

	delete q;
	delete extraObs;
	//Set the VisGraph in the initial state , next time new p  will be added again
	removeDataPointFromVG(initialVisGraph,p);
	return dist_OG;

}

bool ObstructedDistance::checkIntersectionWithSP(float* q,vector< MyShortestPath>& shortestPath_p_qi,double obs[5],
												 VisibilityGraph* initialVisGraph){

	bool intersect=false;
	vector<int> sp;
	//Find shortest path of q
	for(int j=0;j<shortestPath_p_qi.size();j++){
			if(shortestPath_p_qi[j].queryPoints[0]==q[0] && shortestPath_p_qi[j].queryPoints[1]==q[1]){
				sp=shortestPath_p_qi[j].shortestPath;
				break;
			}
	}
	//Read in pair 0 3 1 -> 0 3, 3 1
	for(int i=0;i<sp.size()-1;i++){
		int id=sp[i];
		Point* a = getPointById(initialVisGraph->nodes,id);
		//printf("Point id %d =%d : (%lf,%lf) \n",id,a->id,a->x,a->y);
		id=sp[i+1];
		Point* b = getPointById(initialVisGraph->nodes,id);
		//printf("Point id %d =%d : (%lf,%lf) \n",id,a->id,a->x,a->y);

		tPolygon tPoly=createPolygon(obs);
		tLinestring lineS;
		lineS=createLS(a->x,a->y,b->x,b->y);
		intersect = doesLineAndObstcaleIntersects(lineS,tPoly);
		if(intersect)
			break;
	}
	return intersect;
}

void ObstructedDistance::addOrReplaceSP(float* q,vector< MyShortestPath>& shortestPath_p_qi,vector<int>& shortestPath){

	bool found=false;
	//Replace
	for(int j=0;j<shortestPath_p_qi.size();j++){
			if(shortestPath_p_qi[j].queryPoints[0]==q[0] && shortestPath_p_qi[j].queryPoints[1]==q[1]){
				shortestPath_p_qi[j].shortestPath=shortestPath;
				found=true;
				break;
			}
	}
	if(! found){
		shortestPath_p_qi.push_back(MyShortestPath(shortestPath, q));
	}
}

void ObstructedDistance::replaceObsDist(std::vector < MyStruct >& dist_O_p_qi,float* q,double obsDist){
	for(int j=0;j<dist_O_p_qi.size();j++){
			if(dist_O_p_qi[j].queryPoints[0]==q[0] && dist_O_p_qi[j].queryPoints[1]==q[1]){
				dist_O_p_qi[j].distance=obsDist;
			}
	}
}


int ObstructedDistance::drawAndWriteFileVisEdges(vector<Line*> visEdges) {
	//Remove existing test.txt file
	if (remove("test.txt") != 0)
		perror("Error deleting file");
	else
		puts("File successfully deleted");
	/*FILE *fp;
	fp=fopen("test.txt", "w");
	fclose(fp);
*/
	//Dijkstra algorithm needs to create a vector of size n , where n is the id of the vertex
	int maxNumberOfVertex=-1;
	for (int i = 0; i < visEdges.size(); i++) {
		fileWrite(visEdges[i]->a, visEdges[i]->b);
		if(visEdges[i]->a->id > maxNumberOfVertex){
			maxNumberOfVertex=visEdges[i]->a->id;
		}else
			if(visEdges[i]->b->id > maxNumberOfVertex){
			maxNumberOfVertex=visEdges[i]->b->id;
		}
	}
	return maxNumberOfVertex+1;
}

double ObstructedDistance::computeObstructedDistance(VisibilityGraph* initialVisGraph,float* p, float* q,vector<int>& shortestPath) {
	int maxVertexNum = drawAndWriteFileVisEdges(initialVisGraph->edges);
	double shortestPathDistance = initialVisGraph->findShortestPath(p[0], p[1],
			q[0], q[1],maxVertexNum,shortestPath);
	int i = 0;
	//Print the Shortest Path
	printf("--------The Shortest Path is :");
	while (i<shortestPath.size()) {
		printf("%d ", shortestPath[i]);
		i++;
	}
	//printf("\nObstructed Distance found is %lf\n", shortestPathDistance);
	if(shortestPathDistance==infinity){
		return infinity;
	}

	return shortestPathDistance;
}

void ObstructedDistance::constructInitialVisGraph(VisibilityGraph* initialVisGraph) {

	vector<Obstacle*> obsList;
	ifstream iFile(VISGRAPH_FILE);
	string line;
	Obstacle* obs;
	//ObstacleController* obsController= new ObstacleController();

	/* While there is still a line. */
	while (getline(iFile, line)) {
		/* Printing goes here. */
		//cout << line << endl;
		//Keeping all the polygon string in a vector so that next time we get which obstacles are already added in VisGraph
		obstacleList.push_back(line);
		obs = createObstacle(line);
		obsList.push_back(obs);
	}

	initialVisGraph->setObstacle(obsList);
	vgController = new VisibilityGraphController(initialVisGraph);
	//Construct Vis Graph
	vector<Line*> visEdges = vgController->constructVisGraph();


	iFile.close();

}
void ObstructedDistance::writePointAndQueryPointsInFile(float* p,Point2D queryPoints[],int numOfQueryPoints){
	writePointInFile(p);
	writeQueryPointsInFile(queryPoints,numOfQueryPoints);
	

}
void ObstructedDistance::writeQueryPointsInFile(Point2D queryPoints[],int numOfQueryPoints){

	float *q;
	for(int i=0;i<numOfQueryPoints;i++){
		q = new float[2];
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];
		writePointInFile(q);
		delete q;

	}
}

void ObstructedDistance::writePointInFile(float* p){
	FILE * input = fopen(VISGRAPH_FILE, "a");
	if (input == NULL) {
		printf("Cannot open file %s", VISGRAPH_FILE);
	}
	fprintf(input, "polygon((");
	fprintf(input, "%f %f,%f %f", p[0], p[1], p[0], p[1]);
	fprintf(input, "))\n");

	fclose(input);
}
void ObstructedDistance::writePolygonInFile(double obstacle[5]) {
	
	FILE * input = fopen(VISGRAPH_FILE, "a");
	if (input == NULL) {
		printf("Cannot open file %s", VISGRAPH_FILE);
	}

	char buffer[1024];

	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f,%f %f,%f %f,%f %f%s", "polygon((",
					obstacle[0], obstacle[2], obstacle[1],
					obstacle[2], obstacle[1], obstacle[3],
					obstacle[0], obstacle[3], obstacle[0],
					obstacle[2], "))");

	fprintf(input, "%s\n", buffer);
		
	fclose(input);
}


void shortestPathMatrix(double q_index,vector<int> points,vector<vector<pair<double,double>>>& sp){
	
}
void shortestPathMatrixEuclidean(int q_index,double p_x1,double p_y1,vector<vector<pair<double,double>>>& sp){
	sp[q_index].push_back(make_pair(p_x1,p_y1));
}