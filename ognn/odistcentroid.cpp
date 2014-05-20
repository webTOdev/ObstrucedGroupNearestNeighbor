#include "odistcentroid.h"

using namespace std;

vector<string> obsInRange;
double calculateInitializeThreshold(Point2D centroid[],float* p, Point2D queryPoints[],int numOfQueryPoints);

double calculateInitializeThreshold(Point2D centroid[],float* p, Point2D queryPoints[],int numOfQueryPoints){
	double max=0.0;
	float q[2],c[2];
	c[0]=centroid[0][0];
	c[1]=centroid[0][1];
	for(int i=0;i<numOfQueryPoints;i++){
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];	
		double euclideanDist_c_qi = getDistanceBetweenTwoPoints(c, q);
		double euclideanDist_p_qi = getDistanceBetweenTwoPoints(p, q);	
		if(max<(euclideanDist_c_qi+euclideanDist_p_qi)){
			max=euclideanDist_c_qi+euclideanDist_p_qi;
		}
	}
	return max;
}

double calculateThreshold(Point2D centroid[],float* p, Point2D queryPoints[],int numOfQueryPoints,std::vector<MyStruct>& dist_O_p_qi){
	double max=0.0;
	float q[2],c[2];
	c[0]=centroid[0][0];
	c[1]=centroid[0][1];
	for(int i=0;i<numOfQueryPoints;i++){
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];	
		double euclideanDist_c_qi = getDistanceBetweenTwoPoints(c, q);
		double obsDist_p_qi = dist_O_p_qi[i].distance;	
		if(max<(euclideanDist_c_qi+obsDist_p_qi)){
			max=euclideanDist_c_qi+obsDist_p_qi;
		}
	}
	return max;
}
void ObstructedDistanceCentroid::initialize(Point2D queryPoints[],int numOfQueryPoints){

	L_R.resize(numOfQueryPoints);
	L_C.resize(numOfQueryPoints);
	L_N.resize(numOfQueryPoints);
	float *q;
	//Initialize p as inifinity distance from q_i and L_R to contain q_i
	for(int i=0;i<numOfQueryPoints;i++){
		q = new float[2];
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];	

		L_R[i].push_back(MyStruct(0, q));

	}
}

void initDistanceVector(std::vector<MyStruct>& dist_O_p_qi,Point2D queryPoints[],int numOfQueryPoints){
	float *q;
	//Initialize p as inifinity distance from q_i and L_R to contain q_i
	for(int i=0;i<numOfQueryPoints;i++){
		q = new float[2];
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];	
		
		//printf("\nInitialize Infinity Distance between p %lf,%lf is %lf and q %f,%f\n", p[0],p[1],euclideanDist,q[0],q[1]);
		dist_O_p_qi.push_back(MyStruct(infinty, q));
	}
}

bool pointInListAlready( vector < MyStruct >& L_R_qi,float *p){
	float *q;
	for(int i=0;i<L_R_qi.size();i++){
		q = new float[2];
		
		q[0]=L_R_qi[i].queryPoints[0];
		q[1]=L_R_qi[i].queryPoints[1];
		//printf("p(%lf,%lf) and in list(%lf,%lf)\n",p[0],p[1],q[0],q[1]);
		if(p[0]==q[0] && p[1]==q[1])
			return true;
	}
	return false;

}
bool ObstructedDistanceCentroid::visGraphContainsPoly(char buffer[1024]) {
	
	for(int i=0;i<obstacleList.size();i++){
		if (obstacleList[i].c_str() == buffer)
		return true;
	}

	return false;
}
tPolygon createPolygon(string buffer){
	tPolygon poly;

	bg::read_wkt(buffer.c_str(),poly);
	bg::correct(poly);
	return poly;
}

void ObstructedDistanceCentroid::addNewObstacleInVisGraph(double* obs,VisibilityGraph* initialVisGraph){
	char buffer[1024];
	Clock sw1;
	_snprintf(buffer, sizeof(buffer),
	"%s%f %f,%f %f,%f %f,%f %f,%f %f%s", "polygon((",
	obs[0], obs[2], obs[1],
	obs[2], obs[1], obs[3],
	obs[0], obs[3], obs[0],
	obs[2], "))");
	if (!visGraphContainsPoly(buffer)) {
			Obstacle* newObs = createObstacle(buffer);
			sw1.start();
			initialVisGraph = vgController->addNewObstacleForIncrementalVisGraph(
									initialVisGraph, newObs);
			sw1.stop();
			visGraphConsTime+=sw1.getDiff();
			obstacleList.push_back(buffer);
			obsInRange.push_back(buffer);
	}
}
//Check if any dist_O(qi) has distance greater than threshold
bool findTerminationCondition(std::vector<MyStruct>& dist_O_p_qi,Point2D queryPoints[],int numOfQueryPoints,double t){
	bool condition=false;
	for(int index=0;index<numOfQueryPoints;index++){		
			if(dist_O_p_qi[index].distance>t)
			{
				condition=true;
				return condition;
			}			
	}
	return condition;
}
double ObstructedDistanceCentroid::computeAggObstructedDistance(VisibilityGraph* initialVisGraph,
																float* p, Point2D queryPoints[],int numOfQueryPoints, RTree* rt_obstacle,int function){										
	Clock sw1;
	double dist_OG=-1;
	std::vector < MyStruct > dist_O_p_qi;
	if(globalThreshold==0.0){
		initialize(queryPoints,numOfQueryPoints);
	}
	initDistanceVector(dist_O_p_qi,queryPoints,numOfQueryPoints);
	float centroid[1][2];
	centroidOfQ(queryPoints,numOfQueryPoints,centroid);
	double threshold=calculateInitializeThreshold(centroid,p,queryPoints,numOfQueryPoints);
	//printf("threshold %lf\n",threshold);
	sw1.start();
	addDataPointInVG(initialVisGraph,p);
	sw1.stop();
	visGraphConsTime+=sw1.getDiff();

	//Rectangle_BFN_NNQ will be called only once, so handling it from outside
	double obstacle[5];
	//This call is for the first p
	if(globalThreshold==0.0){
		rt_obstacle->Rectangle_BFN_NNQ(centroid[0], obstacle);
		/*printf("Nearest Obstacle of (%f,%f), is (%f,%f),(%f,%f) dist %lf\n", p[0], p[1],
				obstacle[0], obstacle[2],obstacle[1], obstacle[3],obstacle[4]);*/
		//obsInRange.push_back(obstacle);
		if(obstacle[4]<threshold){		
			addNewObstacleInVisGraph(obstacle,initialVisGraph);
		}
	}
	double* extraObs = new double[5];
	bool firstIteration=true;
	bool condition=false;

	vector<int> shortestPath;
	do{
		//Check whether the obstacle within the new threshold has already been retrieved
		if(threshold>globalThreshold){
			while(1){
				double *nObstacle=new double[5];
				if(! firstIteration){
					if(extraObs[4]<threshold){
						//obsInRange.push_back(extraObs);
						addNewObstacleInVisGraph(extraObs,initialVisGraph);
					}
					//If the last obstacle gave a bigger value than current dmax than no need to check for the next
					else
						break;
				}
				bool treeEmpty = rt_obstacle->retrieve_kth_BFN_Rectangle_NNQ(nObstacle,centroid[0]);
				if(!treeEmpty){
					/*printf("Next Nearest Obstacle is (%f,%f),(%f,%f) dist %lf and threshold %lf\n",
						nObstacle[0], nObstacle[2],nObstacle[1], nObstacle[3],nObstacle[4],threshold);*/
					if(nObstacle[4]<threshold){				
						//obsInRange.push_back(nObstacle);
						addNewObstacleInVisGraph(nObstacle,initialVisGraph);
					}
					else {
						extraObs = new double[5];
						firstIteration=false;
						for(int j=0;j<5;j++)
						{
							extraObs[j]=nObstacle[j];
						}
						break;
					}
				}
				else
					break;
				
			}
		}

		float *q,*v;
		for(int i=0;i<numOfQueryPoints;i++){
		q = new float[2];
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];	
		//Have we already found the real distance between dist_O(q_i,p)?
		if(! pointInListAlready(L_R[i],p)){
			//Take each vertex one by one in L_R[i]
			for(int j=0;j<L_R[i].size();j++){
				v = new float[2];
				v[0]=L_R[i][j].queryPoints[0];
				v[1]=L_R[i][j].queryPoints[1];
				
				if(isVisible(v,p,initialVisGraph)){
					double dist_o_v_q=L_R[i][j].distance;
					relax(v,q,p,initialVisGraph,dist_O_p_qi,dist_o_v_q);
				}
			}
			double dist_obs_p_q=computeObstructedDistance(initialVisGraph,p,q,shortestPath);
			replaceObsDist(dist_O_p_qi,q,dist_obs_p_q);
			//Check we have found already found the real distance between q_i and p
			bool realDistanceFound = isRealDistanceFor_qFound(q,p,dist_O_p_qi,threshold,i);
			if(realDistanceFound)
				continue;
			else{				
				addVerticesOfObsInRangeInLc(obsInRange,initialVisGraph,i);
			}
			Point* vertex;
			double dist_O_v_q=infinty;
			vector < Point* >::iterator it = L_N[i].begin();
			while(it != L_N[i].end()) {
				vertex=(Point*)(*it);
				bool intersect=false;
				vector<Line*> adjacentEdges=initialVisGraph->findEdgesWithThisPoint(vertex);
				for(int k=0;k<obsInRange.size();k++){
					tPolygon poly=createPolygon(obsInRange[k]);
					for(int l=0;l<adjacentEdges.size();l++){
						Line* line=adjacentEdges[l];
						tLinestring lineS=createLS(line->a->x,line->a->y,line->b->x,line->b->y);
						intersect = doesLineAndObstcaleIntersects(lineS,poly);
						if(intersect) break;
					}
					if(intersect) break;

				}
				if(intersect) {
					L_N[i].erase(it);
					L_C[i].push_back(vertex);
				}
				else{
					vector<int> shortestPath;
					float *nVertex=new float[2];
					nVertex[0]=vertex->x;
					nVertex[1]=vertex->y;
					dist_O_v_q=computeObstructedDistance(initialVisGraph,nVertex,q,shortestPath);
					if(dist_O_v_q <= threshold){
						L_N[i].erase(it);
						
						L_R[i].push_back(MyStruct(dist_O_v_q,nVertex));
					}else
					{
						L_N[i].erase(it);
						L_C[i].push_back(vertex);
					}
				}
				++it; // ++i is usually faster than i++. It's a good habit to use it.

			}

			for(int j=0;j<L_C[i].size();j++){
				Point* mPoint=L_C[i][j];	
				float *mVertex=new float[2];
				mVertex[0]=mPoint->x;
				mVertex[1]=mPoint->y;
				dist_O_v_q=computeObstructedDistance(initialVisGraph,mVertex,q,shortestPath);
				if(dist_O_v_q <= threshold){						
					L_R[i].push_back(MyStruct(dist_O_v_q,mVertex));
				}else
				{
					L_N[i].push_back(mPoint);
				}
			}
			L_C[i].clear();
		}
	}

	threshold=calculateThreshold(centroid,p,queryPoints,numOfQueryPoints,dist_O_p_qi);
	condition=findTerminationCondition(dist_O_p_qi,queryPoints,numOfQueryPoints,threshold);
	obsInRange.clear();
	}while(condition);
	//Keeping track of how much radius obstacle has been retrived centering c_Q
	globalThreshold=threshold;
	dist_OG=0.0;
	if(function==0){
		for(int i=0;i<numOfQueryPoints;i++){
			dist_OG+=dist_O_p_qi[i].distance;
		}
	}
	if(function==1){
		std::sort(dist_O_p_qi.begin(), dist_O_p_qi.end(), more_than_key());
		dist_OG=dist_O_p_qi[0].distance;
	}

	//delete q;
	//delete v;
	delete extraObs;
	/*//Set the VisGraph in the initial state remove p if added in graph, next time new p  will be added again
	if(pAdded){
		sw1.start();
		removeDataPointFromVG(initialVisGraph,p);
		sw1.stop();
		visGraphConsTime+=sw1.getDiff();
	}*/
	//The rectangle heap from obs Rtree as it will be recreated for the next p
	//delete rt_obstacle->rectangleNNHeap;
	
	/*float *q;
	for(int i=0;i<L_R.size();i++){
		for(int j=0;j<L_R[i].size();j++){
			q = new float[2];
			q[0]=L_R[i][j].queryPoints[0];
			q[1]=L_R[i][j].queryPoints[1];
			printf("list[%d] (%lf,%lf)\n",i,q[0],q[1]);
		}
	}*/
	return dist_OG;
}

double ObstructedDistanceCentroid::computeObstructedDistance(VisibilityGraph* initialVisGraph,float* p, float* q,vector<int>& shortestPath) {
	int maxVertexNum = drawAndWriteFileVisEdges(initialVisGraph->edges);
	Clock sw1;
	sw1.start();
	double shortestPathDistance = initialVisGraph->findShortestPath(p[0], p[1],
			q[0], q[1],maxVertexNum,shortestPath);
	sw1.stop();
	shortestPathCalcTime+=sw1.getDiff();
	/*int i = 0;
	//Print the Shortest Path
	printf("--------The Shortest Path is :");
	while (i<shortestPath.size()) {
		printf("%d ", shortestPath[i]);
		i++;
	}
	//printf("\nObstructed Distance found is %lf\n", shortestPathDistance);
	*/
	if(shortestPathDistance==infinty){
		return infinty;
	}

	return shortestPathDistance;
}

void ObstructedDistanceCentroid::replaceObsDist(std::vector < MyStruct >& dist_O_p_qi,float* q,double obsDist){
	for(int j=0;j<dist_O_p_qi.size();j++){
			if(dist_O_p_qi[j].queryPoints[0]==q[0] && dist_O_p_qi[j].queryPoints[1]==q[1]){
				dist_O_p_qi[j].distance=obsDist;
			}
	}
}

void ObstructedDistanceCentroid::removeDataPointFromVG(VisibilityGraph* initialVisGraph,float* q){

	char buffer[1024];
	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f%s", "polygon((",
				q[0],q[1], q[0],q[1], "))");
	initialVisGraph = vgController->removeDataPointFromVisGraph(initialVisGraph,initialVisGraph->searchObsWithString(buffer));
}
void ObstructedDistanceCentroid::addDataPointInVG(VisibilityGraph* initialVisGraph,float* q){

	char buffer[1024];
	_snprintf(buffer, sizeof(buffer),
				"%s%f %f,%f %f%s", "polygon((",
				q[0],q[1], q[0],q[1], "))");
	Obstacle* obs = createObstacle(buffer);
	initialVisGraph = vgController->addNewObstacleForIncrementalVisGraph(initialVisGraph, obs);
	obstacleList.push_back(buffer);
}
int ObstructedDistanceCentroid::drawAndWriteFileVisEdges(vector<Line*> visEdges) {
	//Remove existing test.txt file
	if (remove("test.txt") != 0)
		perror("Error deleting file");
	//else
		//puts("File successfully deleted");
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

void ObstructedDistanceCentroid::addVerticesOfObsInRangeInLc(vector<string>& obsInRange,VisibilityGraph* initialVisGraph,int q_index){
	for(int i=0;i<obsInRange.size();i++){
		Obstacle* o = initialVisGraph->searchObsWithString(obsInRange[i]);
		vector<Point*> vertices=o->vertices;
		L_C[q_index].insert(L_C[q_index].end(), vertices.begin(), vertices.end());


	}
}

bool ObstructedDistanceCentroid::isRealDistanceFor_qFound(float* q,float* p,std::vector<MyStruct>& dist_O_p_qi,double threshold,int q_index){
	bool realFound=false;
	double dist_O_p_q;
	int k;
	//Find the current distance between p and q_i 
	for(k=0;k<dist_O_p_qi.size();k++){
		if(dist_O_p_qi[k].queryPoints[0]==q[0] && dist_O_p_qi[k].queryPoints[1]==q[1]){
			dist_O_p_q=dist_O_p_qi[k].distance;
			break;
		}
	}
	//The new distance is less than current distance so update
	if(dist_O_p_q < threshold){
		realFound=true;
		L_R[q_index].push_back(MyStruct(dist_O_p_q, p));
	}

	return realFound;
}

void ObstructedDistanceCentroid::relax(float *v,float* q,float *p,VisibilityGraph* initialVisGraph,std::vector<MyStruct>& dist_O_p_qi,double dist_o_v_q){
	double dist_E_v_p=getDistanceBetweenTwoPoints(v, p);
	double dist_O_p_q;
	int k;
	//Find the current distance between p and q_i 
	for(k=0;k<dist_O_p_qi.size();k++){
		if(dist_O_p_qi[k].queryPoints[0]==q[0] && dist_O_p_qi[k].queryPoints[1]==q[1]){
			dist_O_p_q=dist_O_p_qi[k].distance;
			break;
		}
	}
	//The new distance is less than current distance so update
	if(dist_o_v_q+dist_E_v_p < dist_O_p_q){
		dist_O_p_qi[k].distance=dist_o_v_q+dist_E_v_p;
	}
}
//Check if with any obstacle v-p line intersect
bool ObstructedDistanceCentroid::isVisible(float *v,float* p,VisibilityGraph* initialVisGraph){
	tLinestring lineS;
	bool visible=true;
	lineS=createLS(v[0],v[1],p[0],p[1]);
	vector<Obstacle*> obstacles=initialVisGraph->obstacles;
	for(int i=0;i<obstacles.size();i++){
		bool intersect = doesLineAndObstcaleIntersects(lineS,obstacles[i]->poly);
		if(intersect)
			return false;
	}
	return visible;
}

//Checks whether the line and obstacle intersect
bool ObstructedDistanceCentroid::doesLineAndObstcaleIntersects(tLinestring ls,tPolygon p){
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
tLinestring ObstructedDistanceCentroid::createLS(double x1,double y1, double x2, double y2){
	 tLinestring edge = boost::assign::tuple_list_of(x1, y1)(x2, y2);
	// std::cout << "hull: " << dsv(edge) << std::endl;
	 return edge;
}

void ObstructedDistanceCentroid::writeQueryPointsInFile(Point2D queryPoints[],int numOfQueryPoints){

	if (remove(VISGRAPH_FILE) != 0)
		perror("Error deleting file ");
	float *q;
	for(int i=0;i<numOfQueryPoints;i++){
		q = new float[2];
		q[0] = queryPoints[i][0];
		q[1] = queryPoints[i][1];
		writePointInFile(q);
		delete q;

	}
}

void ObstructedDistanceCentroid::writePointInFile(float* p){
	FILE * input = fopen(VISGRAPH_FILE, "a");
	if (input == NULL) {
		printf("Cannot open file %s", VISGRAPH_FILE);
	}
	fprintf(input, "polygon((");
	fprintf(input, "%f %f,%f %f", p[0], p[1], p[0], p[1]);
	fprintf(input, "))\n");

	fclose(input);
}
void ObstructedDistanceCentroid::constructInitialVisGraph(VisibilityGraph* initialVisGraph) {

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