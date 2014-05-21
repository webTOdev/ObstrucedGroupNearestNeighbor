//***************************************************
//This is the implementation of R*-tree v0.1

//Last revised July 4.
//***************************************************

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>

#include "./rtree/rtree.h"
#include "./rtree/rtnode.h"
#include "./rtree/entry.h"
#include "./blockfile/blk_file.h"
#include "./blockfile/cache.h"
#include "./linlist/linlist.h"main

#include "./rtree/rtree_cmd.h"
#include "./rtree/distance.h"
#include "./his/histogram.h"
#include "ognn/ognn.h"
#include "ognn/odist.h"
#include "ognn/ognngnn.h"
#include "ognn/odistcentroid.h"

//Added by Tanzima
//#include "global.h"
//...............
using namespace std;

//----------------------------------------------------------------------------------------------------------
// SN 21/10/2007 001 <Start>

#define	DIMENSION 2
#define DEFAULT_C 0
const double FLOAT_INFTY = numeric_limits<double>::infinity();
//const double MAXDOUBLE = numeric_limits<double>::infinity();

const double MINX = 0;
const double MINY = 0;
const double MAXX = 10000;
const double MAXY = 10000;
const double THRESHOLD = 0.000000001;

//Experiments Parameters

const int SAMPLE = 1;
//DEFAULT
const int DEFAULT_GRPSIZE = 64;
const int DEFAULT_K = 4;
const double DEFAULT_M_AREAPART = 0.04;
const double DEFAULT_M_RATIO = 1;
const double DEFAULT_R_AREAPART = 0.00005;
const double DEFAULT_R_RATIO = 1;
//MIN
const int MIN_GRPSIZE = 4;
const int MIN_K = 2;
const double MIN_M_AREAPART = 0.02;
const double MIN_M_RATIO = 1;
const double MIN_R_AREAPART = 0.00001;
const double MIN_R_RATIO = 1;
//MAX
const int MAX_GRPSIZE = 256;
const int MAX_K = 32;
const double MAX_M_AREAPART = 0.16;
const double MAX_M_RATIO = 16;
const double MAX_R_AREAPART = 0.0001;
const double MAX_R_RATIO = 16;
//INTERVAL
const int INTERVAL_GRPSIZE = 4;
const int INTERVAL_K = 2;
const double INTERVAL_M_AREAPART = 2;
const double INTERVAL_M_RATIO = 2;
const double INTERVAL_R_AREAPART = 0.00001;
const double INTERVAL_R_RATIO = 2;

//Tanzima
char *DATAFILE = "Datasets/sample.txt";
char *TREEFILE = "Datasets/sample.tree";
char *DATAFILE_MBR = "Datasets/sample_mbr.txt";
char *TREEFILE_MBR = "Datasets/sample_mbr.tree";
char *STREEFILE = "H:/Thesis code/GRPLBSQ-Nusrat/Datasets/cas.tree";
char *DTREEFILE = "H:/Thesis code/GRPLBSQ-Nusrat/Datasets/cad.tree";
char *STREEFILE1 = "H:/Thesis code/GRPLBSQ-Nusrat/Datasets/ca1.tree";
char *STREEFILE2 = "H:/Thesis code/GRPLBSQ-Nusrat/Datasets/ca2.tree";
char *STREEFILE3 = "H:/Thesis code/GRPLBSQ-Nusrat/Datasets/ca3.tree";
char *RESFILE = "H:/Thesis code/GRPLBSQ-Nusrat/Results/result.txt";
char *RECTFILE = "H:/Thesis code/GRPLBSQ-Nusrat/Results/rect.txt";
char *RECTRATIOFILE = "H:/Thesis code/GRPLBSQ-Nusrat/Results/rectratio.txt";
int io_access;
int disktime;
int updatetime;
int counttime;
int kmintime;

const int MAXPOINTS = 10000;
const int MAXTRAJECTORYPOINTS = 5000;
const int MAXTRAJECTORY = 5;
double PI = (double) 3.1415926535;
//typedef double Point2D[DIMENSION];

class Stopwatch {

public:
	void start() {
		c1 = clock();
		t1 = time(0);
	}
	;
	void stop() {
		c2 = clock();
		t2 = time(0);
	}
	;
	int getDiff() {
		return (c2 - c1);
	}
	;

private:
	time_t t1, t2;
	clock_t c1, c2;
};

class Exp_sum_stat {
public:
	//Input parameter
	int k;
	int grpsize;
	double m_area;
	double r_area;

	//Output parameter
	long double stime_sec;
	long double page_faults;
	long double snum_retrievals;

	long double cmintime_sec;
	long double cminmaxtime_sec;
	long double cnum_retrievals[MAX_GRPSIZE];

	Exp_sum_stat() {
		k = DEFAULT_K;
		grpsize = DEFAULT_GRPSIZE;
		m_area = DEFAULT_M_AREAPART;
		r_area = DEFAULT_R_AREAPART;

		stime_sec = 0.0;
		page_faults = 0.0;
		snum_retrievals = 0.0;
		cmintime_sec = 0.0;
		cminmaxtime_sec = 0.0;

		for (int j = 0; j < MAX_GRPSIZE; j++)
			cnum_retrievals[j] = 0.0;
	}
	;
	~Exp_sum_stat() {
	}
	;
};

//Added By Nusrat
class Exp_stat {
public:
	//Input parameter
	int k;
	int grpsize;


	//Output parameter
	long double stime_sec;
	long double io_access;
	long double io_access_obs;
	long double page_faults;
	long float qArea;
	int totalNumberOfPRetrieved;
	long double visGraphConstime_sec;
	long double shortestPathCalctime_sec;



	Exp_stat() {
		k = DEFAULT_K;
		grpsize = DEFAULT_GRPSIZE;
	
		stime_sec = 0.0;
		io_access=0.0;
		io_access_obs=0.0;
		page_faults = 0.0;
		qArea=0.0;
		totalNumberOfPRetrieved=0.0;
		
	}
	;
	~Exp_stat() {
	}
	;
};

// Process
void rect_kGNN_query_max(int n_sample, int g_size, int k, char *s1, char *s2,
		Exp_sum_stat *max_e) {

	Rectangle1 r[MAX_GRPSIZE], m;
	Point2D p[MAX_GRPSIZE];

	//Open input1 File
	FILE *input1, *input2, *dFile;
	char temp[200];

	dFile = fopen("H:/Thesis code/GRPLBSQ-Nusrat/Results/InputFile/debugfilter",
			"a+");
	if (dFile == NULL) {
		printf("Error reading rectdata\n");
	}

	input1 = fopen(s1, "r");
	if (input1 == NULL) {
		printf("Error reading rectdata\n");
	}

	fgets(temp, 200, input1);
	//puts(temp);
	fgets(temp, 200, input1);
	//puts(temp);	

	input2 = fopen(s2, "r");
	if (input2 == NULL) {
		printf("Error reading point location\n");
	}

	fgets(temp, 200, input2);
	//puts(temp);
	fgets(temp, 200, input2);
	//puts(temp);	
	//

	//vector<Pointlocation> _rslt;
	Pointlocation rslt[MAXDATALIMIT];
	int num_of_data = 0;
	int blocksize = 1024;	//4096;

	// R-tree
	Cache *cache = new Cache(DEFAULT_C, blocksize);
	RTree *rt = new RTree(STREEFILE, cache);

	for (int i = 0; i < n_sample; i++) {
		//Experiment
		Stopwatch sw1, sw2, sw3;
		int last_pf = cache->page_faults;
		//..........

		fscanf(input1, "%f%f%f%f", &m.x1, &m.x2, &m.y1, &m.y2);

		for (int j = 0; j < g_size; j++) {
			fscanf(input1, "%f%f%f%f", &r[j].x1, &r[j].x2, &r[j].y1, &r[j].y2);
			fscanf(input2, "%f%f", &p[j][0], &p[j][1]);
		}

		num_of_data = 0;
		sw1.start();
		rt->private_kGNN_max(r, g_size, k, rslt, &num_of_data);
		sw1.stop();
		max_e->stime_sec += sw1.getDiff();
		max_e->snum_retrievals += num_of_data;
		max_e->page_faults += cache->page_faults - last_pf;

		if (num_of_data >= MAXDATALIMIT)
			printf("\nFirstcheck:LIMIT EXCEEDED");

		//For debugging
		FILE * result;
		result = fopen(RESFILE, "w");
		if (result == NULL) {
			printf("Error writing rectfile\n");
		} else {
			// Process & close file
			//fprintf(rFile,"%d\n",k);
			fprintf(result, "%d\n", num_of_data);
			for (int l = 0; l < num_of_data; l++)
				//fprintf(result,"%.5f\t%.5f\t%.5f\n",_rslt[i].x, _rslt[i].y,Dist(_rslt[i],o));
				fprintf(result, "%.5f\t%.5f\t%.5f\t%.5f\n", rslt[l].x,
						rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			fclose(result);
		}
		//.....................................

		//preprocessing to match the data structure
		double tdist;
		Pointlocation trslt[MAXDATALIMIT];
		for (int l = 0; l < num_of_data; l++) {
			trslt[l].x = rslt[l].x;
			trslt[l].y = rslt[l].y;
			trslt[l].dmin = rslt[l].dmin;
			trslt[l].dmax = rslt[l].dmax;
		}

		//min filter

		//For debugging		
		result = fopen(
				"H:/Thesis code/GRPLBSQ-Nusrat/Results/result_minfilter.txt",
				"w");
		if (result == NULL) {
			printf("Error writing rectfile\n");
		} else {
			fprintf(result, "%d\n", num_of_data);
		}
		//.....................................

		double maxdist[MAX_K];

		sw2.start();
		for (int j = 0; j < g_size; j++) {
			//update 
			for (int l = 0; l < num_of_data; l++) {
				tdist = Dist1(trslt[l], p[j]);
				/*
				 if(tdist>trslt[l].dmax)
				 {
				 printf("\nHow?\n");
				 if(MAXRECTDIST1(rslt[l],r[j])<tdist)
				 {
				 printf("\n%f\t%f\t%f\t%f\t%f\t%f\n",r[j].x1,r[j].x2,r[j].y1,r[j].y2,p[j][0],p[j][1]);
				 }
				 }
				 */
				if (tdist > trslt[l].dmin)
					trslt[l].dmin = tdist;
			}
		}

		//sort	
		for (int x = 0; x < k; x++)
			for (int y = x + 1; y < num_of_data; y++) {
				if (trslt[x].dmin > trslt[y].dmin) {
					Pointlocation temprslt;

					temprslt.x = trslt[x].x;
					temprslt.y = trslt[x].y;
					temprslt.dmin = trslt[x].dmin;
					temprslt.dmax = trslt[x].dmax;

					trslt[x].x = trslt[y].x;
					trslt[x].y = trslt[y].y;
					trslt[x].dmin = trslt[y].dmin;
					trslt[x].dmax = trslt[y].dmax;

					trslt[y].x = temprslt.x;
					trslt[y].y = temprslt.y;
					trslt[y].dmin = temprslt.dmin;
					trslt[y].dmax = temprslt.dmax;
				}
			}
		sw2.stop();
		max_e->cmintime_sec += sw2.getDiff();

		//debug

		for (int l = 0; l < num_of_data; l++) {
			fprintf(result, "%f\t%f\t%lf\t%lf\n", trslt[l].x, trslt[l].y,
					trslt[l].dmin, trslt[l].dmax);
		}
		fclose(result);
		//

		//minmax_Filter

		//For debugging		
		result = fopen(
				"H:/Thesis code/GRPLBSQ-Nusrat/Results/result_minmaxfilter.txt",
				"w");
		if (result == NULL) {
			printf("Error writing rectfile\n");
		}
		//.....................................

		sw3.start();

		//update and find maxdist[k]
		for (int l = 0; l < k; l++)
			maxdist[l] = MAXDOUBLE;
		for (int l = 0; l < num_of_data; l++) {
			for (int y = 0; y < k; y++) {
				if (rslt[l].dmax < maxdist[y]) {
					for (int x = k - 1; x > y; x--) {
						maxdist[x] = maxdist[x - 1];
					}
					maxdist[y] = rslt[l].dmax;
					break;
				}
			}
		}

		for (int j = 0; j < g_size; j++) {
			/*
			 //debug
			 fprintf(result,"%d\n",num_of_data);
			 fprintf(result, "%f\n", maxdist[k-1]);
			 for(int l=0; l<num_of_data;l++)
			 {
			 fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			 }
			 fprintf(result, "\n\n\n");
			 //
			 */

			max_e->cnum_retrievals[j] += num_of_data;

			for (int l = 0; l < num_of_data; l++) {
				tdist = Dist1(rslt[l], p[j]);

				if (tdist > rslt[l].dmin)
					rslt[l].dmin = tdist;

			}

			/*
			 //debug
			 fprintf(result,"%d\n",num_of_data);
			 fprintf(result, "%f\n", maxdist[k-1]);
			 for(int l=0; l<num_of_data;l++)
			 {
			 fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			 }

			 //
			 */

			//local pruning
			int x = 0, y = 0;
			if (j < (g_size - 1) && num_of_data > k) {
				for (int l = 0; l < num_of_data; l++) {
					if (rslt[l].dmin <= maxdist[k - 1]) {
						if (l != x) {
							rslt[x].x = rslt[l].x;
							rslt[x].y = rslt[l].y;
							rslt[x].dmin = rslt[l].dmin;
							rslt[x].dmax = rslt[l].dmax;
						}
						x++;
					}
				}

				num_of_data = x;
			}
		}

		//sort	
		for (int x = 0; x < num_of_data; x++)
			for (int y = x + 1; y < num_of_data; y++) {
				if (rslt[x].dmin > rslt[y].dmin) {
					Pointlocation temprslt;

					temprslt.x = rslt[x].x;
					temprslt.y = rslt[x].y;
					temprslt.dmin = rslt[x].dmin;
					temprslt.dmax = rslt[x].dmax;

					rslt[x].x = rslt[y].x;
					rslt[x].y = rslt[y].y;
					rslt[x].dmin = rslt[y].dmin;
					rslt[x].dmax = rslt[y].dmax;

					rslt[y].x = temprslt.x;
					rslt[y].y = temprslt.y;
					rslt[y].dmin = temprslt.dmin;
					rslt[y].dmax = temprslt.dmax;
				}
			}

		sw3.stop();
		max_e->cminmaxtime_sec += sw3.getDiff();

		//debug
		fprintf(result, "%d\n", num_of_data);
		for (int l = 0; l < num_of_data; l++) {
			fprintf(result, "%f\t%f\t%lf\t%lf\n", rslt[l].x, rslt[l].y,
					rslt[l].dmin, rslt[l].dmax);
		}
		fclose(result);
		//

		// debug
		for (int l = 0; l < k; l++)
			if (trslt[l].x != rslt[l].x || trslt[l].y != rslt[l].y)
				fprintf(dFile, "\nFilter result mismatch error at %d", l);

		//

		printf("Iteration %d completed\n", i);

	}
	//debug

	//.....
	fclose(input1);
	fclose(input2);
	fclose(dFile);

	//Rtree
	delete rt;
	delete cache;

}

void rect_kGNN_query_sum(int n_sample, int g_size, int k, char *s1, char *s2,
		Exp_sum_stat *sum_e) {

	Rectangle1 r[MAX_GRPSIZE], m;

	Point2D p[MAX_GRPSIZE];

	//Open input1 File
	FILE *input1, *input2, *dFile;
	char temp[200];

	dFile = fopen("H:/Thesis code/GRPLBSQ-Nusrat/Results/InputFile/debugfilter",
			"a+");
	if (dFile == NULL) {
		printf("Error reading rectdata\n");
	}

	input1 = fopen(s1, "r");
	if (input1 == NULL) {
		printf("Error reading rectdata\n");
	}

	fgets(temp, 200, input1);
	//puts(temp);
	fgets(temp, 200, input1);
	//puts(temp);	

	input2 = fopen(s2, "r");
	if (input2 == NULL) {
		printf("Error reading point location\n");
	}

	fgets(temp, 200, input2);
	//puts(temp);
	fgets(temp, 200, input2);
	//puts(temp);	
	//

	//vector<Pointlocation> _rslt;
	Pointlocation rslt[MAXDATALIMIT];
	int num_of_data = 0;
	int blocksize = 1024;	//4096;

	// R-tree
	Cache *cache = new Cache(DEFAULT_C, blocksize);
	RTree *rt = new RTree(STREEFILE, cache);

	for (int i = 0; i < n_sample; i++) {
		//Experiment
		Stopwatch sw1, sw2, sw3;
		int last_pf = cache->page_faults;
		//..........

		fscanf(input1, "%f%f%f%f", &m.x1, &m.x2, &m.y1, &m.y2);

		for (int j = 0; j < g_size; j++) {
			//add here lines to take input data for source rectangle
			fscanf(input1, "%f%f%f%f", &r[j].x1, &r[j].x2, &r[j].y1, &r[j].y2);
			fscanf(input2, "%f%f", &p[j][0], &p[j][1]);
		}

		num_of_data = 0;
		sw1.start();
		rt->private_kGNN_sum(r, g_size, k, rslt, &num_of_data);
		sw1.stop();
		sum_e->stime_sec += sw1.getDiff();
		sum_e->snum_retrievals += num_of_data;
		sum_e->page_faults += cache->page_faults - last_pf;

		if (num_of_data >= MAXDATALIMIT)
			printf("\nFirstcheck:LIMIT EXCEEDED");

		//For debugging
		FILE * result;
		result = fopen(RESFILE, "w");
		if (result == NULL) {
			printf("Error writing rectfile\n");
		} else {
			// Process & close file
			//fprintf(rFile,"%d\n",k);
			fprintf(result, "%d\n", num_of_data);
			for (int l = 0; l < num_of_data; l++)
				//fprintf(result,"%.5f\t%.5f\t%.5f\n",_rslt[i].x, _rslt[i].y,Dist(_rslt[i],o));
				fprintf(result, "%.5f\t%.5f\t%.5f\t%.5f\n", rslt[l].x,
						rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			fclose(result);
		}
		//.....................................

		//preprocessing to match the data structure

		Pointlocation trslt[MAXDATALIMIT];
		for (int l = 0; l < num_of_data; l++) {
			trslt[l].x = rslt[l].x;
			trslt[l].y = rslt[l].y;
			trslt[l].dmin = rslt[l].dmin;
			trslt[l].dmax = rslt[l].dmax;
		}

		//min filter

		//For debugging		
		result = fopen(
				"H:/Thesis code/GRPLBSQ-Nusrat/Results/result_minfilter.txt",
				"w");
		if (result == NULL) {
			printf("Error writing rectfile\n");
		} else {
			fprintf(result, "%d\n", num_of_data);
		}
		//.....................................

		double maxdist[MAX_K];

		sw2.start();
		for (int j = 0; j < g_size; j++) {
			//update 
			for (int l = 0; l < num_of_data; l++) {
				trslt[l].dmin = trslt[l].dmin - MINRECTDIST1(trslt[l], r[j])
						+ Dist1(trslt[l], p[j]);
			}
		}

		//sort	
		for (int x = 0; x < k; x++)
			for (int y = x + 1; y < num_of_data; y++) {
				if (trslt[x].dmin > trslt[y].dmin) {
					Pointlocation temprslt;

					temprslt.x = trslt[x].x;
					temprslt.y = trslt[x].y;
					temprslt.dmin = trslt[x].dmin;
					temprslt.dmax = trslt[x].dmax;

					trslt[x].x = trslt[y].x;
					trslt[x].y = trslt[y].y;
					trslt[x].dmin = trslt[y].dmin;
					trslt[x].dmax = trslt[y].dmax;

					trslt[y].x = temprslt.x;
					trslt[y].y = temprslt.y;
					trslt[y].dmin = temprslt.dmin;
					trslt[y].dmax = temprslt.dmax;
				}
			}
		sw2.stop();
		sum_e->cmintime_sec += sw2.getDiff();

		//debug

		for (int l = 0; l < num_of_data; l++) {
			fprintf(result, "%f\t%f\t%lf\t%lf\n", trslt[l].x, trslt[l].y,
					trslt[l].dmin, trslt[l].dmax);
		}
		fclose(result);
		//

		//minmax_Filter

		//For debugging		
		result = fopen(
				"H:/Thesis code/GRPLBSQ-Nusrat/Results/result_minmaxfilter.txt",
				"w");
		if (result == NULL) {
			printf("Error writing rectfile\n");
		}
		//.....................................

		sw3.start();
		for (int j = 0; j < g_size; j++) {
			/*
			 //debug
			 fprintf(result,"%d\n",num_of_data);

			 for(int l=0; l<num_of_data;l++)
			 {
			 fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			 }
			 fprintf(result, "\n\n\n");
			 //
			 */

			sum_e->cnum_retrievals[j] += num_of_data;

			//update and find maxdist[k]
			for (int l = 0; l < k; l++)
				maxdist[l] = MAXDOUBLE;

			for (int l = 0; l < num_of_data; l++) {
				rslt[l].dmin = rslt[l].dmin - MINRECTDIST1(rslt[l], r[j])
						+ Dist1(rslt[l], p[j]);
				rslt[l].dmax = rslt[l].dmax - MAXRECTDIST1(rslt[l], r[j])
						+ Dist1(rslt[l], p[j]);

				for (int y = 0; y < k; y++) {
					if (rslt[l].dmax < maxdist[y]) {
						for (int x = k - 1; x > y; x--) {
							maxdist[x] = maxdist[x - 1];
						}
						maxdist[y] = rslt[l].dmax;
						break;
					}
				}
			}

			/*
			 //debug
			 fprintf(result,"%d\n",num_of_data);
			 fprintf(result, "%f\n", maxdist[k-1]);
			 for(int l=0; l<num_of_data;l++)
			 {
			 fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			 }

			 //
			 */

			//local pruning
			int x = 0, y = 0;
			if (j < (g_size - 1) && num_of_data > k) {
				for (int l = 0; l < num_of_data; l++) {
					if (rslt[l].dmin <= maxdist[k - 1]) {
						if (l != x) {
							rslt[x].x = rslt[l].x;
							rslt[x].y = rslt[l].y;
							rslt[x].dmin = rslt[l].dmin;
							rslt[x].dmax = rslt[l].dmax;
						}
						x++;
					}
				}

				num_of_data = x;
			}
		}

		//sort	
		for (int x = 0; x < num_of_data; x++)
			for (int y = x + 1; y < num_of_data; y++) {
				if (rslt[x].dmin > rslt[y].dmin) {
					Pointlocation temprslt;

					temprslt.x = rslt[x].x;
					temprslt.y = rslt[x].y;
					temprslt.dmin = rslt[x].dmin;
					temprslt.dmax = rslt[x].dmax;

					rslt[x].x = rslt[y].x;
					rslt[x].y = rslt[y].y;
					rslt[x].dmin = rslt[y].dmin;
					rslt[x].dmax = rslt[y].dmax;

					rslt[y].x = temprslt.x;
					rslt[y].y = temprslt.y;
					rslt[y].dmin = temprslt.dmin;
					rslt[y].dmax = temprslt.dmax;
				}
			}

		sw3.stop();
		sum_e->cminmaxtime_sec += sw3.getDiff();

		//debug
		fprintf(result, "\n\n%d\n", num_of_data);
		for (int l = 0; l < num_of_data; l++) {
			fprintf(result, "%f\t%f\t%lf\t%lf\n", rslt[l].x, rslt[l].y,
					rslt[l].dmin, rslt[l].dmax);
		}
		fclose(result);
		//

		// debug
		for (int l = 0; l < k; l++)
			if (trslt[l].x != rslt[l].x || trslt[l].y != rslt[l].y)
				fprintf(dFile, "\nFilter result mismatch error at %d", l);

		//

		printf("Iteration %d completed\n", i);

	}
	//debug

	//.....
	fclose(input1);
	fclose(input2);
	fclose(dFile);

	//Rtree
	delete rt;
	delete cache;

}

void rect_GNN_query_sum(int n_sample, int g_size, int k, char *s1,
		Exp_sum_stat *sum_e, int algonum, char res[]) {

	Rectangle1 s[2 * MAX_GRPSIZE], m, m1;
	Rectangle1 d[2 * MAX_GRPSIZE];
	Point2D p[2 * MAX_GRPSIZE];

	//Open input1 File
	FILE *input1, *input2, *dFile;
	char temp[200];

	dFile = fopen("H:/Thesis code/GRPLBSQ-Nusrat/Results/InputFile/debugfilter",
			"a+");
	if (dFile == NULL) {
		printf("Error reading rectdata\n");
	}

	input1 = fopen(s1, "r");
	if (input1 == NULL) {
		printf("Error reading rectdata\n");
	}

	fgets(temp, 200, input1);
	//puts(temp);
	fgets(temp, 200, input1);
	//puts(temp);	

	//vector<Pointlocation> _rslt;
	Pointlocation rslt[2 * kMAX];
	int num_of_data = 0;
	int blocksize = 1024;			//4096;
	Pair A[kMAX];

	// R-tree
	Cache *scache = new Cache(DEFAULT_C, blocksize);
	Cache *dcache = new Cache(DEFAULT_C, blocksize);
	RTree *srt = new RTree(STREEFILE, scache);
	RTree *drt = new RTree(DTREEFILE, dcache);

	Cache *scache1 = new Cache(DEFAULT_C, blocksize);
	Cache *scache2 = new Cache(DEFAULT_C, blocksize);
	Cache *scache3 = new Cache(DEFAULT_C, blocksize);
	RTree *srt1 = new RTree(STREEFILE1, scache1);
	RTree *srt2 = new RTree(STREEFILE2, scache2);
	RTree *srt3 = new RTree(STREEFILE3, scache3);

	for (int i = 0; i < n_sample; i++) {
		//Experiment
		Stopwatch sw1, sw2, sw3;
		int last_pf = scache->page_faults;
		int last_pf1 = dcache->page_faults;
		//..........

		fscanf(input1, "%f%f%f%f", &m.x1, &m.x2, &m.y1, &m.y2);
		//fscanf(input1,"%f%f%f%f",&m1.x1,&m1.x2,&m1.y1,&m1.y2);

		for (int j = 0; j < g_size; j++) {
			fscanf(input1, "%f%f%f%f", &s[j].x1, &s[j].x2, &s[j].y1, &s[j].y2);
			fscanf(input1, "%f%f%f%f", &d[j].x1, &d[j].x2, &d[j].y1, &d[j].y2);	//make changes in the input file to scan destination data
			//fscanf(input2,"%f%f",&p[j][0],&p[j][1]);
		}

		num_of_data = 0;

		sw1.start();
		char c;

		if (algonum == 1) {
			algonum = 1;
			srt->private_GNN_sum(s, d, g_size, rslt, &num_of_data, k, drt, A);
		}
		if (algonum == 2) {

			//srt->twoST_GTP_FA(s,d,g_size, rslt, &num_of_data,k,drt);
			srt2->twoST_GTP_FA(s, d, g_size, rslt, &num_of_data, k, srt3);
			algonum = 2;
		}
		if (algonum == 3) {

			srt1->threeST_GTP_FA(s, d, g_size, rslt, &num_of_data, k, srt2,
					srt3);
			algonum = 3;
		}
		if (algonum == 4) {

			srt1->private_3GNN_sum(s, d, g_size, rslt, &num_of_data, k, srt2,
					srt3, A);
			algonum = 4;
		}

		//srt->private_kGNN_sum(s,g_size,k,rslt, &num_of_data);
		//srt->private_kGNN_sum(s,g_size,k,rslt, &num_of_data,drt);

		sw1.stop();
		sum_e->stime_sec += sw1.getDiff();
		sum_e->snum_retrievals += k; //confusing
		sum_e->page_faults += scache->page_faults - last_pf;
		sum_e->page_faults += dcache->page_faults - last_pf1;

		if (num_of_data >= MAXDATALIMIT)
			printf("\nFirstcheck:LIMIT EXCEEDED");

		//For debugging
		FILE * result;
		result = fopen(res, "a+");
		if (result == NULL) {
			printf("Error writing res rectfile\n");
		} else {
			// Process & close file
			//fprintf(rFile,"%d\n",k);
			fprintf(result, "Sample Size = %d\n", i);
			fprintf(result, "%d\n", k);
			for (int l = 0; l < k; l++) {
				//fprintf(result,"%.5f\t%.5f\t%.5f\n",_rslt[i].x, _rslt[i].y,Dist(_rslt[i],o));
				if (algonum == 1)
					fprintf(result, "%.5f\t%.5f\t%.5f\t%.5f%\t%.5f\n",
							A[l].p1x1, A[l].p1y1, A[l].p2x2, A[l].p2y2,
							A[l].mindist);
				if (algonum == 2)
					fprintf(result, "%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
							rslt[l].x, rslt[l].y, rslt[l].dx, rslt[l].dy,
							rslt[l].dmin, rslt[l].dmax);
				if (algonum == 3)
					fprintf(result,
							"%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
							rslt[l].x, rslt[l].y, rslt[l].dx, rslt[l].dy,
							rslt[l].d1x, rslt[l].d1y, rslt[l].dmin,
							rslt[l].dmax);
			}
			fclose(result);
		}
		//.....................................

		printf("Iteration %d completed\n", i);

	}
	//debug

	//.....
	fclose(input1);
	//fclose(input2);
	fclose(dFile);

	//Rtree
	delete srt;
	delete drt;
	delete scache;
	delete dcache;

}

//Experiments
void summarize_output(char *s1, char *s2, Exp_sum_stat *e, int g) {
	//write in output file

	FILE * outputFile1, *outputFile2;
	outputFile1 = fopen(s1, "a+");
	outputFile2 = fopen(s2, "a+");
	char test;

	if (outputFile1 == NULL) {
		printf("Error writing output\n");
		//char s;
		//scanf("%c",&s);
	}

	fprintf(outputFile1,
			"%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", e->k,
			e->grpsize, e->m_area, e->r_area,
			e->stime_sec / (1.0 * CLOCKS_PER_SEC * SAMPLE),
			e->page_faults / SAMPLE, e->snum_retrievals / SAMPLE,
			e->cmintime_sec / (1.0 * CLOCKS_PER_SEC * SAMPLE),
			e->cminmaxtime_sec / (1.0 * CLOCKS_PER_SEC * SAMPLE));
	//for(int j=0; j<g; j++)
	//fprintf(outputFile2,"%d\t%.5lf\n", j, e->cnum_retrievals[j]/SAMPLE);

	fclose(outputFile1);
	fclose(outputFile2);
}

void exp_vary_k(char *d) {
	char res[100];
	char t_s[100];

	int algonum = 0;
	scanf("%d", &algonum);

	//generate input

	char s1[100], s2[100];

	for (int k = MIN_K; k <= MAX_K; k = k * INTERVAL_K) {
		printf("k = %d\n", k);
		strcpy(res, "C:\\result\\result_k");
		itoa(k, t_s, 10);
		strcat(res, t_s);
		char gap[] = "_";
		strcat(res, gap);
		itoa(algonum, t_s, 10);
		strcat(res, t_s);

		Exp_sum_stat sum_e;
		Exp_sum_stat max_e;

		//initialize input status

		for (int i = 0; i < SAMPLE; i++) {
			sum_e.k = k;
			max_e.k = k;
		}

		//make input filename
		strcpy(s1, "C:/GRPLBSQ/Results/InputFile/Default/fixed_1_default");
		strcpy(s2, "C:/GRPLBSQ/Results/InputFile/Default/fixed_2_default");

		//process
		rect_GNN_query_sum(1, DEFAULT_GRPSIZE, k, s1, &sum_e, algonum, res);
		//rect_kGNN_query_sum(1, DEFAULT_GRPSIZE, k, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &max_e);

		//make outputfile name

		strcpy(s1, "C:/GRPLBSQ/Results/OutputFile/K/sum_fixed_1_k_");
		strcpy(s2, "C:/GRPLBSQ/Results/OutputFile/K/sum_fixed_2_k_");
		itoa(k, t_s, 10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1, s2, &sum_e, DEFAULT_GRPSIZE);

		/*
		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/K/max_fixed_1_k_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/K/max_fixed_2_k_");
		 itoa (k, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		 printf("\n%d:1\n",k);
		 */
		//make input filename
		/*strcpy(s1,"C:/GRPLBSQ/Results/InputFile/Default/variable_1_default");	
		 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/Default/variable_2_default");

		 //process

		 rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &sum_e);
		 //rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &max_e);

		 //make outputfile name

		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/K/sum_variable_1_k_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/K/sum_variable_2_k_");
		 itoa (k, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);


		 summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		 printf("\n%d:2\n",k);
		 */
		/*
		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/K/max_variable_1_k_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/K/max_variable_2_k_");
		 itoa (k, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);


		 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		 printf("\n%d:2\n",k);
		 */
	}
}

void exp_vary_groupsize(char *d) {
	//generate input
	char s1[100], s2[100], t_s[100];
	char res[100];

	int algonum = 0;
	scanf("%d", &algonum);

	for (int g = MIN_GRPSIZE; g <= MAX_GRPSIZE; g = g * INTERVAL_GRPSIZE) {
		strcpy(res, "F:\\result_g");
		itoa(g, t_s, 10);
		strcat(res, t_s);
		char gap[] = "_";
		strcat(res, gap);
		itoa(algonum, t_s, 10);
		strcat(res, t_s);

		Exp_sum_stat sum_e;
		//Exp_sum_stat max_e;

		//initialize input status

		for (int i = 0; i < g; i++) {
			sum_e.grpsize = g;
			//max_e.grpsize=g;
		}

		//make input filename
		strcpy(s1, "C:/GRPLBSQ/Results/InputFile/GroupSize/fixed_1_g_");
		//strcpy(s2,"C:/GRPLBSQ/Results/InputFile/GroupSize/fixed_2_g_");	
		itoa(g, t_s, 10);
		strcat(s1, t_s);
		//strcat(s2, t_s);

		//process
		rect_GNN_query_sum(SAMPLE, g, DEFAULT_K, s1, &sum_e, algonum, res);
		//rect_kGNN_query_sum(SAMPLE, g, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, g, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name

		strcpy(s1, "C:/GRPLBSQ/Results/OutputFile/GroupSize/sum_fixed_1_g_");
		strcpy(s2, "C:/GRPLBSQ/Results/OutputFile/GroupSize/sum_fixed_2_g_");
		itoa(g, t_s, 10);
		//strcat(s2, t_s);

		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1, s2, &sum_e, g);
		//summarize_output(s1,&sum_e,g);

		printf("\n%d:1\n", g);

		/*strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/GroupSize/max_fixed_1_g_");	
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");
		 itoa (g, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);


		 summarize_output(s1,s2,&max_e,g);
		 printf("\n%d:1\n",g);*/

		//make input filename
		/*strcpy(s1,"C:/GRPLBSQ/Results/InputFile/GroupSize/variable_1_g_");	
		 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/GroupSize/variable_2_g_");
		 itoa (g, t_s,10);
		 strcat(s1, t_s);
		 strcat(s2, t_s);


		 //process

		 rect_kGNN_query_sum(SAMPLE, g, DEFAULT_K, s1, s2, &sum_e);
		 //rect_kGNN_query_max(SAMPLE, g, DEFAULT_K, s1, s2, &max_e);

		 //make outputfile name

		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/GroupSize/sum_variable_1_g_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/GroupSize/sum_variable_2_g_");
		 itoa (g, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&sum_e,g);
		 printf("\n%d:2\n",g);


		 /*strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/GroupSize/max_variable_1_g_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/GroupSize/max_variable_2_g_");
		 itoa (g, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&max_e,g);
		 printf("\n%d:2\n",g);*/

	}
}

void exp_vary_M_AREA(char *d) {
	char res[100];
	char t_s[100];

	int algonum = 0;
	scanf("%d", &algonum);

	//generate input
	char s1[100], s2[100];

	for (int m_area = (int) (MIN_M_AREAPART * 100);
			m_area <= (int) (MAX_M_AREAPART * 100);
			m_area = m_area * (INTERVAL_M_AREAPART)) {
		strcpy(res, "F:\\result_m");
		itoa(m_area, t_s, 10);
		strcat(res, t_s);
		char gap[] = "_";
		strcat(res, gap);
		itoa(algonum, t_s, 10);
		strcat(res, t_s);

		Exp_sum_stat sum_e;
		Exp_sum_stat max_e;

		//initialize input status

		for (int i = 0; i < SAMPLE; i++) {
			sum_e.m_area = m_area / 100.0;
			max_e.m_area = m_area / 100.0;
		}

		//make input filename
		strcpy(s1, "C:/GRPLBSQ/Results/InputFile/M_AREA/fixed_1_m_area_");
		strcpy(s2, "C:/GRPLBSQ/Results/InputFile/M_AREA/fixed_2_m_area_");
		itoa(m_area, t_s, 10);
		strcat(s1, t_s);
		strcat(s2, t_s);

		//process
		rect_GNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, &sum_e,
				algonum, res);
		//rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name

		strcpy(s1, "C:/GRPLBSQ/Results/OutputFile/M_AREA/sum_fixed_1_m_area_");
		strcpy(s2, "C:/GRPLBSQ/Results/OutputFile/M_AREA/sum_fixed_2_m_area_");
		itoa(m_area, t_s, 10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		//summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		summarize_output(s1, s2, &sum_e, DEFAULT_GRPSIZE);
		printf("\n%d:1\n", m_area);

		/*
		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/M_AREA/max_fixed_1_m_area_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/M_AREA/max_fixed_2_m_area_");
		 itoa (m_area, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		 printf("\n%d:1\n", m_area);
		 */

		//make input filename
		/*strcpy(s1,"C:/GRPLBSQ/Results/InputFile/M_AREA/variable_1_m_area_");	
		 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/M_AREA/variable_2_m_area_");
		 itoa (m_area, t_s,10);
		 strcat(s1, t_s);
		 strcat(s2, t_s);

		 //process

		 rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		 //rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		 //make outputfile name

		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/M_AREA/sum_variable_1_m_area_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/M_AREA/sum_variable_2_m_area_");
		 itoa (m_area, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		 printf("\n%d:2\n", m_area);

		 /*
		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/M_AREA/max_variable_1_m_area_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/M_AREA/max_variable_2_m_area_");
		 itoa (m_area, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		 printf("\n%d:2\n", m_area);
		 */
	}
}

void exp_vary_R_AREA(char *d) {
	//generate input
	char s1[100], s2[100], t_s[100];

	for (int r_area = (int) (MIN_R_AREAPART * 100000);
			r_area <= (int) (MAX_R_AREAPART * 100000);
			r_area = r_area + (int) (INTERVAL_R_AREAPART * 100000)) {
		Exp_sum_stat sum_e;
		Exp_sum_stat max_e;

		//initialize input status

		for (int i = 0; i < SAMPLE; i++) {
			sum_e.r_area = r_area / 100000.0;
			max_e.r_area = r_area / 100000.0;
		}

		//make input filename
		strcpy(s1, "C:/GRPLBSQ/Results/InputFile/R_AREA/fixed_1_r_area_");
		strcpy(s2, "C:/GRPLBSQ/Results/InputFile/R_AREA/fixed_2_r_area_");
		itoa(r_area, t_s, 10);
		strcat(s1, t_s);
		strcat(s2, t_s);

		//process
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name

		strcpy(s1, "C:/GRPLBSQ/Results/OutputFile/R_AREA/sum_fixed_1_r_area_");
		strcpy(s2, "C:/GRPLBSQ/Results/OutputFile/R_AREA/sum_fixed_2_r_area_");
		itoa(r_area, t_s, 10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		//summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		summarize_output(s1, s2, &sum_e, DEFAULT_GRPSIZE);
		printf("\n%d:1\n", r_area);

		/*
		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/R_AREA/max_fixed_1_r_area_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/R_AREA/max_fixed_2_r_area_");
		 itoa (r_area, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		 printf("\n%d:1\n", r_area);
		 */

		//make input filename
		/*strcpy(s1,"C:/GRPLBSQ/Results/InputFile/R_AREA/variable_1_r_area_");	
		 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/R_AREA/variable_2_r_area_");
		 itoa (r_area, t_s,10);
		 strcat(s1, t_s);
		 strcat(s2, t_s);

		 //process

		 rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		 //rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		 //make outputfile name

		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/R_AREA/sum_variable_1_r_area_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/R_AREA/sum_variable_2_r_area_");
		 itoa (r_area, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		 printf("\n%d:2\n", r_area);

		 /*
		 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/R_AREA/max_variable_1_r_area_");
		 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/R_AREA/max_variable_2_r_area_");
		 itoa (r_area, t_s,10);
		 strcat(s2, t_s);
		 strcat(s1, d);
		 strcat(s2, d);

		 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		 printf("\n%d:2\n", r_area);
		 */
	}
}

void exp_vary_DATASET(char *d) {

	char res[100];
	char t_s[100];

	int algonum = 0;
	scanf("%d", &algonum);

	//generate input
	char s1[100], s2[100];

	strcpy(res, "F:\\result_d");

	char gap[] = "_";
	strcat(res, gap);
	itoa(algonum, t_s, 10);
	strcat(res, t_s);

	Exp_sum_stat sum_e;
	Exp_sum_stat max_e;

	//initialize input status

	//make input filename
	strcpy(s1, "C:/GRPLBSQ/Results/InputFile/Default/fixed_1_default");
	strcpy(s2, "C:/GRPLBSQ/Results/InputFile/Default/fixed_2_default");

	//process
	rect_GNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, &sum_e, algonum,
			res);
	//rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
	//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

	//make outputfile name

	strcpy(s1,
			"C:/GRPLBSQ/Results/OutputFile/DatasetSize/sum_fixed_1_dataset_");
	strcpy(s2,
			"C:/GRPLBSQ/Results/OutputFile/DatasetSize/sum_fixed_2_dataset_");
	strcat(s1, d);
	strcat(s2, d);

	//summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
	summarize_output(s1, s2, &sum_e, DEFAULT_GRPSIZE);
	/*
	 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/DatasetSize/max_fixed_1_dataset_");
	 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/DatasetSize/max_fixed_2_dataset_");
	 strcat(s1, d);
	 strcat(s2, d);

	 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
	 */
	//make input filename
	/*	strcpy(s1,"C:/GRPLBSQ/Results/InputFile/Default/variable_1_default");	
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/Default/variable_2_default");

	 //process

	 rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);*/
	//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);
	//make outputfile name
	/*strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/DatasetSize/sum_variable_1_dataset_");
	 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/DatasetSize/sum_variable_2_dataset_");
	 strcat(s1, d);
	 strcat(s2, d);

	 summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);

	 /*
	 strcpy(s1,"C:/GRPLBSQ/Results/OutputFile/DatasetSize/max_variable_1_dataset_");
	 strcpy(s2,"C:/GRPLBSQ/Results/OutputFile/DatasetSize/max_variable_2_dataset_");
	 strcat(s1, d);
	 strcat(s2, d);

	 summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
	 */
}



//--------------------------GENERATE INPUT RECTANGLES-------------------------
void generate_g_rectangle(int n_sample, int g_size, double m_area_part,
		double m_ratio1, double m_ratio2, char *s1, char *s2) {
	Rectangle1 m[1000], r, d;

	float r_x, r_y, d_x, d_y, random;
	long int M = 10000;
	long double m_area = M * M * m_area_part;
	//long double r_area1 = M*M*r_area_part1;
	//long double r_area2 = M*M*r_area_part2;
	long double r_area, r_ratio, m_ratio;
	long double length, width;

	FILE * inputFile1, *inputFile2, *dFile;
	inputFile1 = fopen(s1, "w");
	inputFile2 = fopen(s2, "w");
	dFile = fopen("C:/GRPLBSQ/Results/InputFile/debug", "a+");

	if (inputFile1 == NULL || inputFile2 == NULL) {
		printf("Error writing rectdata\n");
	}

	fprintf(inputFile1,
			"Sample size:%.d\tMBR Area part:%.5f\tMBR Ratio1:%.5f\tMBR Ratio2:%.5f\n",
			n_sample, m_area, m_ratio1, m_ratio2);
	fprintf(inputFile1, "Number of rect:%.d\n", g_size);
	//fprintf(inputFile2,"Sample size:%.d\tMBR Area part:%.5f\tMBR Ratio:%.5f\tMBR Ratio1:%.5f\n",n_sample,m_area,m_ratio1, m_ratio2);
	//fprintf(inputFile2,"Number of rect:%.d\n",g_size);

	//srand((unsigned)time(0));
	srand(1000);

	//Find n_sample MBR for source
	for (int i = 0; i < n_sample;) {
		m_ratio = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		m_ratio = m_ratio1 + (m_ratio * (m_ratio2 - m_ratio1));

		r_x = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		r_y = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		r_x = (r_x * M);
		r_y = (r_y * M);
		width = 0;
		length = m_ratio * sqrt(m_area / m_ratio);
		if ((r_x - length / 2) < 0 || (r_x + length / 2) > M)
			continue;
		width = m_area / length;
		if ((r_y - width / 2) < 0 || (r_y + width / 2) > M)
			continue;

		m[i].x1 = r_x - length / 2;
		m[i].x2 = r_x + length / 2;
		m[i].y1 = r_y - width / 2;
		m[i].y2 = r_y + width / 2;

		//fprintf(inputFile1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",m[i].x1,m[i].x2,m[i].y1,m[i].y2,(m[i].x2-m[i].x1)*(m[i].y2-m[i].y1),m_area,(m[i].x2-m[i].x1)/(m[i].y2-m[i].y1),m_ratio,length, width);

		i++;
	}

	//Find g_size rectangles for each MBR
	for (int i = 0; i < n_sample; i++) {
		fprintf(inputFile1, "%f\t%f\t%f\t%f\n", m[i].x1, m[i].x2, m[i].y1,
				m[i].y2);

		for (int j = 0; j < g_size;) {

			r_x = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
			r_y = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
			r.x1 = m[i].x1 + (r_x * (m[i].x2 - m[i].x1));
			r.y1 = m[i].y1 + (r_y * (m[i].y2 - m[i].y1));
			r.x2 = m[i].x1 + (r_x * (m[i].x2 - m[i].x1));
			r.y2 = m[i].y1 + (r_y * (m[i].y2 - m[i].y1));

			d_x = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
			d_y = ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
			d.x1 = m[i].x1 + (d_x * (m[i].x2 - m[i].x1));
			d.y1 = m[i].y1 + (d_y * (m[i].y2 - m[i].y1));
			d.x2 = m[i].x1 + (d_x * (m[i].x2 - m[i].x1));
			d.y2 = m[i].y1 + (d_y * (m[i].y2 - m[i].y1));

			if (r.x1 < m[i].x1 || r.x2 > m[i].x2 || r.y1 < m[i].y1
					|| r.y2 > m[i].y2) {
				fprintf(dFile, "\nERROR at %s: %d and %d\n", i, j, s1);
			}
			/*if(r_x < r.x1 || r_x>r.x2 || r_y<r.y1 || r_y>r.y2)
			 {
			 fprintf(dFile,"\nCenter Error at %s: %d and %d",i,j,s1);
			 }*/

			/*r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
			 r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			 r_x = r.x1+(r_x * (r.x2-r.x1));
			 r_y = r.y1+(r_y * (r.y2-r.y1));

			 if(r_x < r.x1 || r_x>r.x2 || r_y<r.y1 || r_y>r.y2)
			 {
			 fprintf(dFile,"\nLocation Error at %s: %d and %d",i,j,s1);
			 }*/

			fprintf(inputFile1, "%f\t%f\t%f\t%f\n", r.x1, r.x2, r.y1, r.y2);
			fprintf(inputFile1, "%f\t%f\t%f\t%f\n", d.x1, d.x2, d.y1, d.y2);
			//fprintf(inputFile2,"%f\t%f\n",r_x,r_y);
			//fprintf(inputFile1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",r.x1,r.x2,r.y1,r.y2,r_x,r_y,(r.x2-r.x1)*(r.y2-r.y1),r_area,(r.x2-r.x1)/(r.y2-r.y1),r_ratio,length, width);
			j++;
		}
	}
	fclose(inputFile1);
	fclose(inputFile2);
	fclose(dFile);
}

void generate_input() {
	char s1[200], s2[200], t_s[100];

	//Vary group size
	for (int g = MIN_GRPSIZE; g <= MAX_GRPSIZE; g = g * INTERVAL_GRPSIZE) {
		strcpy(s1, "C:/GRPLBSQ/Results/InputFile/GroupSize/fixed_1_g_");
		strcpy(s2, "C:/GRPLBSQ/Results/InputFile/GroupSize/fixed_2_g_");
		itoa(g, t_s, 10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, g, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
		generate_g_rectangle(SAMPLE, g, DEFAULT_M_AREAPART, MIN_M_RATIO,
				MAX_M_RATIO, s1, s2);
	}
	/*for(int g=MIN_GRPSIZE; g<=MAX_GRPSIZE; g=g*INTERVAL_GRPSIZE)
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/GroupSize/variable_1_g_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/GroupSize/variable_2_g_");
	 itoa (g, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, g, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	 }*/

	//Vary M_AREA
	for (int m_area = (int) (MIN_M_AREAPART * 100);
			m_area <= (int) (MAX_M_AREAPART * 100);
			m_area = m_area * (INTERVAL_M_AREAPART)) {
		strcpy(s1, "C:/GRPLBSQ/Results/InputFile/M_AREA/fixed_1_m_area_");
		strcpy(s2, "C:/GRPLBSQ/Results/InputFile/M_AREA/fixed_2_m_area_");
		itoa(m_area, t_s, 10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, m_area/100.0,  DEFAULT_R_AREAPART, DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
		generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, m_area / 100.0,
				MIN_M_RATIO, MAX_M_RATIO, s1, s2);
	}
	/*for(int m_area=(int)(MIN_M_AREAPART*100); m_area<=(int)(MAX_M_AREAPART*100); m_area=m_area*(int)(INTERVAL_M_AREAPART))
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/M_AREA/variable_1_m_area_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/M_AREA/variable_2_m_area_");
	 itoa (m_area, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, m_area/100.0,  MIN_R_AREAPART,  MAX_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	 }*/

	/*
	 //Vary M_RATIO
	 for(int m_ratio=MIN_M_RATIO; m_ratio<=MAX_M_RATIO; m_ratio=m_ratio*INTERVAL_M_RATIO)
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/M_RATIO/fixed_1_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/M_RATIO/fixed_2_");
	 itoa (m_ratio, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE,DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, m_ratio, DEFAULT_R_RATIO, DEFAULT_R_RATIO, s1, s2);
	 }
	 for(int m_ratio=MIN_M_RATIO; m_ratio<=MAX_M_RATIO; m_ratio=m_ratio*INTERVAL_M_RATIO)
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/M_RATIO/variable_1_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/M_RATIO/variable_2_");
	 itoa (m_ratio, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, m_ratio, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	 }
	 */
	//Vary R_AREA
	/*for(int r_area=(int)(MIN_R_AREAPART*100000); r_area<=(int)(MAX_R_AREAPART*100000); r_area=r_area+(int)(INTERVAL_R_AREAPART*100000))
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/R_AREA/fixed_1_r_area_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/R_AREA/fixed_2_r_area_");
	 itoa (r_area, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  r_area/100000.0,  r_area/100000.0, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	 }*/
	/*for(int r_area=(int)(MIN_R_AREAPART*100000); r_area<=(int)(MAX_R_AREAPART*100000); r_area=r_area+(int)(INTERVAL_R_AREAPART*100000))
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/R_AREA/variable_1_r_area_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/R_AREA/variable_2_r_area_");
	 itoa (r_area, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  r_area/100000.0,  r_area/100000.0, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	 }*/
	/*
	 //Vary R_RATIO
	 for(int r_ratio=MIN_R_RATIO; r_ratio<=MAX_R_RATIO; r_ratio=r_ratio*INTERVAL_R_RATIO)
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/R_RATIO/fixed_1_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/R_RATIO/fixed_2_");
	 itoa (r_ratio, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE,DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, DEFAULT_M_RATIO, r_ratio, r_ratio, s1, s2);
	 }
	 for(int r_ratio=MIN_R_RATIO; r_ratio<=MAX_R_RATIO; r_ratio=r_ratio*INTERVAL_R_RATIO)
	 {
	 strcpy(s1,"C:/GRPLBSQ/Results/InputFile/R_RATIO/variable_1_");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/R_RATIO/variable_2_");
	 itoa (r_ratio, t_s,10);
	 strcat(s1, t_s);
	 strcat(s2, t_s);
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, DEFAULT_M_RATIO, r_ratio, r_ratio, s1, s2);
	 }
	 */
	//Default
	strcpy(s1,
			"H:/Thesis code/GRPLBSQ-Nusrat/Results/InputFile/Default/fixed_1_default");
	strcpy(s2,
			"H:/Thesis code/GRPLBSQ-Nusrat/Results/InputFile/Default/fixed_2_default");
	puts(s1);
	puts(s2);
	//generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,
			MIN_M_RATIO, MAX_M_RATIO, s1, s2);
	/*strcpy(s1,"C:/GRPLBSQ/Results/InputFile/Default/variable_1_default");
	 strcpy(s2,"C:/GRPLBSQ/Results/InputFile/Default/variable_2_default");
	 puts(s1);
	 puts(s2);
	 generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);*/

}

void createDataPointFile(){
	double temp,x1,x2,y1,y2;
	long id=0;
	char* dataPointFile="Datasets/Germany/hypsogr_0/hypsogr.txt";

	FILE * output = fopen(DATAFILE, "a");
	if (output == NULL) {
			printf("Cannot open file %s", DATAFILE);
	}
	ifstream iFile(dataPointFile);
	string line;

	while (getline(iFile, line)) {
		//cout << line << endl;
		 std::istringstream in(line);
		 in >> temp >> x1 >> y1 >> x2 >> y2;
		 fprintf(output,"%ld %lf %lf %lf %lf\n",id,x1,x1,y1,y1);
		 id++;
		 fprintf(output,"%ld %lf %lf %lf %lf\n",id,x2,x2,y1,y1);
		 id++;
		 fprintf(output,"%ld %lf %lf %lf %lf\n",id,x2,x2,y2,y2);
		 id++;
		 fprintf(output,"%ld %lf %lf %lf %lf\n",id,x1,x1,y2,y2);
		 id++;
	}
	iFile.close();
	fclose(output);


}
bool pointAtBoundary(float* point,float* mbr){
	
	bool atboundary=false;

	if(mbr[0]==point[0] && mbr[2]<=point[1] && mbr[3]>=point[1])
		atboundary=true;
	if(mbr[2]==point[1] && mbr[0]<=point[0] && mbr[1]>=point[0])
		atboundary=true;
	if(mbr[1]==point[0] && mbr[2]<=point[1] && mbr[3]>=point[1])
		atboundary=true;
	if(mbr[3]==point[1] && mbr[0]<=point[0] && mbr[1]>=point[0])
		atboundary=true;

	return atboundary;

	
}
bool intersects(float* mbr,float* bounces){
	
	bool inside;
    bool overlap;

    overlap = TRUE;
    inside = TRUE;

	for (int i = 0; i < 2; i++)
    {
    	//printf("if (mbr[%d] > bounces[%d] ||  mbr[%d] < bounces[%d]) == %f>%f || %f>%f   (overlap=false)\n",2*i,2*i+1,2 * i + 1,2*i,mbr[2*i],bounces[2*i+1],mbr[2 * i + 1],bounces[2*i]);
		if (mbr[2 * i] > bounces[2 * i + 1] ||  mbr[2 * i + 1] < bounces[2 * i])
			overlap = FALSE;
		//printf("if (mbr[%d] > bounces[%d] ||  mbr[%d] < bounces[%d]) == %f>%f || %f>%f   (inside=false)\n",2*i,2*i,2 * i + 1,2*i+1,mbr[2*i],bounces[2*i],mbr[2 * i + 1],bounces[2*i+1]);
		if (mbr[2 * i] < bounces[2 * i] || mbr[2 * i + 1] > bounces[2 * i + 1])
			inside = FALSE;
    }
	if(overlap || inside){
		return true;
	}
}

void generateIntersectFreePoints(){

	float r1[4];
	float r2[4];
	float point[2],mbr[4];

	
	int x,y;
	FILE *input1,*input2,*input3;
	input1 = fopen( "Datasets/sample.txt", "r");
	input3 = fopen( "Datasets/sample_point.txt", "a");
	int index=0;
	 for(int i=0; i<307992; i++)
	 {
		fscanf(input1,"%d%f%f%f%f",&x,&r1[0],&r1[1],&r1[2],&r1[3]);
		//printf("i %f,%f,%f,%f\n",r1[0],r1[1],r1[2],r1[3]);
		input2 = fopen( "Datasets/sample_mbr.txt", "r");
		bool noIntersect=true;
		for(int j=0; j<34334; j++)
		{
			
			fscanf(input2,"%d%f%f%f%f",&y,&r2[0],&r2[1],&r2[2],&r2[3]);
			//printf("j %f,%f,%f,%f",r2[0],r2[1],r2[2],r2[3]);
			//cout<<"\nIntersect? "<<intersects(r1,r2)<<" i= "<<i<<" j= "<<j<<endl;
			if(intersects(r1,r2))
			{
				point[0]=r1[0];point[1]=r1[2];
				if(!(r1[0]==r2[0] && r1[2]==r2[2] && r1[1]==r2[1] && r1[3]==r1[3])){
					if(! pointAtBoundary(point,r2)){
						noIntersect=false;
						break;
					}
					//else
					//	printf("\n at boundary\n");
					//fprintf(input3,"%d\t%d\t\t%f\t%f\t%f\t%f%s%d%s%f\t%f\t%f\t%f\n",x,y,r1[0],r1[2],r1[1],r1[3]," Intersects ",intersects(r1,r2),"  ",r2[0],r2[2],r2[1],r2[3]);
				}
				
			}
		}
		if(noIntersect){
			//printf("\n %d ",x);
			fprintf(input3,"%ld %lf %lf %lf %lf\n",index,r1[0],r1[1],r1[2],r1[3]);
			index++;		
		}
		fclose(input2);

		
	 }
	 fclose(input1);
	 fclose(input3);

}



void createPointMBRFile(){
	double temp,x1,x2,y1,y2;
	int id;

	char* dataPointFile="Datasets/Greece/cities_loc/cities_loc.txt";

	ifstream iFile(dataPointFile);

	FILE * output = fopen(DATAFILE, "a");
	if (output == NULL) {
			printf("Cannot open file %s", DATAFILE);
	}
	string line;

	while (getline(iFile, line)) {
		//cout << line << endl;
		 std::istringstream in(line);
		 in >> id >> x1 >> y1 ;
		 fprintf(output,"%ld %lf %lf %lf %lf\n",id,x1,x1,y1,y1);

	}
	iFile.close();
	fclose(output);




}


void createMBRFile(){
	double temp,x1,x2,y1,y2;
	int id;

	char* dataPointFile="Datasets/Germany/rrlines_0/rrlines.txt";

	ifstream iFile(dataPointFile);

	FILE * output = fopen(DATAFILE_MBR, "a");
	if (output == NULL) {
			printf("Cannot open file %s", DATAFILE);
	}
	string line;

	while (getline(iFile, line)) {
		//cout << line << endl;
		 std::istringstream in(line);
		 in >> id >> x1 >> y1 >> x2 >> y2;
		 fprintf(output,"%ld %lf %lf %lf %lf\n",id,x1,x2,y1,y2);

	}
	iFile.close();
	fclose(output);




}

void range_test(RTree* srt){

	float mbr[4];
	mbr[0]=2345-1000;
	mbr[1]=2345+1000;
	mbr[2]=9823-1000;
	mbr[3]=9823+1000;
	SortedLinList *res_list = new SortedLinList();
	srt -> rangeQuery(mbr, res_list);
	printf("Range Query returned %d entries\n",res_list->get_num());
	res_list->print();

	delete res_list;

}

//Experiments
void print_output(char *s1, Exp_stat *e) {
	//write in output file

	FILE * outputFile1;
	outputFile1 = fopen(s1, "a+");

	char test;

	if (outputFile1 == NULL) {
		printf("Error writing output\n");
		//char s;
		//scanf("%c",&s);
	}

	fprintf(outputFile1,
			"k = %d\tgrp size =%d\tQuery Area = %lf\ttime = %.5lf sec\tio_access_rtree=%.5lf\tio_access_rtree_obs=%.5lf\tp=%d\tVis=%.5lf sec\tSP=%.5lf sec\n",
			e->k,
			e->grpsize, 
			e->qArea,
			e->stime_sec / (1.0 * CLOCKS_PER_SEC * SAMPLE),
			e->io_access,
			e->io_access_obs,
			e->totalNumberOfPRetrieved,
			e->visGraphConstime_sec/ (1.0 * CLOCKS_PER_SEC * SAMPLE),
			e->shortestPathCalctime_sec/ (1.0 * CLOCKS_PER_SEC * SAMPLE)
);
	//for(int j=0; j<g; j++)
	//fprintf(outputFile2,"%d\t%.5lf\n", j, e->cnum_retrievals[j]/SAMPLE);

	fclose(outputFile1);
}

void exp_ognn_sum(float queryPoints[][2],int groupSize,int k,double kNearestNeighbor[][3],RTree *rt_obs,RTree *rt,Cache *cache_obs,Cache *cache){
	Stopwatch sw1;
	Exp_stat *sum_e=new Exp_stat();
	sum_e->k=k;
	sum_e->grpsize=groupSize;
	int last_pf = cache->page_faults;
	//..........

		
	sw1.start();
	OGNN *ognn = new OGNN();
	ognn->ognnMultiPointApproach(queryPoints,groupSize,k,kNearestNeighbor, rt_obs,rt);
	sw1.stop();
	sum_e->stime_sec += sw1.getDiff();
	sum_e->io_access += rt->io_access+rt_obs->io_access;
//	sum_e->page_faults += cache->page_faults - last_pf;


	print_output("result_sum.txt",sum_e);
	delete ognn;
	delete sum_e;
}

void generate_points_inside_rectangle_with_point(float r_x,float r_y,float length,float width){

			//The rectangle considering r_x and r_y as center
			float m[4],r2[4],r1[4];int total=0.0,id;
			m[0] = r_x - length / 2;
			m[1] = r_x + length / 2;
			m[2] = r_y - width / 2;
			m[3] = r_y + width / 2;
			int count=0;
			int group=1;

			while(1){

					float x= fmod(rand(),(m[1]-m[0] + 1)) + m[0];
					float y= fmod(rand(),(m[3]-m[2] + 1)) + m[2];
					r1[0]=x;r1[1]=x;r1[2]=y;r1[3]=y;

					//printf("Rect %lf,%lf,%lf,%lf\n",m[0],m[2],m[1],m[3]);
					//printf("%lf %lf\n",x,y);
					FILE *input1,*input2,*input3;
					input1 = fopen( "Datasets/sample_mbr.txt", "r");
					//input2 = fopen( "Result/Input/queryArea/vary_g_8_k_4_qa_01.txt", "a");
					input2 = fopen( "Result/Input/groupSize/vary_g_2_n_k_4_qa_005.txt", "a");

					
					bool intersect=false;
					for(int j=0; j<34334; j++)
					{				
						fscanf(input1,"%d%f%f%f%f",&id,&r2[0],&r2[1],&r2[2],&r2[3]);
						//printf("j %f,%f,%f,%f",r2[0],r2[1],r2[2],r2[3]);
						//cout<<"\nIntersect? "<<intersects(r1,r2)<<" i= "<<i<<" j= "<<j<<endl;
						if(intersects(r1,r2))
						{
							intersect = true;
							break;
						}					
					
					}
					if(intersect == false){
						fprintf(input2,"%ld %lf %lf\n",total,x,y);
						if(count==group){
							fclose(input2);
							fclose(input1);
							break;
						}
						count++;
						total++;
						
									
					}
					fclose(input2);
					fclose(input1);
			}

}

void generate_nonIntersectingQueryPoints(){
	float defaultQArea=.005,r_x,r_y;float length=MAXX*defaultQArea,width=MAXY*defaultQArea;
	srand(10000);
	int count=0;
	for(int i=0;i<1000;i++){
			//float x = (float)rand()/RAND_MAX*MAXX*0.005;
			//float y = (float)rand()/RAND_MAX*MAXY*.005;
			
			r_x = (float)rand()/RAND_MAX*10000;
			r_y = (float)rand()/RAND_MAX*10000;
				
			if ((r_x - length / 2) < 0 || (r_x + length / 2) > MAXX)
				continue;
			if ((r_y - width / 2) < 0 || (r_y + width / 2) > MAXY)
				continue;
			float r1[4],r2[4];
			r1[0]=r_x;r1[1]=r_x;r1[2]=r_y;r1[3]=r_y;
			bool intersect=false;int id;

			FILE *input1;
			input1 = fopen( "Datasets/sample_mbr.txt", "r");
					
			for(int j=0; j<34334; j++)
			{				
						fscanf(input1,"%d%f%f%f%f",&id,&r2[0],&r2[1],&r2[2],&r2[3]);
						//printf("j %f,%f,%f,%f",r2[0],r2[1],r2[2],r2[3]);
						//cout<<"\nIntersect? "<<intersects(r1,r2)<<" i= "<<i<<" j= "<<j<<endl;
						if(intersects(r1,r2))
						{
							intersect=true;
							break;
						}					
					
			}
			fclose(input1);
			if(intersect)
				continue;
			else{
					generate_points_inside_rectangle_with_point(r_x,r_y,length,width);
					if(count==100)
						break;
					count++;
				}
	}
}

void writeOutputInfile(int k,double kNearestNeighbor[][3],int groupSize,string s){
	FILE * outputFile1;
	outputFile1 = fopen("Result/ognnOutput", "a+");
	fprintf(outputFile1,"\n%s\n",s.c_str());
	fprintf(outputFile1,"\n-------------k=%d , group size=%d -----------------\n",k,groupSize);
		for(int index=0;index<k;index++){
		//printf("\n(%f,%f) has distance %f\n",queryPoints_sorted[j].queryPoints[0],queryPoints_sorted[j].queryPoints[1],queryPoints_sorted[j].distance);
			float* q = new float[3];
			q[0]=kNearestNeighbor[index][0];
			q[1]=kNearestNeighbor[index][1];
			q[2]=kNearestNeighbor[index][2];

			//printf("k=%d, (%f,%f) distance %lf \n",index,q[0], q[1],ognn_sorted[index].distance);
			fprintf(outputFile1,"k=%d, (%f,%f) distance %lf\n",index,q[0], q[1],q[2]);

			delete q;

	}

		fclose(outputFile1);
}

void exp_ognn(float queryPoints[][2],int groupSize,int k,double kNearestNeighbor[][3],RTree *rt_obs,RTree *rt,Cache *cache_obs,Cache *cache,int algoNum,float qArea){
	Stopwatch sw1;
	Exp_stat *sum_e=new Exp_stat();
	sum_e->k=k;
	sum_e->grpsize=groupSize;
	int last_pf = cache->page_faults;
	//..........

		
	sw1.start();
	rt->io_access=0.0;
	rt_obs->io_access=0.0;
	OGNN_GNN *ognn_gnn = new OGNN_GNN();
	char* fileName;
	if(algoNum==1){
	ognn_gnn->ognnUsingEGNN(queryPoints,groupSize,k,kNearestNeighbor, rt_obs,rt,0);
	//fileName="Result/ognnSumUsingEGNN.txt";
	fileName="Result/ognnSumUsingEGNN_obs2.txt";
	writeOutputInfile(k,kNearestNeighbor,groupSize,"\n********************OGNN-GNN-SUM*******************************");
	}
	if(algoNum==2){
	ognn_gnn->ognnSumUsingNN(queryPoints,groupSize,k,kNearestNeighbor, rt_obs,rt,0);
	//fileName="Result/ognnSumUsingNN.txt";
	fileName="Result/ognnSumUsingNN_obs2.txt";
	writeOutputInfile(k,kNearestNeighbor,groupSize,"\n********************OGNN-CENTROID-NN-SUM*******************************\n");
	}
	if(algoNum==3){
	ognn_gnn->ognnUsingEGNN(queryPoints,groupSize,k,kNearestNeighbor, rt_obs,rt,1);
	//fileName="Result/ognnMaxUsingEGNN.txt";
	fileName="Result/ognnMaxUsingEGNN_obs2.txt";
	writeOutputInfile(k,kNearestNeighbor,groupSize,"\n********************OGNN-GNN-MAX*******************************");
	}
	if(algoNum==4){
	ognn_gnn->ognnMaxUsingNN(queryPoints,groupSize,k,kNearestNeighbor, rt_obs,rt,1);
	//fileName="Result/ognnMaxUsingNN.txt";
	fileName="Result/ognnMaxUsingNN_obs2.txt";
	writeOutputInfile(k,kNearestNeighbor,groupSize,"\n********************OGNN-CENTROID-NN-MAX*******************************\n");
	}
	sw1.stop();
	if(rt_obs->io_access == 0)
		return;
	sum_e->stime_sec += sw1.getDiff();
	sum_e->io_access = rt->io_access;
	sum_e->io_access_obs = rt_obs->io_access;
	sum_e->qArea=qArea;
	sum_e->totalNumberOfPRetrieved=ognn_gnn->totalNumberOfPRetrieved,
	sum_e->visGraphConstime_sec=ognn_gnn->visGraphConsTime;
	sum_e->shortestPathCalctime_sec=ognn_gnn->shortestPathCalcTime;
//	sum_e->page_faults += cache->page_faults - last_pf;


	print_output(fileName,sum_e);
	delete ognn_gnn;
	delete sum_e;
}

/*	ognn_gnn->ognnUsingEGNN(queryPoints,3,4,kNearestNeighbor, rt_obs,rt,0);  Algo 1
	ognn_gnn->ognnSumUsingNN(queryPoints,3,4,kNearestNeighbor, rt_obs,rt,0); Algo 2
	ognn_gnn->ognnUsingEGNN(queryPoints,3,4,kNearestNeighbor, rt_obs,rt,1); Algo 3
	ognn_gnn->ognnMaxUsingNN(queryPoints,3,4,kNearestNeighbor, rt_obs,rt,1); Algo 4*/
void exp_ognn_varyk(){
	float queryPoints[8][2];
	int groupSize=8;
	float qArea=0.005;
	
	for (int k = 4; k <= 32; k = k*2 ) {
		for(int algo=2;algo<5;algo++){
			FILE *input1;
			input1 = fopen( "Result/Input/groupSize/vary_g_8_k_4_qa_005.txt", "r");
			for(int sample=0;sample<5;sample++){
				double kNearestNeighbor[64][3]; 
				int blocksize = 1024;			//4096;//1024;//4096;
				int b_length = 1024;
				int dimension = 2;
				
			    int x;
			    for(int i=0; i<groupSize; i++)
			    {
			        fscanf(input1,"%ld%f%f",&x,&queryPoints[i][0],&queryPoints[i][1]);
				}

				Cache *cache = new Cache(0, blocksize);
				Cache *cache_obs = new Cache(0,blocksize);

				RTree *srt = new RTree(TREEFILE,  cache);
				RTree *srt_obs = new RTree(TREEFILE_MBR, cache_obs);
				printf("\n------------  Group Size %d , k = %d , Algo Num = %d  ----------\n",groupSize,k,algo);
				exp_ognn(queryPoints,groupSize,k,kNearestNeighbor, srt_obs,srt,cache_obs,cache,algo,qArea);
				delete cache;
				delete srt;
				delete cache_obs;
				delete srt_obs;
			}
			fclose(input1);
		}
	}
}

void exp_ognn_varyGroupSize(){
	float queryPoints[64][2];
	float qArea=0.005;
	int k=4;
	char* str = "Result/Input/groupSize/vary_g_";
	
	for (int groupSize = 2; groupSize <= 32; groupSize = groupSize*2 ) {
		char dest[120];
		strcpy( dest, str );
		char integer_string[10];
		sprintf(integer_string, "%d", groupSize);
		strcat( dest, integer_string );
		strcat( dest, "_k_4_qa_005.txt" );
		for(int algo=1;algo<5;algo++){
			FILE *input1;
			//input1 = fopen( "Result/Input/groupSize/vary_g_8_k_4_qa_005.txt", "r");
			input1 = fopen( dest, "r");
			for(int sample=0;sample<2;sample++){
				double kNearestNeighbor[64][3]; 
				int blocksize = 1024;			//4096;//1024;//4096;
				int b_length = 1024;
				int dimension = 2;
				
			    int x;
			    for(int i=0; i<groupSize; i++)
			    {
			        fscanf(input1,"%ld%f%f",&x,&queryPoints[i][0],&queryPoints[i][1]);
				}

				Cache *cache = new Cache(0, blocksize);
				Cache *cache_obs = new Cache(0,blocksize);

				RTree *srt = new RTree(TREEFILE,  cache);
				RTree *srt_obs = new RTree(TREEFILE_MBR, cache_obs);
				printf("\n------------  Group Size %d , k = %d , Algo Num = %d  ----------\n",groupSize,k,algo);
				exp_ognn(queryPoints,groupSize,k,kNearestNeighbor, srt_obs,srt,cache_obs,cache,algo,qArea);
				delete cache;
				delete srt;
				delete cache_obs;
				delete srt_obs;
			}
			fclose(input1);
		}
	}
}

void exp_ognn_varyQueryArea(){
	float queryPoints[8][2];
	
	int k=4;
	char* str = "Result/Input/queryArea/vary_qa_";
	int groupSize = 8;
	
	for (int queryArea = 2; queryArea <= 10; queryArea = queryArea+2 ) {
		float qArea=(float)(queryArea)/(float)(1000);
		char dest[120];
		strcpy( dest, str );
		char integer_string[10];
		sprintf(integer_string, "%d", queryArea);
		strcat( dest, integer_string );
		strcat( dest, "_g_8_k_4.txt" );
		for(int algo=1;algo<5;algo++){
			FILE *input1;
			//input1 = fopen( "Result/Input/groupSize/vary_g_8_k_4_qa_005.txt", "r");
			input1 = fopen( dest, "r");
			for(int sample=0;sample<2;sample++){
				double kNearestNeighbor[64][3]; 
				int blocksize = 1024;			//4096;//1024;//4096;
				int b_length = 1024;
				int dimension = 2;
				
			    int x;
			    for(int i=0; i<groupSize; i++)
			    {
			        fscanf(input1,"%ld%f%f",&x,&queryPoints[i][0],&queryPoints[i][1]);
				}

				Cache *cache = new Cache(0, blocksize);
				Cache *cache_obs = new Cache(0,blocksize);

				RTree *srt = new RTree(TREEFILE,  cache);
				RTree *srt_obs = new RTree(TREEFILE_MBR, cache_obs);
				printf("\n------------Query Area %lf,  Group Size %d , k = %d , Algo Num = %d  ----------\n",qArea,groupSize,k,algo);
				exp_ognn(queryPoints,groupSize,k,kNearestNeighbor, srt_obs,srt,cache_obs,cache,algo,qArea);
				delete cache;
				delete srt;
				delete cache_obs;
				delete srt_obs;
			}
			fclose(input1);
		}
	}
}




float val(float oldVal,float min,float max){
	float newVal=(oldVal-min)/(max-min)*MAXX;
	return newVal;
}
//Normalizing real dataset into 0-10000
void normalizeMyRealDataset(){
	 FILE *input1,*input2,*input3,*input4;
	 input1 = fopen( "Datasets/sample_mbr_real.txt", "r");
	// input1 = fopen( "Datasets/sample_real.txt", "r");
	 input2 = fopen( "Datasets/sample_mbr.txt", "a");
	//  input2 = fopen( "Datasets/sample.txt", "a");
	
	 float r1[4];
	 int x;
     /*float max=0.0,min=100000000;
	 for(int i=0; i<307993; i++)
	 //for(int i=0; i<2755; i++)
	 {
		fscanf(input1,"%d%f%f%f%f",&x,&r1[0],&r1[1],&r1[2],&r1[3]);

		if(r1[0]<min) min=r1[0];if(r1[1]<min) min=r1[1];if(r1[2]<min) min=r1[2];if(r1[3]<min) min=r1[3];
		if(r1[0]>max) max=r1[0];if(r1[1]>max) max=r1[1];if(r1[2]>max) max=r1[2];if(r1[3]>max) max=r1[3];


	 }*/
	 float max=403.000000,min=9500.00000;
	 printf("min=%lf max=%lf\n",min,max);
	 fclose(input1);
	// input1 = fopen( "Datasets/sample_real.txt", "r");
	  input1 = fopen( "Datasets/sample_mbr_real.txt", "r");
	 for(int i=0; i<34334; i++)
	 //for(int i=0; i<2755; i++)
	 {
		fscanf(input1,"%d%f%f%f%f",&x,&r1[0],&r1[1],&r1[2],&r1[3]);
		fprintf(input2,"%ld %lf %lf %lf %lf\n",x,val(r1[0],min,max),val(r1[1],min,max),val(r1[2],min,max),val(r1[3],min,max));
		
	 }
	 fclose(input1);
	 fclose(input2);
	
}
//----------------------------------- main -----------------------------------
int main(int argc, char* argv[]) {

	//Create an RTree
	int blocksize = 1024;			//4096;//1024;//4096;
	int b_length = 1024;
	int dimension = 2;

	//createPointMBRFile();

	Cache *cache = new Cache(0, blocksize);
	//normalizeMyRealDataset();
	//generateIntersectFreePoints();
	//generate_nonIntersectingQueryPoints();


	//Create sample input file
	//createDataPointFile();
	//createMBRFile();
	//Uncomment this when you have changed your dataset , otherwise no need to build the rtree everytime
	//RTree *rt = new RTree(DATAFILE, TREEFILE, b_length, cache, dimension);
	//rt->print_tree();
	//delete rt;
	//RTree *srt = new RTree(TREEFILE,  cache);
	//srt->print_tree();

	Cache *cache_obs = new Cache(0,blocksize);
	//Uncomment this when you have changed your dataset , otherwise no need to build the rtree everytime
	//RTree *rt_obs = new RTree(DATAFILE_MBR,TREEFILE_MBR,b_length,cache_obs,dimension);
	//rt_obs->print_tree();
	//range_test(srt_obs);
	//delete rt_obs;

	RTree *srt = new RTree(TREEFILE,  cache);
	RTree *srt_obs = new RTree(TREEFILE_MBR, cache_obs);
	//range_test(srt_obs);
	//srt_obs->print_tree();
	//srt->print_tree();

	/*Point2D queryPoints[1]; int numOfQueryPoints=1;
    double kNearestNeighbor[2][3];

	float centroid[1][2];
	centroid[0][0]=1006.331543;
	centroid[0][1]=2496.926758;

	//printf("\n----------------------------------------Searching for k-NN---------------------------------------\n");
	
	//kNearestNeighbour holds the kGNN Euclidean
	srt_obs->Point_BFN_kGNNQ(centroid, 1, kNearestNeighbor,1,1);
	delete srt_obs;
	printf("\n%lf %lf %lf \n",kNearestNeighbor[0][0],kNearestNeighbor[0][1],kNearestNeighbor[0][2]);
	srt_obs = new RTree(TREEFILE_MBR, cache_obs);
	double obstacle[5];
	float c[2];
	c[0]=1006.331543;
	c[1]=2496.926758;
	srt_obs->Rectangle_BFN_NNQ(c, obstacle);
	printf("\n%lf %lf %lf %lf %lf",obstacle[0],obstacle[1],obstacle[2],obstacle[3],obstacle[4]);
*/

	float m[2];
	/*m[0]=415171 ;
	m[1]=4543907;*/

	m[0]=115.4;
	m[1]=213.2;

	double nearestNeighbor[5];

	//srt->print_tree();
	//Nearest Neighbor query
	/*srt_obs->Rectangle_BFN_NNQ(m, nearestNeighbor);
	printf("Nearest Neighbor of (%f,%f), is (%f,%f),(%f,%f)\n", m[0], m[1],
			nearestNeighbor[0], nearestNeighbor[2],nearestNeighbor[1], nearestNeighbor[3]);
	for(int i=0;i<70;i++){
		if(! srt_obs->retrieve_kth_BFN_Rectangle_NNQ(nearestNeighbor,m))
			printf("kth NN is (%f,%f),(%f,%f)\n",
			nearestNeighbor[0], nearestNeighbor[2],nearestNeighbor[1], nearestNeighbor[3]);
		else break;
		}
*/

//	range_test(rt);

	/*int group_size=3;
	float queryPoints[3][2];
	queryPoints[0][0] = 506364;
	queryPoints[0][1] = 4583290;
	queryPoints[1][0] = 501098;
	queryPoints[1][1] = 4582444;
	queryPoints[2][0] = 501774;
	queryPoints[2][1] = 4581595;*/

	float queryPoints[3][2];
	queryPoints[0][0] = 10;
	queryPoints[0][1] = 10;
	queryPoints[1][0] = 11;
	queryPoints[1][1] = 11;
	queryPoints[2][0] = 8;
	queryPoints[2][1] = 8;

	/*ObstructedDistanceCentroid* oDist = new ObstructedDistanceCentroid();
	VisibilityGraph* initialVisGraph = new VisibilityGraph();
	oDist->writeQueryPointsInFile(queryPoints,3);
	oDist->constructInitialVisGraph(initialVisGraph);
	double obsDist = oDist->computeAggObstructedDistance(initialVisGraph,m,queryPoints,3,srt_obs,1);
	printf("New odist %lf\n",obsDist);
	float m1[2];
	m1[0]=15;
	m1[1]=23;
	obsDist = oDist->computeAggObstructedDistance(initialVisGraph,m1,queryPoints,3,srt_obs,1);
	printf("New odist %lf\n",obsDist);
	delete initialVisGraph;
	delete srt_obs->rectangleNNHeap;

	ObstructedDistance* oDistP = new ObstructedDistance();
	initialVisGraph = new VisibilityGraph();
	oDistP->writeQueryPointsInFile(queryPoints,3);
	oDistP->constructInitialVisGraph(initialVisGraph);
	obsDist = oDistP->computeAggObstructedDistance(initialVisGraph,m1,queryPoints,3,srt_obs,1);
	printf("Old odist %lf\n",obsDist);
*/

	/*ObstructedDistance* oDist = new ObstructedDistance();
	oDist->computeAggObstructedDistance(new VisibilityGraph(),m,queryPoints,3,rt_obs,0);
	double kNearestNeighbor[20][3]; //
	OGNN_GNN *ognn_gnn = new OGNN_GNN();
	ognn_gnn->ognnUsingEGNN(queryPoints,3,4,kNearestNeighbor, srt_obs,srt,0);
	ognn_gnn->ognnSumUsingNN(queryPoints,3,4,kNearestNeighbor, srt_obs,srt,0);
	ognn_gnn->ognnUsingEGNN(queryPoints,3,4,kNearestNeighbor, srt_obs,srt,1);
	ognn_gnn->ognnMaxUsingNN(queryPoints,3,4,kNearestNeighbor, srt_obs,srt,1);
*/

	
	//delete srt->kGNNHeap;

	/*float q[1][2];
	q[0][0] = 12;
	q[0][1] = 12;
	double kNN[20][3]; 
	rt->Point_BFN_kGNNQ(q,5,kNN,1);
	for(int i=0;i<5;i++){
		printf("%d-Group Nearest Neighbor of is ------(%f,%f) is dist %lf------\n", i,
			kNN[i][0], kNN[i][1],kNN[i][2]);
	}

	double kNNP[3];
	for(int i=0;i<6;i++){
		rt->retrieve_kth_BFN_GNNQ(kNNP,q,1);
	//srt->print_tree();
		printf("\n Next %d- Group Nearest Neighbor of (%f,%f) dist %lf",i,kNNP[0],kNNP[1],kNNP[2]);
	}
*/

	//change k
	/*for (int k = 16; k <= 16; k = k + 1) {
		for(int algo=1;algo<5;algo++){
			double kNearestNeighbor[20][3]; //
			printf("\n------------  Group Size %d , k = %d   ----------\n",group_size,k);
			exp_ognn(queryPoints,group_size,k,kNearestNeighbor, srt_obs,srt,cache_obs,cache,algo);
			//delete srt->kGNNHeap;
		}
	}*/

	/*SortedLinList *res_list = new SortedLinList();
	float mbr[4];
	mbr[0]=852.378845-100;
	mbr[1]=852.378845+100;
	mbr[2]=8286.251953-100;
	mbr[3]=8286.251953+100;
	srt_obs->rangeQuery(mbr, res_list);
	res_list->print();
*/
	exp_ognn_varyk();
	exp_ognn_varyGroupSize();
	exp_ognn_varyQueryArea();


	//delete cache;
//	delete srt;
	//delete rt;



	//delete cache_obs;
//	delete srt_obs;
	//delete rt_obs;



	

	

	//delete rt_obs;
//	delete ognn_gnn;

	//generate_input();
	//exp_vary_k("C");
	//exp_vary_groupsize("Z");
	//exp_vary_M_AREA("Z");
	//exp_vary_R_AREA("Z");
	//exp_vary_DATASET("Z");

	//char s1[100], s2[100], ws[100], s3[100], s4[100];
	//Exp_sum_stat max_e;
	//strcpy(s1, "Results/InputFile/input1");
	//strcpy(s2, "Results/InputFile/input2");

	//RTree *rt = new RTree(TREEFILE, b_length, cache, dimension);

	//generate_g_rectangle(5, 5, 0.08, MIN_R_AREAPART, MAX_R_AREAPART, 1,8, 1,8, s1, s2);

	//rect_kGNN_query_max(11, 1024, 8, s1, s2, &max_e);
	//rect_kGNN_query_sum(10, 50, 8, s1, s2, &max_e);
	//check_with_range_query(50, s1, s2);
	//check_validity_rectangle(8);
	//check_validity_filter_and_point(32,100,s1,s2);
	//check_validity_point(8,50,s1,s2);

	//Snapshot queries

	//exp_staticquery();

	/*	int blocksize = 1024;//4096;//1024;//4096;
	 int b_length = 1024;


	 Cache *cache = new Cache(0, blocksize);
	 int dimension=2;

	 RTree *rt = new RTree(DATAFILE, TREEFILE, b_length, cache, dimension);
	 delete rt;
	 delete cache;
	 //RTree *rt = new RTree(TREEFILE, b_length, cache, dimension);



	 /*
	 Rectangle1 r[MAX_GRPSIZE],m;
	 Point2D p[MAX_GRPSIZE];

	 FILE *input1, *input2, *input3, *input4;
	 char temp[200];

	 input3 = fopen( s3, "w");
	 if (input3 == NULL)
	 {
	 printf("Error reading rectdata\n");
	 }

	 input4 = fopen( s4, "w");
	 if (input4== NULL)
	 {
	 printf("Error reading rectdata\n");
	 }

	 input1 = fopen( s1, "r");
	 if (input1 == NULL)
	 {
	 printf("Error reading rectdata\n");
	 }

	 fgets(temp,200,input1);
	 //puts(temp);
	 fgets(temp,200,input1);
	 //puts(temp);

	 input2 = fopen( s2, "r");
	 if (input2 == NULL)
	 {
	 printf("Error reading point location\n");
	 }

	 fgets(temp,200,input2);
	 //puts(temp);
	 fgets(temp,200,input2);
	 //puts(temp);
	 //



	 for(int i=0; i<1000; i++)
	 {
	 fscanf(input1,"%f%f%f%f",&m.x1,&m.x2,&m.y1,&m.y2);
	 for(int j=0; j<1024; j++)
	 {
	 fscanf(input1,"%f%f%f%f",&r[j].x1,&r[j].x2,&r[j].y1,&r[j].y2);
	 fscanf(input2,"%f%f",&p[j][0],&p[j][1]);

	 fprintf(input3,"%f\t%f\t%f\t%f\n",r[j].x1,r[j].x2,r[j].y1,r[j].y2);
	 fprintf(input4,"%f\t%f\n",p[j][0],p[j][1]);
	 if(p[j][0]<r[j].x1 || p[j][0]>r[j].x2 || p[j][1] <r[j].y1 || p[j][1]>r[j].y2)
	 {
	 printf("\nError");
	 printf("\n\nError at %d:%d\n", i, j);
	 printf("\n\nError at %d:%d\n", i, j);
	 }

	 }
	 fprintf(input3,"END\n\n\n");
	 fprintf(input4,"END\n\n\n");
	 }
	 fclose(input1);
	 fclose(input2);
	 fclose(input3);
	 fclose(input4);

	 */

	return 0;
}

