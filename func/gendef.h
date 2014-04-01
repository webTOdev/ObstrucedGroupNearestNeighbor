#ifndef __GENERAL_DEFINITION
#define __GENERAL_DEFINITION

#include <stdio.h>
#include <ctype.h>

//----Valdity region query----------------------------------
#define MAX_REAL 1e+20;


//----BlockFile, CachedBlockFile, Cache---------------------
#define BFHEAD_LENGTH (sizeof(int)*2)    //file header size

#define TRUE 1
#define FALSE 0

#define SEEK_CUR 1
#define SEEK_SET 0
#define SEEK_END 2

typedef char Block[];
//-------------------All-------------------------------------
#define MAXREAL         1e20
#define FLOATZERO       1e-8
#define MAX_DIMENSION   256

#define DIMENSION 2

#define TRUE 1
#define FALSE 0

#define minimum(a, b) (((a) < (b))? (a) : (b)  )
#define maximum(a, b) (((a) > (b))? (a) : (b)  )
//-------------------Class and other data types--------------
class BlockFile;  //for BlockFile definition
class Cache;
class Cacheable   //inherit this class if you wish to design an external
                  //memory structure that can be cached
{
public:
	BlockFile *file;
	Cache *cache;
};
  //==========================================================
class CmdIntrpr  //this is the class of command interpretor.  for a new rtree decescendent
                  //inherit this class to obtain tailored command interpretor
{
public:
	int cnfrm_cmd(char *_msg)
	{ char c = ' ';
	  while (c != 'N' && c != 'Y')
	  { printf("%s (y/n)?", _msg);
	    c = getchar(); 
		char tmp;
		while ((tmp = getchar()) != '\n');
		c = toupper(c); }
	  if (c == 'N') return 0; else return 1; }
  
	void get_cmd(char *_msg, char *_cmd)
	{ printf("%s", _msg);  
	  char *c = _cmd;
	  while (((*c) = getchar()) != '\n')
	    c++;
	  *c = '\0'; } 

	virtual bool build_tree(char *_tree_fname, char *_data_fname, int _b_len, int _dim, int _csize) = 0;
	virtual void free_tree() = 0;
	virtual int qry_sngle(float *_mbr, int *_io_count) = 0;
	/*
	virtual int qry_wrkld(char *_wrkld_fname, int *_io_count, bool _display) = 0;
	*/
	virtual void run() = 0;
	virtual void version() = 0;
};
  //==========================================================
enum SECTION {OVERLAP, INSIDE, S_NONE};
enum R_OVERFLOW {SPLIT, REINSERT, NONE};
enum R_DELETE {NOTFOUND,NORMAL,ERASED};
typedef float *floatptr;
  //==========================================================
struct SortMbr
{
    int dimension;
    float *mbr;
    float *center;
    int index;
};

struct BranchList
{
    int entry_number;
    float mindist;
    float minmaxdist;
    bool section;
};

//Added by Tanzima
typedef float Point2D[DIMENSION];
struct Pointlocation
{
	float x;
	float y;
	float dx;
	float dy;
	float d1x;
	float d1y;
	long double dmin;
	long double dmax;
};

struct DistfromPOI
{
	double d[4];
};

//added by TANZIMA

struct Rectangle1
{
		float x1, y1, x2, y2; // corners' co-ordinates
		//Point2D o; // center
};
//----------------------------------------------------------
//added by TAHRIMA

struct Pair
{
	float p1x1, p1y1, p2x2, p2y2,p3x3, p3y3;
	long double mindist;
};

//------------------------------------------------------------
//-----Global Functions--------------------------------------
void error(char *_errmsg, bool _terminate);

float area(int dimension, float *mbr);
float margin(int dimension, float *mbr);
float overlap(int dimension, float *r1, float *r2);
float* overlapRect(int dimension, float *r1, float *r2);
float objectDIST(float *p1, float *p2, int dim);
float MINMAXDIST(float *Point, float *bounces);
float MINDIST(float *Point, float *bounces, int dim);

bool inside(float &p, float &lb, float &ub);
void enlarge(int dimension, float **mbr, float *r1, float *r2);
bool is_inside(int dimension, float *p, float *mbr);
int pruneBrunchList(float *nearest_distanz, const void *activebrunchList, 
		    int n);
bool section(int dimension, float *mbr1, float *mbr2);
bool section_c(int dimension, float *mbr1, float *center, float radius);

int sort_lower_mbr(const void *d1, const void *d2);
int sort_upper_mbr(const void *d1, const void *d2);
int sort_center_mbr(const void *d1, const void *d2);
int sortmindist(const void *element1, const void *element2);

//---added for validity region---------------------------------
float distinit(float *p, float *e);
float distcont(float *l0, float *l1, float *e);
SECTION section_new(int dimension, float *p, float *q);

// added by Tanzima
float Dist(Pointlocation p, Point2D q);
float Dist1(Pointlocation p, Point2D q);
float Dist(Point2D p, Point2D q);
// added by Tahrima
float Dist(Pointlocation p, Pointlocation q);

//Added by Tanzima for kGNN
float MINRECTDIST(Rectangle1 p, float *bounces);
float MINRECTDIST(float *bounces1, float *bounces);
float MAXRECTDIST(Rectangle1 p, float *bounces);
float MAXRECTDIST(float *bounces1, float *bounces);
float MINRECTDIST1(Pointlocation p, Rectangle1 bounces);
float MAXRECTDIST1(Pointlocation p, Rectangle1 bounces);

#ifdef UNIX
void strupr(char *_msg);
#endif
//-----------------------------------------------------------
#endif