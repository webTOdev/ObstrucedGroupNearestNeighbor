/* rtree.h
   this file defines the class RTree*/

#ifndef __RTREE
#define __RTREE
//------------------------------------------------------------
#include "../func/gendef.h"
#include "../heap/heap.h"

//Added by Tanzima
//#include "../global.h"
//...............

//Added by Tanzima
#include <vector>
using namespace std;

/*#define NOMINMAX
#include <windows.h>
*/
#include <limits>



//added for TP KNN------------------------------------------------
#define SEQUENCE_SENSITIVE false
  //set it true if you want a new partition point every time the order of the NNs change
#define PAST_PAIR 1000

//Added for kGNN by Tanzima


const double MAXDOUBLE = numeric_limits<double>::infinity();
const int kMAX=1000;
const int MAXDATALIMIT=20000;
//------------------------------------------------------------
class LinList;
class SortedLinList;
class Cache;
class RTNode;
class Entry;
//Added by Tanzima
struct Rectangle1;
struct DistfromPOI;
//------------------------------------------------------------
class RTree : public Cacheable
{
public:
//--===on disk===--
	int dimension;                       
	int num_of_data;	                 
    int num_of_dnodes;	                 
    int num_of_inodes;	                 
	int root;                            
	bool root_is_data;                   
//--===others===--
	RTNode *root_ptr;
    bool *re_level;  
    LinList *re_data_cands; 
	LinList *deletelist;

//--===added for TP KNN===--
	int last_pair[PAST_PAIR][2]; //records the pairs of points that produce the minimum inflence time
	int lastcnt; //next last pair to be replaced
	Heap *tpheap;

//--===functions===--
	RTree(char *fname, int _b_length, Cache* c, int _dimension);
    RTree(char *fname, Cache* c);
    RTree(char *inpname, char *fname, int _blength, Cache* c, int _dimension);
    ~RTree();
	void del_root();
	bool delete_entry(Entry *d);
	bool FindLeaf(Entry *e);
    int get_num() { return num_of_data; }
	void insert(Entry *d);
	void load_root();  
	void rangeQuery(float *mbr, SortedLinList *res);
	void read_header(char *buffer);      
	void write_header(char *buffer);
	int update_rslt(Entry *_e, float _dist, Entry *_rslt, 
					 float *_key, int _k);
	//Restored by Nusrat
	void print_tree();

	// This function was added to perform TP-kNN queries by Bobby
	void TPNN_TP(float *_qline, int _k, Entry *_nn, Entry *_rslt, float _max_trvl);
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
	void BFNN(float *_qpt, int _k, Entry *_rslt);
	void BFNNCont(float *qmbr, float *qmbr2, Heap *heap, HeapEntry *e, int k);
	//Added by Tanzima
	float KMin(Point2D m, int c, float cl, int k, Pointlocation _rslt[], int *num_of_data,float d_safe_corner,DistfromPOI _cornertoPOI[]);
	void UpdateCount(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, float r, int *count, DistfromPOI _cornertoPOI[]);
	int UpdateStatus(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, int *count, DistfromPOI _cornertoPOI[]);
	void Rect_kNNQ(Rectangle1 R, float cl, int k, Pointlocation *_rslt, int *num_of_data);
	void Point_BFN_NNQ(Point2D o, double *_rslt);

	//Added by Tanzima for kGNN
	float KMax(int k, Pointlocation _rslt[], int *num_of_data);
	void private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation _rslt[], int *num_of_data);
	void private_kGNN_max(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data);

	void private_GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size, Pointlocation rslt[], int *num_of_data, int k ,RTree *drt,Pair A[]);
	void twoS_GTP_SA(Rectangle1 S[],Rectangle1 D[], int g_size, Pointlocation rslt[], int *num_of_data);
	void private_kGNN_sum(Rectangle1 D[], int g_size, int k, Pointlocation rslt[], int *num_of_data,RTree *drt);
	void twoST_GTP_FA(Rectangle1 S[],Rectangle1 D[], int g_size, Pointlocation rsltp1p2[], int *num_of_data_p1p2,int k,RTree *dtr);
	void threeST_GTP_FA(Rectangle1 S[],Rectangle1 D[], int g_size, Pointlocation rsltp1p2[], int *num_of_data_p1p2,int k,RTree *dtr,RTree *dtr1);
	RTNode* getRTNode(int ref);
	double approximate_mindist(Rectangle1 S[],Rectangle1 D[], int g_size,int k);
	double approximate_mindist1(Rectangle1 S[],Rectangle1 D[], int g_size,int k);
	void private_GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size,Pointlocation rslt[], int *num_of_data_P1,Pointlocation rsltp2[], int *num_of_data_P2,int k,RTree *drt1,Pair A[]);
	void private_3GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size,Pointlocation rslt[], int *num_of_data_P1,int k,RTree *drt,RTree *drt1,Pair A[]);

	//Added by Nusrat for kGNN Euclidean
	void Point_BFN_GNNQ(Point2D o[], double *_rslt,int numOfQueryPoints);
	void Point_BFN_kGNNQ(Point2D o[],int k,double _rslt[][2],int numOfQueryPoints);
	//void retrieve_kth_BFN_GNNQ( double *_rslt,Point2D o[],int numOfQueryPoints);
	void retrieve_kth_BFN_GNNQ( double *_rslt,Point2D o[],int numOfQueryPoints);

	int io_access;
			//These two entries to store the current state of the heap for incrementally retrieving point nn
	Heap *kGNNHeap;
	int latestSon ;

	//These two entries to store the current state of the heap for incrementally retrieving Rectangle nn
	Heap *rectangleNNHeap;
	int latestSonForRectangle ;
	void Rectangle_BFN_NNQ(Point2D o, double *_rslt);
	void retrieve_kth_BFN_Rectangle_NNQ( double *_rslt, Point2D o);
};

#endif // __RTREE
