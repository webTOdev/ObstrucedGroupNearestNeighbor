/*rtree.cpp
  this file implements the RTree class*/
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "rtree.h"
#include "entry.h"
#include "rtnode.h"
#include "distance.h"
#include "../blockfile/cache.h"
#include "../blockfile/blk_file.h"
#include "../linlist/linlist.h"


//Added by Tanzima
extern int io_access;
extern int disktime;
extern int updatetime;
extern int counttime;
extern int kmintime;
//......................
//------------------------------------------------------------
RTree::RTree(char *fname, int _b_length, Cache *c, int _dimension)
  //use this constructor to build a new tree
{
    file = new BlockFile(fname, _b_length);
    cache = c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
	  //note that when a tree is constructed, the root is automatically created
	  //though at this time there is no entry in the root yet.
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr -> block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
RTree::RTree(char *fname, Cache *c)
  //use this constructor to restore a tree from a file
{
    file = new BlockFile(fname, 0);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    file -> read_header(header);
    read_header(header);
	delete [] header;

    root_ptr = NULL;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
RTree::RTree(char *inpname, char *fname, int _b_length, Cache *c, int _dimension)
  // construct new R-tree from a specified input textfile with rectangles
{
    Entry *d;
    FILE *fp;
    file = new BlockFile(fname, _b_length);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);
	delete [] header;

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr->block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;

	int record_count = 0;

    if((fp = fopen(inpname,"r")) == NULL)
    {
      delete this;
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
      while (!feof(fp))
      {
		record_count ++;

		d = new Entry(dimension, NULL);

    	fscanf(fp, "%d", &(d -> son));
		//printf("ID=%d ",d->son);
		//for (int i = 0; i < 2 * dimension; i ++)
		//{
			//fscanf(fp, " %f", &(d -> bounces[i]));
		    fscanf(fp, " %f %f %f %f", &(d->bounces[0]),
		 	&(d->bounces[1]), &(d->bounces[2]), &(d->bounces[3]));
		//}
		//fscanf(fp, "\n");

		//if(record_count==20000 || record_count==20001)	printf(": %f %f %f %f\n", d->bounces[0], d->bounces[2], d->bounces[1], d->bounces[3]);
    	insert(d);

		  //d will be deleted in insert()

		if (record_count % 100 == 0)
		{
			for (int i = 0; i < 79; i ++)  //clear a line
				printf("\b");

			printf("inserting object %d", record_count);
		}
      }
    }
	//TANZIMA
	printf("inserting object %d", record_count);
	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
RTree::~RTree()
{
	char *header = new char[file -> get_blocklength()];
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL)
    {
        delete root_ptr;
        root_ptr = NULL;
    }

	if (cache)
      cache -> flush();

    delete file;

    delete re_data_cands;
	delete deletelist;

    //printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
	//   num_of_inodes, num_of_dnodes, num_of_data);
}
//------------------------------------------------------------
void RTree::del_root()
{
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
bool RTree::delete_entry(Entry *d)
{
	load_root();

	R_DELETE del_ret;
	del_ret=root_ptr->delete_entry(d);

	if (del_ret == NOTFOUND) return false;
	if (del_ret == ERASED) 
		error("RTree::delete_entry--The root has been deleted\n",true);
 
	if (root_ptr -> level > 0 && root_ptr -> num_entries == 1)
		//there is only one entry in the root but the root
		//is not leaf.  in this case, the child of the root is exhalted to root
	{
		root = root_ptr -> entries[0].son;
		delete root_ptr;
		root_ptr = NULL;
		load_root();
		num_of_inodes--;
	}

	//Now will reinsert the entries
	while (deletelist -> get_num() > 0)
	{
		Linkable *e;
		e = deletelist -> get_first();
		Entry *new_e = new Entry(dimension, NULL);
		new_e -> set_from_Linkable(e);
		deletelist -> erase();
		insert(new_e);
	}

	delete root_ptr;
	root_ptr = NULL;

	return true;
}
//------------------------------------------------------------
void RTree::insert(Entry* d)
{
    int i, j;
    RTNode *sn;
    RTNode *nroot_ptr;
    int nroot;
    Entry *de;
    R_OVERFLOW split_root;
    Entry *dc;
    float *nmbr;

    // load root into memory
    load_root();

    // no overflow occured until now
    re_level = new bool[root_ptr -> level + 1];
    for (i = 0; i <= root_ptr -> level; i++)
        re_level[i] = FALSE;

    // insert d into re_data_cands as the first entry to insert
    // make a copy of d because it should be erased later
    Linkable *new_link;
	new_link = d -> gen_Linkable();
	re_data_cands -> insert(new_link);

	delete d;  //we follow the convention that the entry will be deleted when insertion finishes

    j = -1;
    while (re_data_cands -> get_num() > 0)
    {
        // first try to insert data, then directory entries
	    Linkable *d_cand;
		d_cand = re_data_cands -> get_first();
        if (d_cand != NULL)
        {
            // since "erase" deletes the data itself from the
            // list, we should make a copy of the data before
            // erasing it
			dc = new Entry(dimension, NULL);
            dc -> set_from_Linkable(d_cand);
            re_data_cands -> erase();

            // start recursive insert with root
			split_root = root_ptr -> insert(dc, &sn);
        }
        else
	        error("RTree::insert: inconsistent list re_data_cands", TRUE);

    	if (split_root == SPLIT)
    	// insert has lead to split --> new root-page with two sons (i.e. root and sn)
    	{
    	    nroot_ptr = new RTNode(this);
    	    nroot_ptr -> level = root_ptr -> level + 1;
    	    num_of_inodes++;
    	    nroot = nroot_ptr -> block;

    	    de = new Entry(dimension, this);
    	    nmbr = root_ptr -> get_mbr();
    	    memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de->son = root_ptr->block;
    	    de->son_ptr = root_ptr;
    	    nroot_ptr -> enter(de);

    	    de = new Entry(dimension, this);
    	    nmbr = sn -> get_mbr();
    	    memcpy(de -> bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de -> son = sn -> block;
    	    de -> son_ptr = sn;
    	    nroot_ptr->enter(de);

    	    root = nroot;
            root_ptr = nroot_ptr;

            root_is_data = FALSE;
        }
        j++;
    }

    num_of_data++;

    delete [] re_level;

	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
void RTree::load_root()
{
    if (root_ptr == NULL)
        root_ptr = new RTNode(this, root);
}
//------------------------------------------------------------
//mbr has (x_min,x_max,y_min,y_max)
void RTree::rangeQuery(float *mbr, SortedLinList *res)
{
    load_root();
	//Added by Tanzima
	io_access++;
	//..............
    root_ptr -> rangeQuery(mbr,res);

	delete root_ptr;
	root_ptr = NULL;	
}
//------------------------------------------------------------
void RTree::read_header(char *buffer)
{
    int i;

    memcpy(&dimension, buffer, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
int RTree::update_rslt(Entry *_e, float _dist, Entry *_rslt, 
						float *_key, int _k)
{
	for (int i = 0; i < _k; i ++)
	{
		if (_dist < _key[i])
		{
			for (int j = _k - 1; j > i; j --)
			{
				_rslt[j] = _rslt[j - 1];
				_key[j] = _key[j - 1];
			}
			_rslt[i] = *_e;
			_key[i] = _dist;
			return i;
		}
	}
	error("Error in update_rslt\n", true);
	return -1;
}
//------------------------------------------------------------
void RTree::write_header(char *buffer)
{
    int i;

    memcpy(buffer, &dimension, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
//---added for valdity region queries-------------------------
// perform a window query mbr to get the query result into in_objs, put the outer objects retrieved into out_objs_so_far
void RTree::rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far)
{
    load_root();

    root_ptr->rect_win_query(mbr, in_objs, out_objs_so_far);

	delete root_ptr;
	root_ptr = NULL;
}

// perform a window query mbr (excluding the window excr) to get the query result into c_inf_objs
void RTree::rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs)
{
    load_root();

    root_ptr->rect_win_query(mbr, exclmbr, c_inf_objs);

	delete root_ptr;
	root_ptr = NULL;
}

void RTree::BFNN(float *_qpt, int _k, Entry *_rslt)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = (float) MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				if (edist < key[_k - 1])
					update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
			}
			else
			{
				if (edist<key[_k - 1]) 
					//meaning that edist is valid and we insert it to heap
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					heap -> insert(he);
					delete he;
				}
			}
		}
	
		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->key>key[_k - 1]) //the algorithm terminates
					son = -1;
				else if (he->level == 0) 
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
					 else
						son=he->son1;
			}
		}
		delete he;
		//--------------------------------------------------------
	}

	delete [] key;
	delete heap;
}



//RTree::BFNN --- Modified by Tanzima for Rectangle

float RTree::KMin(Point2D m, int c, float cl, int k, Pointlocation _rslt[], int *num_of_data, float d_safe_corner,DistfromPOI _cornertoPOI[])
{
	float *kmindist = new float[k];
	float tempdist;
	for(int i=0; i<k; i++)
	{
			kmindist[i]=(float) MAXREAL;
	}	
	
	for(int i=0; i<*num_of_data;i++)
	{
		if(cl* _cornertoPOI[i].d[c] <= d_safe_corner)
		{
			tempdist=Dist(_rslt[i],m);			

			for(int j=0; j<k; j++)
			{				

				if(kmindist[j]>tempdist)
				{
					for(int l=k-1;l>j;l--)
					{
						kmindist[l]= kmindist[l-1];
					}
					kmindist[j]= tempdist;
					break;
				}
			}
		}
		
	}
	return kmindist[k-1];
}


void RTree::UpdateCount(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, float d_safe_corner, int *count, DistfromPOI _cornertoPOI[])
{
	Point2D c[4];
	c[0][0] = R.x1; c[0][1]=R.y1;
	c[1][0] = R.x1; c[1][1]=R.y2;
	c[2][0] = R.x2; c[2][1]=R.y2;
	c[3][0] = R.x2; c[3][1]=R.y1;
	
		
	_cornertoPOI[*num_of_data-1].d[0]= Dist(_rslt[*num_of_data-1],c[0]);
	_cornertoPOI[*num_of_data-1].d[1]= Dist(_rslt[*num_of_data-1],c[1]);
	_cornertoPOI[*num_of_data-1].d[2]= Dist(_rslt[*num_of_data-1],c[2]);
	_cornertoPOI[*num_of_data-1].d[3]= Dist(_rslt[*num_of_data-1],c[3]);
	
	for(int j=0; j<4; j++)
	{
		if(count[j]<k)
		{
			count[j]=0;
			for(int i=0; i<*num_of_data;i++)
			{
				if(cl*_cornertoPOI[i].d[j] <= d_safe_corner)	count[j]++;
				
			}
		}
	}
	/*
	for(int i=0; i<_rslt.size();i++)
	{
		if(i<_cornertoPOI.size())
		{
			if(cl*_cornertoPOI[i].d[0] <= d_safe_corner)	count[0]++;
			if(cl*_cornertoPOI[i].d[1] <= d_safe_corner)	count[1]++;
			if(cl*_cornertoPOI[i].d[2] <= d_safe_corner)	count[2]++;
			if(cl*_cornertoPOI[i].d[3] <= d_safe_corner)	count[3]++;
		
		}
		else
		{
			Pointlocation p = _rslt[i];
			dist_c_P.d[0]= Dist(p,c[0]);
			if(cl*dist_c_P.d[0] <= d_safe_corner)	count[0]++;
			dist_c_P.d[1]= Dist(p,c[1]);
			if(cl*dist_c_P.d[1] <= d_safe_corner)	count[1]++;
			dist_c_P.d[2]= Dist(p,c[2]);
			if(cl*dist_c_P.d[2] <= d_safe_corner)	count[2]++;
			dist_c_P.d[3]= Dist(p,c[3]);
			if(cl*dist_c_P.d[3] <= d_safe_corner)	count[3]++;
			_cornertoPOI.push_back(dist_c_P);

		}
	}
	*/
}

int RTree::UpdateStatus(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, int *count, DistfromPOI _cornertoPOI[])
{
	Point2D o;
	o[0] = (R.x1+R.x2)/2;
	o[1] = (R.y1+R.y2)/2;
	float r = Dist(_rslt[*num_of_data-1],o);
	
	Point2D c[4];
	c[0][0] = R.x1; c[0][1]=R.y1;
	c[1][0] = R.x1; c[1][1]=R.y2;
	c[2][0] = R.x2; c[2][1]=R.y2;
	c[3][0] = R.x2; c[3][1]=R.y1;

	float d1=Dist(c[0],c[1]);
	float d2=Dist(c[1],c[2]);
	float d_safe_corner = r-0.5*sqrt(d1*d1+d2*d2);

	

	
	UpdateCount(R,cl,k,_rslt,num_of_data,d_safe_corner,count,_cornertoPOI);
	

	for(int i=0; i<4; i++)
	{
		if (count[i]<k)	return 0;
	}

	

	float d_i, d_j,d_max_m,d_max = 0;
	int i,j;
	Point2D m;
	for (i=0; i<4; i++)
	{
		j=(i+1)%4;
		
		m[0]= (c[i][0]+c[j][0])/2;
		m[1]= (c[i][1]+c[j][1])/2;
		d_i= KMin(m,i,cl,k,_rslt,num_of_data,d_safe_corner,_cornertoPOI);
		d_j= KMin(m,j,cl,k,_rslt,num_of_data,d_safe_corner,_cornertoPOI);
		if(d_i>d_j) d_max_m = d_i;
		else d_max_m = d_j;
		if(d_max_m > d_max)
			d_max = d_max_m;
		
	}
	float c_dmax=d1;
	if(d2>d1)	c_dmax=d2;
	float d_safe = r- 0.5 * c_dmax;	

	if(cl*d_max>d_safe)
		return ceil(r+cl*(d_max-d_safe));
	else
		return -1;
}
void RTree::Rect_kNNQ(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data)
{
	int io=0;
	//init status
	int count[4];
	for(int i=0; i<4; i++)	count[i]=0;
	int status=0;
	Point2D o;
	o[0] = (R.x1+R.x2)/2;
	o[1] = (R.y1+R.y2)/2;
	DistfromPOI cornertoPOI[10000];	

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	while (son != -1 && status != -1)
	{
		
		
		RTNode *rtn = new RTNode(this, son);
		//Experiment
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			float o1[2];
			o1[0]=(float)o[0];
			o1[1]=(float)o[1];
			float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;

		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{	
				if (he->level == 0) //p is an object 
				{
					if(*num_of_data==10000)
						//printf("\nGreater than 10000\n");
						error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
					
					_rslt[*num_of_data].x=he->x1;
					_rslt[*num_of_data].y=he->y1;
					*num_of_data=*num_of_data+1;
					//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);
						
					
					
					if(status==0)
					{												
						status = UpdateStatus(R,cl,k,_rslt,num_of_data,count,cornertoPOI);
												
						if(status != -1)	again=true;						
					}
					else if (status < Dist(_rslt[*num_of_data-1],o))
					{							
						status = -1;
					}
					else
					{
						again = true;
					}
				}
				else
				{
					if(status>0)
					{
						if(status<he->key)
							status = -1;
						else
							son=he->son1;
					}
					else
					{
						son=he->son1;
					}
				}
			}
		}
		
		delete he;
	}
	delete heap;
}

//Point Nearest Neighbor query
void RTree::Point_BFN_NNQ(Point2D o, double *_rslt)
{
	
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float o1[2];
			o1[0]=(float)o[0];
			o1[1]=(float)o[1];
			float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);
			printf("%f %f %f %f edist %f\n",rtn->entries[i].bounces[0],rtn->entries[i].bounces[1],rtn->entries[i].bounces[2],rtn->entries[i].bounces[3],edist);

			HeapEntry *he = new HeapEntry();
			he -> key = edist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;

		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->level == 0) //p is an object 
				{
					_rslt[0] = he->x1;
					_rslt[1] = he->y1;
					printf("Point %f,%f Mindist %f\n",he->x1,he->y1,he->key);
					son=-1;
														
				}
				else
				{
					son=he->son1;
					
				}
			}
		}
		
		delete he;
	}
	delete heap;
}
//END

//Added by Nusrat
//Point Group Nearest Neighbor query for numOfQueryPoints
void RTree::Point_BFN_GNNQ(Point2D o[], double *_rslt,int numOfQueryPoints)
{

	for(int j=0;j<numOfQueryPoints;j++){
		printf("%f,%f\n",o[j][0],o[j][1]);
	}
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float o1[2];
			float gnnMinDist=0;
			//Find min Group distance

			for(int j=0;j<numOfQueryPoints;j++){
				o1[0] = (float) o[j][0];
				o1[1] = (float) o[j][1];

				float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);
				gnnMinDist+=edist;
				//printf("%f,%f has mindist %f\n",o1[0],o1[1],edist);
			}

			HeapEntry *he = new HeapEntry();
			he -> key = gnnMinDist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->level == 0) //p is an object
				{
					_rslt[0] = he->x1;
					_rslt[1] = he->y1;
					printf("Point %f,%f Mindist %f\n",he->x1,he->y1,he->key);
					son=-1;

				}
				else
				{
					son=he->son1;

				}
			}
		}

		delete he;
	}
	delete heap;
}
//END

//Added by Nusrat
//Point k Group Nearest Neighbor query for numOfQueryPoints
void RTree::Point_BFN_kGNNQ(Point2D o[], int k,double _rslt[][2],int numOfQueryPoints)
{

	int indexOfGNNRetrieved=0;
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float o1[2];
			float gnnMinDist=0;
			//Find min Group distance

			for(int j=0;j<numOfQueryPoints;j++){
				o1[0] = (float) o[j][0];
				o1[1] = (float) o[j][1];

				float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);
				gnnMinDist+=edist;
				//printf("%f,%f has mindist %f\n",o1[0],o1[1],edist);
			}

			printf("Entry %d : mindist %f\n",i,gnnMinDist);

			HeapEntry *he = new HeapEntry();
			he -> key = gnnMinDist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			if(indexOfGNNRetrieved==k-1)
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->level == 0) //p is an object
				{
					_rslt[indexOfGNNRetrieved][0] = he->x1;
					_rslt[indexOfGNNRetrieved][1] = he->y1;
					printf("Point %f,%f Mindist %f\n",he->x1,he->y1,he->key);
					son=-1;
					indexOfGNNRetrieved++;

				}
				else
				{
					son=he->son1;

				}
			}
		}

		delete he;
	}
	delete heap;
}
//END


// The following code was copied from the implementation of TP-kNN queries from Tony
void RTree::TPNN_TP(float *_qline, int _k, Entry *_nn, Entry *_rslt, float _max_trvl)
{
	float key = _max_trvl; 
	  // the minimum distance that the query point must travel 

//we comment these lines to avoid initing the heap everytime--
//this function is called
//Heap *heap = new Heap();
//heap->init(dimension);
//------------------------------------------------------------
	if (tpheap==NULL)
		error("tpheap is not initialized\n", true);
	tpheap->used=0;
	
	int son = root;
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			//first create an array m2 for e to cal dist----------
			float *m2 = new float [4 * dimension];
			memset(m2, 0, 4 * dimension * sizeof(float));
			memcpy(m2, rtn -> entries[i].bounces, 2 * dimension * sizeof(float));
			//----------------------------------------------------
//testing--------------------------
//if (rtn->entries[i].son==573673 && cnt==84-1)
//	printf("testing...\n");
//---------------------------------

			float edist = (float) MAXREAL;
			if (rtn -> level == 0)
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^leaf node case^^^^^^^^^^^^^^^^^^^^^
			{	
				int eNNsub=-1;
				//if e (i.e., m2) is one of the current NN points-
				//eNNsub stores its subsript in _nn; otherwise, 
				//eNNsub=-1
				for (int j=0; j<_k; j++)
					if (rtn->entries[i].son==_nn[j].son)
					{ eNNsub=j; j=_k;}
				//------------------------------------------------
		
				if (eNNsub==-1 || SEQUENCE_SENSITIVE)
					//namely, if sequence insensitive and e is indeed a NN
					//then we do not handle this entry
				{

					float *m1 = new float [4 * dimension];
					//find the NN that leads to the minimum-------
					//influence time
					int nn_num=-1; //the subsript of the NN to be found
					for (int j = 0; j < _k; j ++)
						//for each of the NN found in the 1st step
					{
						bool yesdo=true; //whether to compute
						  //the inflence time of nn[j] and e
						if (j==eNNsub) yesdo=false;
						
						//check if this pair has produced --------
						//influence time before
						for (int l=0; l<PAST_PAIR; l++)
							if (min(_nn[j].son, rtn->entries[i].son)==min(last_pair[l][0], last_pair[l][1]) &&
								max(_nn[j].son, rtn->entries[i].son)==max(last_pair[l][0], last_pair[l][1]))
							{	yesdo=false; l=PAST_PAIR; }
						//----------------------------------------

						if (yesdo)
						{
//these codes use NNinf===========================================
/*
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							//get the influence time of m2-------- 
							//(for entry e) with respect to(nn[j])
							float this_inf = NNinf(m1, m2, _qline, dimension); 
							//------------------------------------
*/
//================================================================

//these codes use NNinf2==========================================
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							m1[1]=m1[2]; m2[1]=m2[2];
							//create an arry m3 for _qline----------------
							float *m3=new float[2*dimension];
							m3[0]=_qline[0]; m3[1]=_qline[2];
							m3[2]=_qline[4]; m3[3]=_qline[6];
							//--------------------------------------------
							//get the influence time of m2-------- 
							//(for entry e) with respect to(nn[j])
							float this_inf = NNinf2(m1, m2, m3); 
							//------------------------------------
							delete []m3;
//================================================================

							if (this_inf>0 && this_inf<edist)
								//this_inf=0 means that there is another point that has the same distance
								//to the current query position as the current NN. In this implementation,
								//we choose to ignore handling such special cases, which, however, may cause
								//problems for datasets with big cardinality
//							if (this_inf>=0 && this_inf<edist)
							{
								edist=this_inf; nn_num=j;
							}
						}  //END if (yesdo)
					}//END checking all neighbors
					//-------------------------------------------------
					//if (edist<key && edist!=0)
					if (edist<key)
					{
						update_rslt(&(rtn->entries[i]), edist, _rslt, &key, 1);
						_rslt->nn_num=nn_num;
					}
					delete []m1;
				}
			}
//^^^^^^^^^^^^^^^^^^^^^^^^^non-leaf node case^^^^^^^^^^^^^^^^^^^^^
			else
				//Next handle non-leaf node case
			{
				float *m1 = new float [4 * dimension];
				for (int j = 0; j < _k; j ++)
				{
					//first create an array m1 to cal dist--------
					memset(m1, 0, 4 * dimension * sizeof(float));
					memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
					//--------------------------------------------
					float this_mininf = NNmininf(m1, m2, _qline, dimension);
					if (this_mininf < edist)
						edist = this_mininf;
				}
				delete [] m1;

				if (edist < key)
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					tpheap -> insert(he);
					delete he;
				}
			}
			delete [] m2;

		}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	
		delete rtn;

		//get next entry from the heap
		bool again = true;
		while (again)
		{
			again = false;
			HeapEntry *he = new HeapEntry();
			if (!tpheap -> remove(he))  //heap is empty
				son = -1;
			else
				if (he -> key > key)
					//the algorithm can terminate
					son = -1;
				else
					son = he -> son1;
			delete he;
		}
	}
//delete heap;
}



//RTree::private_kGNN --- Created by Tanzima for kGNN

/*
void RTree::private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{	
	printf("private_kGNN_sum -----2nd version\n");
	//printf("k =%d g_size=%d\n",k,g_size);
	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;
	
	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;
	
	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end==0)
	{		
		//printf("outer while loop executing\n");

		RTNode *rtn = new RTNode(this, son);
		
		for (i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			
			edist1 = sqrt(MINRECTDIST(mbr, rtn->entries[i].bounces));	
			if((g_size*edist1)>maxdist[k-1])	continue;

			edist1 = sqrt(MINRECTDIST(R[0], rtn->entries[i].bounces));
			for (j=1; j < g_size; j++)
			{				
				edist1 += sqrt(MINRECTDIST(R[j], rtn->entries[i].bounces));
			}

			if(edist1>maxdist[k-1])	continue;

			edist2 = sqrt(MAXRECTDIST(R[0], rtn->entries[i].bounces));
			for (j=1; j < g_size; j++)
			{	
				edist2 += sqrt(MAXRECTDIST(R[j], rtn->entries[i].bounces));
			}

			//update maxdistk		
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int u=k-1; u>l; u--)
					{
						maxdist[u]=maxdist[u-1];						
					}
					maxdist[l]=edist2;
					break;
				}
			}
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{	
				if (he->level == 0) //p is an object 
				{	
					printf("dequeing--- level = %d\n",he->level);
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;


						*num_of_data=*num_of_data+1;
									
						//get next data  from heap
						again =true;						
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}
		
		delete he;

		//printf("outer while loop terminating\n");
	}
	delete heap;

	printf("2nd version terminated\n");
}
*/

//END


void RTree::private_kGNN_max(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{

	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;
	
	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;
	
	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next

	
	while (son != -1 && end==0)
	{
		RTNode *rtn = new RTNode(this, son);
	
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0, tdist=0;
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);	
			if((edist1)>maxdist[k-1])	continue;

			
			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				tdist = MINRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist1)
					edist1=tdist;			
			}			

			if(edist1>maxdist[k-1])	continue;
			
			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				tdist= MAXRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist2)
					edist2=tdist;			
			}

			

			//update maxdistk		
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];						
					}
					maxdist[l]=edist2;
					break;
				}
			}
		
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;

		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{	

				if (he->level == 0) //p is an object 
				{
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);
						
						//get next data  from heap
						again =true;						
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}
		
		delete he;
	}
	delete heap;
	
}

//added by Tahrima

/*
void RTree::private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data,RTree *drt)
{	
	//printf("inside private_kGNN_sum");
	//RTNode *rtn = new RTNode(this, root);
	//printf("private_kGNN_sum is called");
	drt->private_kGNN_sum(R, g_size, k,rslt,num_of_data);
	printf("Returning from private_kGNN_sum---%d\n",*num_of_data);

}

void RTree::private_GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size,Pointlocation rslt[], int *num_of_data_P1,int k,RTree *drt,Pair A[])
{	
	printf("hereeeeeeeeeee\n");
	//char s;
	//scanf("%c",&s);
	//variable initialization
	int end=0;
	int i, j;
	Pointlocation rsltp2[kMAX];
	int num_of_data_p2 = 0;
	double MinDist = MAXDOUBLE;
	double dist = MAXDOUBLE;
	double currentdist = 0.0;
	Point2D p1,p2;
	
	for(int pair = 0; pair < k;pair++)	
		A[pair].mindist = MAXDOUBLE;

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end == 0)
	{		
		printf("Outer1 loop started\n");

		RTNode *rtn = new RTNode(this, son);
		
		for (i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			
			edist1 = sqrt(MINRECTDIST(S[0], rtn->entries[i].bounces));
			for (j=1; j < g_size; j++)
			{				
				edist1 += sqrt(MINRECTDIST(S[j], rtn->entries[i].bounces));
			}
			edist2 = sqrt(MAXRECTDIST(S[0], rtn->entries[i].bounces));
			for (j=1; j < g_size; j++)
			{	
				edist2 += sqrt(MAXRECTDIST(S[j], rtn->entries[i].bounces));
			}
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			
			//printf("after heap\n");
		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			printf("Outer2 loop started\n");

			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{	
				if (he->level == 0) //p is an object 
				{
					printf("object found\n");
						//enter into result set
						if(*num_of_data_P1==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data_P1].x=he->x1;
						rslt[*num_of_data_P1].y=he->y1;
						rslt[*num_of_data_P1].dmin=he->key;
						rslt[*num_of_data_P1].dmax=he->key1;

						*num_of_data_P1=*num_of_data_P1+1;

						//added p1' to Destination Array
						Rectangle1 *DP = new Rectangle1[g_size+1];
		
						for(i = 0; i < g_size; i++)
						{
							DP[i].x1 = D[i].x1;
							DP[i].y1 = D[i].y1;
							DP[i].x2 = D[i].x2;
							DP[i].y2 = D[i].y2;
						}
						DP[g_size].x1 = rslt[*num_of_data_P1-1].x;
						DP[g_size].y1 = rslt[*num_of_data_P1-1].y;
						DP[g_size].x2 = rslt[*num_of_data_P1-1].x;
						DP[g_size].y2 = rslt[*num_of_data_P1-1].y;
						
						//printf("copy from D to DP\n");

						num_of_data_p2 = 0;						
						
						private_kGNN_sum(DP, g_size+1, k, rsltp2, &num_of_data_p2,drt);//k is always one here
						

						//printf("return from private_kGNN_sum --- %d\n",num_of_data_p2);
						for(int num_p2 = 0 ; num_p2 < num_of_data_p2;num_p2++)
						{
							//printf("calculate the total distance\n");
							//calculate the total distance
							Point2D p;
							p[0] = S[0].x1,p[1] = S[0].y1;
							currentdist = 0.0;

							currentdist += Dist(rslt[*num_of_data_P1-1],p);
							for(i = 1; i < g_size; i++)
							{
								p[0] = S[i].x1;
								p[1] = S[i].y1;
								currentdist += Dist(rslt[*num_of_data_P1-1],p);
							}
							dist = currentdist;
							currentdist += (g_size*Dist(rslt[*num_of_data_P1-1],rsltp2[num_p2]));

							p[0] = D[0].x1,p[1] = D[0].y1;
							currentdist += Dist(rsltp2[num_p2],p);

							for(i = 1; i < g_size; i++)
							{
								p[0] = D[i].x1;
								p[1] = D[i].y1;
								currentdist += Dist(rsltp2[num_p2],p);
							}

							for(int l=0; l<k; l++)
							{
								if(currentdist <= A[l].mindist)
								{
									for(int n=k-1; n>l; n--)
									{
										A[n].mindist=A[n-1].mindist;						
										A[n].p1x1= A[n-1].p1x1;
										A[n].p1y1= A[n-1].p1y1;
										A[n].p2x2= A[n-1].p2x2;
										A[n].p2y2= A[n-1].p2y2;
									}
									A[l].mindist = currentdist;

									A[l].p1x1= rslt[*num_of_data_P1-1].x;
									A[l].p1y1= rslt[*num_of_data_P1-1].y;
									A[l].p2x2= rsltp2[num_p2].x;
									A[l].p2y2= rsltp2[num_p2].y;

									break;
								}
							}

							//printf("calculate the total distance computed\n");
						}
						

						if(dist > A[k-1].mindist)
						{
							end = 1;
							printf("dist is greater than A[k-1].mindist\n");
						}
						else
						{
							//get next data  from heap
							again =true;
						}
														
				}
				else //not leaf node
				{
					son=he->son1;
				}	
			}
			
		}
		
		delete he;
	}
	delete heap;
	
}
*/


/*****************************************************************/

//RTree::private_kGNN --- Created by Tanzima for kGNN


void RTree::private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{	
	printf("private_kGNN_sum -----2nd version\n");
	printf("k =%d g_size=%d\n",k,g_size);
	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;
	
	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;
	
	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end==0)
	{		
		//printf("outer while loop executing\n");

		RTNode *rtn = new RTNode(this, son);
		
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);	
			if((g_size*edist1)>maxdist[k-1])	continue;

			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{				
				edist1 += MINRECTDIST(R[j], rtn->entries[i].bounces);
			}

			if(edist1>maxdist[k-1])	continue;

			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{	
				edist2 += MAXRECTDIST(R[j], rtn->entries[i].bounces);
			}

			//update maxdistk		
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];						
					}
					maxdist[l]=edist2;
					break;
				}
			}
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{	
				if (he->level == 0) //p is an object 
				{	
					//printf("dequeing--- level = %d\n",he->level);
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;


						*num_of_data=*num_of_data+1;
									
						//get next data  from heap
						again =true;						
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}
		
		delete he;

	//	printf("outer while loop terminating\n");
	}
	delete heap;

	printf("2nd version terminated\n");
}


void RTree::private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data,RTree *drt)
{	
	drt->private_kGNN_sum(R, g_size, k,rslt,num_of_data);
	printf("Returning from private_kGNN_sum");
	
}

void RTree::private_GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size,Pointlocation rslt[], int *num_of_data_P1,int k,RTree *drt,Pair A[])
{	
	printf("private_GNN_sum is called");
	//char s;
	//scanf("%c",&s);
	//variable initialization
	int end=0;
	int i, j;
	Pointlocation rsltp2[kMAX];
	int num_of_data_p2 = 0;
	double MinDist = MAXDOUBLE;
	double dist = MAXDOUBLE;
	double currentdist = 0.0;
	Point2D p1,p2;

	for(int pair = 0; pair < k;pair++)	
		A[pair].mindist = MAXDOUBLE;

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end == 0)
	{		
		printf("outer while loop executing\n");

		RTNode *rtn = new RTNode(this, son);
		
		for (i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			
			edist1 = MINRECTDIST(S[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{				
				edist1 += MINRECTDIST(S[j], rtn->entries[i].bounces);
			}
			edist2 = MAXRECTDIST(S[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{	
				edist2 += MAXRECTDIST(S[j], rtn->entries[i].bounces);
			}

			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			printf("inner while loop executing\n");

			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{	
				if (he->level == 0) //p is an object 
				{
					//	printf("After Dequeuing --- level = %d",he->level);
						//enter into result set
						if(*num_of_data_P1==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data_P1].x=he->x1;
						rslt[*num_of_data_P1].y=he->y1;
						rslt[*num_of_data_P1].dmin=he->key;
						rslt[*num_of_data_P1].dmax=he->key1;

						*num_of_data_P1=*num_of_data_P1+1;

						//added p1' to Destination Array
						Rectangle1 DP[kMAX+1];
		
						for(i = 0; i < g_size; i++)
						{
							DP[i].x1 = D[i].x1;
							DP[i].y1 = D[i].y1;
							DP[i].x2 = D[i].x2;
							DP[i].y2 = D[i].y2;
						}
						DP[g_size].x1 = rslt[*num_of_data_P1-1].x;
						DP[g_size].y1 = rslt[*num_of_data_P1-1].y;
						DP[g_size].x2 = rslt[*num_of_data_P1-1].x;
						DP[g_size].y2 = rslt[*num_of_data_P1-1].y;
						
						num_of_data_p2 = 0;
						
						printf("Calling second function");
						private_kGNN_sum(DP, g_size+1, k, rsltp2, &num_of_data_p2,drt);//k is always one here
						printf("Returning from second function");
						for(int num_p2 = 0 ; num_p2 < num_of_data_p2;num_p2++)
						{
							printf("calculate the total distance");
							//calculate the total distance
							Point2D p;
							p[0] = S[0].x1,p[1] = S[0].y1;
							currentdist = 0.0;

							currentdist += Dist(rslt[*num_of_data_P1-1],p);
							for(i = 1; i < g_size; i++)
							{
								p[0] = S[i].x1;
								p[1] = S[i].y1;
								currentdist += Dist(rslt[*num_of_data_P1-1],p);
							}
							dist = currentdist;
							currentdist += (g_size*Dist(rslt[*num_of_data_P1-1],rsltp2[num_p2]));

							p[0] = D[0].x1,p[1] = D[0].y1;
							currentdist += Dist(rsltp2[num_p2],p);

							for(i = 1; i < g_size; i++)
							{
								p[0] = D[i].x1;
								p[1] = D[i].y1;
								currentdist += Dist(rsltp2[num_p2],p);
							}

							for(int l=0; l<k; l++)
							{
								if(currentdist <= A[k-1].mindist)
								{
									for(int n=k-1; n>l; n--)
									{
										A[n].mindist=A[n-1].mindist;						
										A[n].p1x1= A[n-1].p1x1;
										A[n].p1y1= A[n-1].p1y1;
										A[n].p2x2= A[n-1].p2x2;
										A[n].p2y2= A[n-1].p2y2;
									}
									A[l].mindist = currentdist;

									A[l].p1x1= rslt[*num_of_data_P1-1].x;
									A[l].p1y1= rslt[*num_of_data_P1-1].y;
									A[l].p2x2= rsltp2[num_p2].x;
									A[l].p2y2= rsltp2[num_p2].y;

									break;
								}
							}

							printf("calculate the total distance computed");
						}
						

						if(dist > A[k-1].mindist)
						{
							end = 1;
							for(int u=0;u<k;u++){
							printf("\nFinal result---------kth\n");
							printf("\nnumber of p1=%d,u=%d\n",*num_of_data_P1,u);
							printf("p1:%lf %lf p2:%lf %lf mindist:%lf\n\n",A[u].p1x1,A[u].p1y1,A[u].p2x2,A[u].p2y2,A[u].mindist);

							char s;
							scanf("%c",&s);
							}
						}
						else
						{
							//get next data  from heap
							again =true;
						}
														
				}
				else //not leaf node
				{
					son=he->son1;
				}

				
			}
			printf("inner while loop terminated\n");
		}
		
		delete he;
		
	
		printf("outer while loop terminated\n");
	}
	delete heap;
	
}

/****************************************************************/

//another version of group nearest neigbour

void RTree::private_GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size,Pointlocation rslt[], int *num_of_data_P1,Pointlocation rsltp2[], int *num_of_data_P2,int k,RTree *drt1,Pair A[])
{	
	
	//char s;
	//scanf("%c",&s);
	//variable initialization
	int end=0;
	int i, j;
	Pointlocation rsltp3[kMAX];
	int num_of_data_p3 = 0;
	double MinDist = MAXDOUBLE;
	double dist = MAXDOUBLE;
	double currentdist = 0.0;
	Point2D p1,p2;
	
	Rectangle1 *R1 = new Rectangle1[1];
	R1[0].x1 = rslt[*num_of_data_P1-1].x;
	R1[0].y1 = rslt[*num_of_data_P1-1].y;
	R1[0].x2 = rslt[*num_of_data_P1-1].x;
	R1[0].y2 = rslt[*num_of_data_P1-1].y;


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end == 0)
	{		
		RTNode *rtn = new RTNode(this, son);
		
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			
			edist1 = sqrt(MINRECTDIST(R1[0], rtn->entries[i].bounces));			
			edist2 = sqrt(MAXRECTDIST(R1[0], rtn->entries[i].bounces));
	
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{	
				if (he->level == 0) //p is an object 
				{
						
						//enter into result set
						if(*num_of_data_P2==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rsltp2[*num_of_data_P2].x=he->x1;
						rsltp2[*num_of_data_P2].y=he->y1;
						rsltp2[*num_of_data_P2].dmin=he->key;
						rsltp2[*num_of_data_P2].dmax=he->key1;

						*num_of_data_P2=*num_of_data_P2+1;

						//added p1' to Destination Array
						Rectangle1 *DP = new Rectangle1[g_size+1];
		
						for(i = 0; i < g_size; i++)
						{
							DP[i].x1 = D[i].x1;
							DP[i].y1 = D[i].y1;
							DP[i].x2 = D[i].x2;
							DP[i].y2 = D[i].y2;
						}
						DP[g_size].x1 = rsltp2[*num_of_data_P2-1].x;
						DP[g_size].y1 = rsltp2[*num_of_data_P2-1].y;
						DP[g_size].x2 = rsltp2[*num_of_data_P2-1].x;
						DP[g_size].y2 = rsltp2[*num_of_data_P2-1].y;
						
						num_of_data_p3 = 0;						
						
						private_kGNN_sum(DP, g_size+1, k, rsltp3, &num_of_data_p3,drt1);
						
						for(int num_p3 = 0 ; num_p3 < num_of_data_p3;num_p3++)
						{
							//printf("calculate the total distance");
							//calculate the total distance
							Point2D p;
							p[0] = S[0].x1,p[1] = S[0].y1;
							currentdist = 0.0;

							currentdist += Dist(rslt[*num_of_data_P1-1],p);
							for(i = 1; i < g_size; i++)
							{
								p[0] = S[i].x1;
								p[1] = S[i].y1;
								currentdist += Dist(rslt[*num_of_data_P1-1],p);
							}
							
							currentdist += (g_size*Dist(rslt[*num_of_data_P1-1],rsltp2[*num_of_data_P2-1]));
							dist = currentdist;
							currentdist += (g_size*Dist(rsltp2[*num_of_data_P2-1],rsltp3[num_p3]));

							p[0] = D[0].x1,p[1] = D[0].y1;
							currentdist += Dist(rsltp3[num_p3],p);

							for(i = 1; i < g_size; i++)
							{
								p[0] = D[i].x1;
								p[1] = D[i].y1;
								currentdist += Dist(rsltp3[num_p3],p);
							}

							for(int l=0; l<k; l++)
							{
								if(currentdist <= A[l].mindist)
								{
									for(int n=k-1; n>l; n--)
									{
										A[n].mindist=A[n-1].mindist;						
										A[n].p1x1= A[n-1].p1x1;
										A[n].p1y1= A[n-1].p1y1;
										A[n].p2x2= A[n-1].p2x2;
										A[n].p2y2= A[n-1].p2y2;
										A[n].p3x3= A[n-1].p3x3;
										A[n].p3y3= A[n-1].p3y3;
									}
									A[l].mindist = currentdist;

									A[l].p1x1= rslt[*num_of_data_P1-1].x;
									A[l].p1y1= rslt[*num_of_data_P1-1].y;
									A[l].p2x2= rsltp2[*num_of_data_P2-1].x;
									A[l].p2y2= rsltp2[*num_of_data_P2-1].y;
									A[l].p3x3= rsltp3[num_p3].x;
									A[l].p3y3= rsltp3[num_p3].y;

									break;
								}
							}

							//printf("calculate the total distance computed");
						}
						
						if(dist > A[k-1].mindist)
						{
							end = 1;	
						}
						else
						{
							//get next data  from heap
							again =true;
						}
														
				}
				else //not leaf node
				{
					son=he->son1;
				}				
			}
			
		}
		
		delete he;		
	}
	delete heap;
	
}

void RTree::private_3GNN_sum(Rectangle1 S[],Rectangle1 D[], int g_size,Pointlocation rslt[], int *num_of_data_P1,int k,RTree *drt,RTree *drt1,Pair A[])
{	
	
	//char s;
	//scanf("%c",&s);
	//variable initialization
	int end=0;
	int i, j;
	Pointlocation rsltp2[kMAX];
	int num_of_data_p2 = 0;
	double MinDist = MAXDOUBLE;
	double dist = MAXDOUBLE;
	double currentdist = 0.0;
	Point2D p1,p2;
	
	for(int pair = 0; pair < k;pair++)	
		A[pair].mindist = MAXDOUBLE;

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end == 0)
	{		
		RTNode *rtn = new RTNode(this, son);
		
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			
			edist1 = sqrt(MINRECTDIST(S[0], rtn->entries[i].bounces));
			for (j=1; j < g_size; j++)
			{				
				edist1 += sqrt(MINRECTDIST(S[j], rtn->entries[i].bounces));
			}
			edist2 = sqrt(MAXRECTDIST(S[0], rtn->entries[i].bounces));
			for (j=1; j < g_size; j++)
			{	
				edist2 += sqrt(MAXRECTDIST(S[j], rtn->entries[i].bounces));
			}

			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{	
				if (he->level == 0) //p is an object 
				{						
						//enter into result set
						if(*num_of_data_P1==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data_P1].x=he->x1;
						rslt[*num_of_data_P1].y=he->y1;
						rslt[*num_of_data_P1].dmin=he->key;
						rslt[*num_of_data_P1].dmax=he->key1;

						*num_of_data_P1=*num_of_data_P1+1;
						num_of_data_p2 = 0;		

						drt->private_GNN_sum(S,D,g_size,rslt,num_of_data_P1,rsltp2, &num_of_data_p2,k,drt1,A);
						
						Point2D p;
						p[0] = S[0].x1,p[1] = S[0].y1;
						currentdist = 0.0;

						currentdist += Dist(rslt[*num_of_data_P1-1],p);
						for(i = 1; i < g_size; i++)
						{
							p[0] = S[i].x1;
							p[1] = S[i].y1;
							currentdist += Dist(rslt[*num_of_data_P1-1],p);
						}
						dist = currentdist;

						//printf("calculate the total distance computed");
				
						if(dist > A[k-1].mindist)
						{
							end = 1;
						}
						else
						{
							//get next data  from heap
							again =true;
						}														
				}
				else //not leaf node
				{
					son=he->son1;
				}				
			}			
		}		
		delete he;	
	}
	delete heap;
	
}


//END

RTNode* RTree::getRTNode(int ptr)
{
	//printf("inside getRTNode--%d\n",ptr);
	if(ptr == -1)
	{
		ptr = root;
		//printf("Starting from root");
		//char *s;
		//scanf("%c",&s);
	}
	RTNode* rtnode = new RTNode(this,ptr);
	//printf("terminate getRTNode\n");
	return rtnode;
}
void RTree::twoST_GTP_FA(Rectangle1 R[],Rectangle1 D[], int g_size, Pointlocation rsltp1p2[], int *num_of_data_p1p2,int k,RTree *dtr)
{	
	printf("called twoST_GTP_FA\n");
	//variable initialization
	//Variables
	int end=0,ctr_k = 0,i,j,ctr_i,ctr_j,level_r1,level_r2,son=0,dson=0;
	
	double distance = approximate_mindist(R,D,g_size,k);
	
	double ThresholdDist[2*kMAX];
	for(int h = 0; h < k; h++)	
		ThresholdDist[h] = distance;
		//ThresholdDist[h] = MAXDOUBLE;


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son_r1 = root; //this entry is to be visited next
	int son_r2 = -1;
	int empty = 0;

	while (empty != 1 && end==0)
	{	
		//printf("outer while loop1 started\n");

		RTNode *rtn1;
		RTNode *rtn2;
		float r1[4],r2[4];

		if(son != -1){
		//printf("son-r1--%d",son_r1);
		rtn1 = new RTNode(this, son_r1);
		
		ctr_i = rtn1 -> num_entries;
		level_r1 = rtn1->level;
		//printf("-----------------------------------------node of rtn1--level= %d \n",level_r1);
		}else{
			printf("*****************************************heap leaf node of rtn1\n");
			ctr_i = 1;
			level_r1 = 0;
		}
		
		if(dson != -1){		
		rtn2 = dtr->getRTNode(son_r2);
		ctr_j = rtn2 -> num_entries;
		level_r2 = rtn2->level;
		//printf("-----------------------------------------node of rtn2--level= %d\n",level_r2);
		}else{
			printf("*****************************************heap leaf node of rtn2\n");
			ctr_j = 1;
			level_r2 = 0;
		}

		for (i = 0; i < ctr_i; i++)
		{	
			//printf("outer while loop2 started\n");		
			/*char s;
			scanf("%c",&s);*/

			if(son!= -1){
			//copy(r1,rtn1->entries[i].bounces);
				r1[0] = rtn1->entries[i].bounces[0];
				r1[1] = rtn1->entries[i].bounces[1];
				r1[2] = rtn1->entries[i].bounces[2];
				r1[3] = rtn1->entries[i].bounces[3];

				son_r1 = rtn1->entries[i].son;

				//printf("r1 newly initialised.....%d\n",son_r1 );
			}
			
			for(j = 0; j < ctr_j;j++)
			{
				//printf("outer while loop3 started\n");

				/*char s;
				scanf("%c",&s);*/

				if(dson!= -1){
				r2[0] = rtn2->entries[j].bounces[0];
				r2[1] = rtn2->entries[j].bounces[1];
				r2[2] = rtn2->entries[j].bounces[2];
				r2[3] = rtn2->entries[j].bounces[3];

				son_r2 = rtn2->entries[j].son;

				
				}
				
				double edist1=0, edist2=0;

				//dmin(r1,r2) computation

				edist1 = sqrt(MINRECTDIST(R[0], r1));
				for (int m=1; m < g_size; m++)
				{				
					edist1 += sqrt(MINRECTDIST(R[m], r1));
				}

				edist1 += (g_size*sqrt(MINRECTDIST(r1,r2)));

				edist1 += sqrt(MINRECTDIST(D[0], r2));
				for (int m=1; m < g_size; m++)
				{				
					edist1 += sqrt(MINRECTDIST(D[m], r2));
				}

				if(edist1 > ThresholdDist[k-1])	continue;//pruning

				//dmax(r1,r2) computation

				edist2 = sqrt(MAXRECTDIST(R[0], r1));
				for (int m=1; m < g_size; m++)
				{	
					edist2 += sqrt(MAXRECTDIST(R[m], r1));
				}

				edist2 += (g_size*sqrt(MAXRECTDIST(r1,r2)));

				edist2 += sqrt(MAXRECTDIST(D[0], r2));
				for (int m=1; m < g_size; m++)
				{				
					edist2 += sqrt(MAXRECTDIST(D[m], r2));
				}


				//update ThresholdDist
				for(int l=0; l<k; l++)
				{
					if (edist2 < ThresholdDist[l])
					{
						for(int j=k-1; j>l; j--)
						{
							ThresholdDist[j]=ThresholdDist[j-1];						
						}
						ThresholdDist[l]=edist2;
						break;
					}
				}
				
				//Enqueue

				HeapEntry *he = new HeapEntry();
				he -> key = edist1;
				he -> key1 = edist2;
				he -> level = level_r1;	//for R1
				he -> son1 = son_r1;
				he -> dlevel = level_r2;	//for R2
				he -> dson1 = son_r2;
				
				he-> x1 = r1[0]; //for R1
				he-> x2 = r1[1];
				he-> y1 = r1[2];
				he-> y2 = r1[3];
				he-> dx1 = r2[0]; // for R2
				he-> dx2 = r2[1];
				he-> dy1 = r2[2];
				he-> dy2 = r2[3];

				//if(!he -> level || !he -> dlevel)
				//{
					//char *s;
					printf("%d ... %d \n",he -> level,he -> dlevel);
					//scanf("%c",&s);
					//printf("\n");
					//scanf("%c",&s);
					//scanf("%c",&s);
				//}

				heap -> insert(he);
				delete he;			


				//printf("outer while loop3 terminated\n");
			}

			//printf("outer while loop2 terminated\n");
		}
		
		if(son != -1)
			delete rtn1;
		if(dson != -1)
			delete rtn2;
		
		son = 0;
		dson = 0;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{	
			//printf("outer while loop4 started\n");
			/*char s;
			scanf("%c",&s);*/

			again = false;
			if (!heap->remove(he))  //heap is empty
				empty = 1;
			else
			{	
				if (he->level == 0 && he->dlevel == 0) //p is an object 
				{
						//printf("Dequeing-------\n");
						/*char s;
						scanf("%c",&s);*/

						//enter into result set
						if(*num_of_data_p1p2==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rsltp1p2[*num_of_data_p1p2].x=he->x1;
						rsltp1p2[*num_of_data_p1p2].y=he->y1;
						rsltp1p2[*num_of_data_p1p2].dx=he->dx1;
						rsltp1p2[*num_of_data_p1p2].dy=he->dy1;
						rsltp1p2[*num_of_data_p1p2].dmin=he->key;
						rsltp1p2[*num_of_data_p1p2].dmax=he->key1;


						*num_of_data_p1p2=*num_of_data_p1p2+1;
						
						//get next data  from heap
						if(*num_of_data_p1p2 == k)
						{
							end = 1;
							
						}
						else
							again =true;						
						
				}
				else //not leaf node
				{
					//printf("Dequeing-------opposite of level\n");
					if(he->level == 0)
					{
						printf("1st Entry is level-0------son becomes -1\n");
						son = -1;

						r1[0] = he->x1;
						r1[1] = he->x2;
						r1[2] = he->y1;
						r1[3] = he->y2;

					}
					if(he->dlevel == 0)
					{
						printf("2nd Entry is level-0------dson becomes -1\n");
						dson = -1;

						r2[0] = he->dx1;
						r2[1] = he->dx2;
						r2[2] = he->dy1;
						r2[3] = he->dy2;
					}

					son_r1 = he->son1;
					son_r2 = he->dson1;
				}
			}

			//printf("outer while loop4 terminated\n");
			/*char s;
			scanf("%c",&s);*/

		}
		
		delete he;
	}
	delete heap;

	printf("terminated twoST_GTP_FA\n\n\n");
	char s;
			scanf("%c",&s);
}


void RTree::threeST_GTP_FA(Rectangle1 R[],Rectangle1 D[], int g_size, Pointlocation rsltp1p2[], int *num_of_data_p1p2p3,int k_points,RTree *dtr,RTree *dtr1)
{	
	printf("called threeST_GTP_FA\n");
	//variable initialization
	//Variables
	int end=0,i,j,k,ctr_i,ctr_j,ctr_k,level_r1,level_r2,level_r3,son=0,dson=0,dson1=0;
	
	double ThresholdDist[2*kMAX];
	//double distance = approximate_mindist1(R,D,g_size,k_points);

	for(int h = 0; h < k_points; h++)	
		ThresholdDist[h] = MAXDOUBLE;
		//ThresholdDist[h] = distance;

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son_r1 = root; //this entry is to be visited next
	int son_r2 = -1;
	int son_r3 = -1;
	int empty = 0;

	while (empty != 1 && end==0)
	{	
		//printf("++++++++++++++++++++outer while loop1 started\n");

		RTNode *rtn1;
		RTNode *rtn2;
		RTNode *rtn3;
		float r1[4],r2[4],r3[4];

		if(son != -1){
		//printf("son-r1--%d",son_r1);
		rtn1 = new RTNode(this, son_r1);
		ctr_i = rtn1 -> num_entries;
		level_r1 = rtn1->level;
		//printf("-----------------------------------------node of rtn1--level= %d \n",level_r1);
		}else{
			
			ctr_i = 1;
			level_r1 = 0;
		}
		
		if(dson != -1){		
		rtn2 = dtr->getRTNode(son_r2);
		ctr_j = rtn2 -> num_entries;
		level_r2 = rtn2->level;
		//printf("-----------------------------------------node of rtn2--level= %d \n",level_r2);
		}else{
			
			ctr_j = 1;
			level_r2 = 0;
		}

		//Added
		if(dson1 != -1){		
		rtn3 = dtr1->getRTNode(son_r3);
		ctr_k = rtn3 -> num_entries;
		level_r3 = rtn3->level;
		
		}else{
			
			ctr_k = 1;
			level_r3 = 0;
		}
		//////////////////


		for (i = 0; i < ctr_i; i++)
		{	
			

			if(son!= -1){
			
				r1[0] = rtn1->entries[i].bounces[0];
				r1[1] = rtn1->entries[i].bounces[1];
				r1[2] = rtn1->entries[i].bounces[2];
				r1[3] = rtn1->entries[i].bounces[3];

				son_r1 = rtn1->entries[i].son;

				
				
			}
			
			for(j = 0; j < ctr_j;j++)
			{

				if(dson!= -1){

				r2[0] = rtn2->entries[j].bounces[0];
				r2[1] = rtn2->entries[j].bounces[1];
				r2[2] = rtn2->entries[j].bounces[2];
				r2[3] = rtn2->entries[j].bounces[3];

				son_r2 = rtn2->entries[j].son;

				//printf("r2 newly initialised.....%d\n",son_r2 );
				char s;
				//	scanf("%c",&s);
				}

				for(k = 0; k < ctr_k;k++)
				{
					//printf("outer while loop4 started\n");
					if(dson1!= -1){

					r3[0] = rtn3->entries[k].bounces[0];
					r3[1] = rtn3->entries[k].bounces[1];
					r3[2] = rtn3->entries[k].bounces[2];
					r3[3] = rtn3->entries[k].bounces[3];

					son_r3 = rtn3->entries[k].son;

					//printf("r3 newly initialised.....%d\n",son_r3 );
					/*char s;
					scanf("%c",&s);*/
					}

					double edist1=0, edist2=0;

					//dmin(r1,r2) computation

					edist1 = sqrt(MINRECTDIST(R[0], r1));
					for (int m=1; m < g_size; m++)
					{				
						edist1 += sqrt(MINRECTDIST(R[m], r1));
					}

					edist1 += (g_size*sqrt(MINRECTDIST(r1,r2)));
					edist1 += (g_size*sqrt(MINRECTDIST(r2,r3)));
					edist1 += sqrt(MINRECTDIST(D[0], r3));
					for (int m=1; m < g_size; m++)
					{				
						edist1 += sqrt(MINRECTDIST(D[m], r3));
					}

					if(edist1 > ThresholdDist[k_points-1])	
					{
						continue;//pruning
					}

					//dmax(r1,r2) computation

					edist2 = sqrt(MAXRECTDIST(R[0], r1));
					for (int m=1; m < g_size; m++)
					{	
						edist2 += sqrt(MAXRECTDIST(R[m], r1));
					}

					edist2 += (g_size*sqrt(MAXRECTDIST(r1,r2)));
					edist2 += (g_size*sqrt(MAXRECTDIST(r2,r3)));
					edist2 += sqrt(MAXRECTDIST(D[0], r3));
					for (int m=1; m < g_size; m++)
					{				
						edist2 += sqrt(MAXRECTDIST(D[m], r3));
					}


					//update ThresholdDist
					for(int l=0; l<k_points; l++)
					{
						if (edist2 < ThresholdDist[l])
						{
							for(int u=k_points-1; u>l; u--)
							{
								ThresholdDist[u]=ThresholdDist[u-1];						
							}
							ThresholdDist[l]=edist2;
							//printf("updating threshold value\n");
							break;
						}
					}
				
					//Enqueue

					HeapEntry *he = new HeapEntry();
					he -> key = edist1;
					he -> key1 = edist2;
					he -> level = level_r1;	//for R1
					he -> son1 = son_r1;
					he -> dlevel = level_r2;	//for R2
					he -> dson1 = son_r2;
					he -> dlevel1 = level_r3;	//for R2
					he -> dson2 = son_r3;

					he-> x1 = r1[0]; //for R1
					he-> x2 = r1[1];
					he-> y1 = r1[2];
					he-> y2 = r1[3];
					he-> dx1 = r2[0]; // for R2
					he-> dx2 = r2[1];
					he-> dy1 = r2[2];
					he-> dy2 = r2[3];
					he-> d1x1 = r3[0]; // for R3
					he-> d1x2 = r3[1];
					he-> d1y1 = r3[2];
					he-> d1y2 = r3[3];

					//if(!he -> level || !he -> dlevel || !he -> dlevel1)
					//{
						//char *s;
						//printf("%d ... %d ... %d\n",he -> level,he -> dlevel,he -> dlevel1);
						//scanf("%c",&s);
					//}

					heap -> insert(he);
					delete he;			


					//printf("outer while loop4 terminated\n");
				}
				
				//printf("outer while loop3 terminated\n");
			}

			//printf("outer while loop2 terminated\n");
		}
		
		if(son != -1)
			delete rtn1;
		if(dson != -1)
			delete rtn2;
		if(dson1 != -1)
			delete rtn3;

		son = 0;
		dson = 0;
		dson1 = 0;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{	
			//printf("outer while loop5 started&&&&&&&&&&&&&&&&&&&&&&&&\n");
			/*char s;
			scanf("%c",&s);*/

			again = false;
			if (!heap->remove(he))  //heap is empty
			{
				empty = 1;
				//printf("Heap is empty-------*********************************\n");
			}
			else
			{	//printf("Heap is not empty-------*********************************\n");
				if (he->level == 0 && he->dlevel == 0 && he->dlevel1 == 0) //p is an object 
				{
						//printf("Dequeing-------*********************************\nFound point");
						char s;
						//scanf("%c",&s);

						//enter into result set
						if(*num_of_data_p1p2p3==MAXDATALIMIT)
						{
							//printf("\nGreater than 10000\n");
							//char s;
							//scanf("%c",&s);
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						}
						
						rsltp1p2[*num_of_data_p1p2p3].x=he->x1;
						rsltp1p2[*num_of_data_p1p2p3].y=he->y1;
						rsltp1p2[*num_of_data_p1p2p3].dx=he->dx1;
						rsltp1p2[*num_of_data_p1p2p3].dy=he->dy1;
						rsltp1p2[*num_of_data_p1p2p3].d1x=he->d1x1;
						rsltp1p2[*num_of_data_p1p2p3].d1y=he->d1y1;
						rsltp1p2[*num_of_data_p1p2p3].dmin=he->key;
						rsltp1p2[*num_of_data_p1p2p3].dmax=he->key1;


						*num_of_data_p1p2p3=*num_of_data_p1p2p3+1;
						
						if(*num_of_data_p1p2p3 == k_points)
						{
							end = 1;
							
						}
						else
							again =true;						
						
				}
				else //not leaf node
				{
					//printf("Dequeing-------opposite of level\n");
					///printf("heap distance = %d\n",he->key);


					if(he->level == 0)
					{
						//printf("1st Entry is level-0------son becomes -1\n");
						son = -1;

						r1[0] = he->x1;
						r1[1] = he->x2;
						r1[2] = he->y1;
						r1[3] = he->y2;

					}
					if(he->dlevel == 0)
					{
						//printf("2nd Entry is level-0------dson becomes -1\n");
						dson = -1;

						r2[0] = he->dx1;
						r2[1] = he->dx2;
						r2[2] = he->dy1;
						r2[3] = he->dy2;
					}

					if(he->dlevel1 == 0)
					{
						//printf("3rd Entry is level-0------dson1 becomes -1\n");
						dson1 = -1;

						r3[0] = he->d1x1;
						r3[1] = he->d1x2;
						r3[2] = he->d1y1;
						r3[3] = he->d1y2;
					}

					son_r1 = he->son1;
					son_r2 = he->dson1;
					son_r3 = he->dson2;
				}
			}

			//printf("outer while loop5 terminated\n");
			//char s;
			//scanf("%c",&s);

		}
		
		delete he;
	}
	delete heap;

	printf("terminated twoST_GTP_FA\n");
}




double RTree::approximate_mindist(Rectangle1 S[],Rectangle1 D[], int g_size,int k)
{
	double Dist1 = 0.0,Dist2 = 0.0;
	Pointlocation rslt1[2*kMAX+1],rslt2[2*kMAX+1]; 
	int num_of_data_P1 = 0;
	int num_of_data_P2 = 0;
	int i,j;


	private_kGNN_sum(S, g_size, 1, rslt1, &num_of_data_P1);//k is always one here
	
	//added p1' to Destination Array
	Rectangle1 *DP = new Rectangle1[g_size+1];
		
	for(i = 0; i < g_size; i++)
	{
		DP[i].x1 = D[i].x1;
		DP[i].y1 = D[i].y1;
		DP[i].x2 = D[i].x2;
		DP[i].y2 = D[i].y2;
	}
	DP[g_size].x1 = rslt1[num_of_data_P1-1].x;
	DP[g_size].y1 = rslt1[num_of_data_P1-1].y;
	DP[g_size].x2 = rslt1[num_of_data_P1-1].x;
	DP[g_size].y2 = rslt1[num_of_data_P1-1].y;

	private_kGNN_sum(DP, g_size+1, k, rslt2, &num_of_data_P2);

	//printf("calculate the total distance");
	//calculate the total distance
	Point2D p;
	p[0] = S[0].x1,p[1] = S[0].y1;
	
	Dist1 += Dist(rslt1[num_of_data_P1-1],p);
	for(i = 1; i < g_size; i++)
	{
		p[0] = S[i].x1;
		p[1] = S[i].y1;
		Dist1 += Dist(rslt1[num_of_data_P1-1],p);
	}
	Dist1 += (g_size*Dist(rslt1[num_of_data_P1-1],rslt2[num_of_data_P2-1]));

	p[0] = D[0].x1,p[1] = D[0].y1;
	Dist1 += Dist(rslt2[num_of_data_P2-1],p);

	for(i = 1; i < g_size; i++)
	{
		p[0] = D[i].x1;
		p[1] = D[i].y1;
		Dist1 += Dist(rslt2[num_of_data_P2-1],p);
	}


	Rectangle1 *SD = new Rectangle1[g_size+g_size];
		
	for(i = 0; i < g_size; i++)
	{
		SD[i].x1 = S[i].x1;
		SD[i].y1 = S[i].y1;
		SD[i].x2 = S[i].x2;
		SD[i].y2 = S[i].y2;
	}
	for(i = g_size; i < (g_size+g_size); i++)
	{
		SD[i].x1 = D[i].x1;
		SD[i].y1 = D[i].y1;
		SD[i].x2 = D[i].x2;
		SD[i].y2 = D[i].y2;
	}

	private_kGNN_sum(SD, g_size+g_size, 1, rslt1, &num_of_data_P1);//k is always one here
	
		
	for(i = 0; i < g_size; i++)
	{
		DP[i].x1 = D[i].x1;
		DP[i].y1 = D[i].y1;
		DP[i].x2 = D[i].x2;
		DP[i].y2 = D[i].y2;
	}
	DP[g_size].x1 = rslt1[num_of_data_P1-1].x;
	DP[g_size].y1 = rslt1[num_of_data_P1-1].y;
	DP[g_size].x2 = rslt1[num_of_data_P1-1].x;
	DP[g_size].y2 = rslt1[num_of_data_P1-1].y;

	private_kGNN_sum(DP, g_size+1, k, rslt2, &num_of_data_P2);

	//calculate the aggregate distance
	p[0] = S[0].x1,p[1] = S[0].y1;
	
	Dist2 += Dist(rslt1[num_of_data_P1-1],p);
	for(i = 1; i < g_size; i++)
	{
		p[0] = S[i].x1;
		p[1] = S[i].y1;
		Dist1 += Dist(rslt1[num_of_data_P1-1],p);
	}
	Dist2 += (g_size*Dist(rslt1[num_of_data_P1-1],rslt2[num_of_data_P2-1]));

	p[0] = D[0].x1,p[1] = D[0].y1;
	Dist2 += Dist(rslt2[num_of_data_P2-1],p);

	for(i = 1; i < g_size; i++)
	{
		p[0] = D[i].x1;
		p[1] = D[i].y1;
		Dist2 += Dist(rslt2[num_of_data_P2-1],p);
	}

	if(Dist1 < Dist2)	
		return Dist1;
	else				
		return Dist2;


}

double RTree::approximate_mindist1(Rectangle1 S[],Rectangle1 D[], int g_size,int k)
{
	double Dist1 = 0.0,Dist2 = 0.0;
	Pointlocation rsltp1[2*kMAX+1],rsltp2[2*kMAX+1],rsltp3[2*kMAX+1]; 
	int num_of_data_P1 = 0;
	int num_of_data_P2 = 0;
	int num_of_data_P3 = 0;
	int i,j,w,g;
	
	Rectangle1 *P1 = new Rectangle1[2*kMAX+1];
	Rectangle1 *P2 = new Rectangle1[2*kMAX+1];
	Rectangle1 *P3 = new Rectangle1[2*kMAX+1];


	private_kGNN_sum(S, g_size, 1, rsltp1, &num_of_data_P1);//k is always one here
	
	//copy p1 to rectangle
	P1[0].x1 = rsltp1[num_of_data_P1-1].x;
	P1[1].y1 = rsltp1[num_of_data_P1-1].y;
	P1[2].x2 = rsltp1[num_of_data_P1-1].x;
	P1[3].y2 = rsltp1[num_of_data_P1-1].y;

	/**********************************************/

	private_kGNN_sum(P1, num_of_data_P1, 1, rsltp2, &num_of_data_P2);//k is always one here
	P2[0].x1 = rsltp2[num_of_data_P2-1].x;
	P2[1].y1 = rsltp2[num_of_data_P2-1].y;
	P2[2].x2 = rsltp2[num_of_data_P2-1].x;
	P2[3].y2 = rsltp2[num_of_data_P2-1].y;
	

	
	/*********************************************/

	//added p1' to Destination Array
	Rectangle1 *DP = new Rectangle1[g_size+1];
		
	for(i = 0; i < g_size; i++)
	{
		DP[i].x1 = D[i].x1;
		DP[i].y1 = D[i].y1;
		DP[i].x2 = D[i].x2;
		DP[i].y2 = D[i].y2;
	}
	DP[g_size].x1 = rsltp2[num_of_data_P2-1].x;
	DP[g_size].y1 = rsltp2[num_of_data_P2-1].y;
	DP[g_size].x2 = rsltp2[num_of_data_P2-1].x;
	DP[g_size].y2 = rsltp2[num_of_data_P2-1].y;

	private_kGNN_sum(DP, g_size+1, k, rsltp3, &num_of_data_P3);

	//printf("calculate the total distance");
	//calculate the total distance
	Point2D p;
	p[0] = S[0].x1,p[1] = S[0].y1;
	
	Dist1 += Dist(rsltp1[num_of_data_P1-1],p);
	for(i = 1; i < g_size; i++)
	{
		p[0] = S[i].x1;
		p[1] = S[i].y1;
		Dist1 += Dist(rsltp1[num_of_data_P1-1],p);
	}
	Dist1 += (g_size*Dist(rsltp1[num_of_data_P1-1],rsltp2[num_of_data_P2-1]));
	Dist1 += (g_size*Dist(rsltp2[num_of_data_P2-1],rsltp3[num_of_data_P3-1]));

	p[0] = D[0].x1,p[1] = D[0].y1;
	Dist1 += Dist(rsltp3[num_of_data_P3-1],p);

	for(i = 1; i < g_size; i++)
	{
		p[0] = D[i].x1;
		p[1] = D[i].y1;
		Dist1 += Dist(rsltp3[num_of_data_P3-1],p);
	}

	//heuristic-2

	Rectangle1 *SD = new Rectangle1[g_size+g_size];
		
	for(i = 0; i < g_size; i++)
	{
		SD[i].x1 = S[i].x1;
		SD[i].y1 = S[i].y1;
		SD[i].x2 = S[i].x2;
		SD[i].y2 = S[i].y2;
	}
	for(i = 0,g=g_size; i < g_size; i++,g++)
	{
		SD[g].x1 = D[i].x1;
		SD[g].y1 = D[i].y1;
		SD[g].x2 = D[i].x2;
		SD[g].y2 = D[i].y2;
	}

	num_of_data_P1 = 0;
	private_kGNN_sum(SD, g_size+g_size, 1, rsltp1, &num_of_data_P1);//k is always one here
	
		
	for(i = 0; i < g_size; i++)
	{
		DP[i].x1 = D[i].x1;
		DP[i].y1 = D[i].y1;
		DP[i].x2 = D[i].x2;
		DP[i].y2 = D[i].y2;
	}
	DP[g_size].x1 = rsltp1[num_of_data_P1-1].x;
	DP[g_size].y1 = rsltp1[num_of_data_P1-1].y;
	DP[g_size].x2 = rsltp1[num_of_data_P1-1].x;
	DP[g_size].y2 = rsltp1[num_of_data_P1-1].y;
	
	num_of_data_P2=0;
	private_kGNN_sum(DP, g_size+1, 1, rsltp2, &num_of_data_P2);

	DP[g_size].x1 = rsltp2[num_of_data_P2-1].x;
	DP[g_size].y1 = rsltp2[num_of_data_P2-1].y;
	DP[g_size].x2 = rsltp2[num_of_data_P2-1].x;
	DP[g_size].y2 = rsltp2[num_of_data_P2-1].y;
	
	num_of_data_P3=0;
	private_kGNN_sum(DP, g_size+1, k, rsltp3, &num_of_data_P3);



	//calculate the aggregate distance
	p[0] = S[0].x1,p[1] = S[0].y1;
	
	Dist2 += Dist(rsltp1[num_of_data_P1-1],p);
	for(i = 1; i < g_size; i++)
	{
		p[0] = S[i].x1;
		p[1] = S[i].y1;
		Dist1 += Dist(rsltp1[num_of_data_P1-1],p);
	}
	Dist2 += (g_size*Dist(rsltp1[num_of_data_P1-1],rsltp2[num_of_data_P2-1]));
	Dist2 += (g_size*Dist(rsltp2[num_of_data_P2-1],rsltp3[num_of_data_P3-1]));

	p[0] = D[0].x1,p[1] = D[0].y1;
	Dist2 += Dist(rsltp3[num_of_data_P3-1],p);

	for(i = 1; i < g_size; i++)
	{
		p[0] = D[i].x1;
		p[1] = D[i].y1;
		Dist2 += Dist(rsltp3[num_of_data_P3-1],p);
	}

	if(Dist1 < Dist2)	
		return Dist1;
	else				
		return Dist2;


}
