/* this file contains implementation of the class Heap */
#ifndef HEAP_H
#define HEAP_H

#define MAX_HEAP_SIZE 1110000

class HeapEntry;
class Heap;

class HeapEntry
{
public:
	int dim;
	int level;
	int son1;
	int dlevel;
	int dson1;
	int son2;	
	int dlevel1;
	int dson2;
	//int son3;	
	float key;  //entries will be sorted according to this
	float key1;
	//Added for Rect-kNNQ by Tanzima
	float x1;
	float y1;
	float x2;
	float y2;

	//Added for 2S-GTPQ by Tanzima
	float dx1;
	float dy1;
	float dx2;
	float dy2;

	//Added for 3S-GTPQ by Tahrima
	float d1x1;
	float d1y1;
	float d1x2;
	float d1y2;

	//-----functions-----
	HeapEntry();
	~HeapEntry();
	void init_HeapEntry(int _dim);
	void copy(HeapEntry *_he);
};

class Heap
{
public:
	int hsize;        // the heap size
	int used;         // number of used places
	int maxused;
	HeapEntry *cont;  // content of the heap

	//-----functions-----
	Heap();
	~Heap();
	bool check();
	void clean(float _dist);
	void enter(HeapEntry *_he, int _pos);
	void init(int _dim, int _hsize=MAX_HEAP_SIZE);
	void insert(HeapEntry *_he);
	bool remove(HeapEntry *_he);
		void copy(Heap* h);
};

#endif