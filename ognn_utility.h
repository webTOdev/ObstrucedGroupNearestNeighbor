/*
 * ognn_utility.h
 *
 *  Created on: Apr 2, 2014
 *      Author: nut
 */

#ifndef OGNN_UTILITY_H_
#define OGNN_UTILITY_H_

double getDistanceBetweenTwoPoints(Point2D p,Point2D q);
void centroidOfQ(Point2D queryPoints[], int numOfQueryPoints, Point2D centroid[]);
//Use this struct to store the distance between p and q against the q
struct MyStruct
{
    double distance;
    float* queryPoints;

	MyStruct(double d , float* q) { distance = d; queryPoints =  q ;}
};
//Use this struct to sort  q according to distance descending
struct more_than_key
{
    inline bool operator() (const MyStruct& struct1, const MyStruct& struct2)
    {
        return (struct1.distance > struct2.distance);
    }
};

//Use this struct to sort  q according to distance ascending
struct less_than_key
{
    inline bool operator() (const MyStruct& struct1, const MyStruct& struct2)
    {
        return (struct1.distance < struct2.distance);
    }
};

struct MyShortestPath
{
    vector<int> shortestPath;
    float* queryPoints;

	MyShortestPath(vector<int> sp , float* q) { shortestPath = sp; queryPoints =  q ;}
};
class Clock {

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


struct MyVertexStruct
{
    double distance;
    Point* vertex;

	MyVertexStruct(double d , Point* p) { distance = d; vertex =  p ;}
};

#endif /* OGNN_UTILITY_H_ */
