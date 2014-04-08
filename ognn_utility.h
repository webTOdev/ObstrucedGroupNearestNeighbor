/*
 * ognn_utility.h
 *
 *  Created on: Apr 2, 2014
 *      Author: nut
 */

#ifndef OGNN_UTILITY_H_
#define OGNN_UTILITY_H_

double getDistanceBetweenTwoPoints(Point2D p,Point2D q);
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

#endif /* OGNN_UTILITY_H_ */
