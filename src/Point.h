#ifndef POINT_H
#define POINT_H

class Point
{
public:
	Point(double x, double y, double z) : x(x), y(y), z(z) {};
	~Point() = default;
	
	double x;
	double y;
	double z;
};

#endif