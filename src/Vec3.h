#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

class Vec3
{
public:
	Vec3() = default;
	Vec3(double x, double y, double z) : x(x), y(y), z(z) {};
	~Vec3() = default;
	
	double x = 0;
	double y = 0;
	double z = 0;

	double norm() const { return sqrt(x*x + y*y + z*z); }

	void normalize();


};


void Vec3::normalize() {
	double n = norm();
	if (n > 0) {
		x /= n;
		y /= n;
		z /= n;
	}
	else {
		std::cout << "zero vector can't be normalized!"
	}	
}

Vec3 operator-(const Vec3 &v1, const Vec3 &v2) {
	return Vec3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}

Vec3 cross(const Vec3 &v1, const Vec3 &v2) {
	return Vec3(v1.y * v2.z - v1.z * v2.y,
			v1.z * v2.x - v1.x * v2.z,
			v1.x * v2.y - v1.y * v2.x
		);
}




#endif