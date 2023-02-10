//GPISPDS.cuh
#pragma once

#include <device_launch_parameters.h>
#include <math.h>

#define PI 3.14159265

struct devdouble3 {
	union {
		double data[3];
		struct {
			double x;
			double y;
			double z;
		};
	};

	__device__ devdouble3() { x = 0.0; y = 0.0; z = 0.0; }
	__device__ devdouble3(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }

	__device__ devdouble3& operator=(const devdouble3& p)
	{
		x = p.x;    //等于this.x = p.x
		y = p.y;
		z = p.z;
		return *this; //最后把this的值返回
	}

	__device__ devdouble3 operator+(const devdouble3& f)
	{
		devdouble3 temp;
		temp.x = x + f.x;
		temp.y = y + f.y;
		temp.z = z + f.z;
		return temp;
	}

	__device__ devdouble3 operator-(const devdouble3& f)
	{
		devdouble3 temp;
		temp.x = x - f.x;
		temp.y = y - f.y;
		temp.z = z - f.z;
		return temp;
	}

	__device__ devdouble3 operator-()
	{
		devdouble3 temp;
		temp.x = -x;
		temp.y = -y;
		temp.z = -z;
		return temp;
	}

	__device__ devdouble3 operator*(double k)
	{
		devdouble3 temp;
		temp.x = x * k;
		temp.y = y * k;
		temp.z = z * k;
		return temp;
	}

	__device__ devdouble3 operator/(double k)
	{
		double s = (double)1.0 / k;
		devdouble3 temp;
		temp.x = x * s;
		temp.y = y * s;
		temp.z = z * s;
		return temp;
	}

	__device__ void normalize()
	{
		double temp = 1.0 / sqrt(x * x + y * y + z * z);
		x = (double)(x * temp);
		y = (double)(y * temp);
		z = (double)(z * temp);
	}

	__device__ double len()
	{
		return sqrt(x * x + y * y + z * z);
	}

	__device__ double squarelen()
	{
		return x * x + y * y + z * z;
	}

	__device__ devdouble3 cross(devdouble3 v)
	{
		devdouble3 result;
		result.x = this->y * v.z - this->z * v.y;
		result.y = this->z * v.x - this->x * v.z;
		result.z = this->x * v.y - this->y * v.x;
		return result;
	}

	__device__ double dot(devdouble3 v)
	{
		return this->x * v.x + this->y * v.y + this->z * v.z;
	}
};

