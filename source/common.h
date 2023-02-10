//common.h
#pragma once

#include <math.h>
#include <string.h>

struct mydouble3 {
	union {
		double data[3];
		struct {
			double x;
			double y;
			double z;
		};
	};

	mydouble3() { x = 0; y = 0; z = 0; }
	mydouble3(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }

	mydouble3 operator+(const mydouble3& f)
	{
		mydouble3 temp;
		temp.x = x + f.x;
		temp.y = y + f.y;
		temp.z = z + f.z;
		return temp;
	}

	mydouble3 operator-(const mydouble3& f)
	{
		mydouble3 temp;
		temp.x = x - f.x;
		temp.y = y - f.y;
		temp.z = z - f.z;
		return temp;
	}

	mydouble3 operator-()
	{
		mydouble3 temp;
		temp.x = -x;
		temp.y = -y;
		temp.z = -z;
		return temp;
	}

	mydouble3 operator*(double k)
	{
		mydouble3 temp;
		temp.x = x * k;
		temp.y = y * k;
		temp.z = z * k;
		return temp;
	}

	mydouble3 operator/(double k)
	{
		double s = (double)1.0 / k;
		mydouble3 temp;
		temp.x = x * s;
		temp.y = y * s;
		temp.z = z * s;
		return temp;
	}

	void normalize()
	{
		double temp = 1.0 / sqrt(x * x + y * y + z * z);
		x = (double)(x * temp);
		y = (double)(y * temp);
		z = (double)(z * temp);
	}

	double len()
	{
		return sqrt(x * x + y * y + z * z);
	}

	double squarelen()
	{
		return x * x + y * y + z * z;
	}

	mydouble3 cross(mydouble3 v)
	{
		mydouble3 result;
		result.x = this->y * v.z - this->z * v.y;
		result.y = this->z * v.x - this->x * v.z;
		result.z = this->x * v.y - this->y * v.x;
		return result;
	}

	double dot(mydouble3 v)
	{
		return this->x * v.x + this->y * v.y + this->z * v.z;
	}
};

template <typename T>
struct myArray {
	T* elem;
	int curNum;
	int maxNum;

	myArray() { memset(this, 0, sizeof(*this)); }
	myArray(int count)
	{
		curNum = 0;
		maxNum = count;
		elem = (T*)malloc(sizeof(T) * maxNum);
	}
	~myArray()
	{
		if (elem != NULL)
			free(elem);
	}

	void addElem(T data)
	{
		if (elem == NULL)
		{
			curNum = 0;
			maxNum = 8;
			elem = (T*)malloc(sizeof(T) * maxNum);
		}

		if (curNum == maxNum)
		{
			maxNum *= 2;
			elem = (T*)realloc(elem, sizeof(T) * maxNum);
		}

		elem[curNum] = data;
		curNum++;
	}

	bool findElem(T data)
	{
		for (int i = 0; i < this->curNum; i++)
		{
			if (this->elem[i] == data)
				return true;
		}
		return false;
	}

	T operator[](int k)
	{
		return elem[k];
	}

	void clear()
	{
		if (this->elem != NULL)
		{
			free(this->elem);
			memset(this, 0, sizeof(*this));
		}
	}
};

struct icosahedron {
	mydouble3 v[12];
	int f[20][3];
	icosahedron();
	~icosahedron() {};
};

struct AHHGBorder {
	int parent;
	int partner;
	int v[2];
	int subcell;
	int id;
	double n[3];
};

struct AHHGCell {
	int bstart;
	int bnum;
	int rstart;
	int rnum;
	int upcell;
	int subcell;
	int type;
};

struct grid0{
	myArray<mydouble3>* vertexes;
	myArray<AHHGBorder>* borders;
	myArray<AHHGCell>* cells;
	
	void swapCells(int cell0, int cell1);
	void getDiamondPoints(icosahedron* icosa, int id, mydouble3* points);
	mydouble3 getPointOnSphere(mydouble3 v0, mydouble3 v1, double ratio);
	mydouble3 getEdgeNormal(int start, int end);
	void setEdgePartner(int cell0, int edge0, int cell1, int edge1);
	void create();
};

void initgrid0(myArray<mydouble3>* v, myArray<AHHGBorder>* b, myArray<AHHGCell>* c);

