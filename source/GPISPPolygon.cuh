//GPISPPolygon.h
#pragma once

#include <stdio.h>
#include <windows.h>
#include "GPISPDS.cuh"
#include "common.h"

struct LLVertex {
	union {
		double c[2];
		struct {
			double la;
			double lo;
		};
	};

	LLVertex() { la = lo = 0.0; }
	LLVertex(double latitude, double longitude) { la = latitude; lo = longitude; }

	mydouble3  LL2XYZ()
	{
		mydouble3 p;
		double a = la * PI / 180.0;
		double o = lo * PI / 180.0;
		p.x = cos(a) * sin(o);
		p.y = sin(a);
		p.z = cos(a) * cos(o);
		return p;
	}

	LLVertex XYZ2LL(mydouble3 sp)
	{
		LLVertex p;
		//p.lo = atan2(sp.x, sp.z) * 180.0 / PI;
		//p.la = asin(sp.y) * 180.0 / PI;
		p.la = asin(sp.y);
		if (sp.x == 0 && sp.z == 0)
			p.lo = 0;
		else
		{
			double len = sqrt(sp.x * sp.x + sp.z * sp.z);
			p.lo = acos(sp.z / len);
			if (sp.x < 0)
				p.lo = -p.lo;
		}
		p.la = p.la * 180.0 / PI;
		p.lo = p.lo * 180.0 / PI;
		return p;
	}

	void translate(double la_offset, double lo_offset)
	{
		la += la_offset;
		la += la >= 180 ? -360 : la < -180 ? 360 : 0;
		lo += lo_offset;
		lo = lo > 90 ? 180 - lo : lo < -90 ? -180 - lo : lo;
	}
};

struct FakeLLVertex {
	union {
		double c[2];
		struct {
			double fakela;	//latitude
			double lo;	//longitude
		};
	};

	FakeLLVertex() { fakela = lo = 0.0; }

	FakeLLVertex(LLVertex llp) { lo = llp.lo; fakela = llp.la + 90; }

	LLVertex ToLLVertex()
	{
		LLVertex ll;
		ll.lo = this->lo;
		ll.la = this->fakela > 180 ? 270 - this->fakela : this->fakela - 90;
		return ll;
	}

	mydouble3  FakeLL2XYZ()
	{
		double realla = this->fakela > 180 ? 270 - this->fakela : this->fakela - 90;
		mydouble3 p;
		double a = realla * PI / 180.0;
		double o = this->lo * PI / 180.0;
		p.x = cos(a) * sin(o);
		p.y = sin(a);
		p.z = cos(a) * cos(o);
		return p;
	}

	FakeLLVertex XYZ2FakeLL(mydouble3 sp)
	{
		FakeLLVertex fakell;
		fakell.lo = atan2(sp.x, sp.z) * 180.0 / PI;
		fakell.fakela = asin(sp.y) * 180.0 / PI + 90;
		return fakell;
	}

	void translate(double la_offset, double lo_offset)
	{
		this->fakela += la_offset;
		this->fakela += this->fakela >= 360 ? -360 : this->fakela < 0 ? 360 : 0;
		this->lo += lo_offset;
		this->lo += this->lo >= 180 ? -360 : this->lo < -180 ? 360 : 0;
	}
};

struct devFakeLLVertex {
	union {
		double c[2];
		struct {
			double fakela;	//Î³¶È
			double lo;	//¾­¶È
		};
	};

	__device__ devFakeLLVertex() { fakela = lo = 0.0; }

	__device__ devFakeLLVertex(LLVertex llp) { lo = llp.lo; fakela = llp.la + 90; }

	__device__ LLVertex ToLLVertex()
	{
		LLVertex ll;
		ll.lo = this->lo;
		ll.la = this->fakela > 180 ? 270 - this->fakela : this->fakela - 90;
		return ll;
	}

	__device__ devdouble3  FakeLL2XYZ()
	{
		double realla = this->fakela > 180 ? 270 - this->fakela : this->fakela - 90;
		devdouble3 p;
		double a = realla * PI / 180.0;
		double o = this->lo * PI / 180.0;
		p.x = cos(a) * sin(o);
		p.y = sin(a);
		p.z = cos(a) * cos(o);
		return p;
	}

	__device__ devFakeLLVertex XYZ2FakeLL(devdouble3 sp)
	{
		devFakeLLVertex fakell;
		fakell.lo = atan2(sp.x, sp.z) * 180.0 / PI;
		fakell.fakela = asin(sp.y) * 180.0 / PI + 90;
		return fakell;
	}

	__device__ void translate(double la_offset, double lo_offset)
	{
		this->fakela += la_offset;
		this->fakela += this->fakela >= 360 ? -360 : this->fakela < 0 ? 360 : 0;
		this->lo += lo_offset;
		this->lo += this->lo >= 180 ? -360 : this->lo < -180 ? 360 : 0;
	}
};


struct GPISPPolygon {

	int pv_num;
	LLVertex* pllv;
	mydouble3* pv;
	int* pe;
	LLVertex min, max;
	devdouble3* dev_pv;
	int* dev_pe;

	GPISPPolygon() { memset(this, 0, sizeof(*this)); }
	~GPISPPolygon() { clear(); }
	bool isPolyClockwise(LLVertex* vlist, int vnum);
	void loadPolyFreeTempMemory(int* inVNum, LLVertex* pOutVList, LLVertex** pInVList, int inNum, char* buf);
	int loadPolyFromTXT(char* filename);
	void clear();
};

