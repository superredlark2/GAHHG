//GPISPAHHG.cuh
#pragma once

#include "GPISPDS.cuh"
#include "GPISPPolygon.cuh"
#include "GPISPParam.h"
#include "GPISPTP.cuh"
#include "common.h"

#define mymax(a,b) a>b?a:b

struct GAHHGGrid {
	int vstart;
	int vnum;
	int bstart;
	int bnum;
	int cstart;
	int cnum;
	int rnum;
};

struct GAHHG {
	
	GPISPParam* pParam;
	GPISPStat* pStat;
	GPISPPolygon* pSP;
	
	GAHHGGrid* g;
	int gnum;
	myArray<mydouble3>* vertexes;
	myArray<AHHGBorder>* borders;
	myArray<AHHGCell>* cells;

	devdouble3* dev_v;
	AHHGBorder* dev_b;
	AHHGCell* dev_c;
	int* dev_r[3];
	int* dev_se_buf[2];
	int* dev_cell_buf[2];

	GAHHG(GPISPPolygon* poly, GPISPParam* param, GPISPStat* stat)
	{
		memset(this, 0, sizeof(*this)); 
		this->pSP = poly;
		this->pParam = param;
		this->pStat = stat;
	}
	~GAHHG() { clear(); }
	void initDevMemory();
	void copyGrid0HostToDevice();
	void dispatchSELast();
	bool divGrid(int gridid);
	void classifyBottomGrid();
	void checkStorage(int arrayid, int maxnum);
	void clear();
};


struct GPISP {
	GPISPParam* pParam;
	GPISPStat* pStat;
	GPISPPolygon* pPolygon;
	GPISPTP* pTP;
	GAHHG* pAHHG;

	GPISP() { memset(this, 0, sizeof(*this)); }
	~GPISP() { clear(); }
	void createAHHGInit();
	void createAHHG(int maxdepth);
	void classifyAHHG();
	void test();
	void getStat();										
	void exportStat2TXT(char* filename, int desttype);
	void clear();
};
