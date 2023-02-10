//GPISPParam.h
#pragma once

#include <windows.h>

#define TRAVEL_EPSILON 1.0e-5
#define POSITION_EPSILON 1.0e-8

#define CELL_UNDEFINED 0
#define CELL_OUT 1
#define CELL_IN 2

#define INTER_NO 0
#define INTER_X 1
#define INTER_T1 2
#define INTER_T2 3
#define INTER_OVERLAP 4

#define EQU_ZERO(a) (a<POSITION_EPSILON && a>-POSITION_EPSILON)
#define LESS_ZERO(a) (a<-POSITION_EPSILON)
#define GRE_ZERO(a) (a>POSITION_EPSILON)
#define LESS_EQU_ZERO(a) (a<POSITION_EPSILON)
#define GRE_EQU_ZERO(a) (a>-POSITION_EPSILON)

struct myClock {
	LARGE_INTEGER tick1, tick2;
	double dff;
	double time;

	myClock() {
		memset(this, 0, sizeof(*this));
		QueryPerformanceFrequency(&tick1);
		dff = tick1.QuadPart;
	}

	void start() {
		QueryPerformanceCounter(&tick1);
	}

	void end() {
		QueryPerformanceCounter(&tick2);
		__int64 c1, c2;
		c1 = tick1.QuadPart;
		c2 = tick2.QuadPart;
		time = (c2 - c1) / dff;	//unit£ºsecond
	}
};

struct GPISPParam {
	char polyFileName[128];
	char reportFileName[128];
	char resultFileName[128];
	char tpFileName[128];
	int maxDepth;
	
	int blocksize_create;
	int blocksize_classify;
	int blocksize_pisp;
	int max_v_num;
	int max_b_num;
	int max_c_num;
	int max_r_num;
	int max_buf_num;
	bool autoDepth;
	
	GPISPParam() {
		memset(this, 0, sizeof(*this));
		this->maxDepth = 15;
		this->autoDepth = true;
		this->blocksize_create = 256;
		this->blocksize_classify = 256;
		this->blocksize_pisp = 256;
		strcpy(this->polyFileName,"NA0.txt");
		strcpy(this->tpFileName,"NA_query.csv");
		strcpy(this->reportFileName, "report.txt");
		strcpy(this->resultFileName, "result.csv");
	}
};

struct GPISPStat {
	double preTime;
	double createAHHGTime;
	double classifyAHHGTime;
	int AHHG_cpu_space;
	int AHHG_gpu_min_space;
	int AHHG_gpu_vnum;					
	int AHHG_gpu_bnum;
	int AHHG_gpu_cnum;
	int AHHG_gpu_rnum;
	int AHHG_gpu_vspace;
	int AHHG_gpu_bspace;
	int AHHG_gpu_cspace;
	int AHHG_gpu_rspace;
	float gamma;
	double PISPTime;
	double PISPUnitTime;
	int inTpNum;

	GPISPStat() { memset(this, 0, sizeof(*this)); }
};
