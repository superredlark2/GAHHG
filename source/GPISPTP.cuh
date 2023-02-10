//GPISPTP.cuh
#pragma once

#include "GPISPPolygon.cuh"

struct GPISPTP {
	
	int tpNum;
	LLVertex* tpll;
	mydouble3* tp;
	int* result;
	devdouble3* dev_tp;
	int* dev_result;
	
	GPISPTP() { memset(this, 0, sizeof(*this)); }
	~GPISPTP() { clear(); }
	int loadTPFromCSV(char* filename);
	void exportResult2CSV(char* filename);
	void clear();
};
