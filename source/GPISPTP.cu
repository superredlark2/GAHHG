//GPISPTP.cu

#include <cuda_runtime.h>
#include <fstream>
#include <iostream>
#include <strstream>
#include <string.h>
#include "GPISPTP.cuh"
#include "GPISPParam.h"

using namespace std;

int GPISPTP::loadTPFromCSV(char* filename)
{
	printf("Load test points from a CSV file: %s.\n", filename);

	ifstream tpFile(filename, ios::binary);
	if (!tpFile)
	{
		printf("Error. Can not open %s\n.", filename);
		return -1;
	}

	tpFile.seekg(0, ios_base::end);
	streampos size = tpFile.tellg();
	tpFile.seekg(0, ios_base::beg);
	char buf[1024];

	int totalvnum = 0;
	while (!tpFile.eof())
	{
		tpFile.getline(buf, size);
		if (tpFile.fail())
			break;
		double la, lo;
		int ret = sscanf_s(buf, "%lf,%lf", &lo, &la);
		if (ret != 2)
		{
			printf("Error. Line %d of %s. Wrong format.\n", totalvnum+1, filename);
			return -1;
		}
		totalvnum++;
	}

	int multiple = 1;
	this->tpNum = totalvnum * multiple;
	this->tpll = new LLVertex[this->tpNum];
	this->tp = new mydouble3[this->tpNum];
	this->result = new int[this->tpNum];
	 
	tpFile.clear();
	tpFile.seekg(0, ios::beg);
	for (int i = 0; i < totalvnum; i++)
	{
		tpFile.getline(buf, size);
		int ret = sscanf_s(buf, "%lf,%lf", &this->tpll[i].lo, &this->tpll[i].la);
		FakeLLVertex fakelltp(this->tpll[i]);
		this->tp[i] = fakelltp.FakeLL2XYZ();
	}
	tpFile.close();

	if (multiple > 1)
	{
		for (int k = 1; k < multiple; k++)
		{
			for (int i = 0; i < totalvnum; i++)
			{
				this->tpll[multiple * k + i] = this->tpll[i];
				this->tp[multiple * k + i] = this->tp[i];
			}
		}
	}
	
	cudaMalloc((void**)&this->dev_tp, sizeof(devdouble3) * this->tpNum);
	cudaMemcpy(this->dev_tp, this->tp, sizeof(devdouble3) * this->tpNum, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&this->dev_result, sizeof(int) * this->tpNum);

	return 0;
}

void GPISPTP::exportResult2CSV(char* filename)
{
	FILE* fp;
	int ret = fopen_s(&fp, filename, "w");
	if (ret != 0)
	{
		printf("Can not open file %s\n.", filename);
		return;
	}

	for (int i = 0; i < this->tpNum; i++)
	{
		if (this->result[i] == CELL_IN)
			fprintf(fp, "1\n");
		else
			fprintf(fp, "0\n");
	}
	fclose(fp);
}

void GPISPTP::clear()
{
	if (this->tp != NULL)
		delete[] this->tp;
	if (this->tpll != NULL)
		delete[] this->tpll;
	if (this->result != NULL)
		delete[] this->result;
	if (this->dev_tp != NULL)
		cudaFree(this->dev_tp);
	if (this->dev_result != NULL)
		cudaFree(this->dev_result);
	memset(this, 0, sizeof(*this));
}
