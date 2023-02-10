//GPISPPolygon.cu

#include <device_launch_parameters.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <fstream>
#include "GPISPPolygon.cuh"

using namespace std;

void GPISPPolygon::loadPolyFreeTempMemory(int* inVNum, LLVertex* pOutVList, LLVertex** pInVList, int inNum, char* buf)
{
	if (inVNum != NULL)
		delete inVNum;
	if (pOutVList != NULL)
		delete pOutVList;
	if (pInVList != NULL)
	{
		for (int i = 0; i < inNum; i++)
		{
			if (pInVList[i] != NULL)
			{
				delete pInVList[i];
			}
		}
		delete pInVList;
	}
	if (buf != NULL)
		delete buf;
}

int GPISPPolygon::loadPolyFromTXT(char* filename)
{
	printf("Load spherical polygon from %s.\n", filename);

	ifstream polyFile(filename, ios::binary);
	if (!polyFile)
	{
		printf("Error. Can not open %s\n.", filename);
		return -1;
	}

	int outVNum = 0;
	int inNum = 0;
	int temp;
	int* inVNum = NULL;
	int totalvnum = 0;
	int ret=-1;
	int lineindex = 0;
	LLVertex* pOutVList = NULL;
	LLVertex** pInVList = NULL;

	polyFile.seekg(0, ios_base::end);
	streampos size = polyFile.tellg();
	polyFile.seekg(0, ios_base::beg);
	char* buf = new char[size];

	polyFile.getline(buf, size);
	lineindex++;
	ret = sscanf_s(buf, "POLYGON 1 %d", &inNum);
	if (ret != 1)
	{
		loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
		printf("Error. Line %d of %s. Wrong format.\n", lineindex, filename);
		return -1;
	}
	
	if (inNum != 0)
	{
		inVNum = new int[inNum];
		pInVList = new LLVertex * [inNum];
		for (int u = 0; u < inNum; u++)
			pInVList[u] = NULL;
	}
	
	polyFile.getline(buf, size); 
	lineindex++;
	ret = sscanf_s(buf, "0 %d", &outVNum);
	if (ret != 1)
	{
		loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
		printf("Error. Line %d of %s. Wrong format.\n", lineindex, filename);
		return -1;
	}
	else
	{
		if (outVNum < 4)
		{
			loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
			printf("Error. Line %d of %s. Number of edges should be greater than 3. \n", lineindex, filename);
			return -1;
		}
	}

	for (int i = 0; i < outVNum; i++)
	{
		polyFile.getline(buf, size);
		lineindex++;
		double la, lo;
		ret = sscanf_s(buf, "%lf %lf", &lo, &la);
		if (ret != 2)
		{
			loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
			printf("Error. Line %d of %s. Wrong format.\n", lineindex, filename);
			return -1;
		}
	}

	pOutVList = new LLVertex[outVNum];
	totalvnum += outVNum - 1;

	for (int i = 0; i < inNum; i++)
	{
		polyFile.getline(buf, size);
		lineindex++;
		ret = sscanf_s(buf, "%d %d", &temp, &inVNum[i]);
		if (ret != 2)
		{
			loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
			printf("Error. Line %d of %s. Wrong format.\n", lineindex, filename);
			return -1;
		}
		else
		{
			if (inVNum[i] < 4)
			{
				loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
				printf("Error. Line %d of %s. Number of edges should be greater than 3. \n", lineindex, filename);
				return -1;
			}
		}

		for (int j = 0; j < inVNum[i]; j++)
		{
			polyFile.getline(buf, size);
			lineindex++;
			double la, lo;
			ret = sscanf_s(buf, "%lf %lf", &lo, &la);
			if (ret != 2)
			{
				loadPolyFreeTempMemory(inVNum, pOutVList, pInVList, inNum, buf);
				printf("Error. Line %d of %s. Wrong format.\n", lineindex, filename);
				return -1;
			}
		}
		pInVList[i] = new LLVertex[inVNum[i]];
		totalvnum += inVNum[i] - 1;
	}

	this->pv_num = totalvnum;
	this->pllv = new LLVertex[totalvnum];
	this->pv = new mydouble3[totalvnum];
	this->pe = new int[totalvnum * 2];

	polyFile.clear();
	polyFile.seekg(0, ios::beg);
	totalvnum = 0;

	polyFile.getline(buf, size);
	sscanf_s(buf, "POLYGON %d %d", &temp, &inNum);
	polyFile.getline(buf, size);
	sscanf_s(buf, "%d %d", &temp, &outVNum);

	for (int i = 0; i < outVNum - 1; i++)
	{
		polyFile.getline(buf, size);
		sscanf_s(buf, "%lf %lf", &pOutVList[i].lo, &pOutVList[i].la);
		this->pllv[totalvnum] = pOutVList[i];
		this->pv[totalvnum] = pOutVList[i].LL2XYZ();
		totalvnum++;
	}
	polyFile.getline(buf, size);

	for (int i = 0; i < inNum; i++)
	{
		polyFile.getline(buf, size);
		for (int j = 0; j < inVNum[i] - 1; j++)
		{
			polyFile.getline(buf, size);
			sscanf_s(buf, "%lf %lf", &pInVList[i][j].lo, &pInVList[i][j].la);
			pllv[totalvnum] = pInVList[i][j];
			pv[totalvnum] = pInVList[i][j].LL2XYZ();
			totalvnum++;
		}
		polyFile.getline(buf, size);
	}
	polyFile.close();

	this->min.la = this->min.lo = FLT_MAX;
	this->max.la = this->max.lo = -FLT_MAX;
	for (int i = 0; i < outVNum - 1; i++)
	{
		if (pOutVList[i].la < this->min.la) this->min.la = pOutVList[i].la;
		if (pOutVList[i].lo < this->min.lo) this->min.lo = pOutVList[i].lo;
		if (pOutVList[i].la > this->max.la) this->max.la = pOutVList[i].la;
		if (pOutVList[i].lo > this->max.lo) this->max.lo = pOutVList[i].lo;
	}
	for (int i = 0; i < inNum; i++)
	{
		for (int j = 0; j < inVNum[i]; i++)
		{
			if (pInVList[i][j].la < this->min.la) this->min.la = pInVList[i][j].la;
			if (pInVList[i][j].lo < this->min.lo) this->min.lo = pInVList[i][j].lo;
			if (pInVList[i][j].la > this->max.la) this->max.la = pInVList[i][j].la;
			if (pInVList[i][j].lo > this->max.lo) this->max.lo = pInVList[i][j].lo;
		}
	}

	bool* changeFlag = new bool[inNum + 1];
	for (int i = 0; i < inNum + 1; i++)
		changeFlag[i] = false;
	changeFlag[0] = this->isPolyClockwise(pOutVList, outVNum - 1);
	for (int i = 0; i < inNum; i++)
	{
		changeFlag[i + 1] = !(this->isPolyClockwise(pInVList[i], inVNum[i] - 1));
	}

	int index = 0;
	if (!changeFlag[0])
	{
		for (int i = 0; i < outVNum - 1; i++)
		{
			pe[index * 2] = index;
			pe[index * 2 + 1] = index + 1;
			index++;
		}
		pe[index * 2 - 1] = 0;
	}
	else
	{
		for (int i = outVNum - 2; i >= 0; i--)
		{
			pe[index * 2] = i;
			pe[index * 2 + 1] = i - 1;
			index++;
		}
		pe[index * 2 - 1] = outVNum - 2;
	}

	for (int i = 0; i < inNum; i++)
	{
		if (!changeFlag[i + 1])
		{
			int startindex = index;
			for (int j = 0; j < inVNum[i] - 1; j++)
			{
				pe[index * 2] = index;
				pe[index * 2 + 1] = index + 1;
				index++;
			}
			pe[index * 2 - 1] = startindex;
		}
		else
		{
			int startindex = index;
			for (int j = inVNum[i] - 2; j >= 0; j--)
			{
				pe[index * 2] = startindex + j;
				pe[index * 2 + 1] = startindex + j - 1;
				index++;
			}
			pe[index * 2 - 1] = startindex + inVNum[i] - 2;
		}
	}

	cudaMalloc((void**)&this->dev_pv, sizeof(devdouble3) * this->pv_num);
	cudaMemcpy(this->dev_pv, pv, sizeof(devdouble3) * this->pv_num, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&this->dev_pe, sizeof(int) * this->pv_num * 2);
	cudaMemcpy(this->dev_pe, pe, sizeof(int) * this->pv_num * 2, cudaMemcpyHostToDevice);

	delete[] pOutVList;
	for (int i = 0; i < inNum; i++)
		delete[] pInVList[i];
	delete pInVList;
	delete[] buf;

	return 0;
}

bool GPISPPolygon::isPolyClockwise(LLVertex* vlist, int vnum)
{
	double maxvalue = -FLT_MAX;
	int maxindex = 0;
	for (int i = 0; i < vnum; i++)
	{
		if (vlist[i].la > maxvalue)
		{
			maxvalue = vlist[i].la;
			maxindex = i;
		}
	}

	mydouble3 a[2], b;
	if (maxindex == 0)
		a[0] = vlist[vnum - 1].LL2XYZ();
	else
		a[0] = vlist[maxindex - 1].LL2XYZ();
	a[1] = vlist[maxindex].LL2XYZ();
	if (maxindex == vnum - 1)
		b = vlist[0].LL2XYZ();
	else
		b = vlist[maxindex + 1].LL2XYZ();

	mydouble3 e, n;
	e = a[1] - a[0];
	n = a[0].cross(e);
	n.normalize();
	double temp = n.dot(b);
	if (temp > 0) //counterclockwise
		return false;
	else //clockwise
		return true;
}

void GPISPPolygon::clear()
{
	if (this->pv != NULL)
		delete[] this->pv;
	if (this->pllv != NULL)
		delete[] this->pllv;
	if (this->pe != NULL)
		delete[] this->pe;
	if (this->dev_pv != NULL)
		cudaFree(this->dev_pv);
	if (this->dev_pe != NULL)
		cudaFree(this->dev_pe);
	memset(this,0,sizeof(*this));
}
