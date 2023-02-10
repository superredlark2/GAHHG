#include <device_launch_parameters.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <fstream>
#include <iostream>
#include <strstream>
#include <string.h>
#include <cstring>
#include <io.h>
#include <fcntl.h> 
#include "GPISPAHHG.cuh"

using namespace std;

int main(int argc, char* argv[])
{
	int i = 1;
	GPISP gpisp;
	gpisp.pParam = new GPISPParam();

	while (argv[i])
	{
		if (strcmp(argv[i], "/polyfile") == 0)
		{
			i++;
			strcpy_s(gpisp.pParam->polyFileName, argv[i]);
		}
		else if (strcmp(argv[i], "/reportfile") == 0)
		{
			i++;
			strcpy_s(gpisp.pParam->reportFileName, argv[i]);
		}
		else if (strcmp(argv[i], "/resultfile") == 0)
		{
			i++;
			strcpy_s(gpisp.pParam->resultFileName, argv[i]);
		}
		else if (strcmp(argv[i], "/tpfile") == 0)
		{
			i++;
			strcpy_s(gpisp.pParam->tpFileName, argv[i]);
		}
		i++;
	}
	
	int ret;
	gpisp.pPolygon = new GPISPPolygon;
	ret = gpisp.pPolygon->loadPolyFromTXT(gpisp.pParam->polyFileName);
	if (ret == -1)
		return -1;
		
	gpisp.pTP = new GPISPTP;
	ret = gpisp.pTP->loadTPFromCSV(gpisp.pParam->tpFileName);
	if (ret == -1)
		return -1;
	
	gpisp.pStat = new GPISPStat;
	gpisp.pAHHG = new GAHHG(gpisp.pPolygon, gpisp.pParam, gpisp.pStat);
	gpisp.createAHHGInit();
	gpisp.createAHHG(gpisp.pParam->maxDepth);
	gpisp.classifyAHHG();
	
	gpisp.test();
	
	gpisp.getStat();
	gpisp.exportStat2TXT(gpisp.pParam->reportFileName, 1);
	gpisp.pTP->exportResult2CSV(gpisp.pParam->resultFileName);

	return 0;
}


