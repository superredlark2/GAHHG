//GPISPAHHG.cu

#include <device_launch_parameters.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <float.h>
#include <fstream>
#include <iostream>
#include <strstream>
#include <string.h>
#include <cstring>
#include <io.h>
#include "GPISPAHHG.cuh"
#include "scan/scan.cuh"

using namespace std;

__device__ devdouble3 device_getENormal(devdouble3* dev_hcv, int v0, int v1)
{
	devdouble3 e0 = dev_hcv[v0];
	devdouble3 e1 = dev_hcv[v1] - dev_hcv[v0];
	devdouble3 nor = e0.cross(e1);
	nor.normalize();
	return nor;
}

__device__ int device_getUpEdgeForSubcell(
	AHHGCell* dev_c,
	AHHGBorder* dev_b,
	int cellid)
{
	int upcellid = dev_c[cellid].upcell;
	for (int i = 0; i < dev_c[upcellid].bnum; i++)
	{
		if (dev_b[dev_c[upcellid].bstart + i].subcell == cellid)
			return dev_c[upcellid].bstart + i;
	}
	return -1;
}


__device__ void device_setBorderPartner(
	AHHGCell* dev_c,
	AHHGBorder* dev_b,
	int cellindex,
	int edgeindex)
{
	for (int i = 0; i < dev_c[cellindex].bnum; i++)
	{
		if (dev_b[edgeindex].v[0] == dev_b[dev_c[cellindex].bstart + i].v[1] &&
			dev_b[edgeindex].v[1] == dev_b[dev_c[cellindex].bstart + i].v[0] &&
			edgeindex < dev_c[cellindex].bstart + i)
		{
			dev_b[edgeindex].partner = dev_c[cellindex].bstart + i;
			dev_b[dev_c[cellindex].bstart + i].partner = edgeindex;
		}
	}
}


__device__ int device_getCPProp(
	devdouble3* dev_pv,
	int* dev_pe,
	AHHGCell* dev_c,
	int* dev_refs,
	int cell,
	devdouble3 center)
{
	double mindist = FLT_MAX;
	int ne = dev_refs[dev_c[cell].rstart];
	devdouble3 ne_normal;
	int nv_index = -1;

	for (int i = 0; i < dev_c[cell].rnum; i++)
	{
		devdouble3 v[2], p[2], e[2], n[2], inter;
		int se = dev_refs[dev_c[cell].rstart + i];
		double curdist = FLT_MAX;

		v[0] = dev_pv[dev_pe[se * 2]];
		v[1] = dev_pv[dev_pe[se * 2 + 1]];
		e[0] = v[0];
		e[1] = v[1] - v[0];
		n[0] = e[0].cross(e[1]);
		n[0].normalize();
		e[0] = center;
		e[1] = n[0];
		n[1] = e[0].cross(e[1]);
		n[1].normalize();
		inter = n[0].cross(n[1]);
		inter.normalize();
		double temp = inter.dot(v[0]);
		if (LESS_ZERO(temp))
			inter = -inter;

		devdouble3 temp0 = v[0] - inter;
		temp0.normalize();
		devdouble3 temp1 = v[1] - inter;
		temp1.normalize();
		double temp3 = temp0.dot(temp1);
		int cur_nv_index = -1;

		if (LESS_ZERO(temp3))
		{
			curdist = (inter - center).squarelen();
			cur_nv_index = -1;
		}
		else
		{
			curdist = (center - v[0]).squarelen();
			double vdist = (center - v[1]).squarelen();
			if (vdist < curdist)
			{
				curdist = vdist;
				cur_nv_index = dev_pe[se * 2 + 1];
			}
			else
			{
				cur_nv_index = dev_pe[se * 2];
			}
		}

		if (nv_index != -1 && nv_index == cur_nv_index)
		{
			double temp1 = ne_normal.dot(center);
			double temp2 = n[0].dot(center);
			if (LESS_EQU_ZERO(temp1) && GRE_EQU_ZERO(temp2) ||
				GRE_EQU_ZERO(temp1) && LESS_EQU_ZERO(temp2))
			{
				int another_vertex_index;
				if (dev_pe[ne * 2] == nv_index)
					another_vertex_index = dev_pe[ne * 2 + 1];
				else
					another_vertex_index = dev_pe[ne * 2];
				devdouble3 a0 = dev_pv[another_vertex_index] - dev_pv[nv_index];
				if (dev_pe[se * 2] == nv_index)
					another_vertex_index = dev_pe[se * 2 + 1];
				else
					another_vertex_index = dev_pe[se * 2];
				devdouble3 a1 = dev_pv[another_vertex_index] - dev_pv[nv_index];
				devdouble3 b = center - dev_pv[nv_index];
				a0.normalize();
				a1.normalize();
				temp1 = a0.dot(b);
				temp2 = a1.dot(b);
				if (temp2 > temp1)
				{
					mindist = curdist;
					ne = se;
					ne_normal = n[0];
					nv_index = cur_nv_index;
					continue;
				}
			}
		}

		if (curdist < mindist)
		{
			mindist = curdist;
			ne = se;
			ne_normal = n[0];
			nv_index = cur_nv_index;
		}
	}

	double temp = center.dot(ne_normal);
	if (GRE_ZERO(temp))
		return CELL_IN;
	else
		return CELL_OUT;
}


__device__ bool dev_isInterSeg(devdouble3 v0, devdouble3 v1, devdouble3 p0, devdouble3 p1)
{
	devdouble3 v[2];
	v[0] = v0;
	v[1] = v1;

	devdouble3 p[2];
	p[0] = p0;
	p[1] = p1;

	devdouble3 n[2], e[2];
	double ret[2];
	e[0] = v[0];
	e[1] = v[1] - v[0];
	n[0] = e[0].cross(e[1]);
	n[0].normalize();
	ret[0] = n[0].dot(p[0]);
	ret[1] = n[0].dot(p[1]);
	if (LESS_ZERO(ret[0]) && LESS_ZERO(ret[1]) || GRE_ZERO(ret[0]) && GRE_ZERO(ret[1]))
		return false;

	e[0] = p[0];
	e[1] = p[1] - p[0];
	n[1] = e[0].cross(e[1]);
	n[1].normalize();
	ret[0] = n[1].dot(v[0]);
	ret[1] = n[1].dot(v[1]);
	if (LESS_ZERO(ret[0]) && LESS_ZERO(ret[1]) || GRE_ZERO(ret[0]) && GRE_ZERO(ret[1]))
		return false;

	return true;
}


__device__ int device_isPointInsideCell(
	devdouble3 p,
	int cid,
	AHHGCell* dev_c,
	AHHGBorder* dev_b)
{
	int i = 0;
	for (i = 0; i < dev_c[cid].bnum; i++)
	{
		double ret = dev_b[dev_c[cid].bstart + i].n[0] * p.x +
			dev_b[dev_c[cid].bstart + i].n[1] * p.y +
			dev_b[dev_c[cid].bstart + i].n[2] * p.z;
		if (ret < -POSITION_EPSILON)
			return -1;
	}
	return cid;
}

__device__ int device_isEdgeBorderInter(
	devdouble3 p0,
	devdouble3 p1,
	devdouble3 v0,
	devdouble3 v1,
	int& reserved)
{
	devdouble3 n[2], e[4];
	double ret[4];
	e[0] = v0;
	e[1] = v1 - v0;
	n[0] = e[0].cross(e[1]);
	n[0].normalize();
	ret[0] = n[0].dot(p0);
	ret[1] = n[0].dot(p1);

	e[2] = p0;
	e[3] = p1 - p0;
	n[1] = e[2].cross(e[3]);
	n[1].normalize();
	ret[2] = n[1].dot(v0);
	ret[3] = n[1].dot(v1);

	if (EQU_ZERO(ret[0]) && EQU_ZERO(ret[1]))
		return INTER_OVERLAP;

	if (!(LESS_ZERO(ret[2]) && LESS_ZERO(ret[3]) || GRE_ZERO(ret[2]) && GRE_ZERO(ret[3])) && EQU_ZERO(ret[0]))
	{
		reserved = 0;
		return INTER_T1;
	}

	if (!(LESS_ZERO(ret[2]) && LESS_ZERO(ret[3]) || GRE_ZERO(ret[2]) && GRE_ZERO(ret[3])) && EQU_ZERO(ret[1]))
	{
		reserved = 1;
		return INTER_T1;
	}

	if (EQU_ZERO(ret[2]) && EQU_ZERO(ret[3]))
		return INTER_OVERLAP;

	if (!(LESS_ZERO(ret[0]) && LESS_ZERO(ret[1]) || GRE_ZERO(ret[0]) && GRE_ZERO(ret[1])) && EQU_ZERO(ret[2]))
	{
		reserved = 0;
		return INTER_T2;
	}

	if (!(LESS_ZERO(ret[0]) && LESS_ZERO(ret[1]) || GRE_ZERO(ret[0]) && GRE_ZERO(ret[1])) && EQU_ZERO(ret[3]))
	{
		reserved = 1;
		return INTER_T2;
	}

	if (LESS_ZERO(ret[0]) && LESS_ZERO(ret[1]) || GRE_ZERO(ret[0]) && GRE_ZERO(ret[1]) ||
		LESS_ZERO(ret[2]) && LESS_ZERO(ret[3]) || GRE_ZERO(ret[2]) && GRE_ZERO(ret[3]))
		return INTER_NO;

	return INTER_X;
}

__device__ int device_getSegExitBorder(
	devdouble3 p0,
	devdouble3 p1,
	devdouble3* dev_v,
	AHHGCell* dev_c,
	int cid,
	AHHGBorder* dev_b,
	int in,
	int& intertype,
	int& info)
{
	int i = 0;
	for (i = dev_c[cid].bstart; i < dev_c[cid].bstart + dev_c[cid].bnum; i++)
	{
		if (i == in)
			continue;

		int temp;
		int inter = device_isEdgeBorderInter(p0, p1, dev_v[dev_b[i].v[0]], dev_v[dev_b[i].v[1]], temp);

		if (inter == INTER_X)
		{
			intertype = INTER_X;
			return i;
		}

		if (inter == INTER_T1)
		{
			if (temp == 0)
			{
				intertype = INTER_T1;
				info = 0;
				return i;
			}
			else
			{
				intertype = INTER_T1;
				info = 1;
				return -1;
			}
		}

		if (inter == INTER_T2)
		{
			if (in != -1)
			{
				if (dev_b[i].v[temp] == dev_b[in].v[0] || dev_b[i].v[temp] == dev_b[in].v[1])
					continue;
			}
			intertype = INTER_T2;
			info = temp;
			return i;
		}

		if (inter == INTER_OVERLAP)
		{
			intertype = INTER_OVERLAP;
			return i;
		}
	}

	return -1;
}


__device__ bool device_isEEInter(
	devdouble3 p0,
	devdouble3 p1,
	devdouble3 v0,
	devdouble3 v1)
{
	devdouble3 n[2], e[4];
	double ret[4];
	e[0] = v0;
	e[1] = v1 - v0;
	n[0] = e[0].cross(e[1]);
	n[0].normalize();
	ret[0] = n[0].dot(p0);
	ret[1] = n[0].dot(p1);

	e[2] = p0;
	e[3] = p1 - p0;
	n[1] = e[2].cross(e[3]);
	n[1].normalize();
	ret[2] = n[1].dot(v0);
	ret[3] = n[1].dot(v1);

	if (LESS_ZERO(ret[0]) && LESS_ZERO(ret[1]) || GRE_ZERO(ret[0]) && GRE_ZERO(ret[1]) ||
		LESS_ZERO(ret[2]) && LESS_ZERO(ret[3]) || GRE_ZERO(ret[2]) && GRE_ZERO(ret[3]))
		return false;
	else
		return true;
}


__global__ void kernel_DivGrid_calcCenterOffsets(
	int* dev_cell_buf0,
	int* dev_cell_buf1,
	AHHGCell* dev_c,
	int cell_start,
	int cell_num,
	int* dev_se_buf0)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < cell_num)
	{
		if (dev_c[cell_start + tid].rnum == 0)
		{
			dev_cell_buf0[tid] = 0;
			dev_cell_buf1[tid] = 0;
		}
		else
		{
			dev_cell_buf0[tid] = dev_c[cell_start + tid].bnum;
			dev_cell_buf1[tid] = 1;
			atomicAdd(&dev_se_buf0[0], 1);
		}
	}
}

__global__ void kernel_DivGrid_createCenter(
	devdouble3* dev_v,
	int subgrid_hcv_start,
	AHHGBorder* dev_b,
	int subgrid_hce_start,
	AHHGCell* dev_c,
	int cstart,
	int cnum,
	int* dev_cell_buf0,
	int* dev_cell_buf1)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < cnum)
	{
		if (dev_c[cstart + tid].rnum != 0)
		{
			devdouble3 center;
			if (dev_c[cstart + tid].bnum == 4)
				center = (dev_v[dev_b[dev_c[cstart + tid].bstart].v[0]] +
					dev_v[dev_b[dev_c[cstart + tid].bstart].v[1]]) * 0.5;
			else
			{
				for (int i = 0; i < dev_c[cstart + tid].bnum; i++)
					center = center + dev_v[dev_b[dev_c[cstart + tid].bstart + i].v[0]];
				center = center / dev_c[cstart + tid].bnum;
			}
			center.normalize();

			for (int i = 0; i < dev_c[cstart + tid].bnum; i++)
			{
				dev_v[subgrid_hcv_start + dev_cell_buf0[tid] + i] = (dev_v[dev_b[dev_c[cstart + tid].bstart + i].v[0]] + center) * 0.5;
				dev_v[subgrid_hcv_start + dev_cell_buf0[tid] + i].normalize();
			}

			for (int i = 0; i < dev_c[cstart + tid].bnum; i++)
			{
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].id = i;
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].parent = cstart + cnum + dev_cell_buf1[tid];
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].partner = -1;
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].subcell = -1;
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].v[0] = subgrid_hcv_start + dev_cell_buf0[tid] + i;
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].v[1] = subgrid_hcv_start + dev_cell_buf0[tid] + (i + 1) % dev_c[cstart + tid].bnum;
				devdouble3 n = device_getENormal(
					dev_v,
					dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].v[0],
					dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].v[1]);
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].n[0] = n.x;
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].n[1] = n.y;
				dev_b[subgrid_hce_start + dev_cell_buf0[tid] + i].n[2] = n.z;
			}

			dev_c[cstart + cnum + dev_cell_buf1[tid]].bstart = subgrid_hce_start + dev_cell_buf0[tid];
			dev_c[cstart + cnum + dev_cell_buf1[tid]].bnum = dev_c[cstart + tid].bnum;
			dev_c[cstart + cnum + dev_cell_buf1[tid]].rstart = -1;
			dev_c[cstart + cnum + dev_cell_buf1[tid]].rnum = -1;
			dev_c[cstart + cnum + dev_cell_buf1[tid]].upcell = cstart + tid;
			dev_c[cstart + cnum + dev_cell_buf1[tid]].subcell = -1;
			dev_c[cstart + cnum + dev_cell_buf1[tid]].type = CELL_UNDEFINED;
			dev_c[cstart + tid].subcell = cstart + cnum + dev_cell_buf1[tid];

			dev_cell_buf0[tid] += dev_c[cstart + tid].bnum;
			dev_cell_buf1[tid] += 1;
		}
	}
}


__global__ void kernel_DivGrid_calcBorderSCOffsets(
	AHHGBorder* dev_b,
	int bstart,
	int bnum,
	AHHGCell* dev_c,
	int* dev_cell_buf0,
	int* dev_cell_buf1)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < bnum)
	{
		dev_cell_buf0[tid] = 0;
		dev_cell_buf1[tid] = 0;

		if (1)
		{
			AHHGCell* upcell = &dev_c[dev_c[dev_b[bstart + tid].parent].upcell];
			AHHGBorder* upedge = &dev_b[upcell->bstart + dev_b[bstart + tid].id];

			if (upcell->bnum == 4 && dev_b[bstart + tid].id == 0)
				return;

			if (upedge->partner != -1 && dev_c[dev_b[upedge->partner].parent].subcell != -1 &&
				upedge->parent < dev_b[upedge->partner].parent)
			{
				dev_cell_buf0[tid] += 6;
				dev_cell_buf1[tid] += 1;
			}

			if (upedge->partner == -1 || upedge->partner != -1 && dev_c[dev_b[upedge->partner].parent].subcell == -1)
			{
				dev_cell_buf0[tid] += 4;
				dev_cell_buf1[tid] += 1;
			}
		}
	}
}


__global__ void kernel_DivGrid_createBorderSC(
	devdouble3* dev_v,
	AHHGBorder* dev_b,
	int bstart,
	int bnum,
	AHHGCell* dev_c,
	int cstart,
	int cnum,
	int* dev_cell_buf0,
	int* dev_cell_buf1)
{  
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < bnum)
	{
		AHHGCell* upcell = &dev_c[dev_c[dev_b[bstart + tid].parent].upcell];
		AHHGBorder* upedge = &dev_b[upcell->bstart + dev_b[bstart + tid].id];
		int eoffset = bstart + bnum + dev_cell_buf0[tid];
		int coffset = cstart + cnum + dev_cell_buf1[tid];

		if (upcell->bnum == 4 && dev_b[bstart + tid].id == 0)
			return;

		if (upedge->partner != -1 && dev_c[dev_b[upedge->partner].parent].subcell != -1 &&
			upedge->parent < dev_b[upedge->partner].parent)
		{
			AHHGCell* partner_subcell = &dev_c[dev_c[dev_b[upedge->partner].parent].subcell];
			AHHGBorder* partner_subedge = &dev_b[partner_subcell->bstart + dev_b[upedge->partner].id];

			dev_b[eoffset].v[0] = upedge->v[0];
			dev_b[eoffset + 1].v[0] = partner_subedge->v[1];
			dev_b[eoffset + 2].v[0] = partner_subedge->v[0];
			dev_b[eoffset + 3].v[0] = upedge->v[1];
			dev_b[eoffset + 4].v[0] = dev_b[bstart + tid].v[1];
			dev_b[eoffset + 5].v[0] = dev_b[bstart + tid].v[0];

			for (int i = 0; i < 6; i++)
			{
				dev_b[eoffset + i].id = i;
				dev_b[eoffset + i].v[1] = dev_b[eoffset + (i + 1) % 6].v[0];
				dev_b[eoffset + i].parent = coffset;
				dev_b[eoffset + i].partner = -1;
				devdouble3 n = device_getENormal(dev_v, dev_b[eoffset + i].v[0], dev_b[eoffset + i].v[1]);
				dev_b[eoffset + i].n[0] = n.x;
				dev_b[eoffset + i].n[1] = n.y;
				dev_b[eoffset + i].n[2] = n.z;
				dev_b[eoffset + i].subcell = -1;
			}

			dev_c[coffset].bnum = 6;
			dev_c[coffset].bstart = eoffset;
			dev_c[coffset].rstart = -1;
			dev_c[coffset].rnum = -1;
			dev_c[coffset].subcell = -1;
			dev_c[coffset].upcell = upedge->parent;
			dev_c[coffset].type = CELL_UNDEFINED;
			upedge->subcell = coffset;
			dev_b[upedge->partner].subcell = coffset;

			dev_cell_buf0[tid] += 6;
			dev_cell_buf1[tid] ++;
		}

		if (upedge->partner == -1 || upedge->partner != -1 && dev_c[dev_b[upedge->partner].parent].subcell == -1)
		{
			dev_b[eoffset].v[0] = upedge->v[0];
			dev_b[eoffset + 1].v[0] = upedge->v[1];
			dev_b[eoffset + 2].v[0] = dev_b[bstart + tid].v[1];
			dev_b[eoffset + 3].v[0] = dev_b[bstart + tid].v[0];

			for (int i = 0; i < 4; i++)
			{
				dev_b[eoffset + i].id = i;
				dev_b[eoffset + i].v[1] = dev_b[eoffset + (i + 1) % 4].v[0];
				dev_b[eoffset + i].parent = coffset;
				dev_b[eoffset + i].partner = -1;
				devdouble3 n = device_getENormal(dev_v, dev_b[eoffset + i].v[0], dev_b[eoffset + i].v[1]);
				dev_b[eoffset + i].n[0] = n.x;
				dev_b[eoffset + i].n[1] = n.y;
				dev_b[eoffset + i].n[2] = n.z;
				dev_b[eoffset + i].subcell = -1;
			}

			dev_c[coffset].bnum = 4;
			dev_c[coffset].bstart = eoffset;
			dev_c[coffset].rstart = -1;
			dev_c[coffset].rnum = -1;
			dev_c[coffset].subcell = -1;
			dev_c[coffset].upcell = upedge->parent;
			dev_c[coffset].type = CELL_UNDEFINED;
			upedge->subcell = coffset;

			dev_cell_buf0[tid] += 4;
			dev_cell_buf1[tid] ++;
		}
	}
}


__global__ void kernel_DivGrid_setSCPartner(
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int cstart,
	int cnum)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < cnum)
	{
		AHHGCell* curcell = &dev_c[cstart + tid];
		int upedge = device_getUpEdgeForSubcell(dev_c, dev_b, cstart + tid);
		if (upedge == -1)
			return;
		int tempedge;

		if (curcell->bnum == 6)
		{
			int partneredge = dev_b[upedge].partner;
			int partnercell = dev_b[partneredge].parent;
			int u = dev_c[partnercell].bnum;

			tempedge = dev_c[partnercell].bstart + (dev_b[partneredge].id + 1) % u;
			if (dev_b[tempedge].subcell != -1)
				device_setBorderPartner(dev_c, dev_b, dev_b[tempedge].subcell, curcell->bstart);
			else
				dev_b[curcell->bstart].partner = -1;

			tempedge = dev_c[dev_c[partnercell].subcell].bstart + dev_b[partneredge].id;
			dev_b[curcell->bstart + 1].partner = tempedge;
			dev_b[tempedge].partner = curcell->bstart + 1;

			tempedge = dev_c[partnercell].bstart + (dev_b[partneredge].id + u - 1) % u;
			if (dev_b[tempedge].subcell != -1)
				device_setBorderPartner(dev_c, dev_b, dev_b[tempedge].subcell, curcell->bstart + 2);
			else
				dev_b[curcell->bstart + 2].partner = -1;

			tempedge = dev_c[curcell->upcell].bstart + (dev_b[upedge].id + 1) % dev_c[curcell->upcell].bnum;
			if (dev_b[tempedge].subcell != -1)
				device_setBorderPartner(dev_c, dev_b, dev_b[tempedge].subcell, curcell->bstart + 3);
			else
				dev_b[curcell->bstart + 3].partner = -1;

			tempedge = dev_c[dev_c[curcell->upcell].subcell].bstart + dev_b[upedge].id;
			dev_b[curcell->bstart + 4].partner = tempedge;
			dev_b[tempedge].partner = curcell->bstart + 4;

			tempedge = dev_c[curcell->upcell].bstart + (dev_b[upedge].id + dev_c[curcell->upcell].bnum - 1) % dev_c[curcell->upcell].bnum;
			if (dev_b[tempedge].subcell != -1)
				device_setBorderPartner(dev_c, dev_b, dev_b[tempedge].subcell, curcell->bstart + 5);
			else
				dev_b[curcell->bstart + 5].partner = -1;
		}
		else if (curcell->bnum == 4)
		{
			tempedge = dev_c[curcell->upcell].bstart + (dev_b[upedge].id + 1) % dev_c[curcell->upcell].bnum;
			if (dev_b[tempedge].subcell != -1)
				device_setBorderPartner(dev_c, dev_b, dev_b[tempedge].subcell, curcell->bstart + 1);
			else
				dev_b[curcell->bstart + 1].partner = -1;

			tempedge = dev_c[dev_c[curcell->upcell].subcell].bstart + dev_b[upedge].id;
			dev_b[curcell->bstart + 2].partner = tempedge;
			dev_b[tempedge].partner = curcell->bstart + 2;

			tempedge = dev_c[curcell->upcell].bstart + (dev_b[upedge].id + dev_c[curcell->upcell].bnum - 1) % dev_c[curcell->upcell].bnum;
			if (dev_b[tempedge].subcell != -1)
				device_setBorderPartner(dev_c, dev_b, dev_b[tempedge].subcell, curcell->bstart + 3);
			else
				dev_b[curcell->bstart + 3].partner = -1;
		}
	}
}


__global__ void kernel_calcRefNum(
	int* dev_se_buf,
	devdouble3* dev_pv,
	int* dev_pe,
	int pv_num,
	devdouble3* dev_v,
	AHHGCell* dev_c,
	int cnum,
	AHHGBorder* dev_b,
	int* dev_buf,
	int level)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < pv_num)
	{
		
		devdouble3 p0 = dev_pv[dev_pe[tid * 2]];
		devdouble3 p1 = dev_pv[dev_pe[tid * 2 + 1]];
		int curcell = -1;
		
		if (level == 0)
		{
			for (int i = 0; i < cnum; i++)
			{
				curcell = device_isPointInsideCell(p0, i, dev_c, dev_b);
				if (curcell != -1)
					break;
			}
		}
		else
		{
			curcell = device_isPointInsideCell(p0, dev_c[dev_buf[tid]].subcell, dev_c, dev_b);
			
			if (curcell == -1)
			{
				for (int i = 0; i < dev_c[dev_buf[tid]].bnum; i++)
				{
					curcell = dev_b[dev_c[dev_buf[tid]].bstart + i].subcell;
					if (curcell == -1)
						continue;
					curcell = device_isPointInsideCell(p0, curcell, dev_c, dev_b);
					if (curcell != -1)
						break;
				}
			}
			if (curcell == -1)
				return;
		}

		int in = -1, out = -1;

		while (1)
		{
			atomicAdd(&dev_se_buf[tid], 1);

			int intertype, info;
			out = device_getSegExitBorder(p0, p1, dev_v, dev_c, curcell, dev_b, in, intertype, info);

			if (out == -1)
				break;

			if (intertype == INTER_X || intertype == INTER_T1)
			{
				in = dev_b[out].partner;
				curcell = dev_b[in].parent;
			}
			else if (intertype == INTER_T2)
			{
				devdouble3 curpoint = dev_v[dev_b[out].v[info]];
				devdouble3 e = p1 - p0;
				e.normalize();
				curpoint = curpoint + e * TRAVEL_EPSILON;

				int otheredge;
				int u = dev_c[curcell].bnum;

				if (dev_b[out].partner != -1)
				{
					int ret = device_isPointInsideCell(curpoint, dev_b[dev_b[out].partner].parent, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						in = dev_b[out].partner;
						continue;
					}
				}

				if (info == 0)
					otheredge = dev_c[curcell].bstart + (dev_b[out].id + u - 1) % u;
				else
					otheredge = dev_c[curcell].bstart + (dev_b[out].id + 1) % u;

				if (dev_b[otheredge].partner != -1)
				{
					int ret = device_isPointInsideCell(curpoint, dev_b[dev_b[otheredge].partner].parent, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						in = dev_b[otheredge].partner;
						continue;
					}
					else
						return;
				}
				else
				{
					return;
				}
			}
			else if (intertype == INTER_OVERLAP)
			{
				devdouble3 e[2];
				e[0] = p1 - p0;
				e[0].normalize();
				e[1] = dev_v[dev_b[out].v[1]] - dev_v[dev_b[out].v[0]];
				e[1].normalize();
				double temp = e[0].dot(e[1]);
				int u = dev_c[curcell].bnum;
				if (temp > 0)
				{
					int targetedge = dev_c[curcell].bstart + (dev_b[out].id + 1) % u;
					curcell = dev_b[dev_b[targetedge].partner].parent;
					in = dev_b[targetedge].partner;
				}
				else
				{
					int targetedge = dev_c[curcell].bstart + (dev_b[out].id + u - 1) % u;
					curcell = dev_b[dev_b[targetedge].partner].parent;
					in = dev_b[targetedge].partner;
				}
			}
		}
	}
}

__global__ void kernel_writeRef(
	int* dev_refs,
	int* dev_se_buf,
	devdouble3* dev_pv,
	int* dev_pe,
	int pv_num,
	devdouble3* dev_v,
	AHHGCell* dev_c,
	int c_num,
	AHHGBorder* dev_b,
	int* dev_buf,
	int level)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < pv_num)
	{
		devdouble3 p0 = dev_pv[dev_pe[tid * 2]];
		devdouble3 p1 = dev_pv[dev_pe[tid * 2 + 1]];
		int curcell = -1;

		if (level == 0)
		{
			for (int i = 0; i < c_num; i++)
			{
				curcell = device_isPointInsideCell(p0, i, dev_c, dev_b);
				if (curcell != -1)
					break;
			}
		}
		else
		{
			curcell = device_isPointInsideCell(p0, dev_c[dev_buf[tid]].subcell, dev_c, dev_b);
			if (curcell == -1)
			{
				for (int i = 0; i < dev_c[dev_buf[tid]].bnum; i++)
				{
					curcell = dev_b[dev_c[dev_buf[tid]].bstart + i].subcell;
					if (curcell == -1)
						continue;
					curcell = device_isPointInsideCell(p0, curcell, dev_c, dev_b);
					if (curcell != -1)
						break;
				}
			}
			if (curcell == -1)
				return;
			else
				dev_buf[tid] = curcell;
		}

		int in = -1, out = -1;

		while (1)
		{
			int pos = atomicAdd(&dev_se_buf[tid], 1);
			dev_refs[pos] = curcell;

			int intertype, info;
			out = device_getSegExitBorder(p0, p1, dev_v, dev_c, curcell, dev_b, in, intertype, info);

			if (out == -1)
				break;

			if (intertype == INTER_X || intertype == INTER_T1)
			{
				in = dev_b[out].partner;
				curcell = dev_b[in].parent;
			}
			else if (intertype == INTER_T2)
			{
				devdouble3 curpoint = dev_v[dev_b[out].v[info]];
				devdouble3 e = p1 - p0;
				e.normalize();
				curpoint = curpoint + e * TRAVEL_EPSILON;

				int otheredge;
				int u = dev_c[curcell].bnum;

				if (dev_b[out].partner != -1)
				{
					int ret = device_isPointInsideCell(curpoint, dev_b[dev_b[out].partner].parent, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						in = dev_b[out].partner;
						continue;
					}
				}

				if (info == 0)
					otheredge = dev_c[curcell].bstart + (dev_b[out].id + u - 1) % u;
				else
					otheredge = dev_c[curcell].bstart + (dev_b[out].id + 1) % u;

				if (dev_b[otheredge].partner != -1)
				{
					int ret = device_isPointInsideCell(curpoint, dev_b[dev_b[otheredge].partner].parent, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						in = dev_b[otheredge].partner;
						continue;
					}
					else
						return;
				}
				else
				{
					return;
				}
			}
			else if (intertype == INTER_OVERLAP)
			{
				devdouble3 e[2];
				e[0] = p1 - p0;
				e[0].normalize();
				e[1] = dev_v[dev_b[out].v[1]] - dev_v[dev_b[out].v[0]];
				e[1].normalize();
				double temp = e[0].dot(e[1]);
				int u = dev_c[curcell].bnum;
				if (temp > 0)
				{
					int targetedge = dev_c[curcell].bstart + (dev_b[out].id + 1) % u;
					curcell = dev_b[dev_b[targetedge].partner].parent;
					in = dev_b[targetedge].partner;
				}
				else
				{
					int targetedge = dev_c[curcell].bstart + (dev_b[out].id + u - 1) % u;
					curcell = dev_b[dev_b[targetedge].partner].parent;
					in = dev_b[targetedge].partner;
				}
			}
		}
	}
}


__global__ void kernel_su2pre(
	int* count,
	int* dest,
	int len)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < len)
	{
		dest[i] = dest[i] - count[i];
	}
}


__global__ void kernel_EC2CS_setEFlag(
	int* ref_buf,
	int* sedge_buf,
	int sedge_count)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < sedge_count && i != 0)
	{
		ref_buf[sedge_buf[i]] = 1;
	}
}


__global__ void kernel_EC2CS_adjustEFlag(
	int* source,
	int* dest,
	int len)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < len)
	{
		dest[i] = dest[i] + source[i];
	}
}


__global__ void kernel_EC2CS_countEinC(
	int* ref_buf,
	int ref_count,
	int* cell_buf,
	int cell_start)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < ref_count)
	{
		atomicAdd(&cell_buf[ref_buf[i] - cell_start], 1);
	}
}


__global__ void kernel_EC2CS_rearrangeCERef(
	int* ref_source0,
	int* ref_source1,
	int* ref_dest,
	int ref_count,
	int* cell_buf,
	int cell_start)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < ref_count)
	{
		int sedgeindex = ref_source1[i];
		int cellindex = ref_source0[i] - cell_start;
		int pos = atomicAdd(&cell_buf[cellindex], 1);
		ref_dest[pos] = sedgeindex;
	}
}


__global__ void kernel_EC2CS_setCEInfo(
	int* ref_count,
	int* ref_offset,
	AHHGCell* dev_c,
	int cell_start,
	int cell_num)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < cell_num)
	{
		dev_c[cell_start + tid].rnum = ref_count[tid];
		dev_c[cell_start + tid].rstart = ref_offset[tid] - ref_count[tid];
	}
}


__global__ void kernel_Classify_setBottomNECell(
	devdouble3* dev_pv,
	int* dev_pe,
	devdouble3* dev_v,
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int hc_start,
	int hc_num,
	int* dev_refs)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < hc_num)
	{
		AHHGCell* necell = &dev_c[hc_start + tid];
		if (necell->rnum != 0)
		{
			devdouble3 center;
			for (int i = 0; i < necell->bnum; i++)
				center = center + dev_v[dev_b[necell->bstart + i].v[0]];
			center = center / necell->bnum;
			center.normalize();

			necell->type = device_getCPProp(dev_pv, dev_pe, dev_c, dev_refs, hc_start + tid, center);
		}
	}
}


__global__ void kernel_Classify_setBottomDECell(
	devdouble3* dev_pv,
	int* dev_pe,
	devdouble3* dev_v,
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int* dev_refs,
	int hc_start,
	int hc_num)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < hc_num)
	{
		if (dev_c[hc_start + tid].type == CELL_UNDEFINED)
		{
			for (int i = 0; i < dev_c[hc_start + tid].bnum; i++)
			{
				int partneredge = dev_b[dev_c[hc_start + tid].bstart + i].partner;
				if (partneredge != -1 && dev_c[dev_b[partneredge].parent].type != CELL_UNDEFINED)
				{
					devdouble3 v0 = dev_v[dev_b[partneredge].v[0]];
					devdouble3 v1 = dev_v[dev_b[partneredge].v[1]];
					devdouble3 midpoint = (v0 + v1) * 0.5;
					midpoint.normalize();

					int necell = dev_b[partneredge].parent;
					devdouble3 center;
					for (int i = 0; i < dev_c[necell].bnum; i++)
						center = center + dev_v[dev_b[dev_c[necell].bstart + i].v[0]];
					center = center / dev_c[necell].bnum;
					center.normalize();

					int inter_count = 0;
					for (int i = 0; i < dev_c[necell].rnum; i++)
					{
						int se = dev_refs[dev_c[necell].rstart + i];
						devdouble3 u0 = dev_pv[dev_pe[se * 2]];
						devdouble3 u1 = dev_pv[dev_pe[se * 2 + 1]];
						bool ret = dev_isInterSeg(u0, u1, center, midpoint);
						if (ret)
							inter_count++;
					}
					if (inter_count % 2 == 0)
						dev_c[hc_start + tid].type = dev_c[necell].type;
					else
						dev_c[hc_start + tid].type = 3 - dev_c[necell].type;
				}
			}
		}
	}
}


__global__ void kernel_DivGridFast_zeroRefCount(
	AHHGCell* dev_c,
	int cstart,
	int cnum,
	AHHGBorder* dev_b)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < cnum)
	{
		dev_c[cstart + tid].rnum = 0;
		if (cstart == 0)
		{
			dev_c[cstart + tid].rstart = -1;
			dev_c[cstart + tid].subcell = -1;
			dev_c[cstart + tid].upcell = -1;
			dev_c[cstart + tid].type = CELL_UNDEFINED;
			for (int i = 0; i < dev_c[cstart + tid].bnum; i++)
			{
				int bindex = dev_c[cstart + tid].bstart;
				dev_b[bindex + i].subcell = -1;
			}
		}
	}
}


__global__ void kernel_DivGridFast_markNECell(
	devdouble3* dev_pv,
	int* dev_pe,
	int pv_num,
	devdouble3* dev_v,
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int cstart,
	int cnum,
	int level,
	int* dev_buf1,
	int* dev_buf2)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < pv_num)
	{
		devdouble3 p0 = dev_pv[dev_pe[tid * 2]];
		devdouble3 p1 = dev_pv[dev_pe[tid * 2 + 1]];
		int curcell = -1;

		dev_buf1[tid] = dev_buf2[tid];

		if (level == 0)
		{
			for (int i = 0; i < 42; i++)
			{
				curcell = device_isPointInsideCell(p0, i, dev_c, dev_b);
				if (curcell != -1)
					break;
			}
			if (curcell == -1)
				return;
			dev_buf2[tid] = curcell;
		}
		else
		{
			curcell = device_isPointInsideCell(p0, dev_c[dev_buf2[tid]].subcell, dev_c, dev_b);
			if (curcell == -1)
			{
				for (int i = 0; i < dev_c[dev_buf2[tid]].bnum; i++)
				{
					curcell = dev_b[dev_c[dev_buf2[tid]].bstart + i].subcell;
					if (curcell == -1)
						continue;
					curcell = device_isPointInsideCell(p0, curcell, dev_c, dev_b);
					if (curcell != -1)
						break;
				}
			}
			if (curcell == -1)
				return;
			else
				dev_buf2[tid] = curcell;
		}
		
		int in = -1, out = -1;

		while (1)
		{
			dev_c[curcell].rnum = 999;

			int intertype, info;
			out = device_getSegExitBorder(p0, p1, dev_v, dev_c, curcell, dev_b, in, intertype, info);

			if (out == -1)
				break;

			if (intertype == INTER_X || intertype == INTER_T1)
			{
				in = dev_b[out].partner;
				curcell = dev_b[in].parent;
			}
			else if (intertype == INTER_T2)
			{
				devdouble3 curpoint = dev_v[dev_b[out].v[info]];
				devdouble3 e = p1 - p0;
				e.normalize();
				curpoint = curpoint + e * TRAVEL_EPSILON;

				int otheredge;
				int u = dev_c[curcell].bnum;

				if (dev_b[out].partner != -1)
				{
					int ret = device_isPointInsideCell(curpoint, dev_b[dev_b[out].partner].parent, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						in = dev_b[out].partner;
						continue;
					}
				}

				if (info == 0)
					otheredge = dev_c[curcell].bstart + (dev_b[out].id + u - 1) % u;
				else
					otheredge = dev_c[curcell].bstart + (dev_b[out].id + 1) % u;

				if (dev_b[otheredge].partner != -1)
				{
					int ret = device_isPointInsideCell(curpoint, dev_b[dev_b[otheredge].partner].parent, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						in = dev_b[otheredge].partner;
						continue;
					}
					else
						return;
				}
				else
				{
					return;
				}
			}
			else if (intertype == INTER_OVERLAP)
			{
				devdouble3 e[2];
				e[0] = p1 - p0;
				e[0].normalize();
				e[1] = dev_v[dev_b[out].v[1]] - dev_v[dev_b[out].v[0]];
				e[1].normalize();
				double temp = e[0].dot(e[1]);
				int u = dev_c[curcell].bnum;
				if (temp > 0)
				{
					curcell = dev_b[dev_b[dev_c[curcell].bstart + (out - dev_c[curcell].bstart + 1) % u].partner].parent;
					in = dev_b[dev_c[curcell].bstart + (out - dev_c[curcell].bstart + 1) % u].partner;
				}
				else
				{
					curcell = dev_b[dev_b[dev_c[curcell].bstart + (out - dev_c[curcell].bstart + u - 1) % u].partner].parent;
					in = dev_b[dev_c[curcell].bstart + (out - dev_c[curcell].bstart + u - 1) % u].partner;			}
			}
		}
	}
}


__global__ void kernel_Classify_setDECell(
	devdouble3* dev_pv,
	int* dev_pe,
	devdouble3* dev_v,
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int hc_start,
	int hc_num,
	int* dev_refs)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < hc_num)
	{
		if (dev_c[hc_start + tid].rnum == 0 && dev_c[hc_start + tid].type == CELL_UNDEFINED)
		{
			int partneredge = -1;
			int partnercell = -1;
			for (int i = 0; i < dev_c[hc_start + tid].bnum; i++)
			{
				partneredge = dev_b[dev_c[hc_start + tid].bstart + i].partner;
				if (partneredge != -1 && dev_c[dev_b[partneredge].parent].subcell != -1)
				{
					partnercell = dev_b[partneredge].parent;
					break;
				}
			}

			if (partnercell != -1)
			{
				int targetcell = dev_b[partneredge].subcell;
				while (dev_c[targetcell].subcell != -1)
					targetcell = dev_c[targetcell].subcell;

				if (dev_c[targetcell].rnum == 0)
				{
					dev_c[hc_start + tid].type = dev_c[targetcell].type;
					return;
				}
				else
				{
					devdouble3 v0 = dev_v[dev_b[partneredge].v[0]];
					devdouble3 v1 = dev_v[dev_b[partneredge].v[1]];
					devdouble3 midpoint = (v0 + v1) * 0.5;
					midpoint.normalize();

					devdouble3 center;
					for (int i = 0; i < dev_c[targetcell].bnum; i++)
						center = center + dev_v[dev_b[dev_c[targetcell].bstart + i].v[0]];
					center = center / dev_c[targetcell].bnum;
					center.normalize();

					int inter_count = 0;
					for (int i = 0; i < dev_c[targetcell].rnum; i++)
					{
						int se = dev_refs[dev_c[targetcell].rstart + i];
						devdouble3 u0 = dev_pv[dev_pe[se * 2]];
						devdouble3 u1 = dev_pv[dev_pe[se * 2 + 1]];
						bool ret = dev_isInterSeg(u0, u1, center, midpoint);
						if (ret)
							inter_count++;
					}

					if (inter_count % 2 == 0)
						dev_c[hc_start + tid].type = dev_c[targetcell].type;
					else
						dev_c[hc_start + tid].type = 3 - dev_c[targetcell].type;
					return;
				}
			}
		}
	}
}


__global__ void kernel_Classify_setIECell(
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int hc_start,
	int hc_num)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < hc_num)
	{
		if (dev_c[hc_start + tid].type == CELL_UNDEFINED && dev_c[hc_start + tid].subcell == -1)
		{
			for (int i = 0; i < dev_c[hc_start + tid].bnum; i++)
			{
				int partneredge = dev_b[dev_c[hc_start + tid].bstart + i].partner;
				if (partneredge != -1)
				{
					int partnercell = dev_b[partneredge].parent;
					if (dev_c[partnercell].rnum == 0 && dev_c[partnercell].type != CELL_UNDEFINED)
					{
						dev_c[hc_start + tid].type = dev_c[partnercell].type;
						return;
					}
				}
			}
		}
	}
}


__global__ void kernel_Classify_setTopECell(
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int hc_num)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < hc_num)
	{
		for (int j = 0; j < 3; j++)
		{
			if (dev_c[tid].rnum == 0 && dev_c[tid].type == CELL_UNDEFINED)
			{
				for (int i = 0; i < dev_c[tid].bnum; i++)
				{
					int partneredge = dev_b[dev_c[tid].bstart + i].partner;
					if (partneredge != -1 && dev_c[dev_b[partneredge].parent].type != CELL_UNDEFINED)
					{
						dev_c[tid].type = dev_c[dev_b[partneredge].parent].type;
						break;
					}
				}
			}
		}
	}
}


__global__ void kernel_PISP(
	devdouble3* dev_pv,
	int* dev_pe,
	devdouble3* dev_v,
	AHHGBorder* dev_b,
	AHHGCell* dev_c,
	int* dev_refs,
	devdouble3* dev_tp,
	int tpcount,
	int* dev_result)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < tpcount)
	{
		devdouble3 p = dev_tp[tid];

		int curcell = -1;
		for (int i = 0; i < 42; i++)
		{
			curcell = device_isPointInsideCell(p, i, dev_c, dev_b);
			if (curcell != -1)
				break;
		}

		if (curcell == -1)
			return;

		while (dev_c[curcell].subcell != -1)
		{
			bool find = false;
			int ret = device_isPointInsideCell(p, dev_c[curcell].subcell, dev_c, dev_b);
			if (ret != -1)
				curcell = ret;
			else
			{
				for (int i = 0; i < dev_c[curcell].bnum; i++)
				{
					int tempcell = dev_b[dev_c[curcell].bstart + i].subcell;
					if (tempcell == -1)
						continue;
					int ret = device_isPointInsideCell(p, tempcell, dev_c, dev_b);
					if (ret != -1)
					{
						curcell = ret;
						find = true;
						break;
					}
				}
				if (find == false)
					return;
			}
		}

		if (dev_c[curcell].rnum == 0)
		{
			dev_result[tid] = dev_c[curcell].type;
		}
		else
		{
			devdouble3 center;
			for (int i = 0; i < dev_c[curcell].bnum; i++)
				center = center + dev_v[dev_b[dev_c[curcell].bstart + i].v[0]];
			center = center / dev_c[curcell].bnum;
			center.normalize();

			int intercount = 0;
			for (int i = 0; i < dev_c[curcell].rnum; i++)
			{
				int se = dev_refs[dev_c[curcell].rstart + i];
				devdouble3 p0 = dev_pv[dev_pe[se * 2]];
				devdouble3 p1 = dev_pv[dev_pe[se * 2 + 1]];
				int ret = device_isEEInter(p0, p1, center, p);
				if (ret == true)
					intercount++;
			}

			if (intercount % 2 == 0)
				dev_result[tid] = dev_c[curcell].type;
			else
				dev_result[tid] = 3 - dev_c[curcell].type;
		}
	}
}


void GAHHG::initDevMemory()
{
	if (this->dev_v != NULL)
		cudaFree(this->dev_v);
	cudaMalloc((void**)&this->dev_v, sizeof(devdouble3) * this->pParam->max_v_num);

	if (this->dev_b != NULL)
		cudaFree(this->dev_b);
	cudaMalloc((void**)&this->dev_b, sizeof(AHHGBorder) * this->pParam->max_b_num);

	if (this->dev_c != NULL)
		cudaFree(this->dev_c);
	cudaMalloc((void**)&this->dev_c, sizeof(AHHGCell) * this->pParam->max_c_num);

	for (int i = 0; i < 3; i++)
		cudaMalloc((void**)&this->dev_r[i], sizeof(int) * this->pParam->max_r_num);

	cudaMalloc((void**)&this->dev_se_buf[0], sizeof(int) * this->pParam->max_buf_num);
	cudaMalloc((void**)&this->dev_se_buf[1], sizeof(int) * this->pParam->max_buf_num);

	cudaMalloc((void**)&this->dev_cell_buf[0], sizeof(int) * this->pParam->max_buf_num);
	cudaMalloc((void**)&this->dev_cell_buf[1], sizeof(int) * this->pParam->max_buf_num);
}


void GAHHG::copyGrid0HostToDevice()
{
	if (this->vertexes != NULL && this->vertexes->curNum != 0)
	{
		cudaMemcpy(this->dev_v, this->vertexes->elem, sizeof(mydouble3) * this->g[0].vnum, cudaMemcpyHostToDevice);
	}

	if (this->borders != NULL && this->borders->curNum != 0)
	{
		cudaMemcpy(this->dev_b, this->borders->elem, sizeof(AHHGBorder) * this->g[0].bnum, cudaMemcpyHostToDevice);
	}

	if (this->cells != NULL && this->cells->curNum != 0)
	{
		cudaMemcpy(this->dev_c, this->cells->elem, sizeof(AHHGCell) * this->g[0].cnum, cudaMemcpyHostToDevice);
	}

	delete this->cells;
	this->cells = NULL;
	delete this->borders;
	this->borders = NULL;
	delete this->vertexes;
	this->vertexes = NULL;
}


bool GAHHG::divGrid(int gridid)
{
	printf("Subdivide level %d on GPU.\n", gridid);

	int cell_blocks = (int)this->g[gridid].cnum / this->pParam->blocksize_create + 1;
	
	kernel_DivGridFast_zeroRefCount << <cell_blocks, this->pParam->blocksize_create >> > (
		this->dev_c,
		this->g[gridid].cstart,
		this->g[gridid].cnum,
		this->dev_b);

	int se_blocks = (int)this->pSP->pv_num / this->pParam->blocksize_create + 1;
	kernel_DivGridFast_markNECell << <se_blocks, this->pParam->blocksize_create >> > (
			this->pSP->dev_pv,
			this->pSP->dev_pe,
			this->pSP->pv_num,
			this->dev_v,
			this->dev_b,
			this->dev_c,
			this->g[gridid].cstart,
			this->g[gridid].cnum,
			gridid,
			this->dev_r[1],
			this->dev_r[2]);
				
	this->checkStorage(4,this->g[gridid].cnum);

	cudaMemset(this->dev_se_buf[0], 0, sizeof(int));
	kernel_DivGrid_calcCenterOffsets << <cell_blocks, this->pParam->blocksize_create >> > (
		this->dev_cell_buf[0],
		this->dev_cell_buf[1],
		this->dev_c,
		this->g[gridid].cstart,
		this->g[gridid].cnum,
		this->dev_se_buf[0]);

	if (this->pParam->autoDepth == true)
	{
		int U;
		cudaMemcpy(&U, &this->dev_se_buf[0][0], sizeof(int), cudaMemcpyDeviceToHost);
		float k = U * 1.0 / this->pSP->pv_num;
		this->pStat->gamma = k;
		if (k > 0.1 || k < 0.1 && (0.1 - k) < 0.01)
		{
			return false;
		}
	}

	myscan(this->dev_se_buf[0], this->dev_cell_buf[0], this->g[gridid].cnum, false);
	myscan(this->dev_se_buf[1], this->dev_cell_buf[1], this->g[gridid].cnum, false);
	int* temp;
	temp = this->dev_se_buf[0];
	this->dev_se_buf[0] = this->dev_cell_buf[0];
	this->dev_cell_buf[0] = temp;
	temp = this->dev_se_buf[1];
	this->dev_se_buf[1] = this->dev_cell_buf[1];
	this->dev_cell_buf[1] = temp;

	int buf; 
	cudaMemcpy(&buf, &this->dev_cell_buf[0][this->g[gridid].cnum - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->checkStorage(0, this->g[gridid].vstart + this->g[gridid].vnum + buf + 6);

	this->checkStorage(1, this->g[gridid].bstart + this->g[gridid].bnum + (buf + 6) + (buf + 6) * 6);
	int buf2 = buf;

	cudaMemcpy(&buf, &this->dev_cell_buf[1][this->g[gridid].cnum - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->checkStorage(2, this->g[gridid].cstart + this->g[gridid].cnum + (buf + 1) + (buf2 + 6));

	kernel_DivGrid_createCenter << <cell_blocks, this->pParam->blocksize_create >> > (
			this->dev_v,
			this->g[gridid].vstart + this->g[gridid].vnum,
			this->dev_b,
			this->g[gridid].bstart + this->g[gridid].bnum,
			this->dev_c,
			this->g[gridid].cstart,
			this->g[gridid].cnum,
			this->dev_cell_buf[0],
			this->dev_cell_buf[1]);

	cudaMemcpy(&buf, &this->dev_cell_buf[0][this->g[gridid].cnum - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->g[gridid + 1].vstart = this->g[gridid].vstart + this->g[gridid].vnum;
	this->g[gridid + 1].vnum = buf;
	this->g[gridid + 1].bstart = this->g[gridid].bstart + this->g[gridid].bnum;
	this->g[gridid + 1].bnum = buf;
	cudaMemcpy(&buf, &this->dev_cell_buf[1][this->g[gridid].cnum - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->g[gridid + 1].cstart = this->g[gridid].cstart + this->g[gridid].cnum;
	this->g[gridid + 1].cnum = buf;

	this->g[gridid].rnum = buf;

	this->checkStorage(4, this->g[gridid + 1].bnum);
	int sce_blocks = (int)(this->g[gridid + 1].bnum / this->pParam->blocksize_create + 1);
	kernel_DivGrid_calcBorderSCOffsets << <sce_blocks, this->pParam->blocksize_create >> > (
		this->dev_b,
		this->g[gridid + 1].bstart,
		this->g[gridid + 1].bnum,
		this->dev_c,
		this->dev_cell_buf[0],
		this->dev_cell_buf[1]);

	myscan(this->dev_se_buf[0], this->dev_cell_buf[0], this->g[gridid + 1].bnum, false);
	myscan(this->dev_se_buf[1], this->dev_cell_buf[1], this->g[gridid + 1].bnum, false);

	temp = this->dev_se_buf[0];
	this->dev_se_buf[0] = this->dev_cell_buf[0];
	this->dev_cell_buf[0] = temp;
	temp = this->dev_se_buf[1];
	this->dev_se_buf[1] = this->dev_cell_buf[1];
	this->dev_cell_buf[1] = temp;

	kernel_DivGrid_createBorderSC << <sce_blocks, this->pParam->blocksize_create >> > (
		this->dev_v,
		this->dev_b,
		this->g[gridid + 1].bstart,
		this->g[gridid + 1].bnum,
		this->dev_c,
		this->g[gridid + 1].cstart,
		this->g[gridid + 1].cnum,
		this->dev_cell_buf[0],
		this->dev_cell_buf[1]);

	int bsc_start, bsc_num;
	bsc_start = this->g[gridid + 1].cstart + this->g[gridid + 1].cnum;
	cudaMemcpy(&buf, &this->dev_cell_buf[1][this->g[gridid + 1].bnum - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->g[gridid + 1].cnum += buf;
	bsc_num = buf;

	cudaMemcpy(&buf, &this->dev_cell_buf[0][this->g[gridid + 1].bnum - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->g[gridid + 1].bnum += buf;

	int bsc_blocks = (int)(bsc_num / this->pParam->blocksize_create + 1);
	kernel_DivGrid_setSCPartner << <bsc_blocks, this->pParam->blocksize_create >> > (
		dev_b,
		dev_c,
		bsc_start,
		bsc_num);

	return true;
}


void GAHHG::dispatchSELast()
{
	int se_blocks = (int)this->pSP->pv_num / this->pParam->blocksize_create + 1;
	cudaMemset(this->dev_se_buf[0], 0, sizeof(int) * this->pSP->pv_num);
	kernel_calcRefNum << <se_blocks, this->pParam->blocksize_create >> > (
		this->dev_se_buf[0],
		this->pSP->dev_pv,
		this->pSP->dev_pe,
		this->pSP->pv_num,
		this->dev_v,
		this->dev_c,
		this->g[this->gnum - 1].cnum,
		this->dev_b,
		this->dev_r[2],
		this->gnum - 1);

	int lastnum;
	cudaMemcpy(&lastnum, &this->dev_se_buf[0][this->pSP->pv_num - 1], sizeof(int), cudaMemcpyDeviceToHost);

	myscan(this->dev_se_buf[1], this->dev_se_buf[0], this->pSP->pv_num, false);
	
	cudaMemcpy(&this->g[this->gnum - 1].rnum, &this->dev_se_buf[1][this->pSP->pv_num - 1], sizeof(int), cudaMemcpyDeviceToHost);
	this->g[this->gnum - 1].rnum += lastnum;
	
	checkStorage(3, this->g[this->gnum - 1].rnum);

	kernel_writeRef << <se_blocks, this->pParam->blocksize_create >> > (
		this->dev_r[0],
		this->dev_se_buf[1],
		this->pSP->dev_pv,
		this->pSP->dev_pe,
		this->pSP->pv_num,
		this->dev_v,
		this->dev_c,
		this->g[this->gnum - 1].cnum,
		this->dev_b,
		this->dev_r[2],
		this->gnum - 1);

	kernel_su2pre << <se_blocks, this->pParam->blocksize_create >> > (
		this->dev_se_buf[0],
		this->dev_se_buf[1],
		this->pSP->pv_num);

	cudaMemset(this->dev_r[1], 0, sizeof(int) * this->g[this->gnum - 1].rnum);
	kernel_EC2CS_setEFlag << <se_blocks, this->pParam->blocksize_create >> > (
		this->dev_r[1],
		this->dev_se_buf[1],
		this->pSP->pv_num);
	myscan(this->dev_r[2], this->dev_r[1], this->g[this->gnum - 1].rnum, false);

	int ref_blocks = (int)this->g[this->gnum - 1].rnum / this->pParam->blocksize_create + 1;
	kernel_EC2CS_adjustEFlag << <ref_blocks, this->pParam->blocksize_create >> > (
		this->dev_r[1],
		this->dev_r[2],
		this->g[this->gnum - 1].rnum);

	this->checkStorage(4, this->g[this->gnum - 1].cnum);

	cudaMemset(this->dev_cell_buf[0], 0, sizeof(int) * this->g[this->gnum - 1].cnum);
	kernel_EC2CS_countEinC << <ref_blocks, this->pParam->blocksize_create >> > (
		dev_r[0],
		this->g[this->gnum - 1].rnum,
		dev_cell_buf[0],
		this->g[this->gnum - 1].cstart);
		
	myscan(this->dev_cell_buf[1], this->dev_cell_buf[0], this->g[this->gnum - 1].cnum, false);

	kernel_EC2CS_rearrangeCERef << <ref_blocks, this->pParam->blocksize_create >> > (
		this->dev_r[0],
		this->dev_r[2],
		this->dev_r[1],
		this->g[this->gnum - 1].rnum,
		this->dev_cell_buf[1],
		this->g[this->gnum - 1].cstart);

	int cell_blocks = (int)this->g[this->gnum - 1].cnum / this->pParam->blocksize_create + 1;
	kernel_EC2CS_setCEInfo << <cell_blocks, this->pParam->blocksize_create >> > (
		this->dev_cell_buf[0],
		this->dev_cell_buf[1],
		this->dev_c,
		this->g[this->gnum - 1].cstart,
		this->g[this->gnum - 1].cnum);
		
}


void GAHHG::classifyBottomGrid()
{
	printf("Classify level %d cells.\n", this->gnum - 1);

	int cell_blocks = (int)this->g[this->gnum - 1].cnum / this->pParam->blocksize_classify + 1;
	kernel_Classify_setBottomNECell << <cell_blocks, this->pParam->blocksize_classify >> > (
		this->pSP->dev_pv,
		this->pSP->dev_pe,
		this->dev_v,
		this->dev_b,
		this->dev_c,
		this->g[this->gnum - 1].cstart,
		this->g[this->gnum - 1].cnum,
		this->dev_r[0]);

	kernel_Classify_setBottomDECell << <cell_blocks, this->pParam->blocksize_classify >> > (
		this->pSP->dev_pv,
		this->pSP->dev_pe,
		this->dev_v,
		this->dev_b,
		this->dev_c,
		this->dev_r[0],
		this->g[this->gnum - 1].cstart,
		this->g[this->gnum - 1].cnum);

	kernel_Classify_setIECell << <cell_blocks, this->pParam->blocksize_classify >> > (
		this->dev_b,
		this->dev_c,
		this->g[this->gnum - 1].cstart,
		this->g[this->gnum - 1].cnum);
}


void GAHHG::checkStorage(int arrayid, int requirednum)
{
	switch (arrayid)
	{
		case 0:
			if (requirednum > this->pParam->max_v_num)
			{
				int newsize = this->pParam->max_v_num * 2;
				while (newsize < requirednum)
					newsize *= 2;
				devdouble3* newarray;
				cudaMalloc((void**)&newarray, sizeof(devdouble3)*newsize);
				cudaMemcpy(newarray, this->dev_v, sizeof(devdouble3)*this->pParam->max_v_num, cudaMemcpyDeviceToDevice);
				cudaFree(this->dev_v);
				this->dev_v = newarray;
				this->pParam->max_v_num = newsize;
			}
			break;
		case 1:
			if (requirednum > this->pParam->max_b_num)
			{
				int newsize = this->pParam->max_b_num * 2;
				while (newsize < requirednum)
					newsize *= 2;
				AHHGBorder* newarray;
				cudaMalloc((void**)&newarray, sizeof(AHHGBorder) * newsize);
				cudaMemcpy(newarray, this->dev_b, sizeof(AHHGBorder) * this->pParam->max_b_num, cudaMemcpyDeviceToDevice);
				cudaFree(this->dev_b);
				this->dev_b = newarray;
				this->pParam->max_b_num = newsize;
			}
			break;
		case 2:
			if (requirednum > this->pParam->max_c_num)
			{
				int newsize = this->pParam->max_c_num * 2;
				while (newsize < requirednum)
					newsize *= 2;
				AHHGCell* newarray;
				cudaMalloc((void**)&newarray, sizeof(AHHGCell) * newsize);
				cudaMemcpy(newarray, this->dev_c, sizeof(AHHGCell) * this->pParam->max_c_num, cudaMemcpyDeviceToDevice);
				cudaFree(this->dev_c);
				this->dev_c = newarray;
				this->pParam->max_c_num = newsize;
			}
			break;
		case 3:
			if (requirednum > this->pParam->max_r_num)
			{
				int newsize = this->pParam->max_r_num * 2;
				while (newsize < requirednum)
					newsize *= 2;
				int* newarray;
				for (int i = 0; i < 3; i++)
				{
					cudaMalloc((void**)&newarray, sizeof(int) * newsize);
					cudaMemcpy(newarray, this->dev_r[i], sizeof(int) * this->pParam->max_r_num, cudaMemcpyDeviceToDevice);
					cudaFree(this->dev_r[i]);
					this->dev_r[i] = newarray;
				}
				this->pParam->max_r_num = newsize;
			}
			break;
		case 4:
			if (requirednum > this->pParam->max_buf_num)
			{
				int newsize = this->pParam->max_buf_num * 2;
				while (newsize < requirednum)
					newsize *= 2;
				int* newarray;
				for (int i = 0; i < 2; i++)
				{
					cudaMalloc((void**)&newarray, sizeof(int) * newsize);
					cudaMemcpy(newarray, this->dev_cell_buf[i], sizeof(int) * this->pParam->max_buf_num, cudaMemcpyDeviceToDevice);
					cudaFree(this->dev_cell_buf[i]);
					this->dev_cell_buf[i] = newarray;
				}
				for (int i = 0; i < 2; i++)
				{
					cudaMalloc((void**)&newarray, sizeof(int) * newsize);
					cudaMemcpy(newarray, this->dev_se_buf[i], sizeof(int) * this->pParam->max_buf_num, cudaMemcpyDeviceToDevice);
					cudaFree(this->dev_se_buf[i]);
					this->dev_se_buf[i] = newarray;
				}
				this->pParam->max_buf_num = newsize;
			}
			break;
		default:
			break;
	}
}


void GAHHG::clear()
{
	if (this->g != NULL)
		delete this->g;
	if (this->vertexes != NULL)
		delete this->vertexes;
	if (this->borders != NULL)
		delete this->borders;
	if (this->cells != NULL)
		delete this->cells;
	if (this->dev_v != NULL)
		cudaFree(this->dev_v);
	if (this->dev_b != NULL)
		cudaFree(this->dev_b);
	if (this->dev_c != NULL)
		cudaFree(this->dev_c);
	if (this->dev_cell_buf[0] != NULL)
		cudaFree(this->dev_cell_buf[0]);
	if (this->dev_cell_buf[1] != NULL)
		cudaFree(this->dev_cell_buf[1]);
	if (this->dev_se_buf[0] != NULL)
		cudaFree(this->dev_se_buf[0]);
	if (this->dev_se_buf[1] != NULL)
		cudaFree(this->dev_se_buf[1]);
	if (this->dev_r[0] != NULL)
		cudaFree(this->dev_r[0]);
	if (this->dev_r[1] != NULL)
		cudaFree(this->dev_r[0]);
	if (this->dev_r[2] != NULL)
		cudaFree(this->dev_r[0]);
	memset(this,0,sizeof(*this));
}


void GPISP::createAHHGInit()
{
	this->pAHHG->pParam->max_c_num = mymax(this->pPolygon->pv_num, 42*2);
	this->pAHHG->pParam->max_v_num = mymax(this->pPolygon->pv_num, 80*2);
	this->pAHHG->pParam->max_b_num = mymax(this->pPolygon->pv_num*4, 240*2);
	this->pAHHG->pParam->max_r_num = this->pPolygon->pv_num * 2;
	this->pAHHG->pParam->max_buf_num = this->pAHHG->pParam->max_b_num;

	this->pAHHG->gnum = this->pParam->maxDepth;
	this->pAHHG->g = (GAHHGGrid*)malloc(sizeof(GAHHGGrid) * this->pAHHG->gnum);
	memset(this->pAHHG->g, 0, sizeof(GAHHGGrid) * this->pAHHG->gnum);
	this->pAHHG->vertexes = new myArray<mydouble3>(80*2);
	this->pAHHG->borders = new myArray<AHHGBorder>(240*2);
	this->pAHHG->cells = new myArray<AHHGCell>(42*2);
	
	initgrid0(this->pAHHG->vertexes, this->pAHHG->borders, this->pAHHG->cells);
	
	pAHHG->g[0].vstart = 0;
	pAHHG->g[0].vnum = this->pAHHG->vertexes->curNum;
	pAHHG->g[0].bstart = 0;
	pAHHG->g[0].bnum = this->pAHHG->borders->curNum;
	pAHHG->g[0].cstart = 0;
	pAHHG->g[0].cnum = this->pAHHG->cells->curNum;

	this->pAHHG->initDevMemory();
	
	this->pAHHG->copyGrid0HostToDevice();
}

void GPISP::createAHHG(int maxdepth)
{
	myClock clock;
	clock.start();

	this->pAHHG->gnum = maxdepth;
	bool early_terminate = false;
	for (int i = 0; i < this->pAHHG->gnum - 1; i++)
	{
		bool ret = this->pAHHG->divGrid(i);
		if (ret == false)
		{
			early_terminate = true;
			this->pAHHG->gnum = i + 1;
			break;
		}
	}

	if (this->pParam->maxDepth == 1)
	{
		this->pStat->gamma = 42 * 1.0 / this->pPolygon->pv_num;
	}

	if(early_terminate==true && this->pAHHG->gnum != 1)
	{
		int* temp = this->pAHHG->dev_r[1];
		this->pAHHG->dev_r[1] = this->pAHHG->dev_r[2];
		this->pAHHG->dev_r[2] = temp;
	}

	this->pAHHG->dispatchSELast();

	int* tempptr = this->pAHHG->dev_r[0];
	this->pAHHG->dev_r[0] = this->pAHHG->dev_r[1];
	this->pAHHG->dev_r[1] = tempptr;

	cudaDeviceSynchronize();
	clock.end();
	this->pStat->createAHHGTime = clock.time;
}


void GPISP::classifyAHHG()
{
	myClock clock;
	clock.start();

	this->pAHHG->classifyBottomGrid();

	int i = this->pAHHG->gnum - 2;
	while (i >= 0)
	{
		printf("Classify level %d cells.\n", i);

		int cell_blocks = (int)this->pAHHG->g[i].cnum / this->pParam->blocksize_create + 1;
		kernel_Classify_setDECell << <cell_blocks, this->pParam->blocksize_create >> > (
			this->pAHHG->pSP->dev_pv,
			this->pAHHG->pSP->dev_pe,
			this->pAHHG->dev_v,
			this->pAHHG->dev_b,
			this->pAHHG->dev_c,
			this->pAHHG->g[i].cstart,
			this->pAHHG->g[i].cnum,
			this->pAHHG->dev_r[0]);

		kernel_Classify_setIECell << <cell_blocks, this->pParam->blocksize_classify >> > (
			this->pAHHG->dev_b,
			this->pAHHG->dev_c,
			this->pAHHG->g[i].cstart,
			this->pAHHG->g[i].cnum);

		i--;
	};

	kernel_Classify_setTopECell << <1, this->pParam->blocksize_create >> > (
		this->pAHHG->dev_b,
		this->pAHHG->dev_c,
		this->pAHHG->g[0].cnum);

	cudaDeviceSynchronize();
	clock.end();
	this->pStat->classifyAHHGTime = clock.time;
	this->pStat->preTime = this->pStat->createAHHGTime + this->pStat->classifyAHHGTime;
}


void GPISP::test()
{
	myClock timer;
	timer.start();

	int tp_blocks;
	tp_blocks = (int)(this->pTP->tpNum / this->pParam->blocksize_pisp) + 1;
	kernel_PISP << <tp_blocks, this->pParam->blocksize_pisp >> > (
		this->pAHHG->pSP->dev_pv,
		this->pAHHG->pSP->dev_pe,
		this->pAHHG->dev_v,
		this->pAHHG->dev_b,
		this->pAHHG->dev_c,
		this->pAHHG->dev_r[0],
		this->pTP->dev_tp,
		this->pTP->tpNum,
		this->pTP->dev_result);
	cudaDeviceSynchronize();

	timer.end();
	this->pStat->PISPTime = timer.time;

	cudaMemcpy(this->pTP->result, this->pTP->dev_result, sizeof(int)*this->pTP->tpNum,cudaMemcpyDeviceToHost);
}


void GPISP::getStat()
{
	this->pStat->AHHG_gpu_vnum = this->pAHHG->g[this->pAHHG->gnum - 1].vstart + this->pAHHG->g[this->pAHHG->gnum - 1].vnum;
	this->pStat->AHHG_gpu_bnum = this->pAHHG->g[this->pAHHG->gnum - 1].bstart + this->pAHHG->g[this->pAHHG->gnum - 1].bnum;
	this->pStat->AHHG_gpu_cnum = this->pAHHG->g[this->pAHHG->gnum - 1].cstart + this->pAHHG->g[this->pAHHG->gnum - 1].cnum;
	this->pStat->AHHG_gpu_rnum = this->pAHHG->g[this->pAHHG->gnum - 1].rnum;

	this->pStat->AHHG_gpu_vspace = sizeof(devdouble3) * this->pStat->AHHG_gpu_vnum;
	this->pStat->AHHG_gpu_bspace = sizeof(AHHGBorder) * this->pStat->AHHG_gpu_bnum;
	this->pStat->AHHG_gpu_cspace = sizeof(AHHGCell) * this->pStat->AHHG_gpu_cnum;
	this->pStat->AHHG_gpu_rspace = sizeof(int) * this->pStat->AHHG_gpu_rnum;
	
	this->pStat->AHHG_cpu_space = sizeof(GAHHG) + sizeof(GAHHGGrid) * this->pAHHG->gnum + sizeof(GPISPParam) + sizeof(GPISPStat);
	this->pStat->AHHG_gpu_min_space = this->pStat->AHHG_gpu_vspace + 
										this->pStat->AHHG_gpu_bspace + 
										this->pStat->AHHG_gpu_cspace + 
										this->pStat->AHHG_gpu_rspace;
		
	this->pStat->inTpNum = 0;
	for (int i = 0; i < this->pTP->tpNum; i++)
	{
		if (this->pTP->result[i] == CELL_IN)
			this->pStat->inTpNum++;
	}
	this->pStat->PISPUnitTime = this->pStat->PISPTime / this->pTP->tpNum;
}

void GPISP::exportStat2TXT(char* filename, int desttype)
{
	FILE* fp;
	if (desttype == 0)
		fp = stdout;
	else
	{
		int ret = fopen_s(&fp, filename, "a");
		if (ret != 0)
		{
			printf("Can not open file %s\n.", filename);
			return;
		}
	}

	time_t nowtime = time(NULL);
	tm local;
	localtime_s(&local, &nowtime);
	fprintf(fp, "\n----------------------------------------%d/%d/%d %d:%d:%d-----------------------------------------\n",
		local.tm_year + 1900,
		local.tm_mon + 1,
		local.tm_mday,
		local.tm_hour,
		local.tm_min,
		local.tm_sec);

	fprintf(fp, "\nPolygon file: %s\n", this->pParam->polyFileName);
	fprintf(fp, "Number of polygon vertice: %d\n", this->pPolygon->pv_num);

	fprintf(fp, "Query point file: %s\n", this->pParam->tpFileName);
	fprintf(fp, "Number of query points: %d\n", this->pTP->tpNum);

	fprintf(fp, "Depth of GAHHG: %d\n", this->pAHHG->gnum);
	fprintf(fp, "Time of preprocessing: %f seconds (creat GAHHG: %f seconds, classify GAHHG: %f seconds)\n",
		this->pStat->preTime, this->pStat->createAHHGTime, this->pStat->classifyAHHGTime);
	fprintf(fp, "Storage overhead of GAHHG on CPU: %d Bytes\n", this->pStat->AHHG_cpu_space);
	fprintf(fp, "Storage overhead of GAHHG on GPU: %d KB\n",
		(int)(this->pStat->AHHG_gpu_min_space / 1000));

	fprintf(fp, "Total test time for all query points: %f seconds\n", this->pStat->PISPTime);
	fprintf(fp, "Average test time per query point: %.5e seconds\n", this->pStat->PISPUnitTime);
	fprintf(fp, "Test result: inside %d, outside %d\n", this->pStat->inTpNum, this->pTP->tpNum - this->pStat->inTpNum);

	fprintf(fp, "\n---------------------------------------------------------------------------------------------------\n");

	fclose(fp);
}


void GPISP::clear()
{
	if (this->pParam != NULL)
		delete this->pParam;
	if (this->pStat != NULL)
		delete this->pStat;
	if (this->pPolygon != NULL)
		delete this->pPolygon;
	if (this->pTP != NULL)
		delete this->pTP;
	if (this->pAHHG != NULL)
		delete this->pAHHG;
}
