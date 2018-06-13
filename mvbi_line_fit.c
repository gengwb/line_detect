/***************************************************************************************************
* 版权信息：Copyright (c) 2015, 杭州海康威视数字技术股份有限公司
*
* 文件名称：_mvbi_line_fit.c
* 摘    要：直线拟合
*
* 当前版本：V0.1.0
* 作    者：邓志辉
* 日    期：2015-08-05
* 备    注：1.创建   
***************************************************************************************************/

#include <stdlib.h>
#include "hka_defs.h"
#include "mvb_defs.h"
#include "mvbi_line_fit.h"
#include "mvbi_line_fit_type.h"

/***************************************************************************************************
* 功  能：随机选取3个数
* 参  数：*
*         max_num  – I  最大数编号
*         triplet   - O  3个数
* 返回值：HKA_VOID
* 备  注：从[0, max_num - 1]只选取3个数
***************************************************************************************************/
static HKA_VOID MVBI_LineFitRansac_get_random_triplet(HKA_S32 max_num, HKA_S32 triplet[2])
{
	HKA_S32 index   = 0;
	HKA_S32 new_one = HKA_FALSE; 
	HKA_S32 r       = 0;
	HKA_S32 i       = 0;

	while(index < 2)
	{
		new_one = HKA_TRUE;

		r = rand() % max_num;

		for(i = 0; i < index; i++)
		{
			if(r == triplet[i])
			{
				new_one = HKA_FALSE;
				break;
			}
		}

		if(new_one == HKA_TRUE)
		{
			triplet[index] = r;
			index++;
		}
	}
}

/***************************************************************************************************
* 功  能：用两个点计算直线参数
* 参  数：* 
*         p1            -I         点1
*         p2            -I         点2
*         line_coe      -O         直线参数
* 返回值：无
* 备  注：
***************************************************************************************************/
static HKA_S32 MVBI_LineFitRansac_line_coe_32f(HKA_POINT_F p1, HKA_POINT_F p2, HKA_F32 line_coe[3])
{
	HKA_F32 a0   = 0.0f;
	HKA_F32 b0   = 0.0f;
	HKA_F32 c0   = 0.0f;
	HKA_F32 tem  = 0.0f; 
	HKA_F32 temp = 0.0f;

	a0 = 1.0f * (p1.y - p2.y);
	b0 = 1.0f * (p2.x - p1.x);
	c0 = 1.0f * (p1.x * p2.y - p1.y * p2.x);

	temp = a0 * a0 + b0 * b0;
	if (temp < 1e-5f)
	{
		return HKA_FALSE;
	}

	tem = 1.0f / MVB_SQRTF(temp);

	line_coe[0] = a0 * tem;
	line_coe[1] = b0 * tem;
	line_coe[2] = c0 * tem;

	return HKA_TRUE;
}

static HKA_S32 MVBI_LineFitRansac_line_coe_32s(HKA_POINT_I p1, HKA_POINT_I p2, HKA_F32 line_coe[3])
{
	HKA_S32 a0   = 0;
	HKA_S32 b0   = 0;
	HKA_S32 c0   = 0;
	HKA_F32 tem  = 0.0f; 
	HKA_S32 temp = 0;

	a0 = (p1.y - p2.y);
	b0 = (p2.x - p1.x);
	c0 = (p1.x * p2.y - p1.y * p2.x);

	temp = a0 * a0 + b0 * b0;
	if (temp == 0)
	{
		return HKA_FALSE;
	}

	tem = 1.0f / MVB_SQRTF(temp);

	line_coe[0] = a0 * tem;
	line_coe[1] = b0 * tem;
	line_coe[2] = c0 * tem;

	return HKA_TRUE;
}

/***************************************************************************************************
* 功  能：
* 参  数：*
*
* 返回值：状态
* 备  注：
***************************************************************************************************/
static HKA_F32 MVBI_LineFitRansac_line_dist_32f(HKA_F32 line_coef[3], HKA_POINT_F point)
{
	HKA_F32 dist = 0.0f;

	dist = line_coef[0] * point.x + line_coef[1] * point.y + line_coef[2];

	return dist;
}

static HKA_F32 MVBI_LineFitRansac_line_dist_32s(HKA_F32 line_coef[3], HKA_POINT_I point)
{
	HKA_F32 dist = 0.0f;

	dist = line_coef[0] * point.x + line_coef[1] * point.y + line_coef[2];

	return dist;
}

/***************************************************************************************************
* 功  能：计算直线度
* 参  数： 
*         pts         -I  点集
*         pts_num     -I  点的个数
*         mask        -I  点是否参加计算的掩膜
*         line_coe    -I  拟合出的直线系数
* 返回值：
* 备  注：
***************************************************************************************************/
static HKA_F32 MVBI_LineFitRansac_straightness_32f(HKA_POINT_F *pts, HKA_S32 pts_num,  HKA_F32 *line_coe,
											       HKA_S32 kernel_size)
{
	HKA_S32 num  = 0;
	HKA_S32 i    = 0;
	HKA_S32 pr   = 0;
	HKA_F32 dist = 0.0f;
	HKA_F32 msd  = 0.0f;
	HKA_F32 a    = 0.0f;
	HKA_F32 b    = 0.0f;
	HKA_F32 c    = 0.0f;

	a = line_coe[0];
	b = line_coe[1];
	c = line_coe[2];

	for (i = 0; i < pts_num; i++)
	{
		dist = a * pts[i].x + b * pts[i].y + c;
	    msd += (dist * dist);
	}

	if (pts_num == 0)
	{
		msd = 0.0f;
	}
	else
	{
		pr = num * kernel_size * kernel_size;
		msd = (1.0f - msd / pr) * 100.f;
	}

	return  msd;
}

static HKA_F32 MVBI_LineFitRansac_straightness_32s(HKA_POINT_I *pts, HKA_S32 pts_num,  HKA_F32 *line_coe,
												   HKA_S32 kernel_size)
{
	HKA_S32 num  = 0;
	HKA_S32 i    = 0;
	HKA_S32 pr   = 0;
	HKA_F32 dist = 0.0f;
	HKA_F32 msd  = 0.0f;
	HKA_F32 a    = 0.0f;
	HKA_F32 b    = 0.0f;
	HKA_F32 c    = 0.0f;

	a = line_coe[0];
	b = line_coe[1];
	c = line_coe[2];

	for (i = 0; i < pts_num; i++)
	{
		dist = a * pts[i].x + b * pts[i].y + c;
		msd += (dist * dist);
	}

	if (pts_num == 0)
	{
		msd = 0.0f;
	}
	else
	{
		pr = num * kernel_size * kernel_size;
		msd = (1.0f - msd / pr) * 100.f;
	}

	return  msd;
}

/***************************************************************************************************
* 功  能：ranac直线估计处理
* 参  数：*
*         point_data     - I  输入点坐标
*         point_num      - I  输入点个数
*         circle_center  - O  输出圆中心
*         circle_radius  - O  输出圆半径
*         circle_status  - O  输出圆存在与否状态
*         sample_th      - I  采用阈值
*         dist_th        - I  距离阈值
*         max_iter       - I  最大迭代次数
*         consensus      - I  最佳点集
* 返回值：状态
* 备  注：
***************************************************************************************************/
static HKA_S32 MVBI_LineFitRansac_process_32f(HKA_POINT_F *point_data, HKA_S32 point_num,
											  HKA_F32 line_coef[3], HKA_F32 *straight,
											  HKA_F32 sample_ratio, HKA_F32 dist_th, HKA_S32 max_iter, 
										      HKA_POINT_F *consensus)
{
	HKA_S32     trials       = 0;
	HKA_S32     num_iter     = 0;
	HKA_S32     iter         = 0;
	HKA_S32     triplet[2]   = {0};
	HKA_S32     flag         = HKA_FALSE;
	HKA_S32     ninlier      = 0;
	HKA_S32     k            = 0;
	HKA_F32     dist         = 0.0f;
	HKA_F32     dist_dev     = 0.0f;
	HKA_S32     sample_th    = 0;
	HKA_S32     line_sts     = HKA_FALSE;
	HKA_F32     best_coef[3] = {0.0f};
	HKA_F32     est_coef[3]  = {0.0f};
	HKA_S32     sts          = HKA_FALSE;
	HKA_S32     fail_num     = 0;
	HKA_S32     best_num     = 0;

	trials       = (HKA_S32)(MVB_LOGF(1.0f - 0.99f) / MVB_LOGF(1.0f - sample_ratio * sample_ratio * sample_ratio));

	num_iter     = HKA_MIN(max_iter, trials);
	sample_th    = (HKA_S32)MVB_ROUNDF(sample_ratio * point_num);

	for(iter = 0; iter < num_iter;)
	{

		MVBI_LineFitRansac_get_random_triplet(point_num, triplet);

		flag = MVBI_LineFitRansac_line_coe_32f(point_data[triplet[0]], point_data[triplet[1]], est_coef);
		if(flag == HKA_FALSE)
		{
			fail_num++;
			if(fail_num > num_iter)
			{
				sts = HKA_FALSE;
				break;
			}
			continue;
		}

		fail_num = 0;
		ninlier = 0;
		for(k = 0; k < point_num; k++)
		{
			dist = MVBI_LineFitRansac_line_dist_32f(est_coef, point_data[k]);

			if(dist <= dist_th)
			{
				consensus[ninlier].x = point_data[k].x;
				consensus[ninlier].y = point_data[k].y;
				ninlier++;
			}
		}

		if(ninlier >= sample_th)
		{
			// estimate circle
			_MVBI_LineFitLeastSquare_32f(consensus, ninlier, est_coef, &line_sts);
			if(line_sts == HKA_FALSE)
			{
				continue;
			}

			if((ninlier < sample_th) || (ninlier == 0))
			{
				continue;
			}

			if(ninlier > best_num)
			{
				best_num = ninlier;
				MVB_MEMCPY(best_coef, est_coef, 3 * sizeof(HKA_F32));
				sts      = HKA_TRUE;
			}

		}

		iter++;

	}

	if(sts == HKA_TRUE)
	{
		MVB_MEMCPY(line_coef, best_coef, 3 * sizeof(HKA_F32));
		*straight = MVBI_LineFitRansac_straightness_32f(consensus, ninlier,  est_coef, 1);
	}
	else
	{
		MVB_MEMSET(line_coef, 0, 3 * sizeof(HKA_F32));
		*straight = 0.0f;
	}

	return sts;
}

static HKA_S32 MVBI_LineFitRansac_process_32s(HKA_POINT_I *point_data, HKA_S32 point_num,
											  HKA_F32 line_coef[3], HKA_F32 *straight,
											  HKA_F32 sample_ratio, HKA_F32 dist_th, HKA_S32 max_iter, 
											  HKA_POINT_I *consensus)
{
	HKA_S32     trials       = 0;
	HKA_S32     num_iter     = 0;
	HKA_S32     iter         = 0;
	HKA_S32     triplet[2]   = {0};
	HKA_S32     flag         = HKA_FALSE;
	HKA_S32     ninlier      = 0;
	HKA_S32     k            = 0;
	HKA_F32     dist         = 0.0f;
	HKA_F32     dist_dev     = 0.0f;
	HKA_S32     sample_th    = 0;
	HKA_S32     line_sts     = HKA_FALSE;
	HKA_F32     best_coef[3] = {0.0f};
	HKA_F32     est_coef[3]  = {0.0f};
	HKA_S32     sts          = HKA_FALSE;
	HKA_S32     fail_num     = 0;
	HKA_S32     best_num     = 0;

	trials       = (HKA_S32)(MVB_LOGF(1.0f - 0.99f) / MVB_LOGF(1.0f - sample_ratio * sample_ratio * sample_ratio));

	num_iter     = HKA_MIN(max_iter, trials);
	sample_th    = (HKA_S32)MVB_ROUNDF(sample_ratio * point_num);

	for(iter = 0; iter < num_iter;)
	{

		MVBI_LineFitRansac_get_random_triplet(point_num, triplet);

		flag = MVBI_LineFitRansac_line_coe_32s(point_data[triplet[0]], point_data[triplet[1]], est_coef);
		if(flag == HKA_FALSE)
		{
			fail_num++;
			if(fail_num > num_iter)
			{
				sts = HKA_FALSE;
				break;
			}
			continue;
		}

		fail_num = 0;
		ninlier = 0;
		for(k = 0; k < point_num; k++)
		{
			dist = MVBI_LineFitRansac_line_dist_32s(est_coef, point_data[k]);

			if(dist_dev <= dist_th)
			{
				consensus[ninlier].x = point_data[k].x;
				consensus[ninlier].y = point_data[k].y;
				ninlier++;
			}
		}

		if(ninlier >= sample_th)
		{
			// estimate circle
			_MVBI_LineFitLeastSquare_32s(consensus, ninlier, est_coef, &line_sts);
			if(line_sts == HKA_FALSE)
			{
				continue;
			}

			if((ninlier < sample_th) || (ninlier == 0))
			{
				continue;
			}

			if(ninlier > best_num)
			{
				best_num = ninlier;
				MVB_MEMCPY(best_coef, est_coef, 3 * sizeof(HKA_F32));
				sts      = HKA_TRUE;
			}

		}

		iter++;

	}

	if(sts == HKA_TRUE)
	{
		MVB_MEMCPY(line_coef, best_coef, 3 * sizeof(HKA_F32));
		*straight = MVBI_LineFitRansac_straightness_32s(consensus, ninlier,  est_coef, 1);
	}
	else
	{
		MVB_MEMSET(line_coef, 0, 3 * sizeof(HKA_F32));
		*straight = 0.0f;
	}

	return sts;
}

/***************************************************************************************************
* 功  能：分配工作内存
* 参  数：*
*         spec      – IO 算法结构体
*         work_buf   - I  工作内存
*         point_num – I 点个数
*         work_size  - O  工作内存大小
* 返回值：HKA_VOID
* 备注：
***************************************************************************************************/
static HKA_VOID MVBI_LineFitRansac_alloc_mem_32f(MVBI_LF_RANSAC_SPEC_F32 *spec, HKA_U08 *work_buf, 
												 HKA_S32 point_num, HKA_SZT *work_size)
{
	HKA_SZT size        = 0;
	HKA_SZT used_size   = 0;
	HKA_U08 *acc_buf    = HKA_NULL;

	used_size = 0;
	acc_buf   = work_buf;

	size = point_num * sizeof(HKA_POINT_F);
	size = HKA_SIZE_ALIGN_128(size);
	spec->consensus = (HKA_POINT_F*)acc_buf;
	used_size += size;

	*work_size = used_size;
}

static HKA_VOID MVBI_LineFitRansac_alloc_mem_32s(MVBI_LF_RANSAC_SPEC_S32 *spec, HKA_U08 *work_buf, 
												 HKA_S32 point_num, HKA_SZT *work_size)
{
	HKA_SZT size        = 0;
	HKA_SZT used_size   = 0;
	HKA_U08 *acc_buf    = HKA_NULL;

	used_size = 0;
	acc_buf   = work_buf;

	size = point_num * sizeof(HKA_POINT_I);
	size = HKA_SIZE_ALIGN_128(size);
	spec->consensus = (HKA_POINT_I*)acc_buf;
	used_size += size;

	*work_size = used_size;
}

/***************************************************************************************************
* 函数名：MVBI_LineFitLeastSquare_32s MVBI_LineFitLeastSquare_32f
*         _MVBI_LineFitLeastSquare_32s _MVBI_LineFitLeastSquare_32f
* 功  能：直线拟合, 直线方程：c1 *x + c2 * y + c3 = 0
* 参  数：*
*         point_data  - I  输入点坐标
*         point_num   - I  输入点个数
*         line_coef   - O  输出直线系数，[c1, c2, c3]
*         status      - O  拟合状态
* 返回值：HKA_VOID
* 备  注：
***************************************************************************************************/
HKA_VOID _MVBI_LineFitLeastSquare_32s(HKA_POINT_I *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
									  HKA_S32 *status)
{
	HKA_S32 i    = 0;
	HKA_S32 sx   = 0;
	HKA_S32 sy   = 0;
	HKA_S32 sx2  = 0;
	HKA_S32 sy2  = 0;
	HKA_S32 sxy  = 0; 
	HKA_S32 a0   = 0;
	HKA_S32 b0   = 0;
	HKA_S32 a1   = 0;
	HKA_S32 b1   = 0;
	HKA_S32 tem  = 0; 
	HKA_F32 temp = 0.0f; 
	HKA_F32 c1   = 0.0f;
	HKA_F32 c2   = 0.0f;
	HKA_F32 c3   = 0.0f;

	for (i = 0; i < point_num; i++)
	{
		sx  += point_data[i].x;
		sy  += point_data[i].y;
		sx2 += point_data[i].x * point_data[i].x;
		sy2 += point_data[i].y * point_data[i].y;
		sxy += point_data[i].x * point_data[i].y;
	}

	a0 =   point_num * sxy - sx * sy;
	b0 = -(point_num * sx2 - sx * sx);
	a1 =   point_num * sy2 - sy * sy;
	b1 = -a0;

	if( ( HKA_ABS(a0) + HKA_ABS(b0) ) < ( HKA_ABS(a1) + HKA_ABS(b1) ) )
	{
		a0 = a1;
		b0 = b1;
	}

	tem = a0 * a0 + b0 * b0;

	if(tem == 0)
	{
		*status = HKA_FALSE;
		return;
	}

	temp = 1.0f / MVB_SQRTF(1.0f * tem);

	c1 = a0 * temp;
	c2 = b0 * temp;
	c3 = -(c1 * sx + c2 * sy) / point_num; 

	line_coef[0] = c1;
	line_coef[1] = c2;
	line_coef[2] = c3;

	*status = HKA_TRUE;
}

HKA_VOID _MVBI_LineFitLeastSquare_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
									  HKA_S32 *status)
{
	HKA_S32 i   = 0;
	HKA_F32 sx  = 0.0f;
	HKA_F32 sy  = 0.0f;
	HKA_F32 sx2 = 0.0f;
	HKA_F32 sy2 = 0.0f;
	HKA_F32 sxy = 0.0f; 
	HKA_F32 a0  = 0.0f;
	HKA_F32 b0  = 0.0f;
	HKA_F32 a1  = 0.0f;
	HKA_F32 b1  = 0.0f;
	HKA_F32 tem = 0.0f; 
	HKA_F32 c1  = 0.0f;
	HKA_F32 c2  = 0.0f;
	HKA_F32 c3  = 0.0f;

	for (i = 0; i < point_num; i++)
	{
		sx  += point_data[i].x;
		sy  += point_data[i].y;
		sx2 += point_data[i].x * point_data[i].x;
		sy2 += point_data[i].y * point_data[i].y;
		sxy += point_data[i].x * point_data[i].y;
	}

	a0 =   point_num * sxy - sx * sy;
	b0 = -(point_num * sx2 - sx * sx);
	a1 =   point_num * sy2 - sy * sy;
	b1 = -a0;

	if( ( HKA_FABS(a0) + HKA_FABS(b0) ) < ( HKA_FABS(a1) + HKA_FABS(b1) ) )
	{
		a0 = a1;
		b0 = b1;
	}

	tem = a0 * a0 + b0 * b0;

	if(tem < 1e-5f)
	{
		*status = HKA_FALSE;
		return;
	}

	tem = 1.0f / MVB_SQRTF(tem);

	c1 = a0 * tem;
	c2 = b0 * tem;
	c3 = -(c1 * sx + c2 * sy) / point_num; 

	line_coef[0] = c1;
	line_coef[1] = c2;
	line_coef[2] = c3;

	*status = HKA_TRUE;
}

HKA_STATUS MVBI_LineFitLeastSquare_32s(HKA_POINT_I *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                       HKA_S32 *status)
{
	HKA_CHECK_ERROR(HKA_NULL == point_data, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == status, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(point_num < 2, HKA_STS_ERR_NULL_PTR);

	_MVBI_LineFitLeastSquare_32s(point_data, point_num, line_coef, status);
	
	return HKA_STS_OK;
}
HKA_STATUS MVBI_LineFitLeastSquare_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                       HKA_S32 *status)
{
	HKA_CHECK_ERROR(HKA_NULL == point_data, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == status, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(point_num < 2, HKA_STS_ERR_NULL_PTR);

	_MVBI_LineFitLeastSquare_32f(point_data, point_num, line_coef, status);

	return HKA_STS_OK;
}

/***************************************************************************************************
* 函数名：MVBI_LineFitRansacBufferSize_32s MVBI_LineFitRansacBufferSize_32f
*         _MVBI_LineFitRansacBufferSize_32s _MVBI_LineFitRansacBufferSize_32f
* 功  能：获取基于ransac的直线拟合工作内存大小
* 参  数：*
*         max_point_num - I  最大点个数
*         work_size     - O  工作内存大小
* 返回值：HKA_VOID
* 备  注：
***************************************************************************************************/
HKA_VOID _MVBI_LineFitRansacBufferSize_32s(HKA_S32 max_point_num, HKA_SZT *work_size)
{
	MVBI_LF_RANSAC_SPEC_S32 spec     = {0};
	HKA_U08                *work_buf = HKA_NULL;

	work_buf = (HKA_U08*)&spec;

	MVBI_LineFitRansac_alloc_mem_32s(&spec, work_buf, max_point_num, work_size);
}

HKA_VOID _MVBI_LineFitRansacBufferSize_32f(HKA_S32 max_point_num, HKA_SZT *work_size)
{
	MVBI_LF_RANSAC_SPEC_F32 spec     = {0};
	HKA_U08                *work_buf = HKA_NULL;

	work_buf = (HKA_U08*)&spec;

	MVBI_LineFitRansac_alloc_mem_32f(&spec, work_buf, max_point_num, work_size);
}

HKA_STATUS MVBI_LineFitRansacBufferSize_32s(HKA_S32 max_point_num, HKA_SZT *work_size)
{
	HKA_SZT size = 0;

	HKA_CHECK_ERROR(HKA_NULL == work_size, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(max_point_num < 2, HKA_STS_ERR_BAD_ARG);

	_MVBI_LineFitRansacBufferSize_32s(max_point_num, &size);

	*work_size = size;

	return HKA_STS_OK;
}

HKA_STATUS MVBI_LineFitRansacBufferSize_32f(HKA_S32 max_point_num, HKA_SZT *work_size)
{
	HKA_SZT size = 0;

	HKA_CHECK_ERROR(HKA_NULL == work_size, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(max_point_num < 2, HKA_STS_ERR_BAD_ARG);

	_MVBI_LineFitRansacBufferSize_32f(max_point_num, &size);

	*work_size = size;

	return HKA_STS_OK;
}

/***************************************************************************************************
* 函数名：MVBI_LineFitRansac_32s MVBI_LineFitRansac_32f
*         _MVBI_LineFitRansac_32s _MVBI_LineFitRansac_32f
* 功  能：获取基于ransac的直线拟合工作内存大小
* 参  数：*
*         point_data  - I  输入点坐标
*         point_num   - I  输入点个数
*         line_coef   - O  输出直线系数，[c1, c2, c3]
*         status      - O  拟合状态
*         straight    - O  直线度
*         sample_th   - I  采用阈值
*         dist_th     - I  距离阈值
*         max_iter    - I  最大迭代次数
*         work_buf    - I  工作内存
* 返回值：HKA_VOID
* 备  注：
***************************************************************************************************/
HKA_VOID _MVBI_LineFitRansac_32s(HKA_POINT_I *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                 HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                 HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf)
{
	HKA_SZT                 work_size = 0;
	MVBI_LF_RANSAC_SPEC_S32 spec      = {0};
	HKA_POINT_I             *consensus = HKA_NULL;

	MVBI_LineFitRansac_alloc_mem_32s(&spec, work_buf, point_num, &work_size);

	consensus = spec.consensus;

	*status = MVBI_LineFitRansac_process_32s(point_data, point_num, line_coef, straight, sample_ratio, 
		                                     dist_th, max_iter, consensus);


}

HKA_VOID _MVBI_LineFitRansac_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                 HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                 HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf)
{
	HKA_SZT                 work_size  = 0;
	MVBI_LF_RANSAC_SPEC_F32 spec       = {0};
	HKA_POINT_F             *consensus = HKA_NULL;

	MVBI_LineFitRansac_alloc_mem_32f(&spec, work_buf, point_num, &work_size);

	consensus = spec.consensus;

	*status = MVBI_LineFitRansac_process_32f(point_data, point_num, line_coef, straight, sample_ratio, 
		                                     dist_th, max_iter, consensus);

}

HKA_STATUS MVBI_LineFitRansac_32s(HKA_POINT_I *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                  HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                  HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf)
{
	HKA_CHECK_ERROR(HKA_NULL == point_data, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == status, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == straight, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == work_buf, HKA_STS_ERR_NULL_PTR);

	HKA_CHECK_ERROR((sample_ratio < 0.0f) || (sample_ratio > 1.0f), HKA_STS_ERR_BAD_ARG);
	HKA_CHECK_ERROR(dist_th < 0.0f, HKA_STS_ERR_BAD_ARG);
	HKA_CHECK_ERROR(max_iter < 1, HKA_STS_ERR_BAD_ARG);

	_MVBI_LineFitRansac_32s(point_data, point_num, line_coef, status, straight, sample_ratio, 
		                    dist_th, max_iter, work_buf);

	return HKA_STS_OK;
}

HKA_STATUS MVBI_LineFitRansac_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                  HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                  HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf)
{
	HKA_CHECK_ERROR(HKA_NULL == point_data, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == status, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == straight, HKA_STS_ERR_NULL_PTR);
	HKA_CHECK_ERROR(HKA_NULL == work_buf, HKA_STS_ERR_NULL_PTR);

	HKA_CHECK_ERROR((sample_ratio < 0.0f) || (sample_ratio > 1.0f), HKA_STS_ERR_BAD_ARG);
	HKA_CHECK_ERROR(dist_th < 0.0f, HKA_STS_ERR_BAD_ARG);
	HKA_CHECK_ERROR(max_iter < 1, HKA_STS_ERR_BAD_ARG);

	_MVBI_LineFitRansac_32f(point_data, point_num, line_coef, status, straight, sample_ratio, 
		                    dist_th, max_iter, work_buf);

	return HKA_STS_OK;
}

