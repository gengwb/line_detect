/***************************************************************************************************
* 版权信息：Copyright (c) 2016, 杭州海康威视数字技术股份有限公司
* 
* 文件名称：mva_line_detect.c
* 摘    要：直线检测实现
*
* 当前版本：V0.1.4
* 作    者：郑俊君
* 日    期：2016-01-31
* 备    注：1、增加了排除比例参数
*
* 历史版本：V0.1.3
* 作    者：郑俊君
* 日    期：2016-1-28
* 备    注：1、增加了投影宽度
*           2、修改ransac实现方式
*           3、删除了聚类
*           4、更改核的计算方式，改为DOG（different of gauss）
*           5、提升算法的稳定性
*
* 历史版本：V0.1.2
* 作    者：邓志辉
* 日    期：2016-1-12
* 备    注：1、修改算法参数命名，以便更加容易理解
*           2、增加容忍角度
*
* 历史版本：V0.1.1
* 作    者：郑俊君
* 日    期：2015-12-22
* 备    注：1、重构旋转矩形输入形式
*           2、加入搜索傅里叶核，提高检测稳定性
*
* 历史版本：V0.1.0
* 作    者：郑俊君
* 日    期：2015-11-13
* 备    注：1.创建    
***************************************************************************************************/
#include "hka_types.h"
#include "hka_defs.h"
#include "mvb_defs.h"
#include "mvb_types.h"
#include "mvageofind_line_detect.h"
#include "mvageofind_line_detect_type.h"
#include "mvbi_line_fit.h"
#include "_mvbs_filter.h"
#include "_mvbs_sort.h"
#include "_mvbi_data.h"

/***************************************************************************************************
* 功  能：用两个点计算直线参数
* 参  数：* 
*         p1            -I         点1
*         p2            -I         点2
*         line_coe      -O         直线参数
* 返回值：无
* 备  注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_line_coe(HKA_POINT_F p1, HKA_POINT_F p2, MVAGEOFIND_LD_LINE_COE *line_coe)
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
	if (temp == 0.0f)
	{
		temp = 1e-4f;
	}

	tem = 1.0f / MVB_SQRTF(temp);

	line_coe->a = a0 * tem;
	line_coe->b = b0 * tem;
	line_coe->c = c0 * tem;
}

/***************************************************************************************************
* 功  能：带掩膜的用最小二乘法拟合多个点的直线系数
* 参  数： 
*         pts         -I  点集
*         pts_num     -I  点的个数
*         mask        -I  点是否参加计算的掩膜
*         line_coe    -O  拟合出的直线系数
* 返回值：
* 备  注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_line_coe_ls_mask(HKA_POINT_F *pts, HKA_S32 pts_num, HKA_U08 *mask,
												       MVAGEOFIND_LD_LINE_COE *line_coe)
{
	HKA_S32 i       = 0;
	HKA_S32 cur_num = 0;
	HKA_F32 sx      = 0.0f;
	HKA_F32 sy      = 0.0f;
	HKA_F32 sx2     = 0.0f;
	HKA_F32 sy2     = 0.0f;
	HKA_F32 sxy     = 0.0f;
	HKA_F32 a0      = 0.0f;
	HKA_F32 temp    = 0.0f;
	HKA_F32 b0      = 0.0f;
	HKA_F32 a1      = 0.0f;
	HKA_F32 b1      = 0.0f;
	HKA_F32 tem     = 0.0f;

	for (i = 0; i < pts_num; i++)
	{
		if (mask[i])
		{
			sx += pts[i].x;
			sy += pts[i].y;
			sx2 += pts[i].x * pts[i].x;
			sy2 += pts[i].y * pts[i].y;
			sxy += pts[i].x * pts[i].y;
			cur_num++;
		}
	}

	if (cur_num < 2)
	{
		return ;
	}

	a0 = cur_num * sxy - sx * sy;
	b0 = -(cur_num * sx2 - sx * sx);
	a1 = cur_num * sy2 - sy * sy;
	b1 = -a0;

	if (HKA_FABS(a0) + HKA_FABS(b0) < HKA_FABS(a1) + HKA_FABS(b1))
	{
		a0 = a1;
		b0 = b1;
	}

	temp = a0 * a0 + b0 * b0;

	if (0 == temp)
	{
		temp = 1e-4f;
	}

	tem = 1.0f / MVB_SQRTF(temp);

	line_coe->a = a0 * tem;
	line_coe->b = b0 * tem;
	line_coe->c = -(line_coe->a * sx + line_coe->b * sy) / cur_num;
}

/***************************************************************************************************
* 功  能：产生随机数
* 参  数： 
*         nums       -O  随机数
*         limit      -I  随机数的最大值
*         len        -I  随机数个数
*         seed       -I  随机数种子
* 返回值：
* 备  注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_rand_int_nums(HKA_S32 *nums, HKA_S32 limit, 
											        HKA_S32 len, HKA_U32 seed)
{
	HKA_S32 i       = 0;
	HKA_S32 tem     = 0;
	HKA_S32 cur_len = 0;

	//产生随机数种子
	seed += (HKA_U32)time(0);
	MVB_SRAND(seed);

	cur_len = 1;
	nums[0] = MVB_RAND() % limit;
	while (cur_len < len)
	{
		tem = MVB_RAND() % limit;

		for (i = 0; i < cur_len; i++)
		{
			if (nums[i] == tem)
			{
				break;
			}
		}

		if (i >= cur_len)
		{
			nums[cur_len] = tem;
			cur_len++;
		} 
	}
}

/***************************************************************************************************
* 功  能：排除与直线的距离大于设定值的点
* 参  数： 
*         pts          -I  点集
*         pts_num      -I  点的个数
*         dist_thres   -I  距离阈值
*         line_coe     -I  直线系数
*         mask         -O  掩膜，距离小于dist_thres的掩膜值为1；否则为0
* 返回值：
* 备  注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_exclu_distant_pts(HKA_POINT_F *pts, HKA_S32 pts_num, HKA_F32 dist_thres, 
												        HKA_F32 reject_ratio, MVAGEOFIND_LD_LINE_COE *line_coe,
														HKA_U08 *mask, HKA_U08 *buf)								
{
	HKA_S32 i         = 0;
	HKA_S32 count     = 0;
	HKA_F32 dist      = 0.0f;
	HKA_F32 den       = 0.0f;
	HKA_F32 dist_th   = 0.0f;
	HKA_F32 *dist_val = HKA_NULL;
	HKA_F32 *sort_val = HKA_NULL;
	HKA_SIZE_I roi    = {0};


	dist_val = (HKA_F32 *)buf;
	sort_val = (HKA_F32 *)(buf + pts_num * sizeof(HKA_F32));

	den = line_coe->a * line_coe->a + line_coe->b * line_coe->b;
	den = 1.0f / den;

	for (i = 0; i < pts_num; i++)
	{
		dist = HKA_FABS(line_coe->a * pts[i].x + line_coe->b * pts[i].y + line_coe->c);
		dist = (dist * dist * den);

		dist_val[i] = dist;
	}

	_MVBS_SortAscend_32f(dist_val, sort_val, pts_num);

	count = MVB_ROUNDF(pts_num * (1 - reject_ratio));

	dist_th = sort_val[count];

	dist_th = HKA_MIN(dist_th, dist_thres);

	roi.width  = pts_num; 
	roi.height = 1;

	_MVBI_Set_8u_C1R(1, mask, pts_num, roi);

	for (i = 0; i < pts_num; i++ )
	{
		if (dist_val[i] > dist_th )
		{
			mask[i] = 0;
		}
	}
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
static HKA_F32 MVAGEOFIND_LineDetect_compute_straightness(HKA_POINT_F *pts, HKA_S32 pts_num,  
														  HKA_U08 *mask, MVAGEOFIND_LD_LINE_COE *line_coe,
														  HKA_S32 kernel_size)
{
	HKA_S32 num  = 0;
	HKA_S32 i    = 0;
	HKA_S32 PR   = 0;
	HKA_F32 dist = 0.0f;
	HKA_F32 msd  = 0.0f;
	HKA_F32 den  = 0.0f; // 归一化阈值

	if((line_coe->a != 0) || (line_coe->b != 0) )
	{
		den = line_coe->a * line_coe->a + line_coe->b * line_coe->b;
		den = 1.0f / den;
	}

	for (i = 0; i < pts_num; i++)
	{
		if (mask[i])
		{
			dist = line_coe->a * pts[i].x + line_coe->b * pts[i].y + line_coe->c;
			msd += (dist * dist);
			num ++;
		}	
	}
	msd *= den;
	if (num == 0)
	{
		return -1.0f;
	}
	else
	{
		PR = num * kernel_size * kernel_size;
		msd = (1.0f - msd / PR) * 100.f;
		return  msd;
	}
}

/***************************************************************************************************
* 功  能：用ransac方法计算直线系数
* 参  数： 
*         pts           -I  点集
*         pts_num       -I  点的个数
*         max_iter      -I  最大计算次数
*         percent       -I  百分比阈值
*         dist_thres    -I  距离阈值
*         line_coe      -O  直线系数
*         mask          -I  掩膜数据
* 返回值：0：得到满足条件的解；1:没有得到满足条件的解；-1：没有得到解
* 备  注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_ransac_line_coe(HKA_POINT_F *pts, HKA_S32 pts_num, HKA_S32 max_iter,
											         HKA_F32 percent, HKA_F32 dist_thres, HKA_F32 reject_ratio, 
													 HKA_S32 kernel_size, MVAGEOFIND_LD_LINE_COE *line_coe,
											         HKA_F32 *straightness, HKA_U08 *work_buf)										 
{	
	HKA_S32                i           = 0;
	HKA_S32                j           = 0;
	HKA_S32                cur_num     = 0;
	HKA_S32                tar_num     = 0;
	HKA_S32                cur_tar     = 0;
	HKA_F32                den         = 0.0f;
	HKA_F32                cur_dist    = 0.0f;
	HKA_S32                cur_pos[2]  = {0, 0};
	MVAGEOFIND_LD_LINE_COE cur_line    = {0};
	MVAGEOFIND_LD_LINE_COE result_line = {0};
	HKA_U08                *mask       = HKA_NULL;
	HKA_U08                *buf        = HKA_NULL;

	mask = work_buf;
	buf  = (work_buf + pts_num * sizeof(HKA_U08));

	tar_num = MVB_ROUNDF(pts_num * percent);

	for (i = 0; i < max_iter; i++)
	{
		//产生两个点的随机下标，并计算由这两个点所在的直线系数
		j = 0;
		while(j < MVAGEOFIND_LD_RANSAC_RAND_TIMES)
		{
			MVAGEOFIND_LineDetect_rand_int_nums(cur_pos, pts_num, 2, j * 21 + i * 11);

			if (pts[cur_pos[0]].x != pts[cur_pos[1]].x || 
				pts[cur_pos[0]].y != pts[cur_pos[1]].y)
			{
				break;
			}
			j++;
		}
		if (j >= MVAGEOFIND_LD_RANSAC_RAND_TIMES)
		{
			continue;
		}

		MVAGEOFIND_LineDetect_line_coe(pts[cur_pos[0]], pts[cur_pos[1]], &cur_line);

		//逐点判断到计算直线的距离
		cur_num = 0;
		den = cur_line.a * cur_line.a + cur_line.b * cur_line.b;
		den = 1.0f / den;

		for (j = 0; j < pts_num; j++)
		{
			cur_dist = HKA_FABS(cur_line.a * pts[j].x + cur_line.b * pts[j].y + cur_line.c);
			cur_dist = (cur_dist * cur_dist * den);
			if (cur_dist < dist_thres)
			{
				cur_num++;
				mask[j] = 1;
			}
			else
			{
				mask[j] = 0;
			}
		}

		MVAGEOFIND_LineDetect_line_coe_ls_mask(pts, pts_num, mask, line_coe);
		// 		den = line_coe->a * line_coe->a + line_coe->b * line_coe->b;
		// 		den = 1.0f / den;
		// 		sum_dist = 0.0f;
		// 
		// 		for (j = 0; j < pts_num; j++)
		// 		{
		// 			cur_dist = HKA_FABS(line_coe->a * pts[j].x + line_coe->b * pts[j].y + line_coe->c);
		// 			sum_dist += (cur_dist * cur_dist);
		// 			
		// 		}
		// 		sum_dist *= den;
		if (cur_num > cur_tar)
		{
			cur_tar = cur_num;
			//best_dist = sum_dist;
			result_line = *line_coe;
		}
		if (cur_num > tar_num)
		{
			break;
		}
	}


	//判断是否得到满足条件的解
	if (i > max_iter)
	{
		return HKA_FALSE;
	}
	//MVAGEOFIND_LineDetect_line_coe_ls_mask(pts, pts_num, mask, line_coe);


	MVAGEOFIND_LineDetect_exclu_distant_pts(pts, pts_num, dist_thres, reject_ratio, &result_line, mask, buf);

	//再用最小二乘法计算满足条件的点的直线系数
	MVAGEOFIND_LineDetect_line_coe_ls_mask(pts, pts_num, mask, line_coe);

	// 计算直线度
	*straightness = MVAGEOFIND_LineDetect_compute_straightness(pts, pts_num, mask, line_coe, kernel_size);

	//得到满足条件的解
	return HKA_TRUE;
}
/***************************************************************************************************
* 功  能：生成dog
* 参  数：*
*         gauss_kernel  - I  高斯核
*         dog_kernel    - O  dog核
*         kernel_size   - I  核大小

* 返回值：区域像素均值
* 备注：
***************************************************************************************************/
static HKA_VOID  MVAGEOFIND_GenDogKernel_db_32f(HKA_F32 *gauss_kernel, HKA_F32 *dog_kernel, HKA_S32 kernel_size)
{
	HKA_S32 i        = 0;
	HKA_F32 neg_vals = 0.0f;
	HKA_F32 pos_vals = 0.0f;
	HKA_F32 val      = 0.0f;

	if (kernel_size == 1)
	{
		dog_kernel[0] = 0.0f;
	}

	if (kernel_size > 1)
	{
		dog_kernel[0] = gauss_kernel[1] - gauss_kernel[0];
		dog_kernel[kernel_size - 1] = gauss_kernel[kernel_size - 1] - gauss_kernel[kernel_size - 2];
	}

	if (kernel_size > 2)
	{
		for (i = 1; i < kernel_size - 1; i++)
		{
			dog_kernel[i] = 0.5f * (gauss_kernel[i + 1] - gauss_kernel[i - 1]);
		}
	}

	// 归一化
	for (i = 0; i < kernel_size; i++)
	{
		val = dog_kernel[i];

		if (val > 0.0f)
		{
			pos_vals += val;
		}
		if (val < 0.0f)
		{
			neg_vals += val;
		}
	}

	neg_vals = (neg_vals < 0.0f) ? 1.0f / HKA_FABS(neg_vals) : 1.0f;

	if (pos_vals > 0.0f)
	{
		pos_vals = 1.0f / pos_vals;
	}

	for (i = 0; i < kernel_size; i++)
	{
		val = dog_kernel[i];

		if (val > 0.0f)
		{
			dog_kernel[i] *= pos_vals;
		}
		if (val < 0.0f)
		{
			dog_kernel[i] *= neg_vals;
		}
	}

	for (i = 0; i < kernel_size; i ++)
	{
		dog_kernel[i] *= -1.0f;
	}

}

/***************************************************************************************************
* 功  能：生成dog
* 参  数：*
*         gauss_kernel  - I  高斯核
*         dog_kernel    - O  dog核
*         kernel_size   - I  核大小

* 返回值：区域像素均值
* 备注：
***************************************************************************************************/
static HKA_VOID  MVAGEOFIND_GenDogKernel_bd_32f(HKA_F32 *gauss_kernel, HKA_F32 *dog_kernel, HKA_S32 kernel_size)
{
	HKA_S32 i        = 0;
	HKA_F32 neg_vals = 0.0f;
	HKA_F32 pos_vals = 0.0f;
	HKA_F32 val      = 0.0f;

	if (kernel_size == 1)
	{
		dog_kernel[0] = 0.0f;
	}

	if (kernel_size > 1)
	{
		dog_kernel[0] = gauss_kernel[1] - gauss_kernel[0];
		dog_kernel[kernel_size - 1] = gauss_kernel[kernel_size - 1] - gauss_kernel[kernel_size - 2];
	}

	if (kernel_size > 2)
	{
		for (i = 1; i < kernel_size - 1; i++)
		{
			dog_kernel[i] = 0.5f * (gauss_kernel[i + 1] - gauss_kernel[i - 1]);
		}
	}

	// 归一化
	for (i = 0; i < kernel_size; i++)
	{
		val = dog_kernel[i];

		if (val > 0.0f)
		{
			pos_vals += val;
		}
		if (val < 0.0f)
		{
			neg_vals += val;
		}
	}

	neg_vals = (neg_vals < 0.0f) ? 1.0f / HKA_FABS(neg_vals) : 1.0f;

	if (pos_vals > 0.0f)
	{
		pos_vals = 1.0f / pos_vals;
	}

	for (i = 0; i < kernel_size; i++)
	{
		val = dog_kernel[i];

		if (val > 0.0f)
		{
			dog_kernel[i] *= pos_vals;
		}
		if (val < 0.0f)
		{
			dog_kernel[i] *= neg_vals;
		}
	}

}
/***************************************************************************************************
* 功  能：初始化核大小,从黑到白
* 参  数：*
*         kernel_size   - I  核大小
*         kernel         - O  核
* 返回值：区域像素均值
* 备注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_init_kernel_db(HKA_S32 kernel_size, HKA_F32 *kernel)
{
	HKA_F32 gauss_kernel[100] = {0.0f};
	HKA_F32 sigma             = 0.0f;

	sigma = 1.0f * kernel_size / 6.0f;

	_MVBS_GenGaussKernel_32f(gauss_kernel, kernel_size, sigma);

	MVAGEOFIND_GenDogKernel_db_32f(gauss_kernel, kernel, kernel_size);


}
/***************************************************************************************************
* 功  能：初始化核大小,从白到黑
* 参  数：*
*         kernel_size   - I  核大小
*         kernel         - O  核
* 返回值：区域像素均值
* 备注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_init_kernel_bd(HKA_S32 kernel_size, HKA_F32 *kernel)
{
	HKA_F32 gauss_kernel[100] = {0.0f};
	HKA_F32 sigma             = 0.0f;

	sigma = 1.0f * kernel_size / 6.0f;

	_MVBS_GenGaussKernel_32f(gauss_kernel, kernel_size, sigma);

	MVAGEOFIND_GenDogKernel_bd_32f(gauss_kernel, kernel, kernel_size);
	
}

/***************************************************************************************************
* 功  能：计算区域均值
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         point          - I  区域中心点
*         sx_step        - I  搜索方向x的变化步长
*         sy_step        - I  搜索方向y的变化步长
*         bx_step        - I  领域方向x的变化步长
*         bx_step        - I  领域方向y的变化步长
*         kernel_size    - I  核大小
* 返回值：区域像素均值
* 备注：
***************************************************************************************************/
static HKA_F32 MVAGEOFIND_LineDetect_edge_search(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
												 HKA_POINT_F point, HKA_F32 sx_step, HKA_F32 sy_step, 
												 HKA_F32 bx_step, HKA_F32 by_step, HKA_S32 region_width,
												 HKA_S32 edge_width,  HKA_F32 *kernel)									
{
	HKA_S32     i      = 0;
	HKA_S32     j      = 0;
	HKA_S32     x      = 0;
	HKA_S32     y      = 0;
	HKA_S32     width  = 0;
	HKA_S32     height = 0;
	HKA_S32     half1  = 0;
	HKA_S32     half2  = 0;
	HKA_F32     val    = 0.0f;
	HKA_F32     r_val  = 0.0f;
	HKA_F32     scale  = 0.0f;
	HKA_POINT_F p_s    = {0.0f};
	HKA_POINT_F p_d    = {0.0f};

	half1 = edge_width >> 1;
	half2 = region_width >> 1;

	scale = 1.0f / region_width;

	p_s.x = point.x - half1 * sx_step;
	p_s.y = point.y - half1 * sy_step;

	width  = roi_size.width - 1;
	height = roi_size.height - 1;

	for (i = 0; i < edge_width; i ++)
	{	
		r_val = 0.0f;
		for (j = 0; j < region_width; j++)
		{
			p_d.x = p_s.x - half2 * bx_step;
			p_d.y = p_s.y - half2 * by_step;

			x  = HKA_MIN(HKA_MAX(MVB_ROUNDF(p_d.x), 0), width);
			y  = HKA_MIN(HKA_MAX(MVB_ROUNDF(p_d.y), 0), height);
			
			r_val += src[y * src_step + x];

			p_d.x += bx_step;
			p_d.y += by_step;
			
		}
		r_val = r_val * scale;

		val += r_val * kernel[i];

		p_s.x += sx_step;
		p_s.y += sy_step;	
	}

	return val;
}

/***************************************************************************************************
* 搜索沿着宽的方向最优一条边
***************************************************************************************************/
/***************************************************************************************************
* 功  能：搜索最优一条边
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         ray_num        - I  射线的数目
*         edge_strength  - I  边缘的强度
*         edge_points    - O  边缘点集
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_w_db(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i         = 0;
	HKA_F32     angle     = 0;
	HKA_S32     edge_num  = 0;
	HKA_S32     count     = 0;
	HKA_S32     height    = 0;
	HKA_F32     val1      = 0.0f;
	HKA_F32     best_val  = 0.0f;
	HKA_F32     move_step = 0.0f;
	HKA_F32     w_xstep   = 0.0f;
	HKA_F32     w_ystep   = 0.0f;
	HKA_F32     w1_xstep   = 0.0f;
	HKA_F32     w1_ystep   = 0.0f;
	HKA_F32     h_xstep   = 0.0f;
	HKA_F32     h_ystep   = 0.0f;
	HKA_POINT_F point_s   = {0};
	HKA_POINT_F point     = {0};
	HKA_POINT_F point_a   = {0};
	HKA_F32     *kernel   = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break; 
		}
		point = point_s;

		best_val = 1.0f * edge_strength;
		for (count = 0; count < height; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
		           h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);
			if (val1 > best_val)
			{
				best_val = val1;
				point_a  = point;
			}
			point.x += h_xstep;
			point.y += h_ystep;
	
		}

		if (best_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			/*val[edge_num]         = src[MVB_ROUNDF(point_a.y) * src_step + MVB_ROUNDF(point_a.x)];*/
			//val[edge_num]         = best_val;
			edge_num ++;
		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_w_bd(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i         = 0;
	HKA_F32     angle     = 0;
	HKA_S32     edge_num  = 0;
	HKA_S32     count     = 0;
	HKA_S32     height    = 0;
	HKA_F32     val1      = 0.0f;
	HKA_F32     best_val  = 0.0f;
	HKA_F32     move_step = 0.0f;
	HKA_F32     w_xstep   = 0.0f;
	HKA_F32     w_ystep   = 0.0f;
	HKA_F32     w1_xstep   = 0.0f;
	HKA_F32     w1_ystep   = 0.0f;
	HKA_F32     h_xstep   = 0.0f;
	HKA_F32     h_ystep   = 0.0f;
	HKA_POINT_F point_s   = {0};
	HKA_POINT_F point     = {0};
	HKA_POINT_F point_a   = {0};
	HKA_F32     *kernel   = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break; 
		}
		point = point_s;

		best_val = 1.0f * edge_strength;
		for (count = 0; count < height; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);
			if (val1 > best_val)
			{
				best_val = val1;
				point_a  = point;
			}
			point.x += h_xstep;
			point.y += h_ystep;

		}

		if (best_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
/*			val[edge_num]         = src[MVB_ROUNDF(point_a.y) * src_step + MVB_ROUNDF(point_a.x)];*/
			//val[edge_num]         = best_val;
			edge_num ++;
		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}
	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_w_all(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i         = 0;
	HKA_F32     angle     = 0;
	HKA_S32     edge_num  = 0;
	HKA_S32     count     = 0;
	HKA_S32     height    = 0;
	HKA_F32     val1      = 0.0f;
	HKA_F32     best_val  = 0.0f;
	HKA_F32     move_step = 0.0f;
	HKA_F32     w_xstep   = 0.0f;
	HKA_F32     w_ystep   = 0.0f;
	HKA_F32     w1_xstep  = 0.0f;
	HKA_F32     w1_ystep  = 0.0f;
	HKA_F32     h_xstep   = 0.0f;
	HKA_F32     h_ystep   = 0.0f;
	HKA_POINT_F point_s   = {0};
	HKA_POINT_F point     = {0};
	HKA_POINT_F point_a   = {0};
	HKA_F32     *kernel   = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break; 
		}

		point = point_s;

		best_val = 1.0f * edge_strength;

		for (count = 0; count < height; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);
			val1 = HKA_FABS(val1);

			if (val1 > best_val)
			{
				best_val = val1;
				point_a  = point;
			}
			point.x += h_xstep;
			point.y += h_ystep;

		}

		if (best_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;
}

static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_w(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
												        MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
												        MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32     edge_num      = 0;
	HKA_S32     ray_num       = 0;
	HKA_S32     edge_nature   = 0;
	HKA_S32     kernel_size   = 0;
	HKA_S32     edge_strength = 0;
	HKA_S32     region_width  = 0;
	HKA_U08     *buf          = HKA_NULL;
	HKA_POINT_F *edge_points  = HKA_NULL;

	ray_num       = cfg->ray_num;
	// 确保核的大小为奇数
	kernel_size   = ((cfg->kernel_size >> 1) << 1) + 1;
	region_width  = ((cfg->region_width >> 1) << 1) + 1;
	edge_nature   = cfg->edge_polarity;
	edge_strength = cfg->edge_strength;
	edge_points   = spec->edge_points;
	buf           = spec->buf;

	switch (edge_nature)
	{
	case MVBI_EDGE_TYPES_DARK_TO_BRIGHT:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_w_db(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BRIGHT_TO_DARK:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_w_bd(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BOTH:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_w_all(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 搜索沿着宽的方向最后一条边
***************************************************************************************************/
/***************************************************************************************************
* 功  能：搜索最后一条边
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         ray_num        - I  射线的数目
*         edge_strength  - I  边缘的强度
*         edge_points    - O  边缘点集
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_w_db(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)													
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     height      = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     last_val    = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     w1_xstep    = 0.0f;
	HKA_F32     w1_ystep    = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel      = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < height; count ++ )
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);

			if (val1 > edge_strength)
			{
				last_val = val1;
				point_a  = point;
			}
			point.x += h_xstep;
			point.y += h_ystep;
			
		}

		if (last_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_w_bd(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     height      = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     last_val    = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     w1_xstep    = 0.0f;
	HKA_F32     w1_ystep    = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = move_step * MVB_COSF(angle);
	w1_ystep = move_step * MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < height; count ++ )
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);

			if (val1 > edge_strength)
			{
				last_val = val1;
				point_a  = point;
			}

			point.x += h_xstep;
			point.y += h_ystep;

		}

		if (last_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_w_all(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     height      = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     last_val    = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     w1_xstep    = 0.0f;
	HKA_F32     w1_ystep    = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < height; count ++ )
		{

			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);
			val1 = HKA_FABS(val1);
			if (val1 > edge_strength)
			{
				last_val = val1;
				point_a  = point;
			}
			point.x += h_xstep;
			point.y += h_ystep;

		}

		if (last_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_w(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
												        MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
												        MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32     edge_num      = 0;
	HKA_S32     ray_num       = 0;
	HKA_S32     kernel_size   = 0;
	HKA_S32     edge_nature   = 0;
	HKA_S32     edge_strength = 0;
	HKA_S32     region_width  = 0;
	HKA_U08     *buf          = HKA_NULL;
	HKA_POINT_F *edge_points  = HKA_NULL;

	ray_num       = cfg->ray_num;
	kernel_size   = ((cfg->kernel_size >> 1) << 1) + 1;
	region_width  = ((cfg->region_width >> 1) << 1) + 1;
	edge_nature   = cfg->edge_polarity;
	edge_strength = cfg->edge_strength;
	edge_points   = spec->edge_points;
	buf           = spec->buf;

	switch (edge_nature)
	{
	case MVBI_EDGE_TYPES_DARK_TO_BRIGHT:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_w_db(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BRIGHT_TO_DARK:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_w_bd(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BOTH:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_w_all(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 搜索沿着宽的方向第一条边
***************************************************************************************************/
/***************************************************************************************************
* 功  能：搜索第一条边
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         ray_num        - I  射线的数目
*         edge_strength  - I  边缘的强度
*         edge_points    - O  边缘点集
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_w_db(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     height      = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     w1_xstep    = 0.0f;
	HKA_F32     w1_ystep    = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}
		point = point_s;
		
		for (count = 0; count < height; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);
		
			if (val1 >= edge_strength)
			{
				edge_points[edge_num] = point;
				edge_num ++;
				break;
			}

			point.x += h_xstep;
			point.y += h_ystep;
			
		}
		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_w_bd(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     height      = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     w1_xstep    = 0.0f;
	HKA_F32     w1_ystep    = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < height; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);

			if (val1 >= edge_strength)
			{
				edge_points[edge_num] = point;
				edge_num ++;
				break;
			}

			point.x += h_xstep;
			point.y += h_ystep;

		}
		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	

}

static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_w_all(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															 MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															 HKA_S32 edge_strength, HKA_S32 kernel_size, 
															 HKA_S32 region_width, HKA_POINT_F *edge_points,
															 HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     height      = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     w1_xstep    = 0.0f;
	HKA_F32     w1_ystep    = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->width / (ray_num - 1);

	h_xstep = MVB_SINF(-angle);
	h_ystep = MVB_COSF(-angle);

	w1_xstep = MVB_COSF(angle);
	w1_ystep = MVB_SINF(angle);

	w_xstep = move_step * MVB_COSF(angle);
	w_ystep = move_step * MVB_SINF(angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	height = rect->height - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < height; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				h_xstep, h_ystep, w1_xstep, w1_ystep, region_width, kernel_size, kernel);
			val1 = HKA_FABS(val1);

			if (val1 >= edge_strength)
			{
				edge_points[edge_num] = point;
				edge_num ++;
				break;
			}

			point.x += h_xstep;
			point.y += h_ystep;

		}

		point_s.x = point_s.x + w_xstep;
		point_s.y = point_s.y + w_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_w(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
												         MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
												         MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32     edge_num      = 0;
	HKA_S32     ray_num       = 0;
	HKA_S32     edge_nature   = 0;
	HKA_S32     edge_strength = 0;
	HKA_S32     kernel_size   = 0;
	HKA_S32     region_width  = 0;
	HKA_U08     *buf          = HKA_NULL;
	HKA_POINT_F *edge_points  = HKA_NULL;

	ray_num       = cfg->ray_num;
	kernel_size   = ((cfg->kernel_size >> 1)<< 1) + 1;
	region_width  = ((cfg->region_width >> 1)<< 1) + 1;
	edge_nature   = cfg->edge_polarity;
	edge_strength = cfg->edge_strength;
	edge_points   = spec->edge_points;
	buf           = spec->buf;

	switch (edge_nature)
	{
	case MVBI_EDGE_TYPES_DARK_TO_BRIGHT:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_w_db(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BRIGHT_TO_DARK:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_w_bd(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BOTH:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_w_all(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 搜索沿着高的方向最优一条边
***************************************************************************************************/
/***************************************************************************************************
* 功  能：搜索最优一条边
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         ray_num        - I  射线的数目
*         edge_strength  - I  边缘的强度
*         edge_points    - O  边缘点集
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_h_db(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i           = 0;	
	HKA_F32     angle       = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     best_val    = 0.0f;
	HKA_F32     move_step   = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		best_val = 1.0f * edge_strength;
		for (count = 0; count < width; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				   w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);

			if (val1 > best_val)
			{
				best_val = val1;
				point_a = point;
			}

			point.x += w_xstep;
			point.y += w_ystep;	
		}

		if (best_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_h_bd(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i           = 0;	
	HKA_F32     angle       = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     best_val    = 0.0f;
	HKA_F32     move_step   = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		best_val = 1.0f * edge_strength;
		for (count = 0; count < width; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);

			if (val1 > best_val)
			{
				best_val = val1;
				point_a = point;
			}

			point.x += w_xstep;
			point.y += w_ystep;	
		}

		if (best_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
// 			val[edge_num]         = src[MVB_ROUNDF(point_a.y) * src_step + MVB_ROUNDF(point_a.x)];
// 			//val[edge_num]         = best_val;
			edge_num ++;
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	

}

static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_h_all(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;	
	HKA_F32     angle       = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     best_val    = 0.0f;
	HKA_F32     move_step   = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		best_val = 1.0f * edge_strength;
		for (count = 0; count < width; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);

			val1 = HKA_FABS(val1);

			if (val1 > best_val)
			{
				best_val = val1;
				point_a = point;
			}

			point.x += w_xstep;
			point.y += w_ystep;	
		}

		if (best_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_best_line_h(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
												        MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
												        MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32     edge_num      = 0;
	HKA_S32     ray_num       = 0;
	HKA_S32     edge_nature   = 0;
	HKA_S32     edge_strength = 0;
	HKA_S32     kernel_size   = 0;
	HKA_S32     region_width  = 0;
	HKA_U08     *buf          = HKA_NULL;
	HKA_POINT_F *edge_points  = HKA_NULL;

	ray_num       = cfg->ray_num;
	kernel_size   = ((cfg->kernel_size >> 1) << 1) + 1;
	region_width  = ((cfg->region_width >> 1) << 1) + 1;
	edge_nature   = cfg->edge_polarity;
	edge_strength = cfg->edge_strength;
	edge_points   = spec->edge_points;
	buf           = spec->buf;

	switch (edge_nature)
	{
	case MVBI_EDGE_TYPES_DARK_TO_BRIGHT:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_h_db(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BRIGHT_TO_DARK:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_h_bd(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BOTH:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_h_all(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 搜索沿着高的方向最后一条边
***************************************************************************************************/
/***************************************************************************************************
* 功  能：搜索最后一条边
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         ray_num        - I  射线的数目
*         edge_strength  - I  边缘的强度
*         edge_points    - O  边缘点集
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_h_db(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_F32     last_val    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;


		for (count = 0; count < width; count ++ )
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);

			if (val1 >= edge_strength)
			{
				point_a   = point;
				last_val  = val1;
			}

			point.x += w_xstep;
			point.y += w_ystep;
		}

		if (last_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_h_bd(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														   MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
														   HKA_S32 edge_strength, HKA_S32 kernel_size, 
														   HKA_S32 region_width, HKA_POINT_F *edge_points,
														   HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_F32     last_val    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;


		for (count = 0; count < width; count ++ )
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);

			if (val1 >= edge_strength)
			{
				point_a   = point;
				last_val  = val1;
			}

			point.x += w_xstep;
			point.y += w_ystep;
		}

		if (last_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;
}

static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_h_all(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_F32     move_step   = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_F32     last_val    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_POINT_F point_a     = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;


		for (count = 0; count < width; count ++ )
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);
			val1 = HKA_FABS(val1);
			if (val1 >= edge_strength)
			{
				point_a   = point;
				last_val  = val1;
			}

			point.x += w_xstep;
			point.y += w_ystep;
		}

		if (last_val > edge_strength)
		{
			edge_points[edge_num] = point_a;
			edge_num ++;
		}
		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;
}

static HKA_S32 MVAGEOFIND_LineDetect_search_last_line_h(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
												        MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
												        MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32     edge_num      = 0;
	HKA_S32     ray_num       = 0;
	HKA_S32     edge_nature   = 0;
	HKA_S32     kernel_size   = 0;
	HKA_S32     region_width  = 0;
	HKA_S32     edge_strength = 0;
	HKA_U08     *buf          = HKA_NULL;
	HKA_POINT_F *edge_points  = HKA_NULL;

	ray_num       = cfg->ray_num;
	kernel_size   = ((cfg->kernel_size >> 1) <<1) + 1;
	region_width  = ((cfg->region_width >> 1) <<1) + 1;
	edge_nature   = cfg->edge_polarity;
	edge_strength = cfg->edge_strength;
	edge_points   = spec->edge_points;
	buf           = spec->buf;

	switch (edge_nature)
	{
	case MVBI_EDGE_TYPES_DARK_TO_BRIGHT:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_h_db(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BRIGHT_TO_DARK:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_h_bd(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BOTH:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_h_all(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 搜索沿着高的方向第一条边
***************************************************************************************************/
/***************************************************************************************************
* 功  能：搜索第一条边
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         ray_num        - I  射线的数目
*         edge_strength  - I  边缘的强度
*         edge_points    - O  边缘点集
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_h_db(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     move_step   = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_F32     *kernel     = HKA_NULL;
	
	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}
		point = point_s;
		
		for (count = 0; count < width; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);
			point.x += w_xstep;
			point.y += w_ystep;

			if (val1 > edge_strength)
			{
				edge_points[edge_num] = point;
				edge_num ++;
				break;
			}
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_h_bd(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															HKA_S32 edge_strength, HKA_S32 kernel_size, 
															HKA_S32 region_width, HKA_POINT_F *edge_points,
															HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     move_step   = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_bd(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < width; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);
			point.x += w_xstep;
			point.y += w_ystep;

			if (val1 > edge_strength)
			{
				edge_points[edge_num] = point;
				edge_num ++;
				break;
			}
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_h_all(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
															 MVAGEOFIND_LD_RECT *rect, HKA_S32 ray_num, 
															 HKA_S32 edge_strength, HKA_S32 kernel_size, 
															 HKA_S32 region_width, HKA_POINT_F *edge_points,
															 HKA_U08 *buf)
{
	HKA_S32     i           = 0;
	HKA_F32     angle       = 0;
	HKA_S32     edge_num    = 0;
	HKA_S32     count       = 0;
	HKA_S32     width       = 0;
	HKA_F32     val1        = 0.0f;
	HKA_F32     move_step   = 0.0f;
	HKA_F32     w_xstep     = 0.0f;
	HKA_F32     w_ystep     = 0.0f;
	HKA_F32     h_xstep     = 0.0f;
	HKA_F32     h_ystep     = 0.0f;
	HKA_F32     h1_xstep    = 0.0f;
	HKA_F32     h1_ystep    = 0.0f;
	HKA_POINT_F point_s     = {0};
	HKA_POINT_F point       = {0};
	HKA_F32     *kernel     = HKA_NULL;

	kernel = (HKA_F32 *)buf;

	MVAGEOFIND_LineDetect_init_kernel_db(kernel_size, kernel);

	angle = rect->angle;

	move_step = 1.0f * rect->height / (ray_num - 1);

	w_xstep = MVB_COSF(angle);
	w_ystep = MVB_SINF(angle);

	h1_xstep = MVB_SINF(-angle);
	h1_ystep = MVB_COSF(-angle);

	h_xstep = move_step * MVB_SINF(-angle);
	h_ystep = move_step * MVB_COSF(-angle);

	point_s.x = rect->vertex[0].x;
	point_s.y = rect->vertex[0].y;

	width = rect->width - 1;

	for (i = 0; i < ray_num; i++)
	{
		if (   (point_s.x < 0) || (point_s.x >= roi_size.width)
			|| (point_s.y < 0) || (point_s.y >= roi_size.height))
		{
			break;
		}

		point = point_s;

		for (count = 0; count < width; count ++)
		{
			val1 = MVAGEOFIND_LineDetect_edge_search(src, src_step ,roi_size, point,
				w_xstep, w_ystep, h1_xstep, h1_ystep, region_width, kernel_size, kernel);
			val1 = HKA_FABS(val1);
			point.x += w_xstep;
			point.y += w_ystep;

			if (val1 > edge_strength)
			{
				edge_points[edge_num] = point;
				edge_num ++;
				break;
			}
		}

		point_s.x = point_s.x + h_xstep;
		point_s.y = point_s.y + h_ystep;
	}

	return edge_num;	
}

static HKA_S32 MVAGEOFIND_LineDetect_search_first_line_h(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size,
														 MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
														 MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32     edge_num      = 0;
	HKA_S32     ray_num       = 0;
	HKA_S32     edge_nature   = 0;
	HKA_S32     edge_strength = 0;
	HKA_S32     kernel_size   = 0;
	HKA_S32     region_width  = 0;
	HKA_U08     *buf          = HKA_NULL;
	HKA_POINT_F *edge_points  = HKA_NULL;

	ray_num       = cfg->ray_num;
	kernel_size   = ((cfg->kernel_size >> 1) << 1) + 1;
	region_width  = ((cfg->region_width >> 1) << 1) + 1;
	edge_nature   = cfg->edge_polarity;
	edge_strength = cfg->edge_strength;
	edge_points   = spec->edge_points;
	buf           = spec->buf;

	switch (edge_nature)
	{
	case MVBI_EDGE_TYPES_DARK_TO_BRIGHT:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_h_db(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BRIGHT_TO_DARK:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_h_bd(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	case MVBI_EDGE_TYPES_BOTH:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_h_all(src, src_step, roi_size, rect,  ray_num, 
										edge_strength, kernel_size, region_width, edge_points, buf);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 功  能：搜索边缘沿着宽方向的边缘
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_edge_along_width(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
											          MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
											          MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32 edge_type = 0;
	HKA_S32 edge_num  = 0;

	edge_type = cfg->edge_type;

	switch (edge_type)
	{
	case MVBI_LINE_FIND_FIRST:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_w(src, src_step, roi_size, rect, cfg, spec);
		break;
	case MVBI_LINE_FIND_LAST:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_w(src, src_step, roi_size, rect, cfg, spec);
		break;
	case MVBI_LINE_FIND_BEST:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_w(src, src_step, roi_size, rect, cfg, spec);
		break;
	default :
		break;
	}
	return edge_num;	
}

/***************************************************************************************************
* 功  能：搜索边缘沿着高方向的边缘
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         line_w         - I  宽方向直线
*         line_h         - I  高方向直线
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点的数目
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_edge_along_height(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
												       MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg,
												       MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32 edge_type = 0;
	HKA_S32 edge_num  = 0;

	edge_type = cfg->edge_type;

	switch (edge_type)
	{
	case MVBI_LINE_FIND_FIRST:
		edge_num = MVAGEOFIND_LineDetect_search_first_line_h(src, src_step, roi_size, rect, cfg, spec);
		break;
	case MVBI_LINE_FIND_LAST:
		edge_num = MVAGEOFIND_LineDetect_search_last_line_h(src, src_step, roi_size, rect, cfg, spec);
		break;
	case MVBI_LINE_FIND_BEST:
		edge_num = MVAGEOFIND_LineDetect_search_best_line_h(src, src_step, roi_size, rect, cfg, spec);
		break;
	default :
		break;
	}
	return edge_num;	
}

/***************************************************************************************************
* 功  能：计算矩形宽高两个方向的正交直线
* 参  数：*
*         rect           - I  矩形区域
*         line_w         - O  宽方向的直线
*         line_h         - O  高方向的直线
* 返回值：
* 备注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_compute_orth_line_coe(MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LD_LINE_COE *line_w, 
													        MVAGEOFIND_LD_LINE_COE *line_h)
{
	//直线方程为a*x + b*y + c = 0
	if (0.0f == rect->angle)
	{
		line_w->a = 0.0f;
		line_w->b = -1.0f;
		line_w->c = 1.0f * rect->vertex[0].y;

		line_h->a = -1.0f;
		line_h->b = 0.0f;
		line_h->c = 1.0f * rect->vertex[0].x;
	}
	else if (MVB_PI2F == rect->angle)
	{
		line_w->a = -1.0f;
		line_w->b = 0.0f;
		line_w->c = 1.0f * rect->vertex[0].x;
		
		line_h->a = 0.0f;
		line_h->b = -1.0f;
		line_h->c = 1.0f * rect->vertex[0].y;
	}
	else
	{
		line_w->a = MVB_TANF(rect->angle);
		line_w->b = -1.0f;
		line_w->c = - line_w->a * rect->vertex[0].x - line_w->b * rect->vertex[0].y ;

		line_h->a = -1.0f / line_w->a;
		line_h->b = -1.0f;
		line_h->c = - line_h->a * rect->vertex[0].x - line_h->b * rect->vertex[0].y;		
	}
}

/***************************************************************************************************
* 功  能：获取直线的参数
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         flag           - O  搜索方向标志
*         rect           - I  输入的旋转区域
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：边缘点集个数
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_edge_points_detect(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
												        MVAGEOFIND_LD_RECT *rect, MVAGEOFIND_LINE_DETECT_CFG *cfg, 
														MVAGEOFIND_LD_SPEC*spec)
{
	HKA_S32 edge_num       = 0;

	switch (cfg->find_orient)
	{
	case MVAGEOFIND_LINE_FIND_ORIENT_X:
		edge_num = MVAGEOFIND_LineDetect_edge_along_width(src, src_step, roi_size, rect, cfg, spec);
		break;
	case MVAGEOFIND_LINE_FIND_ORIENT_Y:
		edge_num = MVAGEOFIND_LineDetect_edge_along_height(src, src_step, roi_size, rect, cfg, spec);
		break;
	default:
		break;
	}

	return edge_num;
}

/***************************************************************************************************
* 功  能: 获取边缘直线起点以及终点
* 参  数：*
*         line_h         - I  高方向直线
*         line_end       - I  高方向平行直线
*         line_coe       - I  拟合出来的直线
*         line           - O  输出直线系数
* 返回值：状态码
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_compute_start_end_points(MVAGEOFIND_LD_LINE_COE *line_s, 
                                                              MVAGEOFIND_LD_LINE_COE *line_end, 
														      MVAGEOFIND_LD_LINE_COE *line_coe, 
                                                              MVAGEOFIND_LINE_PARA *line)
{
    HKA_F32 c_coe = 0.0f;
    HKA_F32 y_coe = 0.0f;

    //计算直线边缘与上边缘的交点
    c_coe = line_s->a * line_coe->c - line_coe->a * line_s->c;
    y_coe = line_coe->a * line_s->b - line_s->a * line_coe->b;

    if (y_coe == 0.0f)
    {
        return HKA_FALSE;
    }
    //边界水平，查找直线垂直
    if ((HKA_FABS(line_s->a) < 1e-5f) && (HKA_FABS(line_coe->b) < 1e-5f))
    {
        line->s_point.x = -line_coe->c / line_coe->a;
        line->s_point.y = line_s->c;

        line->e_point.x = line->s_point.x;
        line->e_point.y = line_end->c;
        return HKA_TRUE;

    }
    //边界垂直，查找直线水平
    if ((HKA_FABS(line_s->b) < 1e-5f) && (HKA_FABS(line_coe->a) < 1e-5f))
    {
        line->s_point.x = -line_s->c / line_s->a;
        line->s_point.y = line_coe->c;

        line->e_point.x = -line_end->c / line_end->a;
        line->e_point.y = line_coe->c;
        return HKA_TRUE;
    }

	// 交点
	// 	line->s_point.y = (line_s->a * line_coe->c - line_coe->a * line_s->c)
	// 						/ (line_coe->a * line_s->b - line_s->a * line_coe->b);
	// 	line->s_point.x = (line_coe->b * line_s->c - line_s->b * line_coe->c) 
	// 						/ (line_coe->a * line_s->b - line_s->a * line_coe->b);
	// 
	// 	line->e_point.y = (line_end->a * line_coe->c - line_coe->a * line_end->c)
	// 						/ (line_coe->a * line_end->b - line_end->a * line_coe->b);
	// 	line->e_point.x = (line_coe->b * line_end->c - line_end->b * line_coe->c) 
	// 						/ (line_coe->a * line_end->b - line_end->a * line_coe->b);

    line->s_point.y = c_coe / y_coe;
    line->s_point.x = -1.0f *(line_coe->b * line->s_point.y + line_coe->c) / line_coe->a;

    //计算直线边缘与下边缘的交点
    c_coe = line_end->a * line_coe->c - line_coe->a * line_end->c;
    y_coe = line_coe->a * line_end->b - line_end->a * line_coe->b;
    if (y_coe == 0.0f)
    {
        return HKA_FALSE;
    }

    line->e_point.y = c_coe / y_coe;
    line->e_point.x = -1.0f *(line_coe->b * line->e_point.y + line_coe->c) / line_coe->a;

    return HKA_TRUE;
}

/***************************************************************************************************
* 功  能: 获取边缘直线起点以及终点
* 参  数：*
*         rect           - I  输入的旋转区域
*         line_coe       - I  拟合出来的直线系数
*         line           - O  输出的直线参数
*         edge_orient    - I  边缘的方向
* 返回值：状态码
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_confirm_start_end_points(MVAGEOFIND_LD_RECT *rect, 
                                                              MVAGEOFIND_LD_LINE_COE *line_coe,
														      MVAGEOFIND_LINE_PARA *line)
{
	HKA_S32 num    = 0;
	HKA_S32 flag1  = 0;
	HKA_S32 flag2  = 0;
	HKA_POINT_F inter_points[4]     = {0};
	MVAGEOFIND_LD_LINE_COE line_w   = {0};
	MVAGEOFIND_LD_LINE_COE line_h   = {0};
	MVAGEOFIND_LD_LINE_COE line_end = {0};
	MVAGEOFIND_LINE_PARA   line1    = {0};
	MVAGEOFIND_LINE_PARA   line2    = {0};

	//计算正交直线方程
	MVAGEOFIND_LineDetect_compute_orth_line_coe(rect, &line_w, &line_h);

	line_end = line_h;
	line_end.c = -line_end.a * rect->vertex[2].x - line_end.b * rect->vertex[2].y;
	flag1 = MVAGEOFIND_LineDetect_compute_start_end_points(&line_h, &line_end, line_coe, &line1);

	line_end = line_w;
	line_end.c = -line_end.a * rect->vertex[3].x - line_end.b * rect->vertex[3].y;
	flag2 = MVAGEOFIND_LineDetect_compute_start_end_points(&line_w, &line_end, line_coe, &line2);

	if (!(flag1 || flag2))
	{
		return HKA_FALSE;
	}

	if (   ((line1.s_point.x - rect->vertex[0].x) * (line1.s_point.x - rect->vertex[3].x) < 0)
		|| ((line1.s_point.y - rect->vertex[0].y) * (line1.s_point.y - rect->vertex[3].y) < 0))
	{
		inter_points[num++] = line1.s_point;
	}
	if (   ((line2.s_point.x - rect->vertex[1].x) * (line2.s_point.x - rect->vertex[0].x) < 0)
		|| ((line2.s_point.y - rect->vertex[1].y) * (line2.s_point.y - rect->vertex[0].y) < 0))
	{
		inter_points[num++] = line2.s_point;
	}
	if (   ((line1.e_point.x - rect->vertex[2].x) * (line1.e_point.x - rect->vertex[1].x) < 0)
		|| ((line1.e_point.y - rect->vertex[2].y) * (line1.e_point.y - rect->vertex[1].y) < 0))
	{
		inter_points[num++] = line1.e_point;
	}
	if (   ((line2.e_point.x - rect->vertex[3].x) * (line2.e_point.x - rect->vertex[2].x) < 0)
		|| ((line2.e_point.y - rect->vertex[3].y) * (line2.e_point.y - rect->vertex[2].y) < 0))
	{
		inter_points[num++] = line2.e_point;
	}

	if (num != 2)
	{
		return HKA_FALSE;
	}
	line->s_point = inter_points[0];
	line->e_point = inter_points[1];

	return HKA_TRUE;
}

/***************************************************************************************************
* 功  能: 获取矩形框的四个顶点
* 参  数：*
*         center         - I  矩阵中心
*         width          - I  矩阵宽长度
*         height         - I  矩阵高长度
*         angle          - I  矩阵角度
*         vertex         - 0  高方向直线
* 返回值：无
* 备注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_compute_vertex(HKA_POINT_I *center, HKA_S32 width, HKA_S32 height,
													 HKA_F32 angle, HKA_POINT_F *vertex)
{
	HKA_F32 w     = 0.0f;
	HKA_F32 h     = 0.0f;
	HKA_F32 sin_v = 0.0f;
	HKA_F32 cos_v = 0.0f;

	sin_v = MVB_SINF(angle);
	cos_v = MVB_COSF(angle);

	w = 0.5f * width;
	h = 0.5f * height;

	// 宽的左顶点，为起点，顺时针取点
	vertex[0].x = center->x - cos_v * w + sin_v * h;
	vertex[0].y = center->y - sin_v * w - cos_v * h;
	vertex[3].x = center->x - cos_v * w - sin_v * h;
	vertex[3].y = center->y - sin_v * w + cos_v * h;
	vertex[1].x = 2.0f * center->x - vertex[3].x;
	vertex[1].y = 2.0f * center->y - vertex[3].y; 
	vertex[2].x = 2.0f * center->x - vertex[0].x;
	vertex[2].y = 2.0f * center->y - vertex[0].y; 

}
/***************************************************************************************************
* 功  能：拟合直线系数
* 参  数：*
*         edge_points    - I  边缘点集
*         edge_num       - I  边缘点数目
*         line_coe       - O  边缘直线系数
*         buf            - I  工作内存
* 返回值：
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_fit_line_coe(HKA_POINT_F *edge_points, HKA_S32 edge_num,
											      MVAGEOFIND_LD_LINE_COE *line_coe, HKA_F32 *straightness, 
											      HKA_S32 kernel_size, HKA_F32 dist_th, HKA_F32 reject_ratio,
												  HKA_U08 *buf)
{
	HKA_S32                i                                      = 0;
	HKA_S32                sts                                    = 0;
// 	HKA_F32                a                                      = 0.0f;
// 	HKA_F32                b                                      = 0.0f;
// 	HKA_F32                c                                      = 0.0f;
	HKA_F32                s                                      = 0.0f;
	HKA_S32                num                                    = 0;
	HKA_S32                index                                  = 0;
	MVAGEOFIND_LD_LINE_COE coe                                    = {0};
	HKA_F32                est_a[MVAGEOFIND_LD_RANSAC_WORK_TIMES] = {0.0f};
	HKA_F32                est_b[MVAGEOFIND_LD_RANSAC_WORK_TIMES] = {0.0f};
	HKA_F32                est_c[MVAGEOFIND_LD_RANSAC_WORK_TIMES] = {0.0f};
	HKA_F32                est_s[MVAGEOFIND_LD_RANSAC_WORK_TIMES] = {0.0f};

	for (i = 0; i < MVAGEOFIND_LD_RANSAC_WORK_TIMES; i++)
	{
		sts = MVAGEOFIND_LineDetect_ransac_line_coe(edge_points, edge_num, MVAGEOFIND_LD_RANSAC_ITER_TIMES, 
											        0.9f, dist_th, reject_ratio, kernel_size, &coe, &s, buf);
		HKA_CHECK_ERROR(sts == HKA_FALSE, sts);

		est_a[i] = coe.a;
		est_b[i] = coe.b;
		est_c[i] = coe.c;
		est_s[i] = s;
	}
	_MVBS_SortAscend_32f_I(est_a, MVAGEOFIND_LD_RANSAC_WORK_TIMES);
	_MVBS_SortAscend_32f_I(est_b, MVAGEOFIND_LD_RANSAC_WORK_TIMES);
	_MVBS_SortAscend_32f_I(est_c, MVAGEOFIND_LD_RANSAC_WORK_TIMES);
	_MVBS_SortAscend_32f_I(est_s, MVAGEOFIND_LD_RANSAC_WORK_TIMES);

	index = MVAGEOFIND_LD_RANSAC_WORK_TIMES >> 1;

	num = MVAGEOFIND_LD_RANSAC_WORK_TIMES - 1;
	s   = 0.0f;
	for (i = 1; i < num; i++)
	{
// 		a += est_a[i];
// 		b += est_b[i];
// 		c += est_c[i];
		s += est_s[i];
	}
	num = MVAGEOFIND_LD_RANSAC_WORK_TIMES - 2;
	line_coe->a = est_a[index];
	line_coe->b = est_b[index];
	line_coe->c = est_c[index];
// 	line_coe->a = a / num;
// 	line_coe->b = b / num;
// 	line_coe->c = c / num;
	*straightness = s / num;

	return HKA_TRUE;
}

/***************************************************************************************************
* 功  能：聚类
* 参  数：*
*         match_num         - I   聚类点的数量
*         val_data          - I   聚类点像素值
*         points_data       - I   存储点集坐标
*         mark              - I   聚类取大类还是小类标志
* 返回值：HKA_S32
* 备注：
***************************************************************************************************/
// static HKA_S32 MVAGEOFIND_LineDetect_cluster_edge_points(HKA_POINT_F *points, HKA_U08 *val, HKA_S32 match_num)											  										  
// {
// 	HKA_S32     i         = 0;
// 	HKA_S32     count     = 0;
// 	HKA_F32     c1        = 0.0f;
// 	HKA_F32     c1_sum    = 0;
// 	HKA_F32     c2        = 0.0f;
// 	HKA_S32     c2_sum    = 0;
// 	HKA_S32     c1_pre_n  = 0;
// 	HKA_S32     c1_now_n  = 0;
// 	HKA_S32     c2_pre_n  = 0;
// 	HKA_S32     c2_now_n  = 0;
// 	HKA_F32     diff      = 0.0f;
// 	HKA_POINT_F *pre      = HKA_NULL;
// 	HKA_POINT_F *beh      = HKA_NULL;
// 
// 	c1 = val[0];
// 	c2 = val[1];
// 	if (c1 == c2)
// 	{
// 		c1 = c2 + 1;
// 	}
// 
// 	pre = points; 
// 	beh = points;
// 
// 	do 
// 	{
// 		c1_pre_n = c1_now_n;
// 		c2_pre_n = c2_now_n;
// 		c1_now_n = 0;
// 		c2_now_n = 0;
// 		c1_sum   = 0;
// 		c2_sum   = 0;
// 		for (i = 0; i < match_num; i++)
// 		{
// 			if(HKA_FABS(val[i] - c1) < HKA_FABS(val[i] - c2))
// 			{
// 				c1_sum += val[i];
// 				c1_now_n += 1;
// 			}
// 			else
// 			{
// 				c2_sum += val[i];
// 				c2_now_n += 1;
// 			}
// 		}
// 		if ((c1_now_n == 0) || (c2_now_n == 0))
// 		{
// 			break;
// 		}
// 		c1 = 1.0f * c1_sum / c1_now_n;
// 		c2 = 1.0f * c2_sum / c2_now_n;
// 	} while((c1_now_n != c1_pre_n) || (c2_now_n != c2_pre_n));
// 
// // 	diff = 1.0f * HKA_ABS(c2_now_n - c1_now_n) / match_num;
// // 
// // 	if (diff < 0.1f )
// // 	{
// // 		return match_num;
// // 	}
// 
// 	if (c1_now_n >= c2_now_n)
// 	{
// 		for (i = 0; i < match_num; i++)
// 		{
// 			if(HKA_FABS(val[i] - c1) <= HKA_FABS(val[i] - c2))
// 			{
// 				beh[count] = pre[i];
// 				val[count] = val[i];
// 				count ++;
// 			}
// 			if (count == c1_now_n)
// 			{
// 				break;
// 			}
// 		}
// 		return c1_now_n;
// 	}
// 
// 	else
// 	{
// 		for (i = 0; i < match_num; i++)
// 		{
// 			if(HKA_FABS(val[i] - c2) <= HKA_FABS(val[i] - c1))
// 			{
// 				beh[count] = pre[i];
// 				val[count] = val[i];
// 				count ++;
// 			}
// 			if (count == c2_now_n)
// 			{
// 				break;
// 			}
// 		}
// 		return c2_now_n;
// 	}
// }
/***************************************************************************************************
* 功  能 判断角度是否超过容忍角度
* 参  数：*
*         rect_angle     - I  矩形框角度
*         line_angle     - I  直线角度
*         th_angle       - I  容忍角度
* 返回值：状态码
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEONFIND_LineDetect_judge_anlge(HKA_F32 rect_angle, HKA_F32 line_angle, HKA_F32 th_angle)
{
	HKA_F32 delta_angle = 0.0f;

	if (line_angle < -90.f)
	{
		line_angle += 180.f;
	}
	if (line_angle > 90.f)
	{
		line_angle -= 180.f;
	}

	delta_angle = HKA_FABS(line_angle - rect_angle);

	if (delta_angle > th_angle)
	{
		return HKA_FALSE;
	}
	return HKA_TRUE;
}
/***************************************************************************************************
* 功  能：获取直线的参数
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域
*         rect           - I  输入的旋转区域
*         cfg            - I  算法参数
*         spec           - I  内部结构体
* 返回值：状态码
* 备注：
***************************************************************************************************/
static HKA_S32 MVAGEOFIND_LineDetect_process_C1R(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
												  MVAGEOFIND_ROTATED_RECT *rect1, MVAGEOFIND_LINE_PARA *line,
												  MVAGEOFIND_LINE_DETECT_CFG *cfg, MVAGEOFIND_LD_SPEC *spec)										
{
	HKA_S32                edge_num        = 0;
	HKA_S32                sts             = 0;
	HKA_F32                straightness    = 0.0f;
	HKA_U08                *buf            = HKA_NULL;
	HKA_POINT_F            *edge_points    = HKA_NULL;
	MVAGEOFIND_LD_RECT     rect            = {0};
	MVAGEOFIND_LD_LINE_COE line_coe        = {0};
    HKA_F32                angle_tolerance = 0.0f;
    HKA_F32                line_angle      = 0.0f;
    HKA_F32                reject_dist     = 0.0f;
	HKA_F32                reject_ratio    = 0.0f;

	edge_points = spec->edge_points;
	buf         = spec->buf;
	
	rect.angle      = rect1->angle * MVB_PI180F;
	rect.center     = rect1->center;
	rect.width      = rect1->width;
	rect.height     = rect1->height;
    angle_tolerance = 1.0f * cfg->angle_tolerance;
    reject_dist     = 1.0f * cfg->reject_dist;
	reject_ratio    = 1.0f * cfg->reject_ratio;

	// 计算四个顶点
	MVAGEOFIND_LineDetect_compute_vertex(&rect.center, rect.width, rect.height, rect.angle, rect.vertex);

	// 获取点集
	edge_num = MVAGEOFIND_LineDetect_edge_points_detect(src, src_step, roi_size, &rect, cfg, spec);

	if (edge_num < MVAGEOFIND_LD_MIN_FIT_POINTS_NUM)
	{
		return HKA_FALSE;
	}

// 	if (edge_num > 10)
// 	{
// 		edge_num = MVAGEOFIND_LineDetect_cluster_edge_points(edge_points, val, edge_num);
// 	}
	
    sts = MVAGEOFIND_LineDetect_fit_line_coe(edge_points, edge_num, &line_coe, &straightness, 
                                             cfg->kernel_size, reject_dist, reject_ratio, buf);

    HKA_CHECK_ERROR(sts == HKA_FALSE, sts);

	sts = MVAGEOFIND_LineDetect_confirm_start_end_points(&rect, &line_coe, line);
	HKA_CHECK_ERROR(sts == HKA_FALSE, sts);

    line_angle = MVB_ATAN2F(line->e_point.y - line->s_point.y, line->e_point.x - line->s_point.x) * 180 * MVB_RPIF;

	sts = MVAGEONFIND_LineDetect_judge_anlge(rect1->angle, line_angle, angle_tolerance);
	HKA_CHECK_ERROR(HKA_FALSE == sts, sts);

	line->straightness = straightness;
	line->angle        = line_angle;
	return HKA_TRUE;

//     delta_angle = line_angle - rect1->angle;
// 
//     if(HKA_FABS(delta_angle) <= angle_tolerance)
//     {
        
//     }
//     else
//     {
//         return HKA_FALSE;
//     }
	
}

/***************************************************************************************************
* 功  能：检测旋转矩形参数
* 参  数：*
*         rect           - I  旋转矩形区域
*         roi_size       - I  感兴趣区域
* 返回值：状态码
* 备注：
***************************************************************************************************/
static HKA_STATUS MVAGEOFIND_LineDetect_check_rotated_rect(MVAGEOFIND_ROTATED_RECT *rect, HKA_SIZE_I roi_size)
{
	HKA_F32     angle       = 0.0f;
	HKA_POINT_F vertex[4]   = {0.0f};
    HKA_S32     half_width  = 0;
    HKA_S32     half_height = 0;

    half_width  = rect->width >> 1;
    half_height = rect->height >> 1;

	HKA_CHECK_ERROR(   rect->angle < MVAGEOFIND_LINE_DETECT_ROTATE_ANGLE_MIN
					|| rect->angle >= MVAGEOFIND_LINE_DETECT_ROTATE_ANGLE_MAX,
					HKA_STS_ERR_BAD_ARG);

    HKA_CHECK_ERROR(   (rect->center.x < half_width)
                    || (rect->center.y < half_height)
                    || (rect->center.x + half_width >= roi_size.width)
                    || (rect->center.y + half_height >= roi_size.height)
                    || (rect->width < 0)
                    || (rect->height < 0)
                    || (rect->width > roi_size.width)
                    || (rect->height > roi_size.height),
                    HKA_STS_ERR_BAD_ARG);

	// 计算矩形的四个顶点
	angle = rect->angle * MVB_PI180F;

	MVAGEOFIND_LineDetect_compute_vertex(&rect->center, rect->width, rect->height, angle, vertex);

	// 顶点不能越界
	HKA_CHECK_ERROR(   (vertex[0].x < 0)
		            || (vertex[0].y < 0)
					|| (vertex[1].x < 0)
					|| (vertex[1].y < 0)
					|| (vertex[2].x < 0)
					|| (vertex[2].y < 0)
					|| (vertex[3].x < 0)
					|| (vertex[3].y < 0)
					|| (vertex[0].x >= roi_size.width)
					|| (vertex[0].y >= roi_size.height)
					|| (vertex[1].x >= roi_size.width)
					|| (vertex[1].y >= roi_size.height)
					|| (vertex[2].x >= roi_size.width)
					|| (vertex[2].y >= roi_size.height)
					|| (vertex[3].x >= roi_size.width)
					|| (vertex[3].y >= roi_size.height),
					HKA_STS_ERR_BAD_ARG);

	return HKA_STS_OK;
}

/***************************************************************************************************
* 功  能：检测算法参数
* 参  数：*
*         rect           - I  旋转矩形区域
*         roi_size       - I  感兴趣区域
* 返回值：状态码
* 备注：
***************************************************************************************************/
static HKA_STATUS MVAGEOFIND_LineDetect_check_cfg(MVAGEOFIND_LINE_DETECT_CFG *cfg)
{
	HKA_CHECK_ERROR(   (cfg->ray_num < MVAGEOFIND_LINE_DETECT_RAY_NUM_MIN)
					|| (cfg->ray_num > MVAGEOFIND_LINE_DETECT_RAY_NUM_MAX),
					HKA_STS_ERR_BAD_ARG);

	HKA_CHECK_ERROR(   (cfg->edge_polarity < MVBI_EDGE_TYPES_DARK_TO_BRIGHT)
					|| (cfg->edge_polarity > MVBI_EDGE_TYPES_BOTH),
					HKA_STS_ERR_BAD_ARG);

	HKA_CHECK_ERROR(   (cfg->find_orient < MVBI_LINE_FIND_ORIENT_X)
					|| (cfg->find_orient > MVBI_LINE_FIND_ORIENT_Y),
					HKA_STS_ERR_BAD_ARG);

	HKA_CHECK_ERROR(   (cfg->edge_type < MVBI_LINE_FIND_BEST)
					|| (cfg->edge_type > MVBI_LINE_FIND_LAST),
					HKA_STS_ERR_BAD_ARG);

	HKA_CHECK_ERROR(   (cfg->kernel_size < MVAGEOFIND_LINE_DETECT_EDGE_WIDTH_MIN)
		            || (cfg->kernel_size > MVAGEOFIND_LINE_DETECT_EDGE_WIDTH_MAX),
					HKA_STS_ERR_BAD_ARG);

	HKA_CHECK_ERROR(   (cfg->region_width < MVAGEOFIND_LINE_DETECT_EDGE_REGION_WIDTH_MIN)
					|| (cfg->region_width > MVAGEOFIND_LINE_DETECT_EDGE_REGION_WIDTH_MAX),
					HKA_STS_ERR_BAD_ARG);

    HKA_CHECK_ERROR(   (cfg->angle_tolerance < 0)
                    || (cfg->angle_tolerance > 180),
                    HKA_STS_ERR_BAD_ARG);

	HKA_CHECK_ERROR(   (cfg->reject_ratio < 0.0f)
		            || (cfg->reject_ratio > 1.0f),
					HKA_STS_ERR_BAD_ARG);

    HKA_CHECK_ERROR(cfg->reject_dist < 1, HKA_STS_ERR_BAD_ARG);

	return HKA_STS_OK;
}

/***************************************************************************************************
* 功  能：分配工作内存
* 参  数：*
*         spec       - IO  算法结构体
*         work_buf   - I   工作内存
*         roi_size  – I   图像大小
*         work_size  - O   工作内存大小
* 返回值：HKA_VOID
* 备注：
***************************************************************************************************/
static HKA_VOID MVAGEOFIND_LineDetect_alloc_mem(MVAGEOFIND_LD_SPEC *spec, HKA_U08 *work_buf, 
										        HKA_SIZE_I roi_size, HKA_SZT *work_size)
{
	HKA_SZT size       = 0;
	HKA_SZT used_size  = 0;
	HKA_U08 *alloc_buf = HKA_NULL;

	alloc_buf = work_buf;

	//分配拟合点集内存
	spec->edge_points = (HKA_POINT_F *)alloc_buf;
	size       = MVAGEOFIND_LINE_DETECT_RAY_NUM_MAX * sizeof(HKA_POINT_F);
	size       = HKA_SIZE_ALIGN_128(size);
	used_size += size;
	alloc_buf += size;

	//分配拟合点集像素值内存
	spec->val  = (HKA_U08 *)alloc_buf;
	size       = MVAGEOFIND_LINE_DETECT_RAY_NUM_MAX * sizeof(HKA_U08);
	size       = HKA_SIZE_ALIGN_128(size);
	used_size += size;
	alloc_buf += size;

	spec->buf  = (HKA_U08 *)alloc_buf;
	size       = HKA_MAX(roi_size.width, roi_size.height) * sizeof(HKA_U08);
	size       = HKA_SIZE_ALIGN_128(size);
	used_size += size;

	*work_size = used_size;
}

/***************************************************************************************************
* 功  能：获取工作内存的大小
* 参  数：*
*         roi_size       - I  图像尺寸大小
*         work_size      - O  所需工作内存大小
* 返回值：
* 备注：
***************************************************************************************************/
HKA_VOID _MVAGEOFIND_LineDetectBufferSize(HKA_SIZE_I roi_size, HKA_SZT *work_size)
{
	HKA_U08    *work_buf     = HKA_NULL;
	MVAGEOFIND_LD_SPEC spec = {0};

	work_buf = (HKA_U08 *)&spec;

	MVAGEOFIND_LineDetect_alloc_mem(&spec, work_buf, roi_size, work_size);
}

/***************************************************************************************************
* 功  能：获取工作内存的大小
* 参  数：*
*         roi_size       - I  最大图像尺寸
*         work_size      - O  所需工作内存大小
* 返回值：状态码
* 备注：
***************************************************************************************************/
HKA_STATUS MVAGEOFIND_LineDetectBufferSize(HKA_SIZE_I roi_size, HKA_SZT *work_size)
{
	HKA_SZT size = 0;

	HKA_CHECK_ERROR(work_size == HKA_NULL, HKA_STS_ERR_NULL_PTR);

	HKA_CHECK_ERROR(    roi_size.width <= 0
					|| roi_size.height <= 0,
					HKA_STS_ERR_DATA_SIZE);

	_MVAGEOFIND_LineDetectBufferSize(roi_size,  &size);

	HKA_CHECK_ERROR(size > HKA_MAX_MEM_SIZE, HKA_STS_ERR_OVER_MAX_MEM);

	*work_size = size;	

	return HKA_STS_OK;
}

/***************************************************************************************************
* 功  能：直线查找
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域大小
*         rect           - I  输入的旋转矩形
*         find           - O  输出状态
*         line           - O  输出的直线参数
*         cfg            - I  算法检测参数
*         work_buf       - I  工作内存
* 返回值：
* 备注：
***************************************************************************************************/
HKA_VOID _MVAGEOFIND_LineDetect_8u_C1R(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
									   MVAGEOFIND_ROTATED_RECT *rect, HKA_S32 *find,
									   MVAGEOFIND_LINE_PARA *line, MVAGEOFIND_LINE_DETECT_CFG *cfg, 
									   HKA_U08 *work_buf)									 															   
{
	HKA_SZT            work_size =  0;
	MVAGEOFIND_LD_SPEC spec      = {0};

	MVAGEOFIND_LineDetect_alloc_mem(&spec, work_buf, roi_size, &work_size);

	*find = MVAGEOFIND_LineDetect_process_C1R(src, src_step, roi_size, rect, line, cfg, &spec);
}

/***************************************************************************************************
* 功  能：直线查找
* 参  数：*
*         src            - I  输入图像的大小
*         src_step       - I  输入图像的行间距
*         roi_size       - I  感兴趣区域大小
*         rect           - I  输入的旋转矩形
*         find           - O  输出状态
*         line           - O  输出的直线参数
*         cfg            - I  算法检测参数
*         work_buf       - I  工作内存
* 返回值：
* 备注：
***************************************************************************************************/
HKA_STATUS MVAGEOFIND_LineDetect_8u_C1R(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
								        MVAGEOFIND_ROTATED_RECT *rect, HKA_S32 *find,
										MVAGEOFIND_LINE_PARA *line, MVAGEOFIND_LINE_DETECT_CFG *cfg, 
										HKA_U08 *work_buf)									  
{
	HKA_STATUS sts = 0;

	HKA_CHECK_ERROR(src == HKA_NULL, HKA_STS_ERR_NULL_PTR);

	HKA_CHECK_ERROR(   (rect     == HKA_NULL)
					|| (find     == HKA_NULL)
					|| (line     == HKA_NULL)
					|| (cfg      == HKA_NULL)
					|| (work_buf == HKA_NULL),
					HKA_STS_ERR_NULL_PTR);

    HKA_CHECK_ERROR((roi_size.width  <= 0) || (roi_size.height <= 0), HKA_STS_ERR_DATA_SIZE);
    HKA_CHECK_ERROR(src_step < roi_size.width, HKA_STS_ERR_STEP);

	sts = MVAGEOFIND_LineDetect_check_rotated_rect(rect, roi_size);
	HKA_CHECK_ERROR(sts != HKA_STS_OK, sts);

	sts = MVAGEOFIND_LineDetect_check_cfg(cfg);
	HKA_CHECK_ERROR(sts != HKA_STS_OK, sts);

	_MVAGEOFIND_LineDetect_8u_C1R(src, src_step, roi_size, rect, find, line, cfg, work_buf);

	return HKA_STS_OK;
}
