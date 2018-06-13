/***************************************************************************************************
* 版权信息：Copyright (c) 2015, 杭州海康威视数字技术股份有限公司
* 
* 文件名称：mvageofind_line_detect_type.h
* 摘    要：内部函数结构体
*
* 当前版本：V0.1.0
* 作    者：郑俊君
* 日    期：2015-11-13
* 备    注：1.创建    
***************************************************************************************************/
#ifndef _MVAGEOFIND_LINE_DETECT_TYPE_H_
#define _MVAGEOFIND_LINE_DETECT_TYPE_H_

#include "hka_types.h"

#ifdef __cplusplus
extern "C"{
#endif

#define  MVAGEOFIND_LD_RANSAC_RAND_TIMES  (10)   //RANSAC取随机数的次数
#define  MVAGEOFIND_LD_RANSAC_WORK_TIMES  (5)    //RANSAC拟合直线的次数
#define  MVAGEOFIND_LD_RANSAC_ITER_TIMES  (30)
#define  MVAGEOFIND_LD_MIN_FIT_POINTS_NUM (2)

#define  MVBI_LINE_FIND_ORIENT_WIDTH      (1)
#define  MVBI_LINE_FIND_ORIENT_HEIGHT     (2)

/***************************************************************************************************
* 内部数据结构
***************************************************************************************************/
typedef struct _MVAGEOFIND_LD_SPEC_
{
	HKA_POINT_F *edge_points;        //边缘点坐标
	HKA_U08     *val;                //边缘点像素值
	HKA_U08     *buf;                //工作内存
}MVAGEOFIND_LD_SPEC;

// 内部矩形参数
typedef struct _MVAGEOFIND_LD_RECT_
{
	HKA_POINT_I center;
	HKA_POINT_F vertex[4];
	HKA_S32     width;
	HKA_S32     height;
	HKA_F32     angle;
}MVAGEOFIND_LD_RECT;

//直线的系数
typedef struct _MVAGEOFIND_LD_LINE_COE_
{
	HKA_F32 a;
	HKA_F32 b;
	HKA_F32 c;
}MVAGEOFIND_LD_LINE_COE;

#ifdef __cplusplus
}
#endif

#endif //_MVAGEOFIND_LINE_DETECT_TYPE_H_