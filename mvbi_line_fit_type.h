/***************************************************************************************************
* 版权信息：Copyright (c) 2016, 杭州海康威视数字技术股份有限公司
* 
* 文件名称：mvbi_line_fit_type.h
* 摘    要：内部函数结构体
*
* 当前版本：V0.1.0
* 作    者：邓志辉
* 日    期：2016-01-27
* 备    注：1.创建    
***************************************************************************************************/
#ifndef _MVBI_LINE_FIT_TYPE_H_
#define _MVBI_LINE_FIT_TYPE_H_

#include "hka_types.h"

#ifdef __cplusplus
extern "C"{
#endif


/***************************************************************************************************
* 内部数据结构
***************************************************************************************************/
typedef struct _MVBI_LF_RANSAC_SPEC_F32
{
	HKA_POINT_F *consensus;        //边缘点坐标
}MVBI_LF_RANSAC_SPEC_F32;

typedef struct _MVBI_LF_RANSAC_SPEC_S32
{
	HKA_POINT_I *consensus;        //边缘点坐标
}MVBI_LF_RANSAC_SPEC_S32;


#ifdef __cplusplus
}
#endif

#endif //_MVBI_LINE_FIT_TYPE_H_