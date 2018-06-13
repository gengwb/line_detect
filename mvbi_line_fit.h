/***************************************************************************************************
* 
* 版权信息：版权所有 (c) 2016, 杭州海康威视数字技术股份有限公司, 保留所有权利
* 
* 文件名称：mvbi_line_fit.h
* 摘    要：直线拟合算法模块
*
* 当前版本：0.1.0
* 作    者：邓志辉
* 日    期：2016-01-13
* 备    注：创建
***************************************************************************************************/
#ifndef _MVBI_LINE_FIT_H_
#define _MVBI_LINE_FIT_H_

#include "hka_types.h"

#ifdef __cplusplus
extern "C" {
#endif

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
                                      HKA_S32 *status);
HKA_VOID _MVBI_LineFitLeastSquare_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                      HKA_S32 *status);
HKA_STATUS MVBI_LineFitLeastSquare_32s(HKA_POINT_I *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                       HKA_S32 *status);
HKA_STATUS MVBI_LineFitLeastSquare_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                       HKA_S32 *status);

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
HKA_VOID _MVBI_LineFitRansacBufferSize_32s(HKA_S32 max_point_num, HKA_SZT *work_size);
HKA_VOID _MVBI_LineFitRansacBufferSize_32f(HKA_S32 max_point_num, HKA_SZT *work_size);
HKA_STATUS MVBI_LineFitRansacBufferSize_32s(HKA_S32 max_point_num, HKA_SZT *work_size);
HKA_STATUS MVBI_LineFitRansacBufferSize_32f(HKA_S32 max_point_num, HKA_SZT *work_size);

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
                                 HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf);
HKA_VOID _MVBI_LineFitRansac_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                 HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                 HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf);
HKA_STATUS MVBI_LineFitRansac_32s(HKA_POINT_I *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                  HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                  HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf);
HKA_STATUS MVBI_LineFitRansac_32f(HKA_POINT_F *point_data, HKA_S32 point_num, HKA_F32 line_coef[3], 
                                  HKA_S32 *status, HKA_F32 *straight, HKA_F32 sample_ratio, 
                                  HKA_F32 dist_th, HKA_S32 max_iter, HKA_U08 *work_buf);


#ifdef __cplusplus
}
#endif

#endif //_MVBI_LINE_FIT_H_


