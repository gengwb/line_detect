/***************************************************************************************************
* 版权信息：Copyright (c) 2015, 杭州海康威视数字技术股份有限公司
* 
* 文件名称：mvageofind_line_detect.h
* 摘    要：库函数头文件
*
* 当前版本：V0.1.2
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
#ifndef _MVAGEOFIND_LINE_DETECT_H_
#define _MVAGEOFIND_LINE_DETECT_H_

#include "hka_defs.h"
#include "mvb_geofind_lib.h"

#ifdef __cplusplus
extern "C"{
#endif

//直线检测算法参数
/***************************************************************************************************
* 定 义1：矩形的长宽以及顶点定义
* 宽定义：矩形顶点向右的直线
* 高定义：矩形顶点向左的直线
* 顶  点：y坐标最小的顶点
***************************************************************************************************/
typedef struct _MVAGEOFIND_LINE_DETECT_CFG_
{
	HKA_S32 ray_num;         // 射线的数量
	HKA_S32 edge_polarity;   // 边缘属性 1：从黑到白    2：从白到黑   3：最优
	HKA_S32 edge_type;       // 边缘类型 1：最优 2: 第一条边    3：最后一条边 
	HKA_S32 find_orient;     // 边缘方向 1：直线X轴分量多  2：直线Y轴分量多
	HKA_S32 edge_strength;   // 边缘强度
	HKA_S32 kernel_size;     // 边缘宽度
	HKA_S32 region_width;    // 区域宽度
    HKA_S32 angle_tolerance; // 角度容忍度
    HKA_S32 reject_dist;     // 剔除阈值
	HKA_F32 reject_ratio;    // 剔除比例
}MVAGEOFIND_LINE_DETECT_CFG;


//直线输出参数
typedef struct _MVAGEOFIND_LINE_PARA_
{
	HKA_POINT_F s_point;      //起始点
	HKA_POINT_F e_point;      //终点
	HKA_F32     angle;        //角度，角度以图像坐标角度输出
	HKA_F32     straightness; //直线度
}MVAGEOFIND_LINE_PARA;

/***************************************************************************************************
* 功  能：获取工作内存的大小
* 参  数：*
*         roi_size       - I  图像的roi_size
*         work_size      - O  所需工作内存大小
* 返回值：
* 备注：
***************************************************************************************************/
HKA_VOID  _MVAGEOFIND_LineDetectBufferSize(HKA_SIZE_I roi_size, HKA_SZT *work_size);
HKA_STATUS MVAGEOFIND_LineDetectBufferSize(HKA_SIZE_I roi_size, HKA_SZT *work_size);

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
HKA_VOID  _MVAGEOFIND_LineDetect_8u_C1R(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
                                        MVAGEOFIND_ROTATED_RECT *rect, HKA_S32 *find, 
                                        MVAGEOFIND_LINE_PARA *line, MVAGEOFIND_LINE_DETECT_CFG *cfg,
                                        HKA_U08 *work_buf);
HKA_STATUS MVAGEOFIND_LineDetect_8u_C1R(HKA_U08 *src, HKA_S32 src_step, HKA_SIZE_I roi_size, 
								        MVAGEOFIND_ROTATED_RECT *rect, HKA_S32 *find, 
										MVAGEOFIND_LINE_PARA *line, MVAGEOFIND_LINE_DETECT_CFG *cfg,
								        HKA_U08 *work_buf);								  																			  

#ifdef __cplusplus
}
#endif

#endif//_MVAGEOFIND_LINE_DETECT_H_

