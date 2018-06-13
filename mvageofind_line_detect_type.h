/***************************************************************************************************
* ��Ȩ��Ϣ��Copyright (c) 2015, ���ݺ����������ּ����ɷ����޹�˾
* 
* �ļ����ƣ�mvageofind_line_detect_type.h
* ժ    Ҫ���ڲ������ṹ��
*
* ��ǰ�汾��V0.1.0
* ��    �ߣ�֣����
* ��    �ڣ�2015-11-13
* ��    ע��1.����    
***************************************************************************************************/
#ifndef _MVAGEOFIND_LINE_DETECT_TYPE_H_
#define _MVAGEOFIND_LINE_DETECT_TYPE_H_

#include "hka_types.h"

#ifdef __cplusplus
extern "C"{
#endif

#define  MVAGEOFIND_LD_RANSAC_RAND_TIMES  (10)   //RANSACȡ������Ĵ���
#define  MVAGEOFIND_LD_RANSAC_WORK_TIMES  (5)    //RANSAC���ֱ�ߵĴ���
#define  MVAGEOFIND_LD_RANSAC_ITER_TIMES  (30)
#define  MVAGEOFIND_LD_MIN_FIT_POINTS_NUM (2)

#define  MVBI_LINE_FIND_ORIENT_WIDTH      (1)
#define  MVBI_LINE_FIND_ORIENT_HEIGHT     (2)

/***************************************************************************************************
* �ڲ����ݽṹ
***************************************************************************************************/
typedef struct _MVAGEOFIND_LD_SPEC_
{
	HKA_POINT_F *edge_points;        //��Ե������
	HKA_U08     *val;                //��Ե������ֵ
	HKA_U08     *buf;                //�����ڴ�
}MVAGEOFIND_LD_SPEC;

// �ڲ����β���
typedef struct _MVAGEOFIND_LD_RECT_
{
	HKA_POINT_I center;
	HKA_POINT_F vertex[4];
	HKA_S32     width;
	HKA_S32     height;
	HKA_F32     angle;
}MVAGEOFIND_LD_RECT;

//ֱ�ߵ�ϵ��
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