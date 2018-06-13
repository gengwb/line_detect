/***************************************************************************************************
* ��Ȩ��Ϣ��Copyright (c) 2016, ���ݺ����������ּ����ɷ����޹�˾
* 
* �ļ����ƣ�mvbi_line_fit_type.h
* ժ    Ҫ���ڲ������ṹ��
*
* ��ǰ�汾��V0.1.0
* ��    �ߣ���־��
* ��    �ڣ�2016-01-27
* ��    ע��1.����    
***************************************************************************************************/
#ifndef _MVBI_LINE_FIT_TYPE_H_
#define _MVBI_LINE_FIT_TYPE_H_

#include "hka_types.h"

#ifdef __cplusplus
extern "C"{
#endif


/***************************************************************************************************
* �ڲ����ݽṹ
***************************************************************************************************/
typedef struct _MVBI_LF_RANSAC_SPEC_F32
{
	HKA_POINT_F *consensus;        //��Ե������
}MVBI_LF_RANSAC_SPEC_F32;

typedef struct _MVBI_LF_RANSAC_SPEC_S32
{
	HKA_POINT_I *consensus;        //��Ե������
}MVBI_LF_RANSAC_SPEC_S32;


#ifdef __cplusplus
}
#endif

#endif //_MVBI_LINE_FIT_TYPE_H_