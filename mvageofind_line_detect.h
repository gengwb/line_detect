/***************************************************************************************************
* ��Ȩ��Ϣ��Copyright (c) 2015, ���ݺ����������ּ����ɷ����޹�˾
* 
* �ļ����ƣ�mvageofind_line_detect.h
* ժ    Ҫ���⺯��ͷ�ļ�
*
* ��ǰ�汾��V0.1.2
* ��    �ߣ���־��
* ��    �ڣ�2016-1-12
* ��    ע��1���޸��㷨�����������Ա�����������
*           2���������̽Ƕ�
*
* ��ʷ�汾��V0.1.1
* ��    �ߣ�֣����
* ��    �ڣ�2015-12-22
* ��    ע��1���ع���ת����������ʽ
*           2��������������Ҷ�ˣ���߼���ȶ���
*
* ��ʷ�汾��V0.1.0
* ��    �ߣ�֣����
* ��    �ڣ�2015-11-13
* ��    ע��1.����   
***************************************************************************************************/
#ifndef _MVAGEOFIND_LINE_DETECT_H_
#define _MVAGEOFIND_LINE_DETECT_H_

#include "hka_defs.h"
#include "mvb_geofind_lib.h"

#ifdef __cplusplus
extern "C"{
#endif

//ֱ�߼���㷨����
/***************************************************************************************************
* �� ��1�����εĳ����Լ����㶨��
* ���壺���ζ������ҵ�ֱ��
* �߶��壺���ζ��������ֱ��
* ��  �㣺y������С�Ķ���
***************************************************************************************************/
typedef struct _MVAGEOFIND_LINE_DETECT_CFG_
{
	HKA_S32 ray_num;         // ���ߵ�����
	HKA_S32 edge_polarity;   // ��Ե���� 1���Ӻڵ���    2���Ӱ׵���   3������
	HKA_S32 edge_type;       // ��Ե���� 1������ 2: ��һ����    3�����һ���� 
	HKA_S32 find_orient;     // ��Ե���� 1��ֱ��X�������  2��ֱ��Y�������
	HKA_S32 edge_strength;   // ��Եǿ��
	HKA_S32 kernel_size;     // ��Ե���
	HKA_S32 region_width;    // ������
    HKA_S32 angle_tolerance; // �Ƕ����̶�
    HKA_S32 reject_dist;     // �޳���ֵ
	HKA_F32 reject_ratio;    // �޳�����
}MVAGEOFIND_LINE_DETECT_CFG;


//ֱ���������
typedef struct _MVAGEOFIND_LINE_PARA_
{
	HKA_POINT_F s_point;      //��ʼ��
	HKA_POINT_F e_point;      //�յ�
	HKA_F32     angle;        //�Ƕȣ��Ƕ���ͼ������Ƕ����
	HKA_F32     straightness; //ֱ�߶�
}MVAGEOFIND_LINE_PARA;

/***************************************************************************************************
* ��  �ܣ���ȡ�����ڴ�Ĵ�С
* ��  ����*
*         roi_size       - I  ͼ���roi_size
*         work_size      - O  ���蹤���ڴ��С
* ����ֵ��
* ��ע��
***************************************************************************************************/
HKA_VOID  _MVAGEOFIND_LineDetectBufferSize(HKA_SIZE_I roi_size, HKA_SZT *work_size);
HKA_STATUS MVAGEOFIND_LineDetectBufferSize(HKA_SIZE_I roi_size, HKA_SZT *work_size);

/***************************************************************************************************
* ��  �ܣ�ֱ�߲���
* ��  ����*
*         src            - I  ����ͼ��Ĵ�С
*         src_step       - I  ����ͼ����м��
*         roi_size       - I  ����Ȥ�����С
*         rect           - I  �������ת����
*         find           - O  ���״̬
*         line           - O  �����ֱ�߲���
*         cfg            - I  �㷨������
*         work_buf       - I  �����ڴ�
* ����ֵ��
* ��ע��
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

