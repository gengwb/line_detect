/***************************************************************************************************
* 
* ��Ȩ��Ϣ����Ȩ���� (c) 2016, ���ݺ����������ּ����ɷ����޹�˾, ��������Ȩ��
* 
* �ļ����ƣ�mvbi_line_fit.h
* ժ    Ҫ��ֱ������㷨ģ��
*
* ��ǰ�汾��0.1.0
* ��    �ߣ���־��
* ��    �ڣ�2016-01-13
* ��    ע������
***************************************************************************************************/
#ifndef _MVBI_LINE_FIT_H_
#define _MVBI_LINE_FIT_H_

#include "hka_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************************************
* ��������MVBI_LineFitLeastSquare_32s MVBI_LineFitLeastSquare_32f
*         _MVBI_LineFitLeastSquare_32s _MVBI_LineFitLeastSquare_32f
* ��  �ܣ�ֱ�����, ֱ�߷��̣�c1 *x + c2 * y + c3 = 0
* ��  ����*
*         point_data  - I  ���������
*         point_num   - I  ��������
*         line_coef   - O  ���ֱ��ϵ����[c1, c2, c3]
*         status      - O  ���״̬
* ����ֵ��HKA_VOID
* ��  ע��
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
* ��������MVBI_LineFitRansacBufferSize_32s MVBI_LineFitRansacBufferSize_32f
*         _MVBI_LineFitRansacBufferSize_32s _MVBI_LineFitRansacBufferSize_32f
* ��  �ܣ���ȡ����ransac��ֱ����Ϲ����ڴ��С
* ��  ����*
*         max_point_num - I  �������
*         work_size     - O  �����ڴ��С
* ����ֵ��HKA_VOID
* ��  ע��
***************************************************************************************************/
HKA_VOID _MVBI_LineFitRansacBufferSize_32s(HKA_S32 max_point_num, HKA_SZT *work_size);
HKA_VOID _MVBI_LineFitRansacBufferSize_32f(HKA_S32 max_point_num, HKA_SZT *work_size);
HKA_STATUS MVBI_LineFitRansacBufferSize_32s(HKA_S32 max_point_num, HKA_SZT *work_size);
HKA_STATUS MVBI_LineFitRansacBufferSize_32f(HKA_S32 max_point_num, HKA_SZT *work_size);

/***************************************************************************************************
* ��������MVBI_LineFitRansac_32s MVBI_LineFitRansac_32f
*         _MVBI_LineFitRansac_32s _MVBI_LineFitRansac_32f
* ��  �ܣ���ȡ����ransac��ֱ����Ϲ����ڴ��С
* ��  ����*
*         point_data  - I  ���������
*         point_num   - I  ��������
*         line_coef   - O  ���ֱ��ϵ����[c1, c2, c3]
*         status      - O  ���״̬
*         straight    - O  ֱ�߶�
*         sample_th   - I  ������ֵ
*         dist_th     - I  ������ֵ
*         max_iter    - I  ����������
*         work_buf    - I  �����ڴ�
* ����ֵ��HKA_VOID
* ��  ע��
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


