// ��������
#ifndef _BASIC_FUNCTION_H
#define _BASIC_FUNCTION_H

//#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif

// �������
int CrossProduct(float * Y, float * Z, float * YXZ);

// �������
int DotProduct( float * Y, float * Z, float * Dot );

// �������, C = A * B, A: m * nά, B: n * pά, C: m * pά
int MultiplyMat( unsigned int m, unsigned int n, unsigned int p, float * A, float * B, float * C );

// Fast inverse square-root
int invSqrt( float x, float * inv );

// ����ģ
int CalVectorNorm( float * Y, float * VectNorm );

// �������ƺ���
float AmpConstrain( float x, float xMin, float xMax );

// �����䱣���ķ����Һ���
float ArcSinSecure( float x );

// �����䱣���ķ����Һ���
float ArcCosSecure( float x );

// ����Ԫ������Ϊ��λ��Ԫ��
int SetQuaterAsUnit( float * q );

// ��Ԫ��תŷ���ǣ�Z-X-Y˳��
int Quater2EulerZXY( float * q, float * euler);

// ��Ԫ��תŷ���ǣ�Z-Y-X˳��
int Quater2EulerZYX( float * q, float * euler);

// ��ת����ת��Ԫ��(Z-Y-X˳��)
int DCM2QuaterZYX( float * Cnb, float * q );

// ��Ԫ��ת��ת����(Z-Y-X˳��)
int Quter2DCMZYX( float * q, float * Cnb );

// ����Ԫ��src��ÿ��Ԫ��ȡ������ֵ��dst��dst = -src
int OppositeQuater( float * src, float * dst );

// ��Ԫ���Ĺ���( dst = conj(src) )
int QuaterConj( float * src, float * dst );

// ��Ԫ����ģ
int QuaterModulus( float * q, float * QuaterNorm );

// ��Ԫ����һ��
int QuaterNormilize( float * q );

// ��Ԫ���˷�( r = q * p )
int QuaterMultiply( float * q, float * p, float * r );

// ��ԭ����p������Ԫ��q��ת�任�������r,r = q * p * q^(-1)
int CalCoordQuater( float * q, float * p, float * r );

// ��ԭ����p������Ԫ��q��ת��任�������r,r = q^(-1) * p * q
int CalInvCoordQuater( float * q, float * p, float * r );

// ��̬�˶�ѧ����
int AttitudeDynamics( unsigned int n, unsigned int m, float t, const float * x, const float * u, float * dx, void * pData );


#ifdef __cplusplus
}
#endif

#endif     // _BASIC_FUNCTION_H


