// 基础函数
#ifndef _BASIC_FUNCTION_H
#define _BASIC_FUNCTION_H

//#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif

// 向量叉乘
int CrossProduct(float * Y, float * Z, float * YXZ);

// 向量点乘
int DotProduct( float * Y, float * Z, float * Dot );

// 矩阵相乘, C = A * B, A: m * n维, B: n * p维, C: m * p维
int MultiplyMat( unsigned int m, unsigned int n, unsigned int p, float * A, float * B, float * C );

// Fast inverse square-root
int invSqrt( float x, float * inv );

// 向量模
int CalVectorNorm( float * Y, float * VectNorm );

// 区间限制函数
float AmpConstrain( float x, float xMin, float xMax );

// 带区间保护的反正弦函数
float ArcSinSecure( float x );

// 带区间保护的反余弦函数
float ArcCosSecure( float x );

// 将四元数设置为单位四元数
int SetQuaterAsUnit( float * q );

// 四元数转欧拉角（Z-X-Y顺序）
int Quater2EulerZXY( float * q, float * euler);

// 四元数转欧拉角（Z-Y-X顺序）
int Quater2EulerZYX( float * q, float * euler);

// 旋转矩阵转四元数(Z-Y-X顺序)
int DCM2QuaterZYX( float * Cnb, float * q );

// 四元数转旋转矩阵(Z-Y-X顺序)
int Quter2DCMZYX( float * q, float * Cnb );

// 将四元数src的每个元素取反，赋值给dst，dst = -src
int OppositeQuater( float * src, float * dst );

// 四元数的共轭( dst = conj(src) )
int QuaterConj( float * src, float * dst );

// 四元数的模
int QuaterModulus( float * q, float * QuaterNorm );

// 四元数归一化
int QuaterNormilize( float * q );

// 四元数乘法( r = q * p )
int QuaterMultiply( float * q, float * p, float * r );

// 求原向量p经过四元数q旋转变换后的向量r,r = q * p * q^(-1)
int CalCoordQuater( float * q, float * p, float * r );

// 求原向量p经过四元数q旋转逆变换后的向量r,r = q^(-1) * p * q
int CalInvCoordQuater( float * q, float * p, float * r );

// 姿态运动学方程
int AttitudeDynamics( unsigned int n, unsigned int m, float t, const float * x, const float * u, float * dx, void * pData );


#ifdef __cplusplus
}
#endif

#endif     // _BASIC_FUNCTION_H


