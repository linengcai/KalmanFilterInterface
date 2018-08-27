#include <stdio.h>
#include <math.h>
//#include "BasicFunction.h"
//#include "algo_assert.h"
//#include "PDRTrack_Interface.h"

// �������
int CrossProduct(float * Y, float * Z, float * YXZ)
{
	YXZ[0] = (Y[1] * Z[2] - Y[2] * Z[1]);
	YXZ[1] = (Y[2] * Z[0] - Y[0] * Z[2]);
	YXZ[2] = (Y[0] * Z[1] - Y[1] * Z[0]);

	return 0;
}

// �������
int DotProduct( float * Y, float * Z, float * Dot )
{
	float sum = 0.0f;
	for( int i = 0; i != 3; i++ )
	{
		sum += Y[i] * Z[i];
	}
	*Dot = sum;

	return 0;
}

// Fast inverse square-root
int invSqrt( float x, float * inv )
{
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * ( 1.5f - ( halfx * y * y ) );
	y = y * ( 1.5f - ( halfx * y * y ) );
	*inv = y;

	return 0;
}

// �������, C = A * B, A: m * nά, B: n * pά, C: m * pά
int MultiplyMat( unsigned int m, unsigned int n, unsigned int p, float * A, float * B, float * C )
{
	unsigned int i, j, k;

	for( i = 0; i != m; ++i )
	{
		for( j = 0; j != p; ++j )
		{
			C[i*p+j] = 0.0f;
			for( k = 0; k != n; ++k )
			{
				C[i*p+j] += A[i*n+k] * B[k*p+j];
			}
		}
	}

	return 0;
}

// ����ģ
int CalVectorNorm( float * Y, float * VectNorm )
{
	float sum = 0.0f;

	for( unsigned char i = 0; i != 3; i++ )
	{
		sum += Y[i] * Y[i];
	}
	*VectNorm = sqrtf(sum);

	return 0;
}

// �������ƺ���
float AmpConstrain( float x, float xMin, float xMax )
{
	float y = x > xMin ? x : xMin;
	y = y < xMax ? y : xMax;
	return y;
}

// �����䱣���ķ����Һ���
float ArcSinSecure( float x )
{
	const float MaxArg = 0.99999999f;
	return asinf( AmpConstrain( x, -MaxArg, MaxArg ) );
}

// �����䱣���ķ����Һ���
float ArcCosSecure( float x )
{
	const float MaxArg = 0.99999999f;
	return acosf( AmpConstrain( x, -MaxArg, MaxArg ) );
}

// ����Ԫ������Ϊ��λ��Ԫ��
int SetQuaterAsUnit( float * q )
{
	q[0] = 1.0f;
	q[1] = q[2] = q[3] = 0.0f;

	return 0;
}

// ��Ԫ��תŷ���ǣ�Z-X-Y˳��
int Quater2EulerZXY( float * q, float * euler)
{
	float q0q0 = q[0]*q[0];
	float q1q1 = q[1]*q[1];
	float q2q2 = q[2]*q[2];
	float q3q3 = q[3]*q[3];

	euler[0] = ArcSinSecure( 2.0f * ( q[2] * q[3] + q[0] * q[1] ) );
	euler[1] = atan2f( 2.0f * ( q[0] * q[2] - q[1] * q[3] ), q0q0 - q1q1 - q2q2 + q3q3 );
	euler[2] = atan2f( 2.0f * ( q[0] * q[3] - q[1] * q[2] ), q0q0 - q1q1 + q2q2 - q3q3 );

	return 0;
}

// ��Ԫ��תŷ���ǣ�Z-Y-X˳��
int Quater2EulerZYX( float * q, float * euler)
{
	float q0q0 = q[0]*q[0];
	float q1q1 = q[1]*q[1];
	float q2q2 = q[2]*q[2];
	float q3q3 = q[3]*q[3];

	euler[0] = atan2f( 2.0f * ( q[2] * q[3] + q[0] * q[1] ), q0q0 - q1q1 - q2q2 + q3q3 );
	euler[1] = ArcSinSecure( 2.0f * ( q[0] * q[2] - q[1] * q[3] ) );
	euler[2] = atan2f( 2.0f * ( q[1] * q[2] + q[0] * q[3] ), q0q0 + q1q1 - q2q2 - q3q3 );

	return 0;
}

// ��ת����ת��Ԫ��(Z-Y-X˳��)
int DCM2QuaterZYX( float * Cnb, float * q )
{
	float s = sqrtf( Cnb[0] + Cnb[4] + Cnb[8] + 1 );
	q[0] = 0.5f * s;
	s = 0.5f / s;
	q[1] = (Cnb[5] - Cnb[7]) * s;
	q[2] = (Cnb[6] - Cnb[2]) * s;
	q[3] = (Cnb[1] - Cnb[3]) * s;

	return 0;
}

// ��Ԫ��ת��ת����(Z-Y-X˳��)
int Quter2DCMZYX( float * q, float * Cnb )
{
	float q0 = q[0];
	float q1 = q[1];
	float q2 = q[2];
	float q3 = q[3];
	
	float q0q0 = q0 * q0;
	float q0q1 = q0 * q1;
	float q0q2 = q0 * q2;
	float q0q3 = q0 * q3;
	float q1q1 = q1 * q1;
	float q1q2 = q1 * q2;
	float q1q3 = q1 * q3;
	float q2q2 = q2 * q2;
	float q2q3 = q2 * q3;
	float q3q3 = q3 * q3;

	Cnb[0] = q0q0 + q1q1 - q2q2 - q3q3;
	Cnb[1] = 2.0f * ( q1q2 + q0q3 );
	Cnb[2] = 2.0f * ( q1q3 - q0q2 );
	Cnb[3] = 2.0f * ( q1q2 - q0q3 );
	Cnb[4] = q0q0 - q1q1 + q2q2 - q3q3;
	Cnb[5] = 2.0f * ( q2q3 + q0q1 );
	Cnb[6] = 2.0f * ( q0q2 + q1q3 );
	Cnb[7] = 2.0f * ( q2q3 - q0q1 );
	Cnb[8] = q0q0 - q1q1 - q2q2 + q3q3;

	return 0;
}

// ����Ԫ��src��ÿ��Ԫ��ȡ������ֵ��dst��dst = -src
int OppositeQuater( float * src, float * dst )
{
	unsigned char i;

	for( i = 0; i != 4; i++ )
	{
		dst[i] = -src[i];
	}

	return 0;
}

// ��Ԫ���Ĺ���( dst = conj(src) )
int QuaterConj( float * src, float * dst )
{
	dst[0] =  src[0];
	dst[1] = -src[1];
	dst[2] = -src[2];
	dst[3] = -src[3];

	return 0;
}

// ��Ԫ����ģ
int QuaterModulus( float * q, float * QuaterNorm )
{
	float y = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
	*QuaterNorm = sqrtf( y );

	return 0;
}

// ��Ԫ����һ��
int QuaterNormilize( float * q )
{
	int error_check;
	float invQuaterNorm;

	float QuaterNorm = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
	if( QuaterNorm > 1.0e-8f )
	{
		error_check = invSqrt( QuaterNorm, &invQuaterNorm );

		unsigned char i;
		for( i = 0; i != 4; i++ )
		{
			q[i] *= invQuaterNorm;
		}
	}

	return 0;
}

// ��Ԫ���˷�( r = q * p )
int QuaterMultiply( float * q, float * p, float * r )
{

	r[0] = q[0] * p[0] - q[1] * p[1] - q[2] * p[2] - q[3] * p[3];
	r[1] = q[0] * p[1] + q[1] * p[0] + q[2] * p[3] - q[3] * p[2];
	r[2] = q[0] * p[2] + q[2] * p[0] + q[3] * p[1] - q[1] * p[3];
	r[3] = q[0] * p[3] + q[3] * p[0] + q[1] * p[2] - q[2] * p[1];

	return 0;
}

// ��ԭ����p������Ԫ��q��ת�任�������r,r = q * p * q^(-1)
int CalCoordQuater( float * q, float * p, float * r )
{
	float Cnb[9];
	int error_check;

	error_check = Quter2DCMZYX( q, Cnb );

	error_check = MultiplyMat( 3, 3, 1, Cnb, p, r );

	/*float invQ[4], qp[4], a[4];
	unsigned char i;
	QuaterConj( q, invQ );
	a[0] = 0.0f;
	for( i = 0; i < 3; i++ )
	{
		a[i+1] = p[i];
	}
	QuaterMultiply( q, a, qp );
	QuaterMultiply( qp, invQ, a );
	for( i = 0; i != 3; i++ )
	{
		r[i] = a[i+1];
	}*/

	return 0;
}

// ��ԭ����p������Ԫ��q��ת��任�������r,r = q^(-1) * p * q
int CalInvCoordQuater( float * q, float * p, float * r )
{
	float invQ[4];
	int error_check;

	error_check = QuaterConj( q, invQ );

	error_check = CalCoordQuater( invQ, p, r );

	return 0;
}

// ��̬�˶�ѧ����
int AttitudeDynamics( unsigned int n, unsigned int m, float t, const float * x, const float * u, float * dx, void * pData )
{

	float q0 = x[0], q1 = x[1], q2 = x[2], q3 = x[3];
	float wx = u[0], wy = u[1], wz = u[2];

	dx[0] = 0.5f * ( -q1 * wx - q2 * wy - q3 * wz );
	dx[1] = 0.5f * ( q0 * wx - q3 * wy + q2 * wz );
	dx[2] = 0.5f * ( q3 * wx + q0 * wy - q1 * wz );
	dx[3] = 0.5f * ( -q2 * wx + q1 * wy + q0 * wz );

	return 0;
}