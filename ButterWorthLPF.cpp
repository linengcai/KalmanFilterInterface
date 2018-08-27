
#include "ButterWorthLPF.h"

#define MaxStateDim 5
#define STACKSIZE 600
static char ReservedStack[STACKSIZE];
static unsigned int OccupySpaceStack = 0;

static float m_tmp[MaxStateDim], m_dx1[MaxStateDim], m_dx2[MaxStateDim], m_dx3[MaxStateDim], m_dx4[MaxStateDim];

float * FloatAllocate( unsigned int n )
{
	float * p = 0;
	unsigned int occupySize = n * sizeof( float );
	if( OccupySpaceStack + occupySize <= STACKSIZE )
	{
		p = (float *)( ReservedStack + OccupySpaceStack );
		OccupySpaceStack += occupySize;
	}
	return p;
}

int StateEquCLTISISO( unsigned int n, unsigned int m, float t, const float * x, const float * u, 
											float * dx, void * pData )
{
	
	CLTISISO * pCLTI = (CLTISISO *)( pData );
	unsigned int i;
	// 设定 n = 4
	for( i = 0; i < pCLTI->n - 1; ++i )
	{
		// dx[0] = x[1] = 
		// dx[2] = x[3]
		dx[i] = x[i+1];
	}
	// dx[3] = u
	dx[pCLTI->n - 1] = u[0];
	for( i = 0; i != pCLTI->n; ++i )
	{
		// dx[3] = u - (den[3] * x[0] + den[2] * x[1] + den[1] * x[2] + den[0] * x[3])
		dx[pCLTI->n - 1] -= pCLTI->den[pCLTI->n - i] * x[i];
	}
	return 0;
}

int InitCLTISISO( CLTISISO * pCLTI, unsigned int m, const float * num, unsigned int n, const float * den )
{ 
	
	unsigned int i;
	float coef;
	pCLTI->staticGain = 0.0f;
	pCLTI->m = m;
	pCLTI->num = FloatAllocate( m + 1 );
	pCLTI->n = n;
	pCLTI->den = FloatAllocate( n + 1 );
	pCLTI->state = FloatAllocate( n );
	
	// num 分子系数
	for( i = 0; i <= m; ++i )
	{
		pCLTI->num[i] = num[i];
	}
	// den 分母系数
	for( i = 0; i <= n; ++i )
	{
		pCLTI->den[i] = den[i];
	}
	// state N维向量 0
	for( i = 0; i != n; ++i )
	{
		pCLTI->state[i] = 0.0f;
	}

	// den[0]  = 1
	coef = 1.0f / pCLTI->den[0];
	for( i = 0; i <= n; ++i )
	{
		pCLTI->den[i] *= coef;
	}
	for( i = 0; i <= m; ++i )
	{
		pCLTI->num[i] *= coef;
	}

	// N = 1维的情况
	if( n == m )
	{
		pCLTI->staticGain = pCLTI->num[0];
		for( i = 0; i <= m; ++i )
		{
			pCLTI->num[i] -= pCLTI->staticGain * pCLTI->den[i];
		}
		for( i = 0; i < m; ++i )
		{
			pCLTI->num[i] = pCLTI->num[i+1];
		}
		pCLTI->m = pCLTI->m - 1;
	}
	return 0;
}

// n 阶数 omgc -3db的截止频率 计算 分子 分母系数
int InitCLTIAsButtLPF( CLTISISO * pCLTI, unsigned int n, float omgc )
{ 
	
#define MAXFilterOrder 6
	int err_code = 0;
	float omgn = 1.0f;
	unsigned int i;
	float num[1], den[MAXFilterOrder + 1];
	if( n < 1 || n > MAXFilterOrder )
	{
		n = 2; 
	}
	// den 是分子 从系数高到系数低的系数，例如
	// s^4 s^3 s^2 s^1 den[0] 是分子系数
	den[0] = 1.0f;
	switch( n )
	{
	case 1:
		den[1] = 1.0f;
		break;
	case 2:
		den[1] = 1.4142f;
		den[2] = 1.0f;
		break;
	case 3:
		den[1] = 2.0f;
		den[2] = 2.0f;
		den[3] = 1.0f;
		break;
	case 4:
		den[1] = 2.6132f;
		den[2] = 3.4143f;
		den[3] = 2.6132f;
		den[4] = 1.0f;
		break;
	case 5:
		den[1] = 3.2360f;
		den[2] = 5.2359f;
		den[3] = 5.2359f;
		den[4] = 3.2360f;
		den[5] = 1.0f;
		break;
	case 6:
		den[1] = 3.8637f;
		den[2] = 7.4640f;
		den[3] = 9.1415f;
		den[4] = 7.4640f;
		den[5] = 3.8637f;
		den[6] = 1.0f;
		break;
	default:
		den[1] = 1.4142f;
		den[2] = 1.0f;
		n = 2;
		break;
	}
	// 去归一化
	for( i = 0; i <= n; ++i )
	{
		den[i] *= omgn;
		omgn *= omgc;
	}
	num[0] = den[n];
	err_code = InitCLTISISO( pCLTI, 0, num, n, den );
	return 0;
}

// 初始化好了的滤波器 输入 步长 输出
int CalCLTIOutput( CLTISISO * pCLTI, float u, float step, float *response )
{ 
	// 0
	float y = pCLTI->staticGain * u;
	float t = 0.0f;
	unsigned int i;
	int err_code = 0;

	err_code = RK4DifSolver(StateEquCLTISISO, pCLTI->n, 1, t, pCLTI->state, &u, pCLTI, step );

	// m = 0
	for( i = 0; i <= pCLTI->m; ++i )
	{
		// y = num * state
		y += pCLTI->num[pCLTI->m - i] * pCLTI->state[i];
	}
	*response = y;
	return 0;
}

int RK4DifSolver( DifEquType pEquFun, unsigned int n, unsigned int m, float t, float *x, const float *u, void *pData, float step )
{ 
	
	int err_code = 0;
	unsigned int i;
	float h_2 = 0.5f * step;
	
	// dx1 = f( t, x ) 目的输出dx
	err_code = (*pEquFun)( n, m, t, x, u, m_dx1, pData );
	
	// dx2 = f( t + h/2, x + h/2 * dx1 )
	for( i = 0; i != n; ++i )
	{
		m_tmp[i] = x[i] + h_2 * m_dx1[i];
	}
	err_code = (*pEquFun)( n, m, t + h_2, m_tmp, u, m_dx2, pData );
	
	// dx3 = f( t + h/2, x + h/2 * dx2 )
	for( i = 0; i != n; ++i )
	{
		m_tmp[i] = x[i] + h_2 * m_dx2[i];
	}
	err_code = (*pEquFun)( n, m, t + h_2, m_tmp, u, m_dx3, pData );
	
	// dx3 = f( t + h, x + h * dx3 )
	for( i = 0; i != n; ++i )
	{
		m_tmp[i] = x[i] + step * m_dx3[i];
	}
	err_code = (*pEquFun)( n, m, t + step, m_tmp, u, m_dx4, pData );

	for( i = 0; i != n; ++i )
	{
		x[i] += ( m_dx1[i] + 2.0f * ( m_dx2[i] + m_dx3[i] ) + m_dx4[i] ) * step * 0.166666666666667f;
	}
	return 0;
}


