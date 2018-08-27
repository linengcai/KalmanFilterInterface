// 控制器
#ifndef _CONTROLLER_H
#define _CONTROLLER_H

#ifdef __cplusplus
extern "C"
{
#endif
	
// 从预留的堆栈中分配float数组，大小为n
float * FloatAllocate( unsigned int n );

// 微分方程函数类型，n：系统阶次，m：控制量维度，t：自变量
// x：状态量，u：控制量，dx：状态量的微分，pData：可能的额外参数
typedef int (*DifEquType)(unsigned int n, unsigned int m, float t, const float *x, const float *u, float *dx, void *pData);

// 用四阶龙格库塔法更新系统状态和自变量，步长为step
// pEquFun：状态方程函数指针，n：状态量维数，m：控制量维数
// t：自变量，x：状态量，u：控制量，pData：状态方程的额外参数，step：积分步长
int RK4DifSolver( DifEquType pEquFun, unsigned int n, unsigned int m, float t, float *x, const float *u, void *pData, float step );

// 连续时间单入单出线性系统结构体，用传递函数标志
typedef struct
{
	unsigned int n;   // 系统阶次（分母多项式的阶次）
	unsigned int m;   // 分子多项式的阶次
	float *den;   // 分母多项式的系数向量（降幂排列）
	float *num;   // 分子多项式的系数向量（降幂排列）
	float *state;  // 状态向量
	float staticGain;  // 静态增益
}CLTISISO;

// CLTISISO的状态方程，n：系统阶次，m：控制量维度，t：自变量
// x：状态量，u：控制量，dx：状态量的微分，pData：指向CLTISISO对象
int StateEquCLTISISO( unsigned int n, unsigned int m, float t, const float * x, const float * u, 
											float * dx, void * pData );

// 将pCLTI初始化为：tf（num，den），并将状态向量设置为零向量
// m：系统传递函数分子多项式阶次，n：分母多项式阶次
int InitCLTISISO( CLTISISO * pCLTI, unsigned int m, const float * num, unsigned int n, const float * den );

// pCLTI初始化为n阶截止频率为omgc的巴特沃斯低通滤波器
int InitCLTIAsButtLPF( CLTISISO * pCLTI, unsigned int n, float omgc );

// 根据当前状态和输入量u，计算pCLTI下一拍的输出
// pCLTI：SISO系统，u：输入量，step：步长
int CalCLTIOutput( CLTISISO * pCLTI, float u, float step, float * response );


#ifdef __cplusplus
}
#endif

#endif     // _CONTROLLER_H


