// ������
#ifndef _CONTROLLER_H
#define _CONTROLLER_H

#ifdef __cplusplus
extern "C"
{
#endif
	
// ��Ԥ���Ķ�ջ�з���float���飬��СΪn
float * FloatAllocate( unsigned int n );

// ΢�ַ��̺������ͣ�n��ϵͳ�״Σ�m��������ά�ȣ�t���Ա���
// x��״̬����u����������dx��״̬����΢�֣�pData�����ܵĶ������
typedef int (*DifEquType)(unsigned int n, unsigned int m, float t, const float *x, const float *u, float *dx, void *pData);

// ���Ľ��������������ϵͳ״̬���Ա���������Ϊstep
// pEquFun��״̬���̺���ָ�룬n��״̬��ά����m��������ά��
// t���Ա�����x��״̬����u����������pData��״̬���̵Ķ��������step�����ֲ���
int RK4DifSolver( DifEquType pEquFun, unsigned int n, unsigned int m, float t, float *x, const float *u, void *pData, float step );

// ����ʱ�䵥�뵥������ϵͳ�ṹ�壬�ô��ݺ�����־
typedef struct
{
	unsigned int n;   // ϵͳ�״Σ���ĸ����ʽ�Ľ״Σ�
	unsigned int m;   // ���Ӷ���ʽ�Ľ״�
	float *den;   // ��ĸ����ʽ��ϵ���������������У�
	float *num;   // ���Ӷ���ʽ��ϵ���������������У�
	float *state;  // ״̬����
	float staticGain;  // ��̬����
}CLTISISO;

// CLTISISO��״̬���̣�n��ϵͳ�״Σ�m��������ά�ȣ�t���Ա���
// x��״̬����u����������dx��״̬����΢�֣�pData��ָ��CLTISISO����
int StateEquCLTISISO( unsigned int n, unsigned int m, float t, const float * x, const float * u, 
											float * dx, void * pData );

// ��pCLTI��ʼ��Ϊ��tf��num��den��������״̬��������Ϊ������
// m��ϵͳ���ݺ������Ӷ���ʽ�״Σ�n����ĸ����ʽ�״�
int InitCLTISISO( CLTISISO * pCLTI, unsigned int m, const float * num, unsigned int n, const float * den );

// pCLTI��ʼ��Ϊn�׽�ֹƵ��Ϊomgc�İ�����˹��ͨ�˲���
int InitCLTIAsButtLPF( CLTISISO * pCLTI, unsigned int n, float omgc );

// ���ݵ�ǰ״̬��������u������pCLTI��һ�ĵ����
// pCLTI��SISOϵͳ��u����������step������
int CalCLTIOutput( CLTISISO * pCLTI, float u, float step, float * response );


#ifdef __cplusplus
}
#endif

#endif     // _CONTROLLER_H


