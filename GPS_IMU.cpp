// GPS_IMU.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <stdio.h>

#include "BasicFunction.h"
#include "kalman.h"
#include "matrix.h"
#include "ButterWorthLPF.h"

char p_data[1200];

static int low = 0, hdopCount = 0;
static int lowIndex = 0;
static KalmanFilter f;
static char fusionInitFlag = 0; 
extern unsigned int SpaceOccupied;


void update_observe(const KalmanFilter f, double x, double y, double z)
{
	//set_seconds_timestep(f, theta_imu, dis);
	f.observation.data[0][0] = x;
	f.observation.data[1][0] = y;
	f.observation.data[2][0] = z;
	update_kalman(f);
}

// 由 hdop 更新 两个协方差阵
void adjustObservatioNoise(float hdop, KalmanFilter f)
{
	float pos = 0.000001f;
	float pos1 = 0.000001f;
	float mutiply = 0.01f;
	int err_code;

	if (hdop <= 1.0f && hdop > 1.0e-6f)
	{
		if (low) {
			mutiply = 5.0f;
			pos1 = 0.1f;
		}
		else {
			mutiply = 0.02f;
			pos1 = 0.1f;
		}
	}
	else if (hdop > 1.0f && hdop <= 1.5f && !low)
	{
		if (low) {
			mutiply = 10.0f;
		}
		else {
			mutiply = 2.0f;
		}
		pos1 = 0.00001f;
	}
	else
	{
		mutiply = 2000.0f;
		pos1 = 0.00001f;
		low = 1;
	}

	if (hdop >= 1.3f) {
		hdopCount++;
	}
	else {
		hdopCount--;
	}
	if (hdopCount < 0) {
		hdopCount = 0;
	}
	if (hdopCount >= 3) {
		low = 1;
		hdopCount = 0;
	}

	if (low) {
		lowIndex++;
		if (lowIndex >= 10) {
			low = 0;
			lowIndex = 0;
		}
	}
	//	mutiply =0.01f;
	//	pos1 = 10.0f;
	set_identity_matrix(f.observation_noise_covariance);
	scale_matrix(f.observation_noise_covariance, pos * mutiply);
	/*	set_matrix(f.observation_noise_covariance,
			pos * mutiply, 0.0f, 0.0f,
			0.0f, pos * mutiply, 0.0f,
			0.0f, 0.0f, pos * mutiply);*/
	err_code = set_identity_matrix(f.process_noise_covariance);
	scale_matrix(f.process_noise_covariance, pos1);
	/*   set_matrix(f.process_noise_covariance,
			 pos1, 0.0f, 0.0f,
			 0.0f, pos1, 0.0f,
			 0.0f, 0.0f, pos1);*/
}

void set_seconds_timestep(const KalmanFilter f, double theta, float dis)
{
	/* unit_scaler accounts for the relation between position and
	   velocity units */
	   // double unit_scaler = 0.001;
	double steplength;

	steplength = (double)dis / 100.0;

	f.state_transition.data[0][2] = steplength * cos(theta);
	f.state_transition.data[1][2] = -steplength * sin(theta);
}

// 初始化 kalman滤波器，不同的滤波器互不相同
KalmanFilter init_kalman_filter(double noise, double theta, float dis)
{
	/* The state model has four dimensions:
	   x, y, x', y'
	   Each time step we can only observe position, not velocity, so the
	   observation vector has only two dimensions.
	*/
	KalmanFilter f;

	// p_data 清空
	SpaceOccupied = 0;
	// 状态向量维数 量测量维数 给kalman滤波器分配空间
	alloc_filter(p_data, 3, 3, &f);

	/* Assuming the axes are rectilinear does not work well at the
	   poles, but it has the bonus that we don't need to convert between
	   lat/long and more rectangular coordinates. The slight inaccuracy
	   of our physics model is not too important.
	 */
	// 设置 F，也就是状态转移矩阵
	set_identity_matrix(f.state_transition);
	set_seconds_timestep(f, theta, dis);

	/* We observe (x, y) in each time step */
	// 设置 H，量测模型矩阵
	set_identity_matrix(f.observation_model);
	/*
  set_matrix(f.observation_model,
			 1.0, 0.0, 0.0,
			 0.0, 1.0, 0.0,
			 0.0, 0.0, 1.0);
*/
/* Noise in the world. */
	// 设置Q，也就是过程噪声阵
	double pos = 0.0000001;
	set_identity_matrix(f.process_noise_covariance);
	scale_matrix(f.process_noise_covariance, pos);
	/*
  set_matrix(f.process_noise_covariance,
		 pos, 0.0, 0.0,
		 0.0, pos, 0.0,
		 0.0, 0.0, pos);
*/
/* Noise in our observation */
	// 设置R, 也就是量测噪声
	set_identity_matrix(f.observation_noise_covariance);
	scale_matrix(f.observation_noise_covariance, pos * noise);
	/*
  set_matrix(f.observation_noise_covariance,
		 pos * noise, 0.0, 0.0,
		 0.0, pos * noise, 0.0,
		 0.0, 0.0, pos * noise);
*/
/* The start position is totally unknown, so give a high variance */
	// 初始化状态量
	set_zero_matrix(f.state_estimate);
	//  set_matrix(f.state_estimate, 0.0, 0.0, 0.0); 
	// 初始化后验估计后的协方差阵
	set_identity_matrix(f.estimate_covariance);
	double trillion = 50.0f;
	scale_matrix(f.estimate_covariance, trillion);

	return f;
}

void initFusionPara()
{
	low = hdopCount = 0;
	lowIndex = 0;
	fusionInitFlag = 0;
}

int main()
{
    std::cout << "Hello World!\n"; 
	int	err_code; // 代码报错


	// 使用butterworth滤波器
	float data_in = 0.0f;
	float data_out = 0.0f;
	static unsigned char butterInitOK = 0;
	static CLTISISO accLPF;
	if (!butterInitOK) {
		err_code = InitCLTIAsButtLPF(&accLPF, 4, 3.0f * 6.28318530717959f);
		butterInitOK = 1;
	}
	err_code = CalCLTIOutput(&accLPF, data_in, 0.04f, &data_out);

	// 使用kalman滤波器
	if (!fusionInitFlag)
	{
		f = init_kalman_filter(50, 1, 1);
		f.state_estimate.data[2][0] = 0.0f; // 
		f.state_estimate.data[1][0] = 0.0f;
		f.state_estimate.data[0][0] = 0.0f;
		fusionInitFlag = 1;
	}
	
	float hdop = 5;
	adjustObservatioNoise(hdop, f);
	update_observe(f, 1, 1, 1);

}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
