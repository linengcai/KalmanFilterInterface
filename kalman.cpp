/* Kalman filters. */

#include "kalman.h"

int alloc_filter(char *p_data, int state_dimension,
			  int observation_dimension, KalmanFilter * f) {
	
	// 时间戳 状态量维数 观测量维数
  f->timestep = 0;
  f->state_dimension = state_dimension;
  f->observation_dimension = observation_dimension;

  // F 状态转移矩阵 
  f->state_transition = alloc_matrix(p_data, state_dimension,
				    state_dimension);
	// H 观测模型 				
  f->observation_model = alloc_matrix(p_data, observation_dimension,
				     state_dimension);
	// Q 过程噪声				
  f->process_noise_covariance = alloc_matrix(p_data, state_dimension,
					    state_dimension);
	// R 量测噪声
  f->observation_noise_covariance = alloc_matrix(p_data, observation_dimension,
						observation_dimension);
	// Z 量测
  f->observation = alloc_matrix(p_data, observation_dimension, 1);
	// 先验估计状态量
  f->predicted_state = alloc_matrix(p_data, state_dimension, 1);
	// 先验估计后协方差阵
  f->predicted_estimate_covariance = alloc_matrix(p_data, state_dimension,
						 state_dimension);
	// 量测残差					 
  f->innovation = alloc_matrix(p_data, observation_dimension, 1);
	// 残差协方差
  f->innovation_covariance = alloc_matrix(p_data, observation_dimension,
					 observation_dimension);
	// 残差协方差求逆				 
  f->inverse_innovation_covariance = alloc_matrix(p_data, observation_dimension,
						 observation_dimension);
	// kalman增益					 
  f->optimal_gain = alloc_matrix(p_data, state_dimension,
				observation_dimension);
	// 状态量的后验估计结果			
  f->state_estimate = alloc_matrix(p_data, state_dimension, 1);
	// 更新后验估计协方差
  f->estimate_covariance = alloc_matrix(p_data, state_dimension,
				       state_dimension);

  // 中间计算量，无特殊意义
  f->vertical_scratch = alloc_matrix(p_data, state_dimension,
				    observation_dimension);
						
  f->small_square_scratch = alloc_matrix(p_data, observation_dimension,
					observation_dimension);
					
  f->big_square_scratch = alloc_matrix(p_data, state_dimension,
				      state_dimension);
  
  return 0;
}

/* update kalman filter */
int update_kalman(KalmanFilter f) {
	int err_code = 0;
	err_code = predict(f);
	err_code = estimate(f);
	return 0;
}

int predict(KalmanFilter f) {
	int err_code = 0;
  f.timestep++;

  /* Predict the state */
  // x = f * x
  err_code = multiply_matrix(f.state_transition, f.state_estimate,
		  f.predicted_state);

  /* Predict the state estimate covariance */
  // f * P
  err_code = multiply_matrix(f.state_transition, f.estimate_covariance,
		  f.big_square_scratch);
  // f * P * f'
  err_code = multiply_by_transpose_matrix(f.big_square_scratch, f.state_transition,
			       f.predicted_estimate_covariance);
  // p1 = f * P * f' + Q
  err_code = add_matrix(f.predicted_estimate_covariance, f.process_noise_covariance,
	     f.predicted_estimate_covariance);

	
	return 0;
}

int estimate(KalmanFilter f) {
	int err_code = 0;
  /* Calculate innovation */
	// H * X
  err_code = multiply_matrix(f.observation_model, f.predicted_state,
		  f.innovation);
	// Y = z - h * x
  err_code = subtract_matrix(f.observation, f.innovation,
		  f.innovation);


  /* Calculate innovation covariance */
  // P1 * H'
  err_code = multiply_by_transpose_matrix(f.predicted_estimate_covariance,
			       f.observation_model,
			       f.vertical_scratch);

  // Sk = H * P1 * H'
  err_code = multiply_matrix(f.observation_model, f.vertical_scratch,
		  f.innovation_covariance);

  // Sk = H * P1 * H' + R
  err_code = add_matrix(f.innovation_covariance, f.observation_noise_covariance,
	     f.innovation_covariance);


  /* Invert the innovation covariance.
     Note: this destroys the innovation covariance.
     TODO: handle inversion failure intelligently. */
  // Sk 求逆
  err_code = destructive_invert_matrix(f.innovation_covariance,
			    f.inverse_innovation_covariance);

  
  /* Calculate the optimal Kalman gain.
     Note we still have a useful partial product in vertical scratch
     from the innovation covariance. */
  // K = P1 * H * Sk^-1
  err_code = multiply_matrix(f.vertical_scratch, f.inverse_innovation_covariance,
		  f.optimal_gain);


  /* Estimate the state */
  // K * y  y = z - H * x
  err_code = multiply_matrix(f.optimal_gain, f.innovation,
		  f.state_estimate);
  // x = x + K * (z - H * x)
  err_code = add_matrix(f.state_estimate, f.predicted_state,
	     f.state_estimate);


  /* Estimate the state covariance */
  // K * H
  err_code = multiply_matrix(f.optimal_gain, f.observation_model,
		  f.big_square_scratch);
  // I - K * H
  err_code = subtract_from_identity_matrix(f.big_square_scratch);

  // P = (I - K * H) * P1
  err_code = multiply_matrix(f.big_square_scratch, f.predicted_estimate_covariance,
		  f.estimate_covariance);
	
	return 0;
}



