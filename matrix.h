#ifndef __MATRIX_H__
#define __MATRIX_H__

/* 算法类型 */
#define STEP_TYPE 1         // 计步
#define GPS_TYPE 2          // GPS 
#define SWIM_TYPE 3         // 游泳
#define WRIST_TYPE 4        // 翻腕
#define PDR_TYPE 5          // PDR
#define VO2MAX_TYPE 6       // VO2Max
#define BUCKLE_TYPE 7       // 万能扣
#define RADAR_TYPE 8        // 雷达
#define IMAGE_TYPE 9        // 图像处理
#define MAGCALI_TYPE 10     // 磁力计校准
#define UPSTAIR_TYPE 11     // 爬楼

/* 错误码类型 */
#define NULL_STRING -1      // 输入空指针
#define ZERO_DIVISOR -2     // 除零
#define MATRIX_DIME_ERR -3  // 矩阵维数错误，可能在矩阵乘法，初始化单位阵过程中发生
#define MATRIX_SINGULAR -4  // 奇异矩阵 
#define DATA_FLOOD -5       // 数据溢出，可能发生在类型转换或使用反三角函数等类似情况发生 
#define ABNORMAL_INPUT -6   // 输入异常，可能输入了函数不期望的数值

typedef struct {
  /* Dimensions */
  int rows;
  int cols;

  /* Contents of the matrix */
  double** data;
} Matrix;

/* allocate memory for a new array(float) */
float * alloc_float_array(char *p_data, unsigned int n);

/* allocate memory for a new array(double) */
double * alloc_array(char *p_data, unsigned int n);

/* Allocate memory for a new matrix.
   Zeros out the matrix.
   Assert-fails if we are out of memory.
*/
Matrix alloc_matrix(char *p_data, int rows, int cols);

/* Turn m into an zero matrix. */
void set_zero_matrix(Matrix m);

/* Turn m into an identity matrix. */
int set_identity_matrix(Matrix m);

/* Add matrices a and b and put the result in c. */
int add_matrix(Matrix a, Matrix b, Matrix c);

/* Subtract matrices a and b and put the result in c. */
int subtract_matrix(Matrix a, Matrix b, Matrix c);

/* Subtract from the identity matrix in place. */
int subtract_from_identity_matrix(Matrix a);

/* Multiply matrices a and b and put the result in c. */
int multiply_matrix(Matrix a, Matrix b, Matrix c);

/* Multiply matrix a by b-transpose and put the result in c. */
int multiply_by_transpose_matrix(Matrix a, Matrix b, Matrix c);

/* Multiply a matrix by a scalar. */
int scale_matrix(Matrix m, double scalar);

/* Swap rows r1 and r2 of a matrix.
   This is one of the three "elementary row operations". */
int swap_rows(Matrix m, int r1, int r2);

/* Multiply row r of a matrix by a scalar.
   This is one of the three "elementary row operations". */
int scale_row(Matrix m, int r, double scalar);

/* Add a multiple of row r2 to row r1.
   Also known as a "shear" operation.
   This is one of the three "elementary row operations". */
int shear_row(Matrix m, int r1, int r2, double scalar);

/* Invert a square matrix.
   Returns whether the matrix is invertible.
   input is mutated as well by this routine. */
int destructive_invert_matrix(Matrix input, Matrix output);

#endif
