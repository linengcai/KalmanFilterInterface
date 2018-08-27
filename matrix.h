#ifndef __MATRIX_H__
#define __MATRIX_H__

/* �㷨���� */
#define STEP_TYPE 1         // �Ʋ�
#define GPS_TYPE 2          // GPS 
#define SWIM_TYPE 3         // ��Ӿ
#define WRIST_TYPE 4        // ����
#define PDR_TYPE 5          // PDR
#define VO2MAX_TYPE 6       // VO2Max
#define BUCKLE_TYPE 7       // ���ܿ�
#define RADAR_TYPE 8        // �״�
#define IMAGE_TYPE 9        // ͼ����
#define MAGCALI_TYPE 10     // ������У׼
#define UPSTAIR_TYPE 11     // ��¥

/* ���������� */
#define NULL_STRING -1      // �����ָ��
#define ZERO_DIVISOR -2     // ����
#define MATRIX_DIME_ERR -3  // ����ά�����󣬿����ھ���˷�����ʼ����λ������з���
#define MATRIX_SINGULAR -4  // ������� 
#define DATA_FLOOD -5       // ������������ܷ���������ת����ʹ�÷����Ǻ���������������� 
#define ABNORMAL_INPUT -6   // �����쳣�����������˺�������������ֵ

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
