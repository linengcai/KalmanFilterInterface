/* Matrix math. */

//#include "stdafx.h"
#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
#include "matrix.h"

/* total allocated size */
unsigned int SpaceOccupied = 0;

/* allocate memory for an array(float) */
float * alloc_float_array( char *p_data, unsigned int n ) {
	float * p = NULL;
	if( p_data == NULL ) return p;
	unsigned int occupySize = n * sizeof( float );
	p = (float *)( p_data + SpaceOccupied );
	for( unsigned int i = 0; i < n; i++ ) *(p+i) = 0.0f;
	SpaceOccupied += occupySize;
	return p;
}

/* allocate memory for an array(double) */
double * alloc_array( char *p_data, unsigned int n ) {
	double * p = NULL;
	if( p_data == NULL ) return p;
	unsigned int occupySize = n * sizeof( double );
	p = (double *)( p_data + SpaceOccupied );
	for( unsigned int i = 0; i < n; i++ ) *(p+i) = 0.0f;
	SpaceOccupied += occupySize;
	return p;
}

/* This could be reduced to a single malloc if it mattered. */
Matrix alloc_matrix(char *p_data, int rows, int cols) {
  Matrix m;
  m.rows = rows;
  m.cols = cols;
  if( p_data == NULL ) { m.data = NULL; return m; }
  unsigned int occupySize = sizeof(double*) * m.rows;
  m.data = (double**)( p_data + SpaceOccupied );
  SpaceOccupied += occupySize;
  //m.data = (double**) malloc(sizeof(double*) * m.rows);
  for (int i = 0; i < m.rows; ++i) {
    occupySize = sizeof(double) * m.cols;
	m.data[i] = (double*)( p_data + SpaceOccupied );
	SpaceOccupied += occupySize;
	//m.data[i] = (double*) malloc(sizeof(double) * m.cols);
    for (int j = 0; j < m.cols; ++j) {
      m.data[i][j] = 0.0f;
    }
  }
  return m;
}

// 0¾ØÕó
void set_zero_matrix(Matrix m) {
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
			m.data[i][j] = 0.0;
    }
  }
}
// µ¥Î»¾ØÕó
int set_identity_matrix(Matrix m) {
	if( m.rows != m.cols ) return MATRIX_DIME_ERR;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      if (i == j) {
				m.data[i][j] = 1.0;
      } else {
				m.data[i][j] = 0.0;
      }
    }
  }
return 0;
}

int add_matrix(Matrix a, Matrix b, Matrix c) {
	if( a.rows != b.rows || a.rows != c.rows ||
			a.cols != b.cols || a.cols != c.cols ) {
		return MATRIX_DIME_ERR;
	}
  for (int i = 0; i < a.rows; ++i) {
    for (int j = 0; j < a.cols; ++j) {
      c.data[i][j] = a.data[i][j] + b.data[i][j];
    }
  }
return 0;
}

int subtract_matrix(Matrix a, Matrix b, Matrix c) {
	if( a.rows != b.rows || a.rows != c.rows ||
			a.cols != b.cols || a.cols != c.cols ) {
		return MATRIX_DIME_ERR;
	}
  for (int i = 0; i < a.rows; ++i) {
    for (int j = 0; j < a.cols; ++j) {
      c.data[i][j] = a.data[i][j] - b.data[i][j];
    }
  }
return 0;
}

// I - a
int subtract_from_identity_matrix(Matrix a) {
	if( a.rows != a.cols ) return MATRIX_DIME_ERR;
  for (int i = 0; i < a.rows; ++i) {
    for (int j = 0; j < a.cols; ++j) {
      if (i == j) {
				a.data[i][j] = 1.0 - a.data[i][j];
      } else {
				a.data[i][j] = 0.0 - a.data[i][j];
      }
    }
  }
return 0;
}

int multiply_matrix(Matrix a, Matrix b, Matrix c) {
	if( a.cols != b.rows || a.rows != c.rows ||
			b.cols != c.cols ) {
		return MATRIX_DIME_ERR;	
	}
  for (int i = 0; i < c.rows; ++i) {
    for (int j = 0; j < c.cols; ++j) {
      /* Calculate element c.data[i][j] via a dot product of one row of a
	 with one column of b */
      c.data[i][j] = 0.0;
      for (int k = 0; k < a.cols; ++k) {
				c.data[i][j] += a.data[i][k] * b.data[k][j];
      }
    }
  }
return 0;
}

/* This is multiplying a by b-tranpose so it is like multiply_matrix
   but references to b reverse rows and cols. */
int multiply_by_transpose_matrix(Matrix a, Matrix b, Matrix c) {
	if( a.cols != b.cols || a.rows != c.rows || 
		  b.rows != c.cols ) {
		return MATRIX_DIME_ERR; 
	}
  for (int i = 0; i < c.rows; ++i) {
    for (int j = 0; j < c.cols; ++j) {
      /* Calculate element c.data[i][j] via a dot product of one row of a
	 with one row of b */
      c.data[i][j] = 0.0;
      for (int k = 0; k < a.cols; ++k) {
				c.data[i][j] += a.data[i][k] * b.data[j][k];
      }
    }
  }
return 0;
}

int scale_matrix(Matrix m, double scalar) {
	double abs_scalar = scalar >= 0.0 ? scalar : -scalar;
	if( abs_scalar < 1.0e-8 ) return ABNORMAL_INPUT;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      m.data[i][j] *= scalar;
    }
  }
return 0;
}

int swap_rows(Matrix m, int r1, int r2) {
	if( r1 == r2 ) return ABNORMAL_INPUT;
  double* tmp = m.data[r1];
  m.data[r1] = m.data[r2];
  m.data[r2] = tmp;
return 0;
}

int scale_row(Matrix m, int r, double scalar) {
	double abs_scalar = scalar >= 0.0 ? scalar : -scalar;
	if( abs_scalar < 1.0e-8 ) return ABNORMAL_INPUT;
  for (int i = 0; i < m.cols; ++i) {
    m.data[r][i] *= scalar;
  }
return 0;
}

/* Add scalar * row r2 to row r1. */
int shear_row(Matrix m, int r1, int r2, double scalar) {
	if( r1 == r2 ) return ABNORMAL_INPUT;
  for (int i = 0; i < m.cols; ++i) {
    m.data[r1][i] += scalar * m.data[r2][i];
  }
return 0;
}

/* Uses Gauss-Jordan elimination.

   The elimination procedure works by applying elementary row
   operations to our input matrix until the input matrix is reduced to
   the identity matrix.
   Simultaneously, we apply the same elementary row operations to a
   separate identity matrix to produce the inverse matrix.
   If this makes no sense, read wikipedia on Gauss-Jordan elimination.
   
   This is not the fastest way to invert matrices, so this is quite
   possibly the bottleneck. */
// ¾ØÕóÇóÄæ¼ÆËã »áÆÆ»µÊäÈëÁ¿
int destructive_invert_matrix(Matrix input, Matrix output) {
	int err_code = 0;
    err_code = set_identity_matrix(output);

  /* Convert input to the identity matrix via elementary row operations.
     The ith pass through this loop turns the element at i,i to a 1
     and turns all other elements in column i to a 0. */
  for (int i = 0; i < input.rows; ++i) {
    if (input.data[i][i] == 0.0) {
      /* We must swap rows to get a nonzero diagonal element. */
      int r;
      for (r = i + 1; r < input.rows; ++r) {
				if (input.data[r][i] != 0.0) {
					break;
				}
      }
      if (r == input.rows) {
	/* Every remaining element in this column is zero, so this
	   matrix cannot be inverted. */
				return MATRIX_SINGULAR;
      }
      err_code = swap_rows(input, i, r);

      err_code = swap_rows(output, i, r);

    }

    /* Scale this row to ensure a 1 along the diagonal.
       We might need to worry about overflow from a huge scalar here. */

    double scalar = 1.0 / input.data[i][i];
    err_code = scale_row(input, i, scalar);

    err_code = scale_row(output, i, scalar);


    /* Zero out the other elements in this column. */
    for (int j = 0; j < input.rows; ++j) {
      if (i == j) {
	continue;
      }
      double shear_needed = -input.data[j][i];
      err_code = shear_row(input, j, i, shear_needed);

      err_code = shear_row(output, j, i, shear_needed);

    }
  }
  
  return 0;
}
