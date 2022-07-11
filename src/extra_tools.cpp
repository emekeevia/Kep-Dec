#include "extra_tools.h"

Vector matr_vect_mult(Matrix& m, Vector& v){
	Vector temp;
	for(size_t i = 0; i < 3;i++){
		temp[i] = m[i][0]*v[0] + m[i][1]*v[1] + m[i][2]*v[2];
	}
	return temp;
}





