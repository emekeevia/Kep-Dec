#pragma once

#include <iostream>
#include <cmath>
#include <array>
#include <sstream>
#include <vector>

using namespace std;

//using Vector = vector<double>;
//using Matrix = vector<vector<double>>;
using Vector = array<double,3>;
using Matrix = array<array<double,3>,3>;

Vector matr_vect_mult(Matrix& m, Vector& v);

template<typename T>
void grad_to_rad(T& kep){
	kep[0] = kep[0]*M_PI/180;
	kep[1] = kep[1]*M_PI/180;
	kep[4] = kep[4]*M_PI/180;
	kep[5] = kep[5]*M_PI/180;
}


template<typename T>
T round(T element, T round_elem){
	if(abs(element - round_elem) < 0.01){
		return 0.0;
	}
	return element;
}

template<typename T>
T operator+(T v1, T v2){
	if(v1.size() != v2.size()){
		throw("Not equal size!");
	}
	size_t count = v2.size();
	T temp;
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]+v2[i];
	}
	return temp;
}

template<typename T>
vector<T> operator-(vector<T>& v1, vector<T>& v2){
	if(v1.size() != v2.size()){
		throw("Not equal size!");
	}
	size_t count = v1.size();
	T temp(count);
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]-v2[i];
	}
	return temp;
}

template<typename T, int l>
array<T,l> operator-(array<T,l>& v1, array<T,l>& v2){
	size_t count = v1.size();
	array<T,l> temp;
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]-v2[i];
	}
	return temp;
}

template<typename T,typename D>
Vector operator/(const T& v, D a){
	Vector temp;
	if(a != 0){
		for(size_t i = 0; i < v.size();i++){
			temp[i] = v[i]/a;
		}
	}
	return temp;
}

template<typename T, typename D>
T operator*(const T& v, D a){
	T temp;
	for(size_t i = 0; i < 3;i++){
		temp[i] = v[i]*a;
	}
	return temp;
}


template<typename T>
double abs_vect(const T& v){
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


template<typename T>
Vector vect_mult(const T& v1,const T& v2){
	Vector temp;
	temp[0] = v1[1] * v2[2] - v1[2] * v2[1];
	temp[1] = v1[2] * v2[0] - v1[0] * v2[2];
	temp[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return temp;
}

template<typename T>
double scal_mult(const T& v1,const T& v2){
	return static_cast<double>(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
	os << "{";
	bool flag = false;
	for(auto a: v){
		if(!flag){
			os << a;
			flag = true;
		}else{
			os << ", " << a;
		}
	}
	os << "}";
	return os;
}

template<typename T, long long unsigned int l>
ostream& operator<<(ostream& os, const array<T,l>& v){
	os << "[";
	bool flag = false;
	for(auto a: v){
		if(!flag){
			os << a;
			flag = true;
		}else{
			os << ", " << a;
		}

	}
	os << "]";
	return os;
}


