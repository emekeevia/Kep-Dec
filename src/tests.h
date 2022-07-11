#pragma once

#include <iostream>
#include <cassert>
#include <sstream>
#include <exception>
#include <algorithm>
#include "extra_tools.h"
#include "transformation.h"

using namespace std;



class TestRunner { // класс тестирования
	public:
	template <class TestFunc>
	void RunTest(TestFunc func , const string& test_name) {
			try {
				func();
				cerr << test_name << " OK" << endl;
			} catch (runtime_error& e) {
				++fail_count; // увеличиваем счётчик упавших тестов
				cerr << test_name << " fail: " << e.what() << endl;
			}
		}
	~TestRunner() { // деструктор класса TestRunner, в котором анализируем fail_count
		if (fail_count > 0) { //это как раз тот момент, когда
			cerr << fail_count << " unit tests failed. Terminate" << endl;
			exit(1); // завершение программы с кодом возврата 1
		}
	}
		private:
		int fail_count = 0; // счётчик числа упавших тестов
};

template <typename T, long long unsigned int l>
T rel_error(const array<T,l>& t, const array<T,l>& u){
	T er = 0;
	for(size_t i = 0;i < t.size(); i++){
		if(u[i] * t[i] != 0.0){
			er = max(er, abs((t[i]-u[i])/t[i]));
		}else{
			er = max(er,max(abs(t[i]), abs(u[i])));
		}
	}
	return er;
}



//template <typename T>
//void AssertEqual (const vector<T>& t, const vector<T>& u) {
//	// значения двух разных типов для удобства
//	T epsilon = rel_error(t, u);
//	if (epsilon > 0.02) { // если значения не равны, то мы даём знать, что этот assert не сработал
//		ostringstream os;
//		os << "Assertion failed: " << t << "!=" << u;
//		throw runtime_error(os.str()); // бросим исключение с сообщением со значениями t и u
//	}
//}

template <typename T,long long unsigned int l>
void AssertEqual (const array<T,l>& t, const array<T,l>& u) {
	// значения двух разных типов для удобства
	T epsilon = rel_error(t, u);
	if (epsilon > 0.02) { // если значения не равны, то мы даём знать, что этот assert не сработал
		ostringstream os;
		os << "Assertion failed: " << t << "!=" << u << " " << epsilon;
		throw runtime_error(os.str()); // бросим исключение с сообщением со значениями t и u
	}
}

void Test_Dec_to_Kep();
void Test_Kep_to_Dec();
