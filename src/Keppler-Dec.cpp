//============================================================================
// Name        : Keppler-Dec.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================





#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <ctime>

using namespace std;

#include "extra_tools.h"
#include "tests.h"
#include "transformation.h"



void TestAll(){
	TestRunner tr;
	tr.RunTest(Test_Dec_to_Kep, "Test_Dec_to_Kep");
	tr.RunTest(Test_Kep_to_Dec, "Test_Kep_to_Dec");
}


int main() {
	auto start = chrono::high_resolution_clock::now();
	TestAll();
	auto end = chrono::high_resolution_clock::now();
	cout  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" <<endl;
}
