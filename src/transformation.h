#pragma once

#include <array>
#include <iostream>
#include <vector>
#include <cmath>

#include "extra_tools.h"

//using namespace std;

using Vector = array<double,3>;
using Matrix = array<array<double,3>,3>;

struct Kepler_elements;
struct Dec_kordinates;

array<double,6> From_Dec_to_Kep(const vector<double>& r,const vector<double>& v, double mu);

array<double,6> From_Kep_to_Dec(vector<double>& kep, double mu);
