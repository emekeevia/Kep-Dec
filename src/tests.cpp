#include "tests.h"

void Test_Dec_to_Kep(){
	//double mu_E = 398600.4418;//гравитационная постоянная Земли (км^{3}c^{-2})
	double mu_S  = 132712440018.0;//гравитационная постоянная Солнца

	//Тест 1 (настроен вручную)
	vector<double> r = {1.306e8, -1.979e8, 1.837e8};//км
	vector<double> v = {24.004, -6.837, -10.532};//км/c
	AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{60.0, 150.0, 8.794e8, 0.7, 90.0, 45.0});

	//Тест 2 (настроен)
	r = {-5.284e8, -6.854e8, -1.991e8};
	v = {-7.576, -2.527,-9.176};
	AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{120.0, 60.0, 8.794e8, 0.7, 60.0, 135.0});

	//Тест 3 (настроен)
	r = {-3.294e8, -8.003e8, 1.991e8};
	v = {7.916, 7.523, 5.359};
	AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{60.0, 240.0, 8.794e8, 0.7, 150.0, 225.0});

	//Тест 4 (настроен)
	r = {-1.232e8, 2.994e8, 7.446e7};
	v = {-8.851, -9.801, -21.764};
	AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{60.0, 300.0, 8.794e8, 0.7, 225.0, 300.0});

	//Тест 5 (настроен)
	r = {3.884e8, -1.121e8, -1.942e8};
	v = {-19.029, -4.438, -7.687};
	AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{120.0, 180.0, 8.794e8, 0.7, 300.0, 270.0});

	vector<double> kep;
	array<double,6> temp;
	for(double j = 0.0; j < 359.0;){
		for(double i = 0.0; i < 179.0;){
			kep = {i, j, 8.794e8, 0.7, 300.0, 270.0};
			temp = From_Kep_to_Dec(kep, mu_S);
			r[0] = temp[0];
			r[1] = temp[1];
			r[2] = temp[2];
			v[0] = temp[3];
			v[1] = temp[4];
			v[2] = temp[5];
			//From_Dec_to_Kep(r,v, mu_S);
			if(i == 0.0){
				AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{i, 0.0, 8.794e8, 0.7, 0.0, 270.0});
			}else{
				AssertEqual(From_Dec_to_Kep(r,v, mu_S), array<double,6>{i, j, 8.794e8, 0.7, 300.0, 270.0});
			}

			i += 1.0;
		}
		j += 1.0;
	}

}

void Test_Kep_to_Dec(){
	//double mu_E = 398600.4418;//гравитационная постоянная Земли (км^{3}c^{-2})
	double mu_S  = 132712440018.0;//гравитационная постоянная Солнца

	//Тест 1 (настроен вручную)
	vector<double> kep = {60.0, 150.0, 8.794e8, 0.7, 90.0, 45.0};
	array<double,6> r_v = {1.306e8, -1.979e8, 1.837e8 , 24.004, -6.837, -10.532 };
	AssertEqual(From_Kep_to_Dec(kep, mu_S), r_v);

	//Тест 2 (настроен)
	kep = move(vector<double>{120.0, 60.0, 8.794e8, 0.7, 60.0, 135.0});
	r_v = move(array<double,6>{-5.284e8, -6.854e8, -1.991e8, -7.576, -2.527, -9.176});
	AssertEqual(From_Kep_to_Dec(kep, mu_S), r_v);

	//Тест 3 (настроен)
	kep = move(vector<double>{60.0, 240.0, 8.794e8, 0.7, 150.0, 225.0});
	r_v = move(array<double,6>{-3.294e8, -8.003e8, 1.991e8, 7.916, 7.523, 5.359});
	AssertEqual(From_Kep_to_Dec(kep, mu_S), r_v);

	//Тест 4 (настроен)
	kep = move(vector<double>{60.0, 300.0, 8.794e8, 0.7, 225.0, 300.0});
	r_v = move(array<double,6>{-1.232e8, 2.994e8, 7.446e7, -8.851, -9.801, -21.764});
	AssertEqual(From_Kep_to_Dec(kep, mu_S), r_v);

	//Тест 5 (настроен)
	kep = move(vector<double>{120.0, 180.0, 8.794e8, 0.7, 300.0, 270.0});
	r_v = move(array<double,6>{3.884e8, -1.121e8, -1.942e8, -19.029, -4.438, -7.687});
	AssertEqual(From_Kep_to_Dec(kep, mu_S), r_v);

	for(double j = 0.0; j < 359.0;){
		for(double i = 0.0; i < 179.0;){
			kep = move(vector<double>{i, j, 8.794e8, 0.7, 300.0, 270.0});
			From_Kep_to_Dec(kep, mu_S);
			i += 1.0;
		}
		j += 1.0;
	}


}




