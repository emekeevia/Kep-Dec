#include "transformation.h"

array<double,6> From_Dec_to_Kep(const vector<double>& r,const vector<double>& v, double mu){
	array<double,6> kepl;//(i, Omega, a, e, om, nu) =
						   //(наклонение, долгота восходящего угла, большая полуось,
	                       //эксценриситет, аргумент перецентра, истинная аномалия)
	Vector r_vect_mult_v = move(vect_mult(r,v));
	double abs_r_vect_mult_v = abs_vect(r_vect_mult_v);
	Vector k = move(r_vect_mult_v/abs_r_vect_mult_v);
	kepl[0] = acos(k[2]);//наклонение

	if(kepl[0] == 0.0){
		kepl[1] = 0.0;// долгота восходящего узла не определена
	}else{
		kepl[1] = acos(-k[1]/sin(kepl[0]));// долгота восходящего узла
		if(k[0]/sin(kepl[0]) < 0){
			kepl[1] = 2 * M_PI - kepl[1];
		}
	}
	double abs_v = abs_vect(v);
	double abs_r = abs_vect(r);
	double h = abs_v*abs_v - 2 * mu/abs_r; //константа энергии

	Vector c = move(vect_mult(r,v));//вектор константы скорости
	double p = scal_mult(c,c)/mu;//фокальный параметр
	if(h >= 0){
		throw(runtime_error("The orbit is not elliptical"));
	}
	kepl[3] = sqrt(1.0 + h*scal_mult(c,c)/(mu*mu));//эксцентриситет
	kepl[2] = p / (1.0 - kepl[3] * kepl[3]);

	//Найдём истинную аномалию
	kepl[5] = acos((1.0/kepl[3]) * (p/abs_r - 1.0));

	if(scal_mult(r,v) < 0){
		kepl[5] = 2*M_PI - kepl[5];
	}else if(scal_mult(r,v) == 0){
		throw(runtime_error("Undefined anomaly"));
	}

	if(k[2] != 1){
		//Найдем аргумент перецентра
		Vector normir_r = move(r/abs_r);

		Vector I = {1.0, 0.0, 0.0};
		Vector J = {0.0, 1.0, 0.0};

		Vector normir_b = I * cos(kepl[1]) + J * sin(kepl[1]);
		double cos_sum = scal_mult(normir_r, normir_b);
		double sin_sum = scal_mult(k, vect_mult(normir_b, normir_r));
		double omega;

		if(sin_sum >= 0){
			omega = acos(cos_sum) - kepl[5];
		}else{
			omega = 2*M_PI - acos(cos_sum) - kepl[5];
		}

		if(omega < 0){
			omega += 2*M_PI;
		}else if(2*M_PI < omega){
			omega -= 2*M_PI;
		}
		kepl[4] = omega;
	}else{
		kepl[4] = 0.0;
	}

	kepl[0] = kepl[0] * 180/M_PI;
	kepl[1] = kepl[1] * 180/M_PI;
	kepl[4] = kepl[4] * 180/M_PI;
	kepl[5] = kepl[5] * 180/M_PI;

	kepl[0] = round(kepl[0], 360.0);
	kepl[1] = round(kepl[1], 360.0);
	kepl[2] = round(kepl[2], 360.0);
	kepl[3] = round(kepl[3], 360.0);

	return kepl;
}

array<double,6> From_Kep_to_Dec(vector<double>& kep, double mu){
	array<double,6> dec;//первый вектор отвечает за координаты, а второй за скорость
	Matrix B;//матрица перехода от орбитальных координат к декартовым
	Vector b_1, b_2;//первые 2 столбца матрицы B

	if(kep[3] >= 1){
		throw(runtime_error("The orbit is not elliptical"));
	}

	grad_to_rad(kep);

	double E = 2*atan(sqrt((1-kep[3])/(1+kep[3]))*tan(kep[5]/2));//эксцентричная аномалия



	B[0][0] = cos(kep[1]) * cos(kep[4]) - sin(kep[1]) * sin(kep[4]) * cos(kep[0]);
	B[0][1] = -cos(kep[1]) * sin(kep[4]) - sin(kep[1]) * cos(kep[4]) * cos(kep[0]);
	B[0][2] = sin(kep[1]) * sin(kep[0]);

	B[1][0] = sin(kep[1]) * cos(kep[4]) + cos(kep[1]) * sin(kep[4]) * cos(kep[0]);
	B[1][1] = -sin(kep[1]) * sin(kep[4]) + cos(kep[1]) * cos(kep[4]) * cos(kep[0]);
	B[1][2] = -cos(kep[1]) * sin(kep[0]);


	B[2][0] = sin(kep[4]) * sin(kep[0]);
	B[2][1] = cos(kep[4]) * sin(kep[0]);
	B[2][2] = cos(kep[0]);


	b_1 = {B[0][0], B[1][0], B[2][0]};
	b_2 = {B[0][1], B[1][1], B[2][1]};


	Vector r = (b_1*(cos(E) - kep[3]) + b_2*sqrt(1.0 - kep[3] * kep[3])*sin(E))*kep[2];
	for(size_t i = 0; i < 3; i++){
		dec[i] = r[i];
	}


	Vector tau = {-sin(E)/sqrt(1.0 - kep[3]*kep[3]), cos(E), 0.0};
	tau = tau/abs_vect(tau);//направление скорости в системе координат, связанной с плоскостью орбиты


	tau = matr_vect_mult(B, tau);

	double abs_V = sqrt(mu * (2.0/abs_vect(r) - 1.0/kep[2]));
	Vector V = tau * abs_V;

	for(size_t i = 3; i < 6; i++){
		dec[i] = V[i-3];
	}
	kep[0] = kep[0]*180/M_PI;
	kep[1] = kep[1]*180/M_PI;
	kep[4] = kep[4]*180/M_PI;
	kep[5] = kep[5]*180/M_PI;
	return dec;
}




