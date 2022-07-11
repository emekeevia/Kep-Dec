#include "transformation.h"

struct Kepler_elements{
	double i ;
	double Omega;
	double a;
	double e ;
	double omega ;
	double nu ;
	Kepler_elements(double i_ = 0.0, double Omega_ = 0.0, double a_ = 0.0,
					double e_ = 0.0, double omega_ = 0.0, double nu_ = 0.0):i(i_),
					Omega(Omega_), a(a_), e(e_), omega(omega_), nu(nu_){}
};

struct Dec_kordinates{
	Vector r;
	Vector v;
	Dec_kordinates(Vector r_ = {0.0, 0.0, 0.0},Vector v_ = {0.0, 0.0, 0.0}):r(r_),  v(v_){}
};

bool eq(double d1, double d2){
	if(d1 * d2 != 0.0 && abs(d1 - d2/std::max(std::abs(d1), std::abs(d2))) < 1.0e-5){
		return true;
	}else if(d1 * d2 == 0.0 && std::max(std::abs(d1), std::abs(d2)) < 1.0e-5){
		return true;
	}
	return false;
}
bool operator==(Kepler_elements& k1, Kepler_elements& k2){
	return eq(k1.i, k2.i) && eq(k1.Omega, k2.Omega) && eq(k1.a, k2.a) &&
			eq(k1.e, k2.e) && eq(k1.omega, k2.omega) && eq(k1.nu, k2.nu);
}

Kepler_elements From_Dec_to_Kep(const Dec_kordinates& dec, double mu){
	Kepler_elements kepl;//(i, Omega, a, e, om, nu) =
						   //(наклонение, долгота восходящего угла, большая полуось,
	                       //эксценриситет, аргумент перецентра, истинная аномалия)
	Vector c = std::move(vect_mult(dec.r, dec.v));//вектор константы скорости
	double abs_r_vect_mult_v = abs_vect(c);
	Vector k = std::move(c/abs_r_vect_mult_v);
	kepl.i = acos(k[2]);//наклонение

	if(kepl.i > 0.5e-5){
		kepl.Omega = std::acos(-k[1]/sin(kepl.i));// долгота восходящего узла
		if(k[0]/std::sin(kepl.i) < 0){
			kepl.Omega = 2 * M_PI - kepl.Omega;
		}
	}

	double abs_v = abs_vect(dec.v);
	double abs_r = abs_vect(dec.r);
	double h = abs_v*abs_v - 2 * mu/abs_r; //константа энергии

	double c_in_sq = scal_mult(c,c);
	double p = c_in_sq/mu;//фокальный параметр
	if(h >= 0){
		throw(runtime_error("The orbit is not elliptical"));
	}
	kepl.e = std::sqrt(1.0 + h*c_in_sq/(mu*mu));//эксцентриситет
	kepl.a = p / (1.0 - kepl[3] * kepl[3]);

	//Найдём истинную аномалию
	kepl.nu = std::acos((1.0/kepl.e) * (p/abs_r - 1.0));

	if(scal_mult(dec.r,dec.v) < 0.0){
		kepl.nu = 2*M_PI - kepl.nu;
	}else if(scal_mult(dec.r,dec.v) > 0.0 && scal_mult(dec.r,dec.v) < 0.5e-5 ){
		throw(runtime_error("Undefined anomaly"));
	}

	if(k[2] != 1.0){
		//Найдем аргумент перецентра
		Vector normir_r = std::move(dec.r/abs_r);

		Vector I = {1.0, 0.0, 0.0};
		Vector J = {0.0, 1.0, 0.0};

		Vector normir_b = I * std::cos(kepl.Omega) + J * std::sin(kepl.Omega);
		double cos_sum = scal_mult(normir_r, normir_b);
		double sin_sum = scal_mult(k, vect_mult(normir_b, normir_r));

		if(sin_sum >= 0){
			kepl.omega = std::acos(cos_sum) - kepl.nu;
		}else{
			kepl.omega = 2*M_PI - std::acos(cos_sum) - kepl.nu;
		}

		if(kepl.omega < 0){
			kepl.omega += 2*M_PI;
		}else if(2*M_PI < kepl.omega){
			kepl.omega -= 2*M_PI;
		}

	}

//	kepl[0] = kepl[0] * 180/M_PI;
//	kepl[1] = kepl[1] * 180/M_PI;
//	kepl[4] = kepl[4] * 180/M_PI;
//	kepl[5] = kepl[5] * 180/M_PI;
//
//	kepl[0] = round(kepl[0], 360.0);
//	kepl[1] = round(kepl[1], 360.0);
//	kepl[2] = round(kepl[2], 360.0);
//	kepl[3] = round(kepl[3], 360.0);

	return kepl;
}

Dec_kordinates From_Kep_to_Dec(const Kepler_elements& kep, double mu){
	Dec_kordinates dec;//первый вектор отвечает за координаты, а второй за скорость
	Matrix B;//матрица перехода от орбитальных координат к декартовым
	Vector b_1, b_2;//первые 2 столбца матрицы B

	if(kep.e >= 1.0){
		throw(runtime_error("The orbit is not elliptical"));
	}

	//grad_to_rad(kep);

	double E = 2*std::atan(std::sqrt((1-kep.e)/(1+kep.e))*std::tan(kep.nu/2));//эксцентричная аномалия

	double cos_zero = std::cos(kep.i);
	double sin_zero = std::sin(kep.i);
	double cos_first = std::cos(kep.Omega);
	double sin_first = std::sin(kep.Omega);
	double cos_four = std::cos(kep.omega);
	double sin_four = std::sin(kep.omega);


	B[0][0] = cos_first * cos_four - sin_first * sin_four * cos_zero;
	B[0][1] = -cos_first * sin_four - sin_first * cos_four * cos_zero;
	B[0][2] = sin_first * sin_zero;

	B[1][0] = sin_first * cos_four + cos_first * sin_four * cos_zero;
	B[1][1] = -sin_first * sin_four + cos_first * cos_four * cos_zero;
	B[1][2] = -cos_first * sin_zero;


	B[2][0] = sin_four * sin_zero;
	B[2][1] = cos_four * sin_zero;
	B[2][2] = cos_zero;


	b_1 = {B[0][0], B[1][0], B[2][0]};
	b_2 = {B[0][1], B[1][1], B[2][1]};


	dec.r = (b_1*(std::cos(E) - kep.e) + b_2*std::sqrt(1.0 - kep.e * kep.e)*std::sin(E))*kep.a;



	Vector tau = {-std::sin(E)/std::sqrt(1.0 - kep.e*kep.e), std::cos(E), 0.0};
	tau = tau/abs_vect(tau);//направление скорости в системе координат, связанной с плоскостью орбиты


	tau = matr_vect_mult(B, tau);

	double abs_V = std::sqrt(mu * (2.0/abs_vect(dec.r) - 1.0/kep.a));
	dec.v = tau * abs_V;

//	kep[0] = kep[0]*180/M_PI;
//	kep[1] = kep[1]*180/M_PI;
//	kep[4] = kep[4]*180/M_PI;
//	kep[5] = kep[5]*180/M_PI;
	return dec;
}




