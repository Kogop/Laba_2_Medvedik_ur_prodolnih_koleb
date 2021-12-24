#pragma once
#include <iostream>
#include <fstream>
using namespace std;
const int N = 50, Nm = (N - 2) * (N - 2);
double L = 1, h = L / (N - 1.0);//квадрат L x L

double A1[N][N], SLAE[Nm][Nm + 1];

//alpha - коэффициент перед производной в левом условии
//betha - коэффициент перед функцией в левом условии
double alpha = 1, betha = 1;
//gamma - коэффициент перед производной в правом условии
//delta - коэффициент перед функцией в правом условии
double gamma = 1, delta = 1;
/*		Equation:
*	 U_xx+U_yy=f(x,y)
*/
double f(double x, double y) {
	if ((0.48 < x < 0.52)&&(0.48 < y < 0.52)) { return -1000; }
	else {
		return 50;
	}
	//return -x;
	//return 0;
	//return 5+cos(x) + y * y;
}
double U(int j) {
	return 0;
}
//Граничные условия
//Нижний участок y=0
double phi1(double x) {
	return 5;
}
//Верхний участок y=L
double phi2(double x) {
	return 2;
}
//Левый участок x=0
double phi3(double y) {
	return 6;
}
//Правый участок x=L
double phi4(double y) {
	return 1;
}

void Gauss(int k, double Matrix[Nm][Nm + 1]) {
	if (Matrix[k][k] != 1) {
		double T = Matrix[k][k];
		for (int j = k; j < Nm + 1; j++) {//нормирование строки
			Matrix[k][j] = Matrix[k][j] / T;
		}
	}
	for (int i = 0; i < Nm; i++) { //проходим по столбцу
		if ((Matrix[i][k] != 0) & (i != k)) {
			double T = Matrix[i][k];
			Matrix[i][k] = 0;
			for (int j = k + 1; j < Nm + 1; j++) { //проходим по двум строкам и вычитаем их
				Matrix[i][j] -= Matrix[k][j] * T;
			}
		}
	}
	if (k < Nm - 1) {
		Gauss(k + 1, Matrix);
	}
}

void vyvod(double Matr[Nm][Nm + 1]) {
	for (int i = 0; i < Nm; i++) {
		for (int j = 0; j < Nm + 1; j++) {
			cout << Matr[i][j] << " ";
		}
		cout << endl;
	}
}

void vyvod(double Matr[N][N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%4f ", Matr[i][j]);
		}
		cout << endl;
	}
}

void Method() {
	//Заносим граничные условия 
	for (int i = 0; i < N; i++){
		A1[i][0] = phi3(i * h);
		A1[0][i] = phi2(i * h);
		A1[i][N - 1] = phi4(i * h);
		A1[N - 1][i] = phi1(i * h);
	}
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			double d = f(h * j, h * i);
			int _i = i - 1, _j = j - 1;
			int Str = _i * (N - 2) + _j;

			SLAE[Str][Str] = -4;

			if (i - 1 == 0) { d -= A1[i - 1][j]; }
			else { SLAE[Str][Str - (N-2)] = 1; }

			if (i + 1 == N - 1) { d -= A1[i + 1][j]; }
			else { SLAE[Str][Str + (N - 2)] = 1; }

			if (j - 1 == 0) { d -= A1[i][j - 1]; }
			else { SLAE[Str][Str - 1] = 1; }

			if (j + 1 == N - 1) { d -= A1[i][j + 1]; }
			else { SLAE[Str][Str + 1] = 1; }

			SLAE[Str][Nm] = d;
			//cout << _i * Nm + _j << endl;
		}
	}
	Gauss(0, SLAE);
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			int Str = (i-1) * (N - 2) + (j-1);
			A1[i][j] = SLAE[Str][Nm];
		}
	}
}
void File() {
	ofstream f1("eliptich.txt", ios_base::trunc);
	f1 << "X Y Z Mod\n";
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f1 << i * h << " " << j * h << " 0 " << A1[i][j] << endl;
		}
	}
	f1.close();
}