#pragma once
#include <iostream>
using namespace std;
const int N = 200, K = 10;
const int Nm = N - 1;
double L = 5, T = 5, a = 1,
h = L / N, teta = T / K;
double U_k[N + 1], U_kp1[N + 1];   //матрицы для явной схемы
double A[Nm][Nm + 1];
double sigma = a * a * teta / h / h;

//alpha - коэффициент перед производной в левом условии
//betha - коэффициент перед функцией в левом условии
double alpha = 1, betha = 1;
//gamma - коэффициент перед производной в правом условии
//delta - коэффициент перед функцией в правом условии
double gamma = 1, delta = 1;

double U(int j) {
	return teta * a * a * (U_k[j + 1] - 2 * U_k[j] + U_k[j - 1]) / h / h + U_k[j];
}

//начальное условие
double begU(double x) {
	return 5;
}
//правое граничное условие
double RightU(double t) {
	return 1;
}
//левое граничное условие
double LeftU(double t) {
	//return alpha*(U_kp1[1]-U_kp1[0]
	return 11;
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
//явная схема
void yavnaya() {
	cout << "sigma: " << sigma << endl;
	for (int i = 0; i <= N; i++) {
		U_k[i] = begU(i * h);
	}
	U_k[0] = LeftU(0);
	U_k[N] = RightU(0);
	for (int k = 1; k <= K; k++) {
		for (int j = 1; j <= N - 1; j++) {
			U_kp1[j] = U(j);
		}
		//U_kp1[0]=
		for (int j = 1; j <= N - 1; j++) {
			U_k[j] = U_kp1[j];
		}
	}
	cout << "yavnaya: " << endl;
	ofstream f1("parabolich_yavn.txt", ios_base::trunc);
	for (int i = 0; i <= N; i++) {
		f1 << i * h << " " << U_k[i] << endl;
		cout << i * h << " " << U_k[i] << endl;
	}
	f1.close();
}

void neyavnaya() {
	//неявная схема
	//при к=0
	for (int i = 0; i <= N; i++) {
		U_k[i] = begU(i * h);
		if (i == 0) {
			U_k[i] = LeftU(0);
		}
		if (i == N) {
			U_k[i] = RightU(0);
		}
	}

	double b = -(1 + 2 * sigma);
	for (int k = 1; k <= K; k++) {
		for (int i = 1; i <= N - 1; i++) {
			double a, c, d;
			if (i == 1) {
				d = -(U_k[i] + sigma * LeftU(k * teta)); c = sigma;
				A[0][0] = b;
				A[0][1] = c;
				A[0][Nm] = d;
			}
			else {
				if (i == Nm) {
					a = sigma;
					d = -(U_k[i] + sigma * RightU(k * teta));
					A[i - 1][i - 2] = a;
					A[i - 1][i - 1] = b;
					A[i - 1][Nm] = d;
				}
				else {
					a = sigma;
					c = sigma;
					d = -U_k[i];
					A[i - 1][i - 2] = a;
					A[i - 1][i - 1] = b;
					A[i - 1][i] = c;
					A[i - 1][Nm] = d;
				}
			}
		}
		Gauss(0, A);
		for (int i = 1; i <= N - 1; i++) {
			U_k[i] = A[i - 1][Nm];
		}
	}
	cout << "Neyavnaya: " << endl;
	for (int i = 0; i <= N; i++) {
		//f1 << i * h << " " << U_k[i] << endl;
		cout << i * h << " " << U_k[i] << endl;
	}
}