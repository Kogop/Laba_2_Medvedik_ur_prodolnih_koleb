#include <iostream>
#include <fstream>
//#include "eliptich.h"
//#include "parabol.h"
//#include "hiperbol.h"

using namespace std;

//  parabol
//
//using namespace std;
//const int N = 6, K = 10;
//const int Nm = N - 1;
//double L = 5, T = 5, a = 0.1,
//h = L / N, teta = T / K;
//double U_k[N + 1], U_kp1[N + 1];   //матрицы для явной схемы
//double A[Nm][Nm + 1];
//double sigma = a * a * teta / h / h;
//
////alpha - коэффициент перед производной в левом условии
////betha - коэффициент перед функцией в левом условии
//double alpha = 1, betha = 1;
////gamma - коэффициент перед производной в правом условии
////delta - коэффициент перед функцией в правом условии
//double gamma = 1, delta = 1;
//
//double U(int j) {
//	return teta * a * a * (U_k[j + 1] - 2 * U_k[j] + U_k[j - 1]) / h / h + U_k[j];
//}
//
////начальное условие
//double begU(double x) {
//	return 5;
//}
////правое граничное условие
//double RightU(double t) {
//	return 1;
//}
////левое граничное условие
//double LeftU(double t) {
//	//return alpha*(U_kp1[1]-U_kp1[0]
//	return 11;
//}
//double AnalitRew(double x,double t) {
//
//	return exp(-a * t) * sin(x);
// }
//void Gauss(int k, double Matrix[Nm][Nm + 1]) {
//	if (Matrix[k][k] != 1) {
//		double T = Matrix[k][k];
//		for (int j = k; j < Nm + 1; j++) {//нормирование строки
//			Matrix[k][j] = Matrix[k][j] / T;
//		}
//	}
//	for (int i = 0; i < Nm; i++) { //проходим по столбцу
//		if ((Matrix[i][k] != 0) & (i != k)) {
//			double T = Matrix[i][k];
//			Matrix[i][k] = 0;
//			for (int j = k + 1; j < Nm + 1; j++) { //проходим по двум строкам и вычитаем их
//				Matrix[i][j] -= Matrix[k][j] * T;
//			}
//		}
//	}
//	if (k < Nm - 1) {
//		Gauss(k + 1, Matrix);
//	}
//}
//

//
//void yavnaya() {
//	cout << "sigma: " << sigma << endl;
//	for (int i = 0; i <= N; i++) {
//		U_k[i] = begU(i * h);
//	}
//	U_k[0] = LeftU(0);
//	U_k[N] = RightU(0);
//	for (int k = 1; k <= K; k++) {
//		for (int j = 1; j <= N - 1; j++) {
//			U_k[j] = U(j);
//		}
//	}
//	cout << "yavnaya: " << endl;
//	ofstream f1("parabolich_yavn.txt", ios_base::trunc);
//	for (int i = 0; i <= N; i++) {
//		f1 << i * h << " " << U_k[i] << endl;
//		cout << i * h << " " << U_k[i] << endl;
//	}
//	f1.close();
//}
//
//
//

//void neyavnaya() {
//	//неявная схема
//	//при к=0
//	for (int i = 0; i <= N; i++) {
//		U_k[i] = begU(i * h);
//		if (i == 0) {
//			U_k[i] = LeftU(0);
//		}
//		if (i == N) {
//			U_k[i] = RightU(0);
//		}
//	}
//
//	double b = -(1 + 2 * sigma);
//	for (int k = 1; k <= K; k++) {
//		for (int i = 1; i <= N - 1; i++) {
//			double a, c, d;
//			if (i == 1) {
//				d = -(U_k[i] + sigma * LeftU(k * teta)); c = sigma;
//				A[0][0] = b;
//				A[0][1] = c;
//				A[0][Nm] = d;
//			}
//			else {
//				if (i == Nm) {
//					a = sigma;
//					d = -(U_k[i] + sigma * RightU(k * teta));
//					A[i - 1][i - 2] = a;
//					A[i - 1][i - 1] = b;
//					A[i - 1][Nm] = d;
//				}
//				else {
//					a = sigma;
//					c = sigma;
//					d = -U_k[i];
//					A[i - 1][i - 2] = a;
//					A[i - 1][i - 1] = b;
//					A[i - 1][i] = c;
//					A[i - 1][Nm] = d;
//				}
//			}
//		}
//		Gauss(0, A);
//		for (int i = 1; i <= N - 1; i++) {
//			U_k[i] = A[i - 1][Nm];
//		}
//	}
//	cout << "Neyavnaya: " << endl;
//	for (int i = 0; i <= N; i++) {
//		//f1 << i * h << " " << U_k[i] << endl;
//		cout << i * h << " " << U_k[i] << endl;
//	}
//}



////////////////////////////////////////////////////////////////////////////

// hiperbol 

using namespace std;
const int N = 40, K = 100;
const int Nm = N - 1;
double L = 5, T = 5,
h = L / N, teta = T / K;
double a = 1, b = 1, c = 1;
double U_k[N + 1], U_kp1[N + 1], U_km1[N + 1];   //матрицы для явной схемы
double A[Nm][Nm + 1];
double sigma = a * a * teta * teta / h / h;

/*
U_tt=a^2*U_xx+b*U_x+c*U+f(x,t)
*/

//alpha - коэффициент перед производной в левом условии
//betha - коэффициент перед функцией в левом условии
double alpha = 1, betha = 1;
//gamma - коэффициент перед производной в правом условии
//delta - коэффициент перед функцией в правом условии
double gamma = 1, delta = 1;

double f(double x, double t) {
	//return 0;
	return 2 * t;
}

double U(int j, int k) {
	return U_k[j + 1] * (sigma + b * teta * teta / 2 / h) + U_k[j] * (-2 * sigma + 2 + c * teta * teta) +
		U_k[j - 1] * (sigma - b * teta * teta / 2 / h) - U_km1[j] + teta * teta * f(j * h, teta * k);
}

//первое начальное условие
double begU(double x) {
	return 0;
}
//второе начальное условие
double dbegU(double x) {
	return 0.01 * cos(x);
}
//правое граничное условие
double RightU(double t) {
	return 0;
}
//левое граничное условие
double LeftU(double t) {
	return 0;
}

//void Gauss(int k, double Matrix[Nm][Nm + 1]) {
//	if (Matrix[k][k] != 1) {
//		double T = Matrix[k][k];
//		for (int j = k; j < Nm + 1; j++) {//нормирование строки
//			Matrix[k][j] = Matrix[k][j] / T;
//		}
//	}
//	for (int i = 0; i < Nm; i++) { //проходим по столбцу
//		if ((Matrix[i][k] != 0) & (i != k)) {
//			double T = Matrix[i][k];
//			Matrix[i][k] = 0;
//			for (int j = k + 1; j < Nm + 1; j++) { //проходим по двум строкам и вычитаем их
//				Matrix[i][j] -= Matrix[k][j] * T;
//			}
//		}
//	}
//	if (k < Nm - 1) {
//		Gauss(k + 1, Matrix);
//	}
//}
//
//void vyvod(double Matr[Nm][Nm + 1]) {
//	for (int i = 0; i < Nm; i++) {
//		for (int j = 0; j < Nm + 1; j++) {
//			cout << Matr[i][j] << " ";
//		}
//		cout << endl;
//	}
//}

//явная схема
void yavnaya() {
	cout << "sigma: " << sigma << endl;
	//k=0
	for (int i = 0; i <= N; i++) {
		U_km1[i] = begU(i * h);
	}
	U_k[0] = LeftU(0);
	U_k[N] = RightU(0);
	//k=1
	for (int i = 0; i <= N; i++) {
		U_k[i] = begU(i * h) + dbegU(i * h) * teta;
	}
	U_k[0] = LeftU(teta);
	U_k[N] = RightU(teta);
	for (int k = 2; k <= K; k++) {
		for (int j = 1; j <= N - 1; j++) {
			U_kp1[j] = U(j, k);
		}
		for (int j = 1; j <= N - 1; j++) {
			U_km1[j] = U_k[j];
			U_k[j] = U_kp1[j];
		}
	}
	cout << "yavnaya: " << endl;
	ofstream f1("hiperbol_yavn.txt", ios_base::trunc);
	for (int i = 0; i <= N; i++) {
		f1 << i * h << " " << U_k[i] << endl;
		cout << i * h << " " << U_k[i] << endl;
	}
	f1.close();
}
//
//void neyavnaya() {
//	//неявная схема
//	//k=0
//	for (int i = 0; i <= N; i++) {
//		U_km1[i] = begU(i * h);
//	}
//	U_k[0] = LeftU(0);
//	U_k[N] = RightU(0);
//	//k=1
//	for (int i = 0; i <= N; i++) {
//		U_k[i] = begU(i * h) + dbegU(i * h) * teta;
//	}
//	U_k[0] = LeftU(teta);
//	U_k[N] = RightU(teta);
//
//	double b1 = -(1 + 2 * sigma);
//	for (int k = 2; k <= K; k++) {
//		for (int i = 1; i <= N - 1; i++) {
//			double a1, c1, d1;
//			if (i == 1) {
//				d1 = U_km1[i] - 2 * U_k[i] - sigma * LeftU(teta * k); c1 = sigma;
//				A[0][0] = b1;
//				A[0][1] = c1;
//				A[0][Nm] = d1;
//			}
//			else {
//				if (i == Nm) {
//					a1 = sigma;
//					d1 = U_km1[i] - 2 * U_k[i] - sigma * RightU(teta * k);
//					A[i - 1][i - 2] = a1;
//					A[i - 1][i - 1] = b1;
//					A[i - 1][Nm] = d1;
//				}
//				else {
//					a1 = sigma;
//					c1 = sigma;
//					d1 = U_km1[i] - 2 * U_k[i];
//					A[i - 1][i - 2] = a1;
//					A[i - 1][i - 1] = b1;
//					A[i - 1][i] = c1;
//					A[i - 1][Nm] = d1;
//				}
//			}
//		}
//		//vyvod(A);
//		//system("pause");
//		Gauss(0, A);
//		for (int j = 1; j <= N - 1; j++) {
//			U_km1[j] = U_k[j];
//			U_k[j] = A[j - 1][Nm];
//		}
//	}
//	cout << "Neyavnaya: " << endl;
//	ofstream f1("giperbol_neyavn.txt", ios_base::trunc);
//	for (int i = 0; i <= N; i++) {
//		f1 << i * h << " " << U_k[i] << endl;
//		cout << i * h << " " << U_k[i] << endl;
//	}
//	f1.close();
//}
///////////////////////////////////////


int main() {
	
	yavnaya();

	//neyavnaya();

	return 1;
}
