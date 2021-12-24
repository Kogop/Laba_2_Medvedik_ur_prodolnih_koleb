#include <iostream>
#include <fstream>
//#include "eliptich.h"
//#include "parabol.h"
//#include "hiperbol.h"

using namespace std;

//
//
//using namespace std;
//const int N = 6, K = 10;
//const int Nm = N - 1;
//double L = 5, T = 5, a = 0.1,
//h = L / N, teta = T / K;
//double U_k[N + 1], U_kp1[N + 1];   //������� ��� ����� �����
//double A[Nm][Nm + 1];
//double sigma = a * a * teta / h / h;
//
////alpha - ����������� ����� ����������� � ����� �������
////betha - ����������� ����� �������� � ����� �������
//double alpha = 1, betha = 1;
////gamma - ����������� ����� ����������� � ������ �������
////delta - ����������� ����� �������� � ������ �������
//double gamma = 1, delta = 1;
//
//double U(int j) {
//	return teta * a * a * (U_k[j + 1] - 2 * U_k[j] + U_k[j - 1]) / h / h + U_k[j];
//}
//
////��������� �������
//double begU(double x) {
//	return 5;
//}
////������ ��������� �������
//double RightU(double t) {
//	return 1;
//}
////����� ��������� �������
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
//		for (int j = k; j < Nm + 1; j++) {//������������ ������
//			Matrix[k][j] = Matrix[k][j] / T;
//		}
//	}
//	for (int i = 0; i < Nm; i++) { //�������� �� �������
//		if ((Matrix[i][k] != 0) & (i != k)) {
//			double T = Matrix[i][k];
//			Matrix[i][k] = 0;
//			for (int j = k + 1; j < Nm + 1; j++) { //�������� �� ���� ������� � �������� ��
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

////////////////////////////////////////////////////////////////////////////

using namespace std;
const int N = 40, K = 100;
const int Nm = N - 1;
double L = 5, T = 5,
h = L / N, teta = T / K;
double a = 1, b = 1, c = 1;
double U_k[N + 1], U_kp1[N + 1], U_km1[N + 1];   //������� ��� ����� �����
double A[Nm][Nm + 1];
double sigma = a * a * teta * teta / h / h;

/*
U_tt=a^2*U_xx+b*U_x+c*U+f(x,t)
*/

//alpha - ����������� ����� ����������� � ����� �������
//betha - ����������� ����� �������� � ����� �������
double alpha = 1, betha = 1;
//gamma - ����������� ����� ����������� � ������ �������
//delta - ����������� ����� �������� � ������ �������
double gamma = 1, delta = 1;

double f(double x, double t) {
	//return 0;
	return 2 * t;
}

double U(int j, int k) {
	return U_k[j + 1] * (sigma + b * teta * teta / 2 / h) + U_k[j] * (-2 * sigma + 2 + c * teta * teta) +
		U_k[j - 1] * (sigma - b * teta * teta / 2 / h) - U_km1[j] + teta * teta * f(j * h, teta * k);
}

//������ ��������� �������
double begU(double x) {
	return 0;
}
//������ ��������� �������
double dbegU(double x) {
	return 0.01 * cos(x);
}
//������ ��������� �������
double RightU(double t) {
	return 0;
}
//����� ��������� �������
double LeftU(double t) {
	return 0;
}

//void Gauss(int k, double Matrix[Nm][Nm + 1]) {
//	if (Matrix[k][k] != 1) {
//		double T = Matrix[k][k];
//		for (int j = k; j < Nm + 1; j++) {//������������ ������
//			Matrix[k][j] = Matrix[k][j] / T;
//		}
//	}
//	for (int i = 0; i < Nm; i++) { //�������� �� �������
//		if ((Matrix[i][k] != 0) & (i != k)) {
//			double T = Matrix[i][k];
//			Matrix[i][k] = 0;
//			for (int j = k + 1; j < Nm + 1; j++) { //�������� �� ���� ������� � �������� ��
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

//����� �����
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


///////////////////////////////////////


int main() {
	
	yavnaya();



	return 1;
}
