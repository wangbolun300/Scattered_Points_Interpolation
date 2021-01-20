#include"basis.h"
#include<iostream>
double Ni0(const int i, const double u, const std::vector<double> &U) {
	if (u >= U[i] && u < U[i + 1]) return 1.0;
	return 0;
}
double handle_division(const double a, const double b) {
	if (b == 0) return 1;
	else return a / b;
}
double Nip(const int i, const int p, const double u, const std::vector<double> &U) {
	if (p == 0) {
		return Ni0(i, u, U);
	}
	
	double result1 = handle_division(u - U[i], U[i + p] - U[i])*Nip(i, p - 1, u, U);
	double result2 = handle_division(U[i + p + 1] - u, U[i + p + 1] - U[i + 1])*Nip(i + 1, p - 1, u, U);
	return result1 + result2;
}


