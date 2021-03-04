#ifdef SPARSE_INTERP_WITH_GMP
#include<exact_calculation.h>
#include <Eigen/Dense>
#include<iostream>

// the B-spline basis function N_{i,0}(u). U is the knot vector
double Ni0_(const int i, const double u, const std::vector<double> &U) {
	if (u >= U[i] && u < U[i + 1]) return 1.0;

	return 0.0;
}

Rational handle_division(const Rational a, const Rational b) {

	if (b == 0) {
		// if the denominator is 0, then this term is 0
		return 0;
	}
	else return a / b;
}
// the B-spline basis function N_{i,p}(u)
Rational Nip_(const int i, const int p, const double u, const std::vector<double> &U) {
	if (p == 0) {
		return Ni0_(i, u, U);
	}
	if (u == U.back()) {
		if (i == U.size() - 2 - p) {
			return 1;
		}
		else {
			return 0;
		}
	}
	Rational r1, r2;
	r1 = Rational(u) - Rational(U[i]);
	r2 = Rational(U[i + p]) - Rational(U[i]);
	Rational result1 = handle_division(r1, r2)*Rational(Nip_(i, p - 1, u, U));

	r1 = Rational(U[i + p + 1]) - Rational(u);
	r2 = Rational(U[i + p + 1]) - Rational(U[i + 1]);
	Rational result2 = handle_division(r1, r2)*Rational(Nip_(i + 1, p - 1, u, U));
	return result1 + result2;
}
Rational Nip_Rational(const int i, const int p, const double u, const std::vector<double> &U) {
	return Nip_(i, p, u, U);
}

void print_rational_matrix(const MatrixXs& matrix) {
	for (int i = 0; i < matrix.rows(); i++) {
		for (int j = 0; j < matrix.cols(); j++) {
			Rational r = matrix(i, j);
			double v = r.to_double();
			std::cout << v << ", ";
		}
		std::cout << std::endl;
	}
}
int rank(const MatrixXs& matrix) {
	const auto m= matrix.fullPivLu();
	/*int r = m.rank();
	return r;*/
	print_rational_matrix(m);
	return 1;
}

//bool equation_has_solution(const Eigen::MatrixXd& A,
//	const Eigen::VectorXd& b) {
//	Eigen::MatrixXd  Ab;
//	int rankA = rank(A);
//	Ab.resize(A.rows(), A.cols() + 1);
//	Ab << A, b;
//	int rankAb = rank(Ab);
//	//std::cout << "Ab,\n" << Ab << std::endl;
//	if (rankA == rankAb) {
//		return true;
//	}
//	return false;
//}

#endif