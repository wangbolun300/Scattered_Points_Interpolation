#ifdef SPARSE_INTERP_WITH_GMP
#include<exact_calculation.h>
#include <Eigen/Dense>
#include<iostream>
#include <Types.hpp>

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
	//const auto m= matrix.fullPivLu();
	Eigen::FullPivLU<MatrixXs> lu(matrix);
	MatrixXs u = lu.matrixLU().triangularView<Eigen::Upper>();

	int rk = 0;
	int rmin = std::min(u.rows(), u.cols());
	
	for (int i = 0; i < rmin; i++) {
		
		if (u(i, i) == 0) {
			break;
		}
		else {
			rk += 1;
		}
	}
	//std::cout << "u is \n";
	//print_rational_matrix(u);
	std::cout << "rank is " << rk << std::endl;
	return rk;
}

MatrixXs build_matrix_A_Rational(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const std::vector<int> row_id) {
	assert(paras.cols() == 2);

	int nu = U.size() - 2 - degree1;// n + 1 = number of u control points
	int nv = V.size() - 2 - degree2;// n + 1 = number of v control points
	int m = paras.rows() - 1;// m + 1 = the number of data points
	MatrixXs result;
	result.resize(row_id.size(), (nu + 1)*(nv + 1));

	for (int l = 0; l < row_id.size(); l++) {
		int i = row_id[l];

		double u = paras(i, 0);
		double v = paras(i, 1);

		for (int j = 0; j < result.cols(); j++) {
			// get the indices of N_r(u) and N_q(v)
			int r = j / (nv + 1);
			int q = j - r * (nv + 1);
			Rational N1 = Nip_Rational(r, degree1, u, U);
			Rational N2 = Nip_Rational(q, degree2, v, V);


			result(l, j) = N1 * N2;
		}
	}

	return result;
}

// build vector b for interpolation problem (Ax=b)

Eigen::VectorXd build_Vector_b_rational(const Eigen::MatrixXd& points, const int dimension,
	const std::vector<int> row_id) {
	Eigen::VectorXd result;
	result.resize(row_id.size());
	for (int l = 0; l < row_id.size(); l++) {
		int i = row_id[l];
		result(l) = points(i, dimension);
	}
	return result;
}
bool equation_has_solution_rational(const MatrixXs& A,
	const Eigen::VectorXd& b) {
	MatrixXs  Ab;
	int rankA = rank(A);
	std::cout << "rank A = " << rankA << std::endl;
	Ab.resize(A.rows(), A.cols() + 1);
	MatrixXs br(b.size(), 1);
	for (int i = 0; i < br.rows(); i++) {
		br(i, 0) = b[i];
	}


	Ab << A, br;
	int rankAb = rank(Ab);
	//std::cout << "Ab,\n" << Ab << std::endl;
	if (rankA == rankAb) {
		return true;
	}
	return false;
}
bool selected_rows_have_solution_rational(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	const std::vector<int> &row_id, const int dimension) {

	MatrixXs A = build_matrix_A_Rational(degree1, degree2, U, V, paras, row_id);

	Eigen::VectorXd b = build_Vector_b_rational(points, dimension, row_id);

	return equation_has_solution_rational(A, b);
}



#endif