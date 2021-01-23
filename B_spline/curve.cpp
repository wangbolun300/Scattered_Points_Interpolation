#include"curve.h"
#include <Eigen/Dense>
#include<iostream>
// return a point position of a given curve
Vector3d BsplinePoint(const int degree, const std::vector<double>& U, const double para,
	const std::vector<Vector3d> &pts) {
	Eigen::Vector3d result = Eigen::Vector3d(0, 0, 0);
	for (int i = 0; i < pts.size(); i++) {
		double base = Nip(i, degree, para, U);
		//std::cout << "base " << base << std::endl;
		result += base * pts[i];
	}
	return result;
}
// return a point position of a given curve
Vector3d BsplinePoint(const int degree, const std::vector<double>& U, const double para,
	const Eigen::MatrixXd& pts) {
	Eigen::Vector3d result = Eigen::Vector3d(0, 0, 0);
	for (int i = 0; i < pts.rows(); i++) {
		double base = Nip(i, degree, para, U);
		//std::cout << "base " << base << std::endl;
		result += base * pts.row(i);
	}
	return result;
}

//TODO this error is not what we want
Eigen::MatrixXd slove_linear_system(const Eigen::MatrixXd& A, const Eigen::MatrixXd &b,
	const bool check_error, double &relative_error) {
	Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
	//Eigen::VectorXd x = A.fullPivLu().solve(b);

	if (check_error) {
		relative_error = (A*x - b).norm() / b.norm();
	}
	return x;
}
Eigen::VectorXd slove_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd &b,
	const bool check_error, double &relative_error) {
	Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
	//Eigen::VectorXd x = A.fullPivLu().solve(b);

	if (check_error) {
		relative_error = (A*x - b).norm() / b.norm();
	}
	return x;
}
Eigen::MatrixXd build_matrix_N(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras) {
	int n = U.size() - 2 - degree; // there are n+1 control points;
	int m = paras.size() - 1;  // there are m+1 points to fit
	Eigen::MatrixXd result(m - 1, n - 1);
	for (int i = 0; i < m - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			result(i, j) = Nip(j + 1, degree, paras[i + 1], U);
		}
	}
	return result;
}

// use NTN*D = R to solve the n-1 control points D
Eigen::MatrixXd build_matrix_R(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras, const std::vector<Vector3d>& points) {
	int n = U.size() - 2 - degree; // there are n+1 control points;
	int m = paras.size() - 1;  // there are m+1 points to fit
	Eigen::MatrixXd result(n - 1, 3);
	std::vector<Vector3d> ri(m);
	for (int i = 1; i < m; i++) {// from 1 to m-1
		ri[i] = 
			points[i] - points[0] * Nip(0, degree, paras[i], U) - points[m] * Nip(n, degree, paras[i], U);
	}
	
	for (int i = 0; i < n - 1; i++) {
		Vector3d row(0, 0, 0);
		for (int j = 1; j < m; j++) { // j from 1 to m-1
			row += Nip(i + 1, degree, paras[j], U)*ri[j];
		}
		result.row(i) = row;
	}
	return result; 
}

Eigen::MatrixXd solve_curve_control_points(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras, const std::vector<Vector3d>& points) {
	int npoints = points.size();
	int n = U.size() - 2 - degree; // there are n+1 control points;
	Eigen::MatrixXd result(n + 1, 3);
	assert(npoints == paras.size());
	Eigen::MatrixXd N = build_matrix_N(degree, U, paras);
	Eigen::MatrixXd NTN = N.transpose()*N;
	Eigen::MatrixXd R = build_matrix_R(degree, U, paras, points);// NTN*D=R

	bool check_error = true;
	double error;
	Eigen::MatrixXd interior = slove_linear_system(NTN, R, check_error, error);
	std::cout << "error, " << error << std::endl;
	
	result.row(0) = points[0];
	result.row(n) = points[npoints - 1];
	result.middleRows(1,n-1) = interior;
	return result;
}