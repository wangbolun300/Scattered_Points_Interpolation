#pragma once
// This file contains the benchmark surfaces.
//
//
//
//
#include<cmath>
#include <Types.hpp>
double peak_function(const double x, const double y) {
	double r = 3 * pow(1 - x, 2)*exp(-x * x - (y + 1)*(y + 1)) - 10 * (0.2*x - pow(x, 3) - pow(y, 5))*
		exp(-x * x - y * y) - 1 / 3 * exp(-pow(x + 1, 2) - y * y);
	if (fabs(r) < SCALAR_ZERO) {
		r = 0;
	}
	r /= 2;
	return r;
}

Eigen::MatrixXd get_peak_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = 6 * para3d[0];
		double y = 6 * para3d[1];
		ver.row(i - skip) << x, y, peak_function(x, y);
	}
	ver.row(nbr - skip - 4) << -6, -6, peak_function(-6, -6);
	ver.row(nbr - skip - 3) << -6, 6, peak_function(-6, 6);
	ver.row(nbr - skip - 2) << 6, 6, peak_function(6, 6);
	ver.row(nbr - skip - 1) << 6, -6, peak_function(6, -6);
	return ver;
}

double contour_function(const double x, const double y) {
	return 4 * x*exp(-x * x - y * y);
}
Eigen::MatrixXd get_contour_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	double s = 2;
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = s * para3d[0];
		double y = s * para3d[1];
		ver.row(i - skip) << x, y, contour_function(x, y);
	}
	ver.row(nbr - skip - 4) << -s, -s, contour_function(-s, -s);
	ver.row(nbr - skip - 3) << -s, s, contour_function(-s, s);
	ver.row(nbr - skip - 2) << s, s, contour_function(s, s);
	ver.row(nbr - skip - 1) << s, -s, contour_function(s, -s);
	return ver;
}
double hyperbolic_function(const double x, const double y) {
	return x *x - y * y;
}
Eigen::MatrixXd get_hyperbolic_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	double s = 1;
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = s * para3d[0];
		double y = s * para3d[1];
		ver.row(i - skip) << x, y, hyperbolic_function(x, y);
	}
	ver.row(nbr - skip - 4) << -s, -s, hyperbolic_function(-s, -s);
	ver.row(nbr - skip - 3) << -s, s, hyperbolic_function(-s, s);
	ver.row(nbr - skip - 2) << s, s, hyperbolic_function(s, s);
	ver.row(nbr - skip - 1) << s, -s, hyperbolic_function(s, -s);
	return ver;
}
double sinus_function(const double x, const double y) {
	return sin(3 * 3.1415926*(x*x + y * y)) / 10;
}
Eigen::MatrixXd get_sinus_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	double s = 1;
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = s * para3d[0];
		double y = s * para3d[1];
		ver.row(i - skip) << x, y, sinus_function(x, y);
	}
	ver.row(nbr - skip - 4) << -s, -s, sinus_function(-s, -s);
	ver.row(nbr - skip - 3) << -s, s, sinus_function(-s, s);
	ver.row(nbr - skip - 2) << s, s, sinus_function(s, s);
	ver.row(nbr - skip - 1) << s, -s, sinus_function(s, -s);
	return ver;
}

Vector3d bilinear_function(const Vector3d& v0s, const Vector3d&v0e, const Vector3d& v1s, const Vector3d& v1e,
	const double u, const double v) {
	Vector3d m0 = (v0e - v0s)*u + v0s;
	Vector3d m1 = (v1e - v1s)*u + v1s;
	Vector3d m = (m1 - m0)*v + m0;
	return m;
}
Eigen::MatrixXd get_bilinear_sample_points(const int nbr, const int skip, Eigen::MatrixXd& param) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	param.resize(nbr - skip, 2);
	double s = 1;// domain scale is [-s, s]x[-s, s]
	Vector3d v0s(0, 0, 1);
	Vector3d v0e(1, 1, 0);
	Vector3d v1s(0, 1, 1);
	Vector3d v1e(0, 0, 0);
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double u = (1 + para3d[0]) / 2;
		double v = (1 + para3d[1]) / 2;
		param.row(i - skip) << u, v;
		ver.row(i - skip) = bilinear_function(v0s, v0e, v1s, v1e, u, v);
	}

	param.row(nbr - skip - 4) << 0, 0;
	param.row(nbr - skip - 3) << 0, 1;
	param.row(nbr - skip - 2) << 1, 1;
	param.row(nbr - skip - 1) << 1, 0;
	
	ver.row(nbr - skip - 4) = bilinear_function(v0s, v0e, v1s, v1e, 0, 0);
	ver.row(nbr - skip - 3) = bilinear_function(v0s, v0e, v1s, v1e, 0, 1);
	ver.row(nbr - skip - 2) = bilinear_function(v0s, v0e, v1s, v1e, 1, 1);
	ver.row(nbr - skip - 1) = bilinear_function(v0s, v0e, v1s, v1e, 1, 0);
	
	return ver;
}


// u0 v0 are in [0,1].
// the function u in [1,2], v in [-pi/2, pi/2]
Vector3d snail_function(double u0, double v0) {
	const double pi = 3.1415926;
	double u = 1+ u0;
	double v = -pi/2 + pi * v0;
	double x = u * cos(v)*sin(u);
	double y = u * cos(u)*cos(v);
	double z = -u * sin(v);
	return Vector3d(x, y, z);
}