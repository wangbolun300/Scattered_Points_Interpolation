#include"test.h"
#include"curve.h"
#include<iostream>
#include<surface.h>
#include<cmath>

Eigen::MatrixXd vector_to_matrix_3d(const std::vector<Vector3d>& v) {
	Eigen::MatrixXd result(v.size(), 3);
	for (int i = 0; i < v.size(); i++) {
		result.row(i) = v[i];
	}
	return result;
}
void vertices_to_edges(const Eigen::MatrixXd& pts, Eigen::MatrixXi &edges) {
	edges.resize(pts.rows() - 1, 2);
	for (int i = 0; i < edges.rows(); i++) {
		edges(i, 0) = i; edges(i, 1) = i + 1;
	}
}
void test_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color) {
	int nbr_curve_pts = 500;
	// 8 control points
	std::vector<double> U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };

	int Qnbr = 12;
	int degree = 3;
	std::vector<Vector3d> pts(Qnbr);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	pts[8] = Vector3d(0, 2, 3);
	pts[9] = Vector3d(0, 3, 0);
	pts[10] = Vector3d(0, 2, 2);
	pts[11] = Vector3d(0, 4, 0);
	std::vector<double> paras = Centripetal_parameterization(pts);
	for (int i = 0; i < paras.size(); i++) {
		std::cout << paras[i] << std::endl;
	}
	Eigen::MatrixXd Control = solve_curve_control_points(degree, U, paras, pts);
	std::cout << "control points:" << std::endl << Control << std::endl;

	Vector3d color1(0, 0, 0), color2(0.5, 0.5, 0.5), color3(1, 0, 0);
	//////////////////////
	// target points and color
	target_pts = vector_to_matrix_3d(pts);

	target_pts_color.resize(Qnbr, 3);
	for (int i = 0; i < Qnbr; i++) {
		target_pts_color.row(i) = color3;
	}
	// target points and color
	////////////////////////
	//////////////////////
	// control points and color
	control_pts = Control;

	control_pts_color.resize(Control.rows(), 3);
	for (int i = 0; i < Control.rows(); i++) {
		control_pts_color.row(i) = color1;
	}
	// control points and color
	////////////////////////
	////////////////////////
	//set curve points and color
	curve_pts.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		double temp_para = i / double(nbr_curve_pts);
		curve_pts.row(i) = BsplinePoint(degree, U, temp_para, control_pts);
	}

	curve_pts_color.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		curve_pts_color.row(i) = color2;
	}
	//set curve points and color
	/////////////////////////


}
void visual_curve_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color) {
	int nbr_curve_pts = 500;
	// 8 control points
	//std::vector<double> U_init = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };
	std::vector<double> U_init = { {0,0,0,0,0.1,1,1,1,1} };
	int Qnbr = 12;
	int degree = 3;
	std::vector<Vector3d> pts(Qnbr);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	pts[8] = Vector3d(0, 2, 3);
	pts[9] = Vector3d(0, 3, 0);
	pts[10] = Vector3d(0, 2, 2);
	pts[11] = Vector3d(0, 4, 0);

	std::vector<double> paras = Centripetal_parameterization(pts);
	std::cout << "paras " << std::endl;
	for (int i = 0; i < paras.size(); i++) {
		std::cout << paras[i] << std::endl << std::endl;
	}
	std::vector<double> result_vector = fix_knot_vector_to_interpolate_curve(degree, U_init, paras, pts);
	std::cout << "fixed " << std::endl;
	for (int i = 0; i < result_vector.size(); i++) {
		std::cout << result_vector[i] << std::endl << std::endl;
	}
	Eigen::MatrixXd Control = solve_curve_control_points(degree, result_vector, paras, pts);
	std::cout << "control points:" << std::endl << Control << std::endl;

	Vector3d color1(0, 0, 0), color2(0.5, 0.5, 0.5), color3(1, 0, 0);
	//////////////////////
	// target points and color
	target_pts = vector_to_matrix_3d(pts);

	target_pts_color.resize(Qnbr, 3);
	for (int i = 0; i < Qnbr; i++) {
		target_pts_color.row(i) = color3;
	}
	// target points and color
	////////////////////////
	//////////////////////
	// control points and color
	control_pts = Control;

	control_pts_color.resize(Control.rows(), 3);
	for (int i = 0; i < Control.rows(); i++) {
		control_pts_color.row(i) = color1;
	}
	// control points and color
	////////////////////////
	////////////////////////
	//set curve points and color
	curve_pts.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		double temp_para = i / double(nbr_curve_pts);
		curve_pts.row(i) = BsplinePoint(degree, result_vector, temp_para, control_pts);
	}

	curve_pts_color.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		curve_pts_color.row(i) = color2;
	}
	//set curve points and color
	/////////////////////////

}
void test_curve_knot_fixing() {
	/*Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
		Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color*/
	int nbr_curve_pts = 500;
	// 8 control points
	std::vector<double> U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };

	int Qnbr = 12;
	int degree = 3;
	std::vector<Vector3d> pts(Qnbr);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	pts[8] = Vector3d(0, 2, 3);
	pts[9] = Vector3d(0, 3, 0);
	pts[10] = Vector3d(0, 2, 2);
	pts[11] = Vector3d(0, 4, 0);

	std::vector<double> paras = Centripetal_parameterization(pts);
	std::cout << "paras " << std::endl;
	for (int i = 0; i < paras.size(); i++) {
		std::cout << paras[i] << std::endl << std::endl;
	}
	std::vector<double> result_vector = fix_knot_vector_to_interpolate_curve(degree, U, paras, pts);
	std::cout << "fixed " << std::endl;
	for (int i = 0; i < result_vector.size(); i++) {
		std::cout << result_vector[i] << std::endl << std::endl;
	}
}
double peak_function(const double x, const double y) {
	double r = 3 * pow(1 - x, 2)*exp(-x * x - (y + 1)*(y + 1)) - 10 * (0.2*x - pow(x, 3) - pow(y, 5))*
		exp(-x * x - y * y) - 1 / 3 * exp(-pow(x + 1, 2) - y * y);
	return r;
}

void B_spline_surface_to_mesh(const Bsurface &surface, const int pnbr, Eigen::MatrixXd &ver, Eigen::MatrixXi& faces) {
	std::vector<std::vector<Vector3d>> pts;
	ver.resize(pnbr*pnbr, 3);
	int verline = 0;
	for (int i = 0; i < pnbr; i++) {
		for (int j = 0; j < pnbr; j++) {
			double upara = double(i) / (pnbr - 1);
			double vpara = double(j) / (pnbr - 1);
			ver.row(verline) = BSplineSurfacePoint(surface, upara, vpara);
			verline++;
		}
	}
	faces.resize(2 * (pnbr - 1)*(pnbr - 1), 3);
	int fline = 0;
	for (int i = 0; i < pnbr - 1; i++) {
		for (int j = 0; j < pnbr - 1; j++) {
			faces.row(fline) = Vector3i(i + pnbr*j, i + pnbr*(j + 1), i + pnbr*(1 + j) + 1);
			faces.row(fline + 1) = Vector3i(i + pnbr * j, i + pnbr * (1 + j) + 1, i + pnbr * j + 1);
			fline += 2;
		}
	}
}
void test_surface_visual(Eigen::MatrixXd &ver, Eigen::MatrixXi& faces) {
	Bsurface surface;
	surface.degree1 = 3;
	surface.degree2 = 3;
	surface.U = { {0,0,0,0,0.1,0.3,0.5,0.7,0.9,1,1,1,1} };// nu=13-2-3=8
	surface.V = { {0,0,0,0,0.1,0.2,0.3,0.5,0.7,0.8,0.9,1,1,1,1} };// nv=10
	int nu = surface.nu();
	int nv = surface.nv();
	std::vector<std::vector<Vector3d>> control(nu + 1);
	for (int i = 0; i < control.size(); i++) {
		control[i].resize(nv + 1);
		for (int j = 0; j < control[i].size(); j++) {
			double x = -3 + double(i) * 6 / nu;
			double y = -3 + double(j) * 6 / nv;
			double z = peak_function(x, y);
			control[i][j] = Vector3d(x, y, z);
		}
	}
	surface.control_points = control;

	B_spline_surface_to_mesh(surface, 50, ver, faces);
}

