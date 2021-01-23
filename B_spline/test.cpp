#include"test.h"
#include"curve.h"
#include<iostream>

Eigen::MatrixXd vector_to_matrix_3d(const std::vector<Vector3d>& v) {
	Eigen::MatrixXd result(v.size(),3);
	for (int i = 0; i < v.size(); i++) {
		result.row(i) = v[i];
	}
	return result;
}

void test_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color, 
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color) {
	int nbr_curve_pts = 100;
	// 8 control points
	std::vector<double> U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };

	int Qnbr = 8;
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
	for (int i = 0; i < Control.rows(); i++) {
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
