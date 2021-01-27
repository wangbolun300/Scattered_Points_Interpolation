#include"basis.h"
#include"curve.h"
#include <igl/opengl/glfw/Viewer.h>
#include<iostream>
#include<array>
#include"test.h"
void show_basis(const std::vector<double>&b_vec) {
	for (int i = 0; i < b_vec.size(); i++) {
		std::cout << b_vec[i] << " ";
	}
	std::cout << std::endl;
}
void test_opengl() {
	// Inline mesh of a cube
	const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 1.0, 0.0,
		0.0, 1.0, 1.0,
		1.0, 0.0, 0.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 0.0,
		1.0, 1.0, 1.0).finished();
	const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) <<
		1, 7, 5,
		1, 3, 7,
		1, 4, 3,
		1, 2, 4,
		3, 8, 7,
		3, 4, 8,
		5, 7, 8,
		5, 8, 6,
		1, 5, 6,
		1, 6, 2,
		2, 6, 8,
		2, 8, 4).finished().array() - 1;

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.launch();
}


void draw_a_line() {
	Eigen::MatrixXd pts(3,3),Color(3,3);
	//
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Color(i, j) = 0;
			pts(i, j) = 0;
		}
	}
	pts(0, 0) = 1; pts(0, 1) = 2; pts(0, 2) = 3;
	std::cout << pts << std::endl;
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_points(pts,Color);
	viewer.launch();
}
void draw_a_curve() {
	int p = 3;
	// 8 control points
	std::vector<double>U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };
	std::array<double, 8> us = { {0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9} };
	std::vector<Vector3d> pts(8);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	igl::opengl::glfw::Viewer viewer;
	Eigen::MatrixXd  Color(8, 3), Hcolor(100, 3), OneColor(1, 3), PTS(8, 3), Curve(100, 3);
	OneColor(0, 0) = 0; OneColor(0, 1) = 0; OneColor(0, 2) = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 3; j++) {
			Color(i, j) = 0;
			
		}
		PTS.row(i) = pts[i];
		
	}
	for (int i = 0; i < 100; i++) {
		Curve.row(i) = BsplinePoint(3, U, i / double(100), pts);
		Hcolor(i, 0) = 0.5; Hcolor(i, 1) = 0.5; Hcolor(i, 2) = 0.5;
	}
	std::cout << Curve << std::endl;
	viewer.data().set_points(PTS, Color);
	viewer.data().add_points(Curve, Hcolor);
	viewer.launch();
	//for (int i = 0; i < 8; i++) {
	//	
	//}
}
void plot_fitting_result() {
	Eigen::MatrixXd control_pts;
	Eigen::MatrixXd control_pts_color;
	Eigen::MatrixXd curve_pts; Eigen::MatrixXd curve_pts_color;
	Eigen::MatrixXd target_pts; Eigen::MatrixXd target_pts_color;
	test_fitting(control_pts, control_pts_color,
			curve_pts, curve_pts_color, target_pts, target_pts_color);
	Eigen::MatrixXi edges;
	vertices_to_edges(curve_pts, edges);
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_points(control_pts, control_pts_color);
	//viewer.data().add_points(curve_pts, curve_pts_color);
	viewer.data().add_points(target_pts, target_pts_color);

	viewer.data().set_edges(curve_pts, edges, curve_pts_color.row(0));
	viewer.launch();
}
int main() {
	//test_opengl();
	//int p = 3;
	//std::vector<double>U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };
	//draw_a_line();
	//draw_a_curve();
	//test_fitting();
	//plot_fitting_result();
	test_curve_knot_fixing();
	return 0;
}