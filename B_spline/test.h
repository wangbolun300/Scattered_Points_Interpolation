#pragma once
#include"curve.h"
void test_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color);
void vertices_to_edges(const Eigen::MatrixXd& pts, Eigen::MatrixXi &edges);
void test_curve_knot_fixing();
void visual_curve_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color);
void test_surface_visual(Eigen::MatrixXd &ver, Eigen::MatrixXi& faces);
void parameter_grid_to_mesh(const Eigen::MatrixXd& uv, Eigen::MatrixXd &ver, Eigen::MatrixXi& edges);
void test_surface_knot_preprocessing(Eigen::MatrixXd &points, Eigen::MatrixXd& knotP, Eigen::MatrixXi& knotE);
void test_knot_fixing(Eigen::MatrixXd &points, Eigen::MatrixXd& knotP, Eigen::MatrixXi& knotE);