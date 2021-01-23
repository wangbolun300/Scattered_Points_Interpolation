#pragma once
#include"curve.h"
void test_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color);
void vertices_to_edges(const Eigen::MatrixXd& pts, Eigen::MatrixXi &edges);