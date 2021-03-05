#pragma once
#include<Eigen/Core>
#include<array>
#include<iostream>
#include<vector>
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector2i Vector2i;
static const int STAIR_FORWARD = 0;
static const int STAIR_BACKWARD = 1;
static const int STAIR_WHOLE = 2;
static const int STAIR_HIGHEST = 3;
void print_vector(const std::vector<double>& input);
void print_vector(const std::vector<int>& input);

#ifdef SPARSE_INTERP_WITH_GMP
#include<Rational.hpp>
typedef Eigen::Matrix<Rational, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
bool selected_rows_have_solution_rational(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	const std::vector<int> &row_id, const int dimension);
#endif