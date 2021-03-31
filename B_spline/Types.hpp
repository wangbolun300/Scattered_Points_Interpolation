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

static const double SCALAR_ZERO = 1e-8;
void print_vector(const std::vector<double>& input);
void print_vector(const std::vector<int>& input);

// type converters
void vertices_to_edges(const Eigen::MatrixXd& pts, Eigen::MatrixXi &edges);
Eigen::MatrixXd vector_to_matrix_3d(const std::vector<Vector3d>& v);
std::vector<Vector3d> matrix3d_to_vector(const Eigen::MatrixXd& v);
Eigen::MatrixXd list_to_matrix_3d(const std::vector<std::vector<double>>& v);
Eigen::MatrixXd list_to_matrix_3d(const std::vector<Vector3d>& v, const std::vector<int>& selected);

#ifdef SPARSE_INTERP_WITH_GMP
#include<Rational.hpp>
typedef Eigen::Matrix<Rational, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
bool selected_rows_have_solution_rational(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	const std::vector<int> &row_id, const int dimension);
#endif