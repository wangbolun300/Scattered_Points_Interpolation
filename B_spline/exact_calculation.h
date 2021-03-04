#pragma once
#ifdef SPARSE_INTERP_WITH_GMP
#include<Rational.hpp>
#include<Eigen/Core>
#include<vector>
typedef Eigen::Matrix<Rational, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

Rational Nip_Rational(const int i, const int p, const double u, const std::vector<double> &U);
int rank(const MatrixXs& matrix);


#endif