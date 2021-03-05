#pragma once
#ifdef SPARSE_INTERP_WITH_GMP
#include<Rational.hpp>
#include<Eigen/Core>
#include<vector>
#include<Types.hpp>


Rational Nip_Rational(const int i, const int p, const double u, const std::vector<double> &U);
int rank(const MatrixXs& matrix);


#endif