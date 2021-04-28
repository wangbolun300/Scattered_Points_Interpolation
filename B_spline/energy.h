#pragma once
#include <basis.h>
std::vector<double> polynomial_simplify(const std::vector<double>& poly);
std::vector<double> polynomial_add(const std::vector<double>& poly1, const std::vector<double>& poly2);
std::vector<double> polynomial_times(const std::vector<double>& poly1, const std::vector<double>& poly2);
std::vector<double> polynomial_times(const std::vector<double>& poly1, const double& nbr);
std::vector<double> Nip_func(const int i, const int p, const double u, const std::vector<double> &U);
double polynomial_value(const std::vector<double>& poly, const double para);
std::vector<double> polynomial_integration(const std::vector<double>& poly);
double polynomial_integration(const std::vector<double>& poly, const double lower, const double upper);

// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
// the integration domain is [u1, u2]
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2);