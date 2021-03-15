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