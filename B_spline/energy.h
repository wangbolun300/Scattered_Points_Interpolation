#pragma once
#include <sparse_interp/basis.h>
#include<sparse_interp/surface.h>
namespace SIBSplines{


// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
// the integration domain is [u1, u2]
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2,
	PolynomialBasis& basis, const bool uv);
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2);


void push_control_point_list_into_surface(Bsurface& surface, const std::vector<Vector3d>& cps);
void push_p_lambda_vector_to_control_points(const Eigen::MatrixXd &pl,
												const int dimension, std::vector<Vector3d> &control_points);
}