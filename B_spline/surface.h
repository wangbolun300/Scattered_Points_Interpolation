#pragma once
#include<sparse_interp/basis.h>
#include <array>
namespace SIBSplines{
Vector3d BSplineSurfacePoint(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V, const double upara,
	const double vpara, const std::vector<std::vector<Vector3d>>& control);



std::vector<double> get_iso_line_parameters(const int degree1, const int degree2, const bool v_direction, const int line_id,
	const std::vector<double>& Ugrid, const std::vector<double>& Vgrid, const Eigen::MatrixXi& grid_map);


void output_timing();
}
