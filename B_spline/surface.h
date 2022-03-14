#pragma once
#include<basis.h>
#include <array>
Vector3d BSplineSurfacePoint(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V, const double upara,
	const double vpara, const std::vector<std::vector<Vector3d>>& control);
struct Bsurface {
	int degree1;
	int degree2;
	std::vector<double> U;  
	std::vector<double> V;  
	double upara;
	double vpara; 
	std::vector<std::vector<Vector3d>> control_points;
	int nu();// nu + 1 is the number of control points in u direction
	int nv();
};

struct per_too_large {
	bool flag;
};
Vector3d BSplineSurfacePoint(const Bsurface& surface, const double upara, const double vpara);

std::vector<double> get_iso_line_parameters(const int degree1, const int degree2, const bool v_direction, const int line_id,
	const std::vector<double>& Ugrid, const std::vector<double>& Vgrid, const Eigen::MatrixXi& grid_map);
void generate_interpolation_knot_vectors(int degree1, int degree2,
	std::vector<double>& Uknot, std::vector<double>& Vknot,
	const Eigen::MatrixXd& param_original,
	double &per_ours, const double per, const int target_steps, const bool enable_max_fix_nbr);

void output_timing();

double max_interpolation_err(const Eigen::MatrixXd&ver, const Eigen::MatrixXd& param, Bsurface& surface);
Eigen::MatrixXd interpolation_err_for_apprximation(const Eigen::MatrixXd&ver,
	const Eigen::MatrixXd& param, Bsurface& surface, double &max_err);
