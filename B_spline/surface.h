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
	int nu();
	int nv();
};
Vector3d BSplineSurfacePoint(const Bsurface& surface, const double upara, const double vpara);


// pre-processing to fix the problem that too many points corresponding to one interval of 
// u-v domain
void fix_surface_grid_parameter_too_many(const int degree1, const std::vector<double>& Uin,
	const int degree2, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras,
	std::vector<double>& Uout, std::vector<double>& Vout);

// check if there is solution(s) for interpolation problem
bool equation_has_solution(const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b);

bool equation_has_solution(const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b, int& rank_diff);