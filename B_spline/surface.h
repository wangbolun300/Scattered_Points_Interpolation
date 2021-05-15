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
std::vector<double> knot_vector_insert_one_value(const std::vector<double>& U, const double value);

void fix_knot_vector_to_interpolate_surface(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, 
	std::vector<double>& Uout, std::vector<double>& Vout);

bool selected_rows_have_solution(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	const std::vector<int> &row_id, const int dimension);

void easist_way_to_fix_knot_vector_to_interpolate_surface(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	std::vector<double>& Uout, std::vector<double>& Vout);

// get feasible control point matrix. if checking v direction control points (v_direction=true), make sure that the Uin knot vector
// is already fixed. for the explanation of Ugrid, Vgrid and UVmap, see function generate_UV_grid() in 'mesh_processing.h'
Eigen::MatrixXi get_feasible_control_point_matrix(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin, const bool v_direction,
	const std::vector<double>& Ugrid, const std::vector<double>&Vgrid, const Eigen::MatrixXi& UVmap,
	const int nbr_para, std::vector<std::vector<std::array<int, 2>>>&para_to_feasible);

// calculate weights and select ACP according to the weight
Eigen::MatrixXi calculate_active_control_points_from_feasible_control_points(const Eigen::MatrixXi& fcp, const bool v_direction,
	const std::vector<double> &Uknot, const std::vector<double> &Vknot,
	const Eigen::MatrixXd& paras, const int degree1, const int degree2, 
	std::vector<std::vector<std::array<int, 2>>> &para_to_feasible);

std::vector<double> get_iso_line_parameters_from_ACP(const Eigen::MatrixXi&ACP, const int id, const Eigen::MatrixXd& paras, const bool v_direction);

std::vector<double> get_iso_line_parameters(const int degree1, const int degree2, const bool v_direction, const int line_id,
	const std::vector<double>& Ugrid, const std::vector<double>& Vgrid, const Eigen::MatrixXi& grid_map);
void generate_interpolation_knot_vectors(const bool start_from_v_direction, int degree1, int degree2,
	std::vector<double>& Uknot, std::vector<double>& Vknot,
	const Eigen::MatrixXd& param_original, Eigen::MatrixXd& param_perturbed, const Eigen::MatrixXi& F, const int mesh_perturbation_level);