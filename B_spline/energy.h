#pragma once
#include <basis.h>
#include<surface.h>
class PolynomialBasis {
public:
	PolynomialBasis(Bsurface& surface);
	PolynomialBasis();
	void init(Bsurface& surface);
	std::vector<double> poly(const int id, const double value, const bool UVknot);
	void clear();
	std::vector<std::vector<std::vector<double>>> Ubasis;
	std::vector<std::vector<std::vector<double>>> Vbasis;
private:
	std::vector<double> Uknot;
	std::vector<double> Vknot;
	int degree1;
	int degree2;
	int nu;
	int nv;
	int inited = false;
	std::vector<std::vector<std::vector<double>>> calculate(const bool uorv);// 0 checking u; 1 checking v
};
class PartialBasis {
public:
	//PartialBasis(PolynomialBasis& basis, Bsurface& surface);
	PartialBasis(Bsurface& surface);
	PartialBasis();
	void init(Bsurface& surface);
	std::vector<double> poly(const int id, const double value, const bool UVknot, int partial);
	std::vector<double> Uknot;
	std::vector<double> Vknot;
	int degree1;
	int degree2;
	void clear();
private:
	std::vector<std::vector<std::vector<double>>> Ubasis;
	std::vector<std::vector<std::vector<double>>> Vbasis;
	std::vector<std::vector<std::vector<double>>> Ubasis_1;
	std::vector<std::vector<std::vector<double>>> Vbasis_1;
	std::vector<std::vector<std::vector<double>>> Ubasis_2;
	std::vector<std::vector<std::vector<double>>> Vbasis_2;
	std::vector<std::vector<std::vector<double>>> do_partial(const
		std::vector<std::vector<std::vector<double>>>&basis);
};
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
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2,
	PolynomialBasis& basis, const bool uv);
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2);

void solve_control_points_for_fairing_surface(Bsurface& surface, const Eigen::MatrixXd& paras,
	const Eigen::MatrixXd & points, PartialBasis& basis);

Eigen::MatrixXd surface_energy_calculation(Bsurface& surface, PartialBasis& basis,
	const int discrete, Eigen::MatrixXd &energy_uu, Eigen::MatrixXd &energy_vv, Eigen::MatrixXd& energy_uv);

// [U[which],U[which+1]) is the problematic one
void detect_max_energy_interval(Bsurface& surface, const Eigen::MatrixXd& energy, const Eigen::MatrixXd &energy_uu,
	const Eigen::MatrixXd & energy_vv, bool& uorv, int &which, double &em);
void iteratively_approximate_method(int degree1, int degree2,
	std::vector<double>& Uknot, std::vector<double>& Vknot,
	const Eigen::MatrixXd& param, const Eigen::MatrixXd& ver,
	const double tolerance,
	std::vector<Bsurface> &surfaces, const double per);