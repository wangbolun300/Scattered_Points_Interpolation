#pragma once
#include <Eigen/Core>
#include <array>
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

//#define NO_SELECTING_ACP
//#define NAIVE_SELECTING_ACP
//#define WEIGHT_NAIVE

namespace SIBSplines
{
	typedef Eigen::SparseMatrix<double> SparseMatrixXd;
	typedef Eigen::Triplet<double> Trip;
	typedef Eigen::Vector3d Vector3d;
	typedef Eigen::Vector3i Vector3i;
	typedef Eigen::Vector2d Vector2d;
	typedef Eigen::Vector2i Vector2i;

	static const double SCALAR_ZERO = 1e-8;
	
	struct per_too_large
	{
		bool flag;
	};

	namespace ply_operations
	{
		std::vector<double> polynomial_simplify(const std::vector<double> &poly);
		std::vector<double> polynomial_add(const std::vector<double> &poly1, const std::vector<double> &poly2);
		std::vector<double> polynomial_times(const std::vector<double> &poly1, const std::vector<double> &poly2);
		std::vector<double> polynomial_times(const std::vector<double> &poly1, const double &nbr);
		double polynomial_value(const std::vector<double> &poly, const double para);
		std::vector<double> polynomial_integration(const std::vector<double> &poly);
		double polynomial_integration(const std::vector<double> &poly, const double lower, const double upper);
	}
	
	class PartialBasis;
	class Bsurface
	{
	public:
	Bsurface(){};
		int degree1;
		int degree2;
		std::vector<double> U;
		std::vector<double> V;
		double upara;
		double vpara;
		std::vector<std::vector<Vector3d>> control_points;
		int nu(); // nu + 1 is the number of control points in u direction
		int nv();
		void generate_interpolation_knot_vectors(int degree1, int degree2,
												 std::vector<double> &Uknot, std::vector<double> &Vknot,
												 const Eigen::MatrixXd &param_original,
												 double &per_ours, const double per, const int target_steps, const bool enable_max_fix_nbr);
		void generate_interpolation_knot_vectors(int degree1, int degree2,
												 std::vector<double> &Uknot, std::vector<double> &Vknot,
												 const Eigen::MatrixXd &param_original,
												 const double per_ours, const double per, const int target_steps, const bool enable_max_fix_nbr, per_too_large &per_flag);
		void solve_control_points_for_fairing_surface(Bsurface &surface, const Eigen::MatrixXd &paras,
													  const Eigen::MatrixXd &points, PartialBasis &basis);
		// calculate thin-plate-energy in region [Ui, U(i+1)]x[Vj, V(j+1)]
		Eigen::MatrixXd surface_energy_calculation(Bsurface &surface, PartialBasis &basis,
												   const int discrete, Eigen::MatrixXd &energy_uu, Eigen::MatrixXd &energy_vv, Eigen::MatrixXd &energy_uv);

		// [U[which],U[which+1]) is the problematic one
		void detect_max_energy_interval(Bsurface &surface, const Eigen::MatrixXd &energy, const Eigen::MatrixXd &energy_uu,
										const Eigen::MatrixXd &energy_vv, bool &uorv, int &which, double &em);

		Eigen::MatrixXd interpolation_err_for_apprximation(const Eigen::MatrixXd &ver,
														   const Eigen::MatrixXd &param, Bsurface &surface, double &max_err);
		double max_interpolation_err(const Eigen::MatrixXd &ver, const Eigen::MatrixXd &param, Bsurface &surface);
		void surface_visulization(Bsurface &surface, const int nbr, Eigen::MatrixXd &v, Eigen::MatrixXi &f);
		Vector3d BSplineSurfacePoint(const Bsurface &surface, const double upara, const double vpara);
		Vector3d BSplineSurfacePoint(const std::vector<std::vector<std::vector<double>>> &upolys, const std::vector<std::vector<std::vector<double>>> &vpolys,
									 const std::vector<std::vector<std::vector<double>>> &tpolys, double u, double v, double t, const std::vector<double> &U, const std::vector<double> &V,
									 const std::vector<double> &T, const std::vector<std::vector<std::vector<Vector3d>>> &CPs, const int udegree,
									 const int vdegree, const int tdegree);
		// insert in total n knots into U and V. 
		void RefineKnots(int nbr);

		void constructRegularF(const int vnbr, const int rnbr, Eigen::MatrixXi &F)
		{

			int cnbr = vnbr / rnbr;
			F.resize((cnbr - 1) * (rnbr - 1), 4);
			int fnbr = 0;
			for (int i = 0; i < cnbr - 1; i++)
			{
				for (int j = 0; j < rnbr - 1; j++)
				{
					int v0 = i * rnbr + j;
					int v1 = i * rnbr + j + 1;
					int v2 = (i + 1) * rnbr + j + 1;
					int v3 = (i + 1) * rnbr + j;
					F(fnbr, 0) = v0;
					F(fnbr, 1) = v1;
					F(fnbr, 2) = v2;
					F(fnbr, 3) = v3;
					fnbr++;
				}
			}
			// std::cout<<"check F, \n"<<F<<"\n";
		}
	};
	class Bcurve
	{
	public:
	Bcurve(){};
		int degree;
		std::vector<double> U;
		// double upara;
		std::vector<Vector3d> control_points;
		int nu(); // nu + 1 is the number of control points in u direction
		bool curve_can_be_interpolated(const std::vector<double> &U, const int degree, const std::vector<double> &paras,
									   int &prob_id);
		bool curve_can_be_interpolated(const std::vector<double> &U, const int degree, const Eigen::VectorXd &paras, int &prob_id);
		std::vector<double> fix_knot_vector_to_interpolate_curve(const int degree, const std::vector<double> &init_vec,
																		 const std::vector<double> &paras, const double per, bool &fully_fixed, const int fix_nbr);
		void solve_control_points_for_fairing_curve(Bcurve &curve, const std::vector<double> &paras,
													const std::vector<Vector3d> &pts, const double a, const double b);
		// the output is the curve.control_points
		// trying to find a curve minimizing the energy, while interpolating the points whose parameters are paras.
		void solve_control_points_for_fairing_curve(Bcurve &curve, const std::vector<double> &paras,
													const Eigen::MatrixXd &points, const double a, const double b);

		// least square method solve the control points
		Eigen::MatrixXd solve_curve_control_points(const int degree, const std::vector<double> &U,
												   const std::vector<double> &paras, const std::vector<Vector3d> &points);
		Vector3d BsplinePoint(const int degree, const std::vector<double> &U, const double para,
							  const std::vector<Vector3d> &pts);

		Vector3d BsplinePoint(const int degree, const std::vector<double> &U, const double para,
							  const Eigen::MatrixXd &pts);
	};
	class PolynomialBasis
	{
	public:
		PolynomialBasis(Bsurface &surface);// get all the bases of the B-spline surface into polynomial format
		PolynomialBasis();
		void init(Bsurface &surface);
		std::vector<double> poly(const int id, const double value, const bool UVknot);
		void clear();
		// the (i,j) th element of the basis is the basis function defined on {U_i,U_{i+1}}, the j ranges from 0 to degree.
		std::vector<std::vector<std::vector<double>>> Ubasis;
		std::vector<std::vector<std::vector<double>>> Vbasis;
		std::vector<double> Uknot;
		std::vector<double> Vknot;
		int degree1;
		int degree2;
		std::vector<std::vector<std::vector<double>>> calculate_single(const int degree, const std::vector<double> &knotVector);
		
	private:
		int nu;
		int nv;
		int inited = false;
		std::vector<std::vector<std::vector<double>>> calculate(const bool uorv); // 0 checking u; 1 checking v
	};

	// get the first or second partial differential polynomial of U or V knot vectors
	class PartialBasis
	{
	public:
		// PartialBasis(PolynomialBasis& basis, Bsurface& surface);
		PartialBasis(Bsurface &surface);
		PartialBasis();
		void init(Bsurface &surface);
		void init(PolynomialBasis &pb);
		// return a certain basis (in 0, 1 or 2 order partial differential) of N_{id}(value)
		// UVknot==1 -> check v
		std::vector<double> poly(const int id, const double value, const bool UVknot, int partial);
		std::vector<double> Uknot;
		std::vector<double> Vknot;
		int degree1;
		int degree2;
		void clear();

	// private:
		std::vector<std::vector<std::vector<double>>> Ubasis;
		std::vector<std::vector<std::vector<double>>> Vbasis;
		std::vector<std::vector<std::vector<double>>> Ubasis_1;
		std::vector<std::vector<std::vector<double>>> Vbasis_1;
		std::vector<std::vector<std::vector<double>>> Ubasis_2;
		std::vector<std::vector<std::vector<double>>> Vbasis_2;
		// do 1 order differential to a list (with different value as indices) of a list (with the id 'i' in N_{i,p}(u)) of polynomials
		std::vector<std::vector<std::vector<double>>> do_partial(const std::vector<std::vector<std::vector<double>>> &basis);
	};
	// type converters
	void vertices_to_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXi &edges);
	Eigen::MatrixXd vector_to_matrix_3d(const std::vector<Vector3d> &v);

	Eigen::MatrixXd slove_linear_system(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b,
										const bool check_error, double &relative_error);
	void write_points(const std::string& file, const Eigen::MatrixXd& ver);
	void write_csv(const std::string &file, const std::vector<std::string> titles, const std::vector<double> data);
	bool write_triangle_mesh(const std::string filename, const Eigen::MatrixXd &ver, Eigen::MatrixXi &f);
	// parametrize the mesh onto [0, 1]x[0, 1]
	void mesh_parameterization(
		const std::string &meshfile, Eigen::MatrixXd &V, Eigen::MatrixXd &param, Eigen::MatrixXi &F);
	void print_vector(const std::vector<double> &input);
	void print_vector(const std::vector<int> &input);

	std::vector<double> Nip_func(const int i, const int p, const double u, const std::vector<double> &U);

	std::vector<double> knot_vector_insert_one_value(const std::vector<double> &U, const double value);

	std::vector<double> fix_knot_vector_to_interpolate_curve_WKW(const int degree, const std::vector<double> &init_vec,
																 const std::vector<double> &paras, const double per, bool &fully_fixed, const int fix_nbr = -1);
	std::vector<double> basisValues(const int whichItv, const int degree, const std::vector<std::vector<std::vector<double>>> &basis, const double param);
}	