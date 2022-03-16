#include <sparse_interp/Types.hpp>

// This is an example of using our code to generating a B-spline
//  surface by interpolating the vertices of a triangle mesh

using namespace SIBSplines;
std::string example_root_path(SI_MESH_DIR);
void run_mesh_reconstruction()
{
	double precision = 0;
	Eigen::MatrixXd ver;
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, paramout;
	std::string modelname = "tiger.obj";
	std::string meshfile = example_root_path + modelname;
	std::cout << "reading mesh model: " << meshfile << std::endl;
	// mesh parametrization, and print out the parametrization result as a obj mesh.
	mesh_parameterization(meshfile, ver, param, F);
	paramout.resize(param.rows(), 3);
	Eigen::VectorXd param_zero = Eigen::VectorXd::Zero(param.rows());
	paramout << param, param_zero;
	write_triangle_mesh(example_root_path + "param_" + modelname, paramout, F);

	// construct the surface object
	Bsurface surface;
	// set up the initial parameters.
	int nbr = param.rows();					// the number of data points
	surface.degree1 = 3;					// degree of u direction
	surface.degree2 = 3;					// degree of v direction
	surface.U = {{0, 0, 0, 0, 1, 1, 1, 1}}; // the initial U knot vector
	surface.V = surface.U;					// the initial V knot vector
	int target_steps = 10;					// the number of iterations for constructing lists $L$.
	bool enable_max_fix_nbr = true;			// progressively update the knot vectors to make the two knot vectors balanced in length.
	double delta = 0.9;						// the parameter to improve the solving stability
	double per = 0.5;						// the parameter inherited from [Wen-Ke Wang et al, 2008, CAD]
	// generate knot vectors to make sure the data points can be interpolated
	surface.generate_interpolation_knot_vectors(surface.degree1, surface.degree2, surface.U, surface.V, param, delta, per, target_steps, enable_max_fix_nbr);
	std::cout << "knot vectors generated" << std::endl;

	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	int visual_nbr = 200; // the discretization scale for the output surface. The mesh will be 200x200

	// basis contains all the basis functions and their 1 and 2 order diffenrential form.
	PartialBasis basis(surface);

	// 	solve the control points to obtain the surface.
	surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
	std::cout << "surface solved" << std::endl;

	// convert B-spline surface into a triangle mesh
	surface.surface_visulization(surface, visual_nbr, SPs, SFs);

	precision = surface.max_interpolation_err(ver, param, surface);
	std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;

	write_points(example_root_path + "pts" + std::to_string(nbr) + "_m_" + modelname, ver);
	write_triangle_mesh(example_root_path + "intp_" + "p" + std::to_string(nbr) + "_m_" + modelname, SPs, SFs);
}

int main()
{

	run_mesh_reconstruction();

	std::cout << "done !" << std::endl;
	return 0;
}
