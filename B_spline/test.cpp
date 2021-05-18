#include"test.h"
#include"curve.h"
#include<iostream>
#include<surface.h>
#include<cmath>
#include<mesh_processing.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/harmonic.h>
#include<energy.h>
#include <igl/write_triangle_mesh.h>
void test_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color) {
	int nbr_curve_pts = 500;
	// 8 control points
	std::vector<double> U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };

	int Qnbr = 12;
	int degree = 3;
	std::vector<Vector3d> pts(Qnbr);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	pts[8] = Vector3d(0, 2, 3);
	pts[9] = Vector3d(0, 3, 0);
	pts[10] = Vector3d(0, 2, 2);
	pts[11] = Vector3d(0, 4, 0);
	std::vector<double> paras = Centripetal_parameterization(pts);
	for (int i = 0; i < paras.size(); i++) {
		std::cout << paras[i] << std::endl;
	}
	Eigen::MatrixXd Control = solve_curve_control_points(degree, U, paras, pts);
	std::cout << "control points:" << std::endl << Control << std::endl;

	Vector3d color1(0, 0, 0), color2(0.5, 0.5, 0.5), color3(1, 0, 0);
	//////////////////////
	// target points and color
	target_pts = vector_to_matrix_3d(pts);

	target_pts_color.resize(Qnbr, 3);
	for (int i = 0; i < Qnbr; i++) {
		target_pts_color.row(i) = color3;
	}
	// target points and color
	////////////////////////
	//////////////////////
	// control points and color
	control_pts = Control;

	control_pts_color.resize(Control.rows(), 3);
	for (int i = 0; i < Control.rows(); i++) {
		control_pts_color.row(i) = color1;
	}
	// control points and color
	////////////////////////
	////////////////////////
	//set curve points and color
	curve_pts.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		double temp_para = i / double(nbr_curve_pts);
		curve_pts.row(i) = BsplinePoint(degree, U, temp_para, control_pts);
	}

	curve_pts_color.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		curve_pts_color.row(i) = color2;
	}
	//set curve points and color
	/////////////////////////


}
void visual_curve_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color) {
	int nbr_curve_pts = 500;
	// 8 control points
	//std::vector<double> U_init = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };
	std::vector<double> U_init = { {0,0,0,0,1,1,1,1} };
	int Qnbr = 12;
	int degree = 3;
	std::vector<Vector3d> pts(Qnbr);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	pts[8] = Vector3d(0, 2, 3);
	pts[9] = Vector3d(0, 3, 0);
	pts[10] = Vector3d(0, 2, 2);
	pts[11] = Vector3d(0, 4, 0);

	std::vector<double> paras = Centripetal_parameterization(pts);
	std::cout << "paras " << std::endl;
	for (int i = 0; i < paras.size(); i++) {
		std::cout << paras[i] << std::endl << std::endl;
	}
	std::vector<double> result_vector =
	//{ {0,0,0,0,0.1,1,1,1,1} };
		//fix_knot_vector_to_interpolate_curve_boolean(degree, U_init, paras);
		fix_knot_vector_to_interpolate_curve(degree, U_init, paras, pts);
	std::cout << "fixed " << std::endl;
	for (int i = 0; i < result_vector.size(); i++) {
		std::cout << result_vector[i] << std::endl << std::endl;
	}
	Eigen::MatrixXd Control = solve_curve_control_points(degree, result_vector, paras, pts);
	std::cout << "control points:" << std::endl << Control << std::endl;

	Vector3d color1(0, 0, 0), color2(0.5, 0.5, 0.5), color3(1, 0, 0);
	//////////////////////
	// target points and color
	target_pts = vector_to_matrix_3d(pts);

	target_pts_color.resize(Qnbr, 3);
	for (int i = 0; i < Qnbr; i++) {
		target_pts_color.row(i) = color3;
	}
	// target points and color
	////////////////////////
	//////////////////////
	// control points and color
	control_pts = Control;

	control_pts_color.resize(Control.rows(), 3);
	for (int i = 0; i < Control.rows(); i++) {
		control_pts_color.row(i) = color1;
	}
	// control points and color
	////////////////////////
	////////////////////////
	//set curve points and color
	curve_pts.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		double temp_para = i / double(nbr_curve_pts);
		curve_pts.row(i) = BsplinePoint(degree, result_vector, temp_para, control_pts);
	}

	curve_pts_color.resize(nbr_curve_pts, 3);
	for (int i = 0; i < nbr_curve_pts; i++) {
		curve_pts_color.row(i) = color2;
	}
	//set curve points and color
	/////////////////////////

}
void test_curve_knot_fixing() {
	/*Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
		Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color*/
	int nbr_curve_pts = 500;
	// 8 control points
	std::vector<double> U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };

	int Qnbr = 12;
	int degree = 3;
	std::vector<Vector3d> pts(Qnbr);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	pts[8] = Vector3d(0, 2, 3);
	pts[9] = Vector3d(0, 3, 0);
	pts[10] = Vector3d(0, 2, 2);
	pts[11] = Vector3d(0, 4, 0);

	std::vector<double> paras = Centripetal_parameterization(pts);
	std::cout << "paras " << std::endl;
	for (int i = 0; i < paras.size(); i++) {
		std::cout << paras[i] << std::endl << std::endl;
	}
	std::vector<double> result_vector = fix_knot_vector_to_interpolate_curve(degree, U, paras, pts);
	std::cout << "fixed " << std::endl;
	for (int i = 0; i < result_vector.size(); i++) {
		std::cout << result_vector[i] << std::endl << std::endl;
	}
}
double peak_function(const double x, const double y) {
	double r = 3 * pow(1 - x, 2)*exp(-x * x - (y + 1)*(y + 1)) - 10 * (0.2*x - pow(x, 3) - pow(y, 5))*
		exp(-x * x - y * y) - 1 / 3 * exp(-pow(x + 1, 2) - y * y);
	if (fabs(r) < SCALAR_ZERO) {
		r = 0;
	}
	return r;
}

void B_spline_surface_to_mesh(const Bsurface &surface, const int pnbr, Eigen::MatrixXd &ver, Eigen::MatrixXi& faces) {
	std::vector<std::vector<Vector3d>> pts;
	ver.resize(pnbr*pnbr, 3);
	int verline = 0;
	for (int i = 0; i < pnbr; i++) {
		for (int j = 0; j < pnbr; j++) {
			double upara = double(i) / (pnbr - 1);
			double vpara = double(j) / (pnbr - 1);
			ver.row(verline) = BSplineSurfacePoint(surface, upara, vpara);
			verline++;
		}
	}
	faces.resize(2 * (pnbr - 1)*(pnbr - 1), 3);
	int fline = 0;
	for (int i = 0; i < pnbr - 1; i++) {
		for (int j = 0; j < pnbr - 1; j++) {
			faces.row(fline) = Vector3i( i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
			faces.row(fline + 1) = Vector3i( i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
			fline += 2;
		}
	}
}

void parameter_grid_to_mesh(const Eigen::MatrixXd& uv, Eigen::MatrixXd &ver, Eigen::MatrixXi& edges) {
	assert(uv.cols() == 2);
	ver.resize(uv.rows() * 4, 3);
	edges.resize(uv.rows() * 2, 2);
	for (int i = 0; i < uv.rows(); i++) {
		double u = uv(i, 0);
		double v = uv(i, 1);
		ver.row(4 * i) = Vector3d(0, v, 0);
		ver.row(4 * i + 1) = Vector3d(1, v, 0);
		ver.row(4 * i + 2) = Vector3d(u, 0, 0);
		ver.row(4 * i + 3) = Vector3d(u, 1, 0);
	}
	for (int i = 0; i < uv.rows(); i++) {
		edges.row(2 * i) = Vector2i(4 * i, 4 * i + 1);
		edges.row(2 * i + 1) = Vector2i(4 * i + 2, 4 * i + 3);
	}

}
void knot_intervals_to_mesh(const int degree1, const int degree2, 
	const std::vector<double>& U, const std::vector<double>& V, 
	Eigen::MatrixXd &ver, Eigen::MatrixXi& edges) {
	int udiff = U.size() - 2 * degree1, vdiff = V.size() - 2 * degree2;
	Eigen::MatrixXd uv(udiff + vdiff, 2);
	for (int i = 0; i < uv.rows(); i++) {
		if (i < udiff) {
			uv.row(i) << U[i + degree1], 0;
		}
		else {
			uv.row(i) << 0, V[i - udiff + degree2];
		}
	}
	parameter_grid_to_mesh(uv, ver, edges);
}
void test_surface_visual(Eigen::MatrixXd &ver, Eigen::MatrixXi& faces) {
	Bsurface surface;
	surface.degree1 = 3;
	surface.degree2 = 3;
	surface.U = { {0,0,0,0,0.1,0.3,0.5,0.7,0.9,1,1,1,1} };// nu=13-2-3=8
	surface.V = { {0,0,0,0,0.1,0.2,0.3,0.5,0.7,0.8,0.9,1,1,1,1} };// nv=10
	int nu = surface.nu();
	int nv = surface.nv();
	std::vector<std::vector<Vector3d>> control(nu + 1);
	for (int i = 0; i < control.size(); i++) {
		control[i].resize(nv + 1);
		for (int j = 0; j < control[i].size(); j++) {
			double x = -3 + double(i) * 6 / nu;
			double y = -3 + double(j) * 6 / nv;
			double z = peak_function(x, y);
			control[i][j] = Vector3d(x, y, z);
		}
	}
	surface.control_points = control;

	B_spline_surface_to_mesh(surface, 50, ver, faces);
}

void test_surface_knot_preprocessing(Eigen::MatrixXd &points, Eigen::MatrixXd& knotP, Eigen::MatrixXi& knotE) {
	const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
	const std::string filename = path + "camel_smallest.obj";

	Eigen::MatrixXd V; Eigen::MatrixXi F;
	Eigen::MatrixXd  param;
	mesh_parameterization(filename, V, param, F);

	int degree1 = 3, degree2 = 3;
	std::vector<double> vecU = { {0,0,0,0,0.1,1,1,1,1} };
	std::vector<double> vecV = { {0,0,0,0,0.1,1,1,1,1} };
	std::vector<double> Uout, Vout;
	fix_surface_grid_parameter_too_many(degree1, vecU, degree2, vecV, param, Uout, Vout);
	std::cout << "fixed U is \n"; print_vector(Uout);
	std::cout << "fixed V is \n"; print_vector(Vout);
	points.resize(param.rows(), 3);
	for (int i = 0; i < points.rows(); i++) {
		points.row(i) << param.row(i), 0;
	}
	knot_intervals_to_mesh(degree1, degree2, Uout, Vout, knotP, knotE);
} 

// here the points are the parameters or the 3d positions
void test_knot_fixing(Eigen::MatrixXd &points, Eigen::MatrixXd& knotP, Eigen::MatrixXi& knotE) {
	const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
	const std::string filename = path + "camel_smallest.obj";
	// camel_small_open.obj is the problematic one
	Eigen::MatrixXd Ver; Eigen::MatrixXi F;
	Eigen::MatrixXd  param;
	mesh_parameterization(filename, Ver, param, F);
	points.resize(param.rows(), 3);
	for (int i = 0; i < points.rows(); i++) {
		points.row(i) << param.row(i), 0;
	}
	int degree1 = 3, degree2 = 3;
	std::vector<double> vecU = { {0,0,0,0,0.1,1,1,1,1} };
	std::vector<double> vecV = { {0,0,0,0,0.1,1,1,1,1} };
	std::vector<double> Uout, Vout;


	//knot_intervals_to_mesh(degree1, degree2, Uout, Vout, knotP, knotE);
	std::cout << "before fixing" << std::endl;
	std::cout << "para size " << param.rows() << std::endl;
	fix_knot_vector_to_interpolate_surface(degree1, degree2, vecU, vecV, param, Ver, Uout, Vout);
	//easist_way_to_fix_knot_vector_to_interpolate_surface(degree1, degree2, vecU, vecV, param, Ver, Uout, Vout);
	std::cout << "fixed U" << std::endl; print_vector(Uout);
	std::cout << "fixed V" << std::endl; print_vector(Vout);
	std::vector<int> list;
	for (int i = 0; i < param.rows(); i++) {
		list.push_back(i);
	}
	for (int i = 0; i < 3; i++) {
		bool solvable = selected_rows_have_solution_rational(degree1, degree2, Uout, Vout, param, Ver, list, i);
		std::cout << "solvable test, i=" << i << ", solvable= " << solvable << std::endl;
	}
	knot_intervals_to_mesh(degree1, degree2, Uout, Vout, knotP, knotE);
	//std::cout << "param\n" << param << std::endl;
}
Eigen::MatrixXi inverse_faces(const Eigen::MatrixXi& F) {
	Eigen::MatrixXi result(F.rows(), 3);
	for (int i = 0; i < F.rows(); i++) {
		result.row(i) << F(i, 0), F(i, 2), F(i, 1);
	}
	return result;
}

// 
int map_ij_to_list_id(const int i, const int j,const int isize, const int jsize) {
	
	return i * jsize + j;
}
void fit_border_with_curve(const Eigen::MatrixXd& paras, const Eigen::MatrixXd &ver, const Eigen::MatrixXi& map,
	const bool is_v, const int line_id, Bcurve& curve, const double a, const double b, std::vector<Vector3d> &points) {
	std::vector<double> paralist;
	
	if (is_v) {
		//curve.U = V;
		for (int i = 0; i < map.rows(); i++) {// iso-v
			
			if (map(i, line_id) >= 0) {
				double value = paras(map(i, line_id), 0);
				Vector3d vertex = ver.row(map(i, line_id));
				paralist.push_back(value);
				points.push_back(vertex);
			}
			
		}
	}
	else {
		//curve.U = U;
		
		for (int i = 0; i < map.cols(); i++) {// iso-u
			if (map(line_id, i) >= 0) {
				double value = paras(map(line_id, i), 1);
				paralist.push_back(value);
				Vector3d vertex = ver.row(map(line_id, i));
				points.push_back(vertex);
			}
		}
	}
	std::cout << "paras and points" << std::endl;
	print_vector(paralist);
	curve.degree = 3;
	std::vector<double> kv = { {0,0,0,0,0.1,0.2,0.4,0.5,0.6,0.7,0.8,1,1,1,1} };
	
	std::vector<double> kv_new = fix_knot_vector_to_interpolate_curve(curve.degree, kv, paralist, points);
	std::cout << "knot vector get fixed" << std::endl;

	

	curve.U = kv_new;
	print_vector(kv_new);
	solve_control_points_for_fairing_curve(curve, paralist, points, a, b);
}
void curve_visulization(const Bcurve&curve, const int nbr, Eigen::MatrixXd&e0, Eigen::MatrixXd &e1) {
	e0.resize(nbr, 3);
	e1.resize(nbr, 3);
	for (int i = 0; i < nbr; i++) {
		double para0 = double(i) / (nbr);
		double para1 = double(i + 1) / (nbr);
		e0.row(i) = BsplinePoint(curve.degree, curve.U, para0, curve.control_points);
		e1.row(i) = BsplinePoint(curve.degree, curve.U, para1, curve.control_points);
	}
	return;
}

void surface_visulization(Bsurface& surface, const int nbr, Eigen::MatrixXd & v, Eigen::MatrixXi &f) {
	B_spline_surface_to_mesh(surface, nbr, v, f);
	return;
}

// if no boundary is specified, this code will project vertices on x-y plane, and find the convex hull of the 2d vertices,
// then use the boundary of convex hull as boundary of u-v domain, then do parametrization in [0,1]x[0,1]
void find_boundar_and_parametrization(const Eigen::MatrixXd& ver, Eigen::MatrixXd &param, Eigen::MatrixXi &F) {
	
	Eigen::MatrixXd paras(ver.rows(), 2); //the 2d projection of the points 
	paras << ver.col(0), ver.col(1);
	Eigen::VectorXi loop;
	std::cout << "start find border" << std::endl;
	// find the convex hull on 2d,
	find_border_loop(paras, loop);
	if (0) {// if user is interested in this part
		Eigen::MatrixXd bdver(loop.size(), 3);
		for (int i = 0; i < loop.size(); i++) {
			bdver.row(i) = ver.row(loop[i]);
		}
		std::cout << "finish find border" << std::endl;
		Eigen::MatrixXi edges;
		vertices_to_edges(bdver, edges);
		
	}
	Eigen::MatrixXd boundary_uv;
	map_vertices_to_square(ver, loop, boundary_uv);
	std::cout << "boundary uv\n" << boundary_uv << std::endl;
	// triangulation
	constrained_delaunay_triangulation(paras, loop, F);
	igl::harmonic(ver, F, loop, boundary_uv, 1, param);// parametrization finished
	std::cout << "ver\n" << ver << std::endl;
	//std::cout << "paras\n" << paras << std::endl;
	std::cout << "loop\n" << loop << std::endl;
	std::cout << "F\n" << F << std::endl;
	
}

// direct project x y as u v, and rescale the scene into [0, 1]x[0, 1] 
void direct_project_x_y_and_parametrization(const Eigen::MatrixXd& ver, Eigen::MatrixXd &param, Eigen::MatrixXi &F) {

	Eigen::MatrixXd paras(ver.rows(), 2), para_new; //the 2d projection of the points 
	paras << ver.col(0), ver.col(1);
	para_new = paras;
	Eigen::VectorXi loop;
	std::cout << "start find border" << std::endl;
	// find the convex hull on 2d,
	find_border_loop(paras, loop);
	if (0) {// if user is interested in this part
		Eigen::MatrixXd bdver(loop.size(), 3);
		for (int i = 0; i < loop.size(); i++) {
			bdver.row(i) = ver.row(loop[i]);
		}
		std::cout << "finish find border" << std::endl;
		Eigen::MatrixXi edges;
		vertices_to_edges(bdver, edges);

	}
	
	double xmin = paras(0, 0);
	double xmax = paras(0, 0);
	double ymin = paras(0, 1);
	double ymax = paras(0, 1);
	for (int i = 0; i < paras.rows(); i++) {
		double x = paras(i, 0);
		double y = paras(i, 1);
		if (xmin > x) {
			xmin = x;
		}
		if (xmax < x) {
			xmax = x;
		}
		if (ymin > y) {
			ymin = y;
		}
		if (ymax < y) {
			ymax = y;
		}
	}
	std::cout << "bounding box: (" << xmin << ", " << ymin << "), (" << xmax << ", " << ymax << ")" << std::endl;
	
	double xsize = xmax - xmin;
	double ysize = ymax - ymin;
	for (int i = 0; i < paras.rows(); i++) {
		double x = paras(i, 0);
		double y = paras(i, 1);
		x = x - xmin;
		y = y - ymin;
		x = x / xsize;
		y = y / ysize;
		para_new(i, 0) = x;
		para_new(i, 1) = y;
	}
	param = para_new;
	for (int i = 0; i < param.rows(); i++) {
		assert(param(i, 0) <= 1 && param(i, 0) >= 0 && param(i, 1) <= 1 && param(i, 1) >= 0);
	}
	// triangulation
	constrained_delaunay_triangulation(para_new, loop, F);
	std::cout << "check directly parametrization\n" << param << std::endl;
}
Eigen::MatrixXd get_peak_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = 6 * para3d[0];
		double y = 6 * para3d[1];
		ver.row(i - skip) << x, y, peak_function(x, y);
	}
	ver.row(nbr - skip - 4) << -6, -6, peak_function(-6, -6);
	ver.row(nbr - skip - 3) << -6, 6, peak_function(-6, 6);
	ver.row(nbr - skip - 2) << 6, 6, peak_function(6, 6);
	ver.row(nbr - skip - 1) << 6, -6, peak_function(6, -6);
	return ver;
}
double contour_function(const double x, const double y) {
	return 4 * x*exp(-x * x - y * y);
}
Eigen::MatrixXd get_contour_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	double s = 2;
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = s * para3d[0];
		double y = s * para3d[1];
		ver.row(i - skip) << x, y, contour_function(x, y);
	}
	ver.row(nbr - skip - 4) << -s, -s, contour_function(-s, -s);
	ver.row(nbr - skip - 3) << -s, s, contour_function(-s, s);
	ver.row(nbr - skip - 2) << s, s, contour_function(s, s);
	ver.row(nbr - skip - 1) << s, -s, contour_function(s, -s);
	return ver;
}
double hyperbolic_function(const double x, const double y) {
	return x *x - y * y;
}
Eigen::MatrixXd get_hyperbolic_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	double s = 1;
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = s * para3d[0];
		double y = s * para3d[1];
		ver.row(i - skip) << x, y, hyperbolic_function(x, y);
	}
	ver.row(nbr - skip - 4) << -s, -s, hyperbolic_function(-s, -s);
	ver.row(nbr - skip - 3) << -s, s, hyperbolic_function(-s, s);
	ver.row(nbr - skip - 2) << s, s, hyperbolic_function(s, s);
	ver.row(nbr - skip - 1) << s, -s, hyperbolic_function(s, -s);
	return ver;
}
double sinus_function(const double x, const double y) {
	return sin(3 * 3.1415926*(x*x + y * y)) / 10;
}
Eigen::MatrixXd get_sinus_sample_points(const int nbr, const int skip) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	double s = 1;
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double x = s * para3d[0];
		double y = s * para3d[1];
		ver.row(i - skip) << x, y, sinus_function(x, y);
	}
	ver.row(nbr - skip - 4) << -s, -s, sinus_function(-s, -s);
	ver.row(nbr - skip - 3) << -s, s, sinus_function(-s, s);
	ver.row(nbr - skip - 2) << s, s, sinus_function(s, s);
	ver.row(nbr - skip - 1) << s, -s, sinus_function(s, -s);
	return ver;
}

Vector3d bilinear_function(const Vector3d& v0s, const Vector3d&v0e, const Vector3d& v1s, const Vector3d& v1e,
	const double u, const double v) {
	Vector3d m0 = (v0e - v0s)*u + v0s;
	Vector3d m1 = (v1e - v1s)*u + v1s;
	Vector3d m = (m1 - m0)*v + m0;
	return m;
}
Eigen::MatrixXd get_bilinear_sample_points(const int nbr, const int skip, Eigen::MatrixXd& param) {
	Eigen::MatrixXd ver;
	ver.resize(nbr - skip, 3);
	param.resize(nbr - skip, 2);
	double s = 1;// domain scale is [-s, s]x[-s, s]
	Vector3d v0s(0, 0, 1);
	Vector3d v0e(1, 1, 0);
	Vector3d v1s(0, 1, 1);
	Vector3d v1e(0, 0, 0);
	for (int i = 0; i < nbr - 4; i++) {
		Vector3d para3d = Vector3d::Random();
		if (i < skip) {
			continue;
		}
		double u = (1 + para3d[0]) / 2;
		double v = (1 + para3d[1]) / 2;
		param.row(i - skip) << u, v;
		ver.row(i - skip) = bilinear_function(v0s, v0e, v1s, v1e, u, v);
	}

	param.row(nbr - skip - 4) << 0, 0;
	param.row(nbr - skip - 3) << 0, 1;
	param.row(nbr - skip - 2) << 1, 1;
	param.row(nbr - skip - 1) << 1, 0;
	
	ver.row(nbr - skip - 4) = bilinear_function(v0s, v0e, v1s, v1e, 0, 0);
	ver.row(nbr - skip - 3) = bilinear_function(v0s, v0e, v1s, v1e, 0, 1);
	ver.row(nbr - skip - 2) = bilinear_function(v0s, v0e, v1s, v1e, 1, 1);
	ver.row(nbr - skip - 1) = bilinear_function(v0s, v0e, v1s, v1e, 1, 0);
	
	return ver;
}


void get_mesh_vertices_and_parametrization(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
	Eigen::MatrixXd  &param) {
	const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
	const std::string filename = path + "camel_small_open.obj";
	mesh_parameterization(filename, V, param, F);
}

// method = 0, peak; 
// method = 1, contour; 
// method = 2, hyperbolic; 
// method = 3, sinus;
// method = 4, bilinear
void get_function_vertices_and_parametrization(const int nbr, const int skip, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
	Eigen::MatrixXd  &param, const int method) {
	//V = get_peak_sample_points(nbr, skip);
	bool need_param = true;
	switch (method) {
	case 0:
		srand(15);
		V = get_peak_sample_points(nbr, skip);
	case 1:
		srand(7);
		V = get_contour_sample_points(nbr, skip);
	case 2:
		srand(19);
		V = get_hyperbolic_sample_points(nbr, skip);
	case 3:
		srand(24);
		V = get_sinus_sample_points(nbr, skip);
	case 4:
		srand(5);
		
		V = get_bilinear_sample_points(nbr, skip, param);
		Eigen::VectorXi loop;
		
		find_border_loop(param, loop);
		
		constrained_delaunay_triangulation(param, loop, F);
		
		return;
	}


	if (need_param) {
		direct_project_x_y_and_parametrization(V, param, F);
	}
	
}

void make_peak_exmple() {
	Eigen::MatrixXd ver;
	int nbr = 200;// nbr of points
	int skip = 0;
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, param_perturbed;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	
	get_function_vertices_and_parametrization(nbr, skip, ver, F, param, 4);
	
	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int perturb_itr = 5;
	double per_ours = 0.2;
	double per = 0.1;
	int target_steps = 10; 
	bool enable_max_fix_nbr = true;

	mesh_parameter_perturbation(param, F, param_perturbed, perturb_itr);
	std::cout << "param_perturbed\n" << param_perturbed << std::endl;
	igl::opengl::glfw::Viewer viewer;
	/*viewer.data().set_mesh(param_perturbed, F);
	viewer.launch();*/
	generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, param_perturbed, F, perturb_itr,per_ours,per, target_steps, enable_max_fix_nbr);
	

	Eigen::MatrixXd paramout(ver.rows(), 3), zeros(ver.rows(),1);
	zeros = Eigen::MatrixXd::Constant(ver.rows(), 1, 0);
	
	paramout << param_perturbed, zeros;

	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout<<"before solving control points"<<std::endl;
		solve_control_points_for_fairing_surface(surface, param_perturbed, ver);
		
		surface_visulization(surface, 100, SPs, SFs);
	}

	bool write_file = true;
	if (write_file)
	{
		const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
		igl::write_triangle_mesh(path + "0517_missing.obj", SPs, SFs);
		//igl::write_triangle_mesh(path + "0517_missing_param.obj", paramout, F);
	}
	output_timing();
	

	
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;

	/*
	viewer.data().set_edges(bdver, edges, fcolor);*/
	//viewer.data().add_points(vector_to_matrix_3d(inter_pts), ecolor);
	
	// see the linear interpolated surface
	//viewer.data().set_mesh(param, F);
	//viewer.data().add_edges(edge0, edge1, pcolor);

	//viewer.data().set_mesh(param_perturbed, F);
	viewer.data().clear();
	viewer.data().set_mesh(SPs, SFs);
	viewer.data().add_points(ver, ecolor);
	viewer.launch();
}