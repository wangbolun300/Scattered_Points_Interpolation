#include"test.h"
#include"curve.h"
#include<iostream>
#include<surface.h>
#include<cmath>
#include<mesh_processing.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/harmonic.h>
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
	std::vector<double> U_init = { {0,0,0,0,0.1,1,1,1,1} };
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
			faces.row(fline) = Vector3i(i + pnbr * j, i + pnbr * (j + 1), i + pnbr * (1 + j) + 1);
			faces.row(fline + 1) = Vector3i(i + pnbr * j, i + pnbr * (1 + j) + 1, i + pnbr * j + 1);
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

void make_peak_exmple() {
	Eigen::MatrixXd ver; Eigen::MatrixXi faces;
	int nbr = 100;// nbr of points
	ver.resize(nbr, 3);
	for (int i = 0; i < nbr; i++) {
		Vector3d para3d = Vector3d::Random();
		double x = -3 + 6 * para3d[0];
		double y = -3 + 6 * para3d[1];
		ver.row(i) << x, y, peak_function(x, y);
	}
	
	Eigen::MatrixXd paras(ver.rows(), 2); //the parameters of the points 
	paras << ver.col(0), ver.col(1);
	Eigen::VectorXi loop;

	std::cout << "start find border" << std::endl;

	find_border_loop(paras, loop);
	Eigen::MatrixXd bdver(loop.size(), 3);
	for (int i = 0; i < loop.size(); i++) {
		bdver.row(i) = ver.row(loop[i]);
	}
	std::cout << "finish find border" << std::endl;
	Eigen::MatrixXi edges;
	vertices_to_edges(bdver, edges);
	Eigen::MatrixXi F;
	constrained_delaunay_triangulation(paras, loop, F);
	Eigen::MatrixXd boundary_uv;
	map_vertices_to_square(ver, loop, boundary_uv);
	Eigen::MatrixXd param;
	igl::harmonic(ver,F,loop,boundary_uv,1,param);

	igl::opengl::glfw::Viewer viewer;
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;
	/*viewer.data().add_points(ver, ecolor);
	viewer.data().add_points(bdver, fcolor);
	viewer.data().set_edges(bdver, edges, fcolor);*/
	viewer.data().set_mesh(param, F);
	viewer.launch();
}