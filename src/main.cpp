#include"basis.h"
#include"curve.h"
#include <igl/opengl/glfw/Viewer.h>
#include<iostream>
#include<array>
#include"test.h"
#include<mesh_processing.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/write_triangle_mesh.h>
#include<exact_calculation.h>
#include<energy.h>
#include<integration.h>
igl::opengl::glfw::Viewer global_viewer;
void show_basis(const std::vector<double>&b_vec) {
	for (int i = 0; i < b_vec.size(); i++) {
		std::cout << b_vec[i] << " ";
	}
	std::cout << std::endl;
}
void test_opengl() {
	// Inline mesh of a cube
	const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 1.0, 0.0,
		0.0, 1.0, 1.0,
		1.0, 0.0, 0.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 0.0,
		1.0, 1.0, 1.0).finished();
	const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) <<
		1, 7, 5,
		1, 3, 7,
		1, 4, 3,
		1, 2, 4,
		3, 8, 7,
		3, 4, 8,
		5, 7, 8,
		5, 8, 6,
		1, 5, 6,
		1, 6, 2,
		2, 6, 8,
		2, 8, 4).finished().array() - 1;

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.launch();
}


void draw_a_line() {
	Eigen::MatrixXd pts(3,3),Color(3,3);
	//
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Color(i, j) = 0;
			pts(i, j) = 0;
		}
	}
	pts(0, 0) = 1; pts(0, 1) = 2; pts(0, 2) = 3;
	std::cout << pts << std::endl;
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_points(pts,Color);
	viewer.launch();
}
void draw_a_curve() {
	int p = 3;
	// 8 control points
	std::vector<double>U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };
	std::array<double, 8> us = { {0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9} };
	std::vector<Vector3d> pts(8);
	pts[0] = Vector3d(0, 0, 0);
	pts[1] = Vector3d(0, 1, 0);
	pts[2] = Vector3d(0, 1, 2);
	pts[3] = Vector3d(0, 1, 1);
	pts[4] = Vector3d(0, 2, 1);
	pts[5] = Vector3d(0, 2, 0);
	pts[6] = Vector3d(0, 3, 1);
	pts[7] = Vector3d(0, 3, 3);
	igl::opengl::glfw::Viewer viewer;
	Eigen::MatrixXd  Color(8, 3), Hcolor(100, 3), OneColor(1, 3), PTS(8, 3), Curve(100, 3);
	OneColor(0, 0) = 0; OneColor(0, 1) = 0; OneColor(0, 2) = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 3; j++) {
			Color(i, j) = 0;
			
		}
		PTS.row(i) = pts[i];
		
	}
	for (int i = 0; i < 100; i++) {
		Curve.row(i) = BsplinePoint(3, U, i / double(100), pts);
		Hcolor(i, 0) = 0.5; Hcolor(i, 1) = 0.5; Hcolor(i, 2) = 0.5;
	}
	std::cout << Curve << std::endl;
	viewer.data().set_points(PTS, Color);
	viewer.data().add_points(Curve, Hcolor);
	viewer.launch();
	//for (int i = 0; i < 8; i++) {
	//	
	//}
}
void plot_fitting_result() {
	Eigen::MatrixXd control_pts;
	Eigen::MatrixXd control_pts_color;
	Eigen::MatrixXd curve_pts; Eigen::MatrixXd curve_pts_color;
	Eigen::MatrixXd target_pts; Eigen::MatrixXd target_pts_color;
	visual_curve_fitting(control_pts, control_pts_color,
			curve_pts, curve_pts_color, target_pts, target_pts_color);
	Eigen::MatrixXi edges;
	vertices_to_edges(curve_pts, edges);
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_points(control_pts, control_pts_color);
	//viewer.data().add_points(curve_pts, curve_pts_color);
	viewer.data().add_points(target_pts, target_pts_color);
	std::cout << "target pts nbr, " << target_pts.rows() << std::endl;
	viewer.data().set_edges(curve_pts, edges, curve_pts_color.row(0));
	viewer.launch();
}

// scale is the length of the axis arrows
void draw_axis(const double scale) {
	Eigen::MatrixXd red(1, 3), green(1, 3), blue(1, 3);
	Eigen::MatrixXd origin(1, 3), xaxis(1, 3), yaxis(1, 3),zaxis(1, 3);
	origin.row(0) = Vector3d(0, 0, 0);
	red.row(0) = Vector3d(1, 0, 0);
	green.row(0) = Vector3d(0, 1, 0);
	blue.row(0) = Vector3d(0, 0, 1);
	Eigen::MatrixXd V1,V2, C;
	V1.resize(3, 3); V2.resize(3, 3);
	V1 << origin, origin, origin;
	V2 << red, green, blue;

	C.resize(3, 3);
	C << red, green, blue;
	V1 = V1 * scale;
	V2 = V2 * scale;
	//Eigen::MatrixXi E(3,2);
	//E.row(0) << 0, 1;
	//E.row(1) << 0, 2;
	//E.row(2) << 0, 3;
	global_viewer.data().add_edges(V1, V2, C);
	global_viewer.data().line_width = 1;
}





void visual_and_chop_mesh(const bool write_mesh) {
	const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\simple\\";
	const std::string filename = path + "blendspheres.obj";
	Eigen::MatrixXd V, Vclean; Eigen::MatrixXi F;
	read_and_visual_mesh(filename, V, F);
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	fcolor.row(0) = Vector3d(0, 0.5, 0.5); ecolor.row(0) = Vector3d(0, 0, 0);

	///////////////////////////////////////////
	//// parameterization part
	//Eigen::VectorXi bnd;
	//igl::boundary_loop(F, bnd);// boundary vertices detection
	//Eigen::MatrixXd bnd_uv, param;
	///*igl::map_vertices_to_circle(V, bnd, bnd_uv);*/
	//map_vertices_to_square(V, bnd, bnd_uv);
	////exit(0);
	//igl::harmonic(V, F, bnd, bnd_uv, 1, param);
	///////////////////////////////////////////////
	bool remove_faces = false;
	Eigen::MatrixXi newF,Fclean;
	if (remove_faces){
		remove_some_faces(2, 38, false, V, F, newF);
		remove_some_faces(1, 32, false, V, newF, newF);
	}
	else {
		newF = F;
	}
	
	generate_clean_mesh_data_for_parametrization(V, newF, Vclean, Fclean);
	if (write_mesh)
		igl::write_triangle_mesh(path + "camel_smallest.obj", Vclean, Fclean);
	Eigen::VectorXi bnd;
	igl::boundary_loop(Fclean, bnd);
	Eigen::MatrixXd visual_bnd(bnd.size(), 3);
	for (int i = 0; i < bnd.size(); i++) {
		visual_bnd.row(i) = Vclean.row(bnd[i]);
	}
	//std::cout << "boundary ids,\n" << bnd << std::endl;
	global_viewer.data().set_mesh(Vclean, Fclean);
	global_viewer.data().add_points(visual_bnd, fcolor);
	draw_axis(10);
	global_viewer.launch();


}


void test_mesh_parameterization() {
	const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
	const std::string filename = path + "camel_smallest.obj";

	Eigen::MatrixXd V; Eigen::MatrixXi F;
	Eigen::MatrixXd  param, para_perturbed;
	mesh_parameterization(filename, V, param, F);
	
	mesh_parameter_perturbation(param, F, para_perturbed, 2);
	/////////////////////////////////////////////
	Eigen::MatrixXd grid_ver; Eigen::MatrixXi grid_edges;
	parameter_grid_to_mesh(param, grid_ver, grid_edges);

	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	fcolor.resize(1, 3);
	fcolor << 1, 0, 0;

	bool see_girds = true;
	bool see_triangles = false;
	if (see_girds) {
		global_viewer.data().set_edges(grid_ver, grid_edges, fcolor);
	}
	//
	std::cout << "original paras\n" << param << std::endl;
	std::cout << "perturbed paras\n" << para_perturbed << std::endl;
	if (see_triangles) {
		global_viewer.data().set_mesh(para_perturbed, F);
	}
	//
	draw_axis(10);

	
	global_viewer.launch();


}
void visual_surface() {
	Eigen::MatrixXd ver;
	Eigen::MatrixXi faces;
	test_surface_visual(ver, faces);
	global_viewer.data().set_mesh(ver, faces);
	global_viewer.launch();
}
void visual_surface_processing() {
	Eigen::MatrixXd points; 
	Eigen::MatrixXd knotP; 
	Eigen::MatrixXi knotE;
	test_surface_knot_preprocessing(points, knotP, knotE);
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.5, 0.5, 0.5;
	global_viewer.data().set_edges(knotP, knotE, fcolor);
	global_viewer.data().add_points(points, ecolor);
	draw_axis(2);
	global_viewer.launch();
}
void visual_surface_knot_fixing() {
	Eigen::MatrixXd points;
	Eigen::MatrixXd knotP;
	Eigen::MatrixXi knotE;
	test_knot_fixing(points, knotP, knotE);
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.5, 0.5, 0.5;
	global_viewer.data().set_edges(knotP, knotE, fcolor);
	global_viewer.data().add_points(points, ecolor);
	draw_axis(2);
	global_viewer.launch();
}

//void test_rational() {
//	std::vector<double> U = { {0,0,0,0,0.5,1,1,1,1} };
//	Rational re = Nip_Rational(1, 3, 0.3, U);
//	std::cout << "rational number is, " << re.to_double() << std::endl;
//	Rational a = 5;
//	Rational b = -a;
//	std::cout << "-a= " << (-a).to_double()<<" "<<
//		b.to_double() << std::endl;
//	b = 1;
//	b /= a;
//	std::cout << "b " << b <<" a "<<a<< std::endl<<std::endl;
//
//	MatrixXs mt;
//	mt.resize(3, 2);
//	mt << 1, 1.5,
//		2, 2,
//		3, 1;
//	rank(mt);
//}

void test_poly() {
	std::vector<double> poly1 = { {0,1,2,0} };
	std::vector<double> poly2 = { {3,2} };
	std::vector<double> result = polynomial_times(poly1, poly2);

	/*std::vector<double> U={ {0, 0, 0, 0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 1, 1, 1, 1} };
	double u = 0.88;
	std::cout << "Nip " << Nip(7, 3, u, U) << std::endl;
	std::cout << "Nip " << polynomial_value(Nip_func(7, 3, u, U),u) << std::endl;
	std::cout << "poly operation result " << std::endl;
	print_vector(result);*/
	std::cout << "poly intergration, " << std::endl;
	print_vector(polynomial_integration(poly1));
}
void test_orientation() {
	Vector2d p(0, 0), q(0, 1), r(1, 0);
	std::cout << "ori " << orient_2d(p, q, r) << std::endl;
}
void test_interpolation() {
	int Fid = 0;
	double u = 1;
	double v = 1;
	Eigen::MatrixXd vertices(3, 3);
	Vector3d v0(0, 0, 0), v1(4, 0, 2), v2(0, 3, 1);
	vertices.row(0) = v0;
	vertices.row(1) = v1;
	vertices.row(2) = v2;
	Eigen::MatrixXd param(3, 2);
	Vector2d pa0(0, 0), pa1(4, 0), pa2(0, 3);
	
	param.row(0) = pa0;
	param.row(1) = pa1;
	param.row(2) = pa2;
	std::cout << "param \n" << param << std::endl;
	Eigen::MatrixXi F(1,3);
	F << 0, 1, 2;
	Vector3d point =linear_interpolation( Fid,  u,  v,
		 vertices, param,
		F);
	std::cout << "point\n" << point << std::endl;
}

void test1() {
	//int degree = 3;
	//std::vector<double> init_vec = { {0,0,0,0,0.1,0.13,0.15,0.2,1,1,1,1} };
	//std::vector<double> paras = { {0,0.1, 0.2, 0.3, 0.4 ,0.5 , 0.6, 0.7, 0.8} };
	//fix_knot_vector_to_interpolate_curve_WKW(degree, init_vec, paras);
}

// method = 0, peak; 
// method = 1, contour; 
// method = 2, hyperbolic; 
// method = 3, sinus;
// method = 4, bilinear

int main() {
	const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
	int model = 4;
	int nbr = 50;
	double par = 0.5;// ours
	double per = 0.9;
	std::string tail = "";
	run_ours(model, nbr, par, path, tail, per);
	//test_opengl();
	//int p = 3;
	//std::vector<double>U = { {0,0,0,0,0.1,0.4,0.7,0.9,1,1,1,1} };
	//draw_a_line();
	//draw_a_curve();
	//test_fitting();
	//test1();
//plot_fitting_result();
	//test_curve_knot_fixing();
	//visual_mesh();
	//test_mesh_parameterization();
	//visual_and_chop_mesh(false);
	//visual_surface();
	//visual_surface_processing();
	//visual_surface_knot_fixing();
	//test_rational();
	//test_poly();
	//go_through_temperature_interpolation();
//make_peak_exmple();
	//run_Seungyong();
	return 0;
}