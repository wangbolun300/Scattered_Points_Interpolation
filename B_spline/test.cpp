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
#include<igl/Timer.h>

double peak_function(const double x, const double y) {
	double r = 3 * pow(1 - x, 2)*exp(-x * x - (y + 1)*(y + 1)) - 10 * (0.2*x - pow(x, 3) - pow(y, 5))*
		exp(-x * x - y * y) - 1 / 3 * exp(-pow(x + 1, 2) - y * y);
	if (fabs(r) < SCALAR_ZERO) {
		r = 0;
	}
	r /= 2;
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
void surface_visulization(std::vector<Bsurface>& surfaces, const int pnbr, Eigen::MatrixXd & ver, 
	Eigen::MatrixXi &faces, const int without) {
	std::vector<std::vector<Vector3d>> pts;
	ver.resize(pnbr*pnbr, 3);
	int verline = 0;
	for (int i = 0; i < pnbr; i++) {
		for (int j = 0; j < pnbr; j++) {
			double upara = double(i) / (pnbr - 1);
			double vpara = double(j) / (pnbr - 1);
			ver.row(verline) = Vector3d(0, 0, 0);
			for (int si = 0; si < surfaces.size()- without; si++) {
				ver.row(verline) += BSplineSurfacePoint(surfaces[si], upara, vpara);
			}
			
			verline++;
		}
	}
	faces.resize(2 * (pnbr - 1)*(pnbr - 1), 3);
	int fline = 0;
	for (int i = 0; i < pnbr - 1; i++) {
		for (int j = 0; j < pnbr - 1; j++) {
			faces.row(fline) = Vector3i(i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
			faces.row(fline + 1) = Vector3i(i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
			fline += 2;
		}
	}
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
	//std::cout << "check directly parametrization\n" << param << std::endl;
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
// u0 v0 are in [0,1].
// the function u in [1,2], v in [-pi/2, pi/2]
Vector3d snail_function(double u0, double v0) {
	const double pi = 3.1415926;
	double u = 1+ u0;
	double v = -pi/2 + pi * v0;
	double x = u * cos(v)*sin(u);
	double y = u * cos(u)*cos(v);
	double z = -u * sin(v);
	return Vector3d(x, y, z);
}

//#define B_SPLINE_WITH_NOISE
double noise_scale = 0.05;
// method = 0, peak; 
// method = 1, contour; 
// method = 2, hyperbolic; 
// method = 3, sinus;
// method = 4, bilinear
void input_parameters_get_model_sample_points( Eigen::MatrixXd &V, Eigen::MatrixXi &F,
	const Eigen::MatrixXd  &param, const int method) {
	const int nbr_models = 6;

	const auto mfunction_value = [](const int method, const double x, const double y) {
		Vector3d v0s(0, 0, 1);
		Vector3d v0e(1, 1, 0);
		Vector3d v1s(0, 1, 1);
		Vector3d v1e(0, 0, 0);
		switch (method) {
		case 0:
			return Vector3d(x, y, peak_function(x, y));
		case 1:
#ifdef B_SPLINE_WITH_NOISE
		{
			Vector3d norm(4 * exp(-x * x - y * y) - 8 * x*x*exp(-x * x - y * y), -8 * x*y*exp(-x * x - y * y), -1);
			norm = norm.normalized();
			Vector3d noise = norm * noise_scale*Vector3d::Random()[0];
			return Vector3d(x + noise[0], y + noise[1], contour_function(x, y) + noise[2]);
		}
#else
			return Vector3d(x, y, contour_function(x, y));
#endif
		case 2:
#ifdef B_SPLINE_WITH_NOISE
		{
			Vector3d norm(2 * x, -2 * y, -1);
			norm = norm.normalized();
			Vector3d noise = norm * noise_scale*Vector3d::Random()[0];
			return Vector3d(x + noise[0], y + noise[1], hyperbolic_function(x, y) + noise[2]);
		}
#else
			return Vector3d(x, y, hyperbolic_function(x, y));
#endif
		case 3:
			return Vector3d(x, y, sinus_function(x, y));
		case 4:
			return bilinear_function(v0s, v0e, v1s, v1e, x, y);
		case 5:
			return snail_function(x, y);
		}
	};
	std::vector<double> scale= { {3,2,1,1,1} };
	std::vector<int> seed = { {3,7,33,10,5,0} };
	Eigen::MatrixXd ver;
	ver.resize(param.rows(), 3);
	
	double s = scale[method];// domain scale is [-s, s]x[-s, s]
	int end = param.rows();
	srand(seed[method]);
	bool need_param = true;
	if (method == 4 || method == 5) {
		for (int i = 0; i < end; i++) {
			ver.row(i) = mfunction_value(method, param(i,0), param(i,1));
		}
		Eigen::VectorXi loop;

		find_border_loop(param, loop);
		constrained_delaunay_triangulation(param, loop, F);
		
	}
	else {
		for (int i = 0; i < end; i++) {
			double x = s * (param(i, 0) * 2 - 1);
			double y = s * (param(i, 1) * 2 - 1);
			ver.row(i) = mfunction_value(method, x, y);
		}
	}
	V = ver;
	if (1) {
		double xmin = INFINITE, ymin = INFINITE, zmin = INFINITE;
		double xmax = -INFINITE, ymax = -INFINITE, zmax = -INFINITE;
		int pnbr = 100;
		Eigen::MatrixXi faces;
		std::vector<std::vector<Vector3d>> pts;
		ver.resize(pnbr*pnbr, 3);
		int verline = 0;
		for (int i = 0; i < pnbr; i++) {
			for (int j = 0; j < pnbr; j++) {
				double upara = double(i) / (pnbr - 1);
				double vpara = double(j) / (pnbr - 1);
				double u, v;
				if (method == 4 || method == 5) {
					u = upara;
					v = vpara;
				}
				else {
					u = 2 * s * upara - s;
					v = 2 * s * vpara - s;
				}
				ver.row(verline) = mfunction_value(method, u, v);
				if (xmin > ver(verline, 0)) {
					xmin = ver(verline, 0);
				}
				if (ymin > ver(verline, 1)) {
					ymin = ver(verline, 1);
				}
				if (zmin > ver(verline, 2)) {
					zmin = ver(verline, 2);
				}
				if (xmax < ver(verline, 0)) {
					xmax = ver(verline, 0);
				}
				if (ymax < ver(verline, 1)) {
					ymax = ver(verline, 1);
				}
				if (zmax < ver(verline, 2)) {
					zmax = ver(verline, 2);
				}
				verline++;
			}
		}
		std::cout << "model bounding box \n" << xmin << ", " << xmax << "\n" << ymin << ", " << ymax << "\n" << zmin << ", " << zmax << std::endl;
		faces.resize(2 * (pnbr - 1)*(pnbr - 1), 3);
		int fline = 0;
		for (int i = 0; i < pnbr - 1; i++) {
			for (int j = 0; j < pnbr - 1; j++) {
				faces.row(fline) = Vector3i(i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
				faces.row(fline + 1) = Vector3i(i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
				fline += 2;
			}
		}
		igl::write_triangle_mesh("./model_" + std::to_string(method) + ".obj", ver, faces);
		std::cout << "write original model" << std::endl;
	}

}



void get_model_sample_points(const int nbr, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
	Eigen::MatrixXd  &param, const int method, bool corners,std::string path) {
	const int nbr_models = 6;

	const auto mfunction_value = [](const int method, const double x, const double y) {
		Vector3d v0s(0, 0, 1);
		Vector3d v0e(1, 1, 0);
		Vector3d v1s(0, 1, 1);
		Vector3d v1e(0, 0, 0);
		switch (method) {
		case 0:
			return Vector3d(x, y, peak_function(x, y));
		case 1:
#ifdef B_SPLINE_WITH_NOISE
		{
			Vector3d norm(4 * exp(-x * x - y * y) - 8 * x*x*exp(-x * x - y * y), -8 * x*y*exp(-x * x - y * y), -1);
			norm = norm.normalized();
			Vector3d noise = norm * noise_scale*Vector3d::Random()[0];
			return Vector3d(x+noise[0], y+noise[1], contour_function(x, y)+noise[2]);
		}
#else
			return Vector3d(x, y, contour_function(x, y));
#endif
		case 2:
#ifdef B_SPLINE_WITH_NOISE
		{
			Vector3d norm(2*x, -2*y, -1);
			norm = norm.normalized();
			Vector3d noise = norm * noise_scale*Vector3d::Random()[0];
			return Vector3d(x + noise[0], y + noise[1], hyperbolic_function(x, y) + noise[2]);
		}
#else
			return Vector3d(x, y, hyperbolic_function(x, y));
#endif
		case 3:
			return Vector3d(x, y, sinus_function(x, y));
		case 4:
			return bilinear_function(v0s, v0e, v1s, v1e, x, y);
		case 5:
			return snail_function(x, y);
		}
	};
	std::vector<double> scale = { {3,2,1,1,1} };
	std::vector<int> seed = { {3,7,33,10,5,0} };
	Eigen::MatrixXd ver;
	ver.resize(nbr, 3);
	param.resize(nbr, 2);
	double s = scale[method];// domain scale is [-s, s]x[-s, s]
	int end = corners ? nbr - 4 : nbr;
	srand(seed[method]);
	bool need_param = true;
	if (method == 4 || method == 5) {
		for (int i = 0; i < end; i++) {
			Vector3d para3d = Vector3d::Random();
			double u = (1 + para3d[0]) / 2;
			double v = (1 + para3d[1]) / 2;
			param.row(i) << u, v;
			ver.row(i) = mfunction_value(method, u, v);
		}
		if (corners) {
			param.row(nbr - 4) << 0, 0;
			param.row(nbr - 3) << 0, 1;
			param.row(nbr - 2) << 1, 1;
			param.row(nbr - 1) << 1, 0;

			ver.row(nbr - 4) = mfunction_value(method, 0, 0);
			ver.row(nbr - 3) = mfunction_value(method, 0, 1);
			ver.row(nbr - 2) = mfunction_value(method, 1, 1);
			ver.row(nbr - 1) = mfunction_value(method, 1, 0);
		}
		Eigen::VectorXi loop;

		find_border_loop(param, loop);
		constrained_delaunay_triangulation(param, loop, F);
		need_param = false;

	}
	else {
		for (int i = 0; i < end; i++) {
			Vector3d para3d = Vector3d::Random();
			double x = s * para3d[0];
			double y = s * para3d[1];
			ver.row(i) = mfunction_value(method, x, y);
		}
		if (corners) {
			ver.row(nbr - 4) =  mfunction_value(method, -s, -s);
			ver.row(nbr - 3) =  mfunction_value(method, -s, s);
			ver.row(nbr - 2) =  mfunction_value(method, s, s);
			ver.row(nbr - 1) = mfunction_value(method, s, -s);
		}

	}
	V = ver;
	if (need_param) {
		direct_project_x_y_and_parametrization(V, param, F);
	}
	if (1) {
		double xmin = INFINITE, ymin = INFINITE, zmin = INFINITE;
		double xmax = -INFINITE, ymax = -INFINITE, zmax = -INFINITE;
		int pnbr = 100;
		Eigen::MatrixXi faces;
		std::vector<std::vector<Vector3d>> pts;
		ver.resize(pnbr*pnbr, 3);
		int verline = 0;
		for (int i = 0; i < pnbr; i++) {
			for (int j = 0; j < pnbr; j++) {
				double upara = double(i) / (pnbr - 1);
				double vpara = double(j) / (pnbr - 1);
				double u, v;
				if (method == 4 || method == 5) {
					u = upara;
					v = vpara;
				}
				else {
					u = 2 * s * upara - s;
					v = 2 * s * vpara - s;
				}
				ver.row(verline) = mfunction_value(method, u, v);
				if (xmin > ver(verline, 0)) {
					xmin = ver(verline, 0);
				}
				if (ymin > ver(verline, 1)) {
					ymin = ver(verline, 1);
				}
				if (zmin > ver(verline, 2)) {
					zmin = ver(verline, 2);
				}
				if (xmax < ver(verline, 0)) {
					xmax = ver(verline, 0);
				}
				if (ymax < ver(verline, 1)) {
					ymax = ver(verline, 1);
				}
				if (zmax < ver(verline, 2)) {
					zmax = ver(verline, 2);
				}
				verline++;
			}
		}
		std::cout << "model bounding box \n" << xmin << ", " << xmax << "\n" << ymin << ", " << ymax << "\n" << zmin << ", " << zmax << std::endl;
		faces.resize(2 * (pnbr - 1)*(pnbr - 1), 3);
		int fline = 0;
		for (int i = 0; i < pnbr - 1; i++) {
			for (int j = 0; j < pnbr - 1; j++) {
				faces.row(fline) = Vector3i(i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
				faces.row(fline + 1) = Vector3i(i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
				fline += 2;
			}
		}
		igl::write_triangle_mesh(path + "model_" + std::to_string(method) + ".obj", ver, faces);
		std::cout << "write original model" << std::endl;
	}


}
// method = 0, peak; 
// method = 1, contour; 
// method = 2, hyperbolic; 
// method = 3, sinus;
// method = 4, bilinear
//void get_function_vertices_and_parametrization(const int nbr, const int skip, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
//	Eigen::MatrixXd  &param, const int method) {
//	//V = get_peak_sample_points(nbr, skip);
//	bool need_param = true;
//	switch (method) {
//	case 0:
//		srand(15);
//		V = get_peak_sample_points(nbr, skip);
//		break;
//	case 1:
//		srand(7);
//		V = get_contour_sample_points(nbr, skip);
//		break;
//	case 2:
//		srand(33);
//		V = get_hyperbolic_sample_points(nbr, skip);
//		break;
//	case 3:
//		srand(24);
//		V = get_sinus_sample_points(nbr, skip);
//		break;
//	case 4:
//		srand(5);
//		
//		V = get_bilinear_sample_points(nbr, skip, param);
//		Eigen::VectorXi loop;
//		
//		find_border_loop(param, loop);
//		
//		constrained_delaunay_triangulation(param, loop, F);
//		need_param = false;
//		break;
//		
//	}
//
//
//	if (need_param) {
//		direct_project_x_y_and_parametrization(V, param, F);
//	}
//	/*for (int i = 0; i < V.rows(); i++) {
//		for (int j = 0; j < 3; j++) {
//			if (fabs(V(i, j)) < SCALAR_ZERO) {
//				(V(i, j) = 0;
//			}
//		}
//	}*/
//	return;
//}
double perturbed_distance(const Eigen::MatrixXd&p1, const Eigen::MatrixXd&p2) {
	double dis = 0;
	for (int i = 0; i < p1.rows(); i++) {
		double new_dis = (p1.row(i) - p2.row(i)).norm();
		if (new_dis > dis) {
			dis = new_dis;
		}
	}
	return dis;
}


void write_points(const std::string& file, const Eigen::MatrixXd& ver) {
	std::ofstream fout;
	fout.open(file);
	for (int i = 0; i < ver.rows(); i++) {
		fout 
			//<< std::setprecision(17) 
			<< "v " << ver(i, 0) << " " << ver(i, 1) << " " << ver(i, 2) << std::endl;
	}
	fout.close();
}
void write_control_pts(std::vector<std::vector<Vector3d>>& cps, std::string file) {
	std::ofstream fout;
	fout.open(file);
	for (int i = 0; i < cps.size(); i++) {
		for (int j = 0; j < cps[i].size(); j++) {
			fout << "v " << cps[i][j][0] << " " << cps[i][j][1] << " " << cps[i][j][2] << std::endl;
		}
	}
	fout.close();
}
void write_csv(const std::string& file, const std::vector<std::string> titles,const std::vector<double> data) {
	std::ofstream fout;
	fout.open(file);
	for (int i = 0; i < titles.size()-1; i++) {
		fout << titles[i] << ",";
	}
	fout << titles.back() << std::endl;
	for (int i = 0; i < data.size() - 1; i++) {
		fout << data[i] << ",";
	}
	fout << data.back() << std::endl;
	fout.close();
}
void write_svg_pts(const std::string& file, const Eigen::MatrixXd &param) {
	std::ofstream fout;
	fout.open(file);
	double scale = 1000;
	double dot = 10;
	double displace = 30;
	std::string color = "#ED6E46";
	std::string black = "#000000";
	fout << "<svg>" << std::endl;
	for (int i = 0; i < param.rows(); i++) {
		double x = displace + param(i, 0)*scale;
		double y = displace + param(i, 1)*scale;
		fout << " <circle cx=\"" << x << "\" cy = \"" << y << "\" r = \"" << dot << "\" fill = \"" << color << "\" />" << std::endl;
	}
	fout << "<polyline points=\"";
	
	std::vector<double> xlist = { {0,0,1,1} };
	std::vector<double> ylist = { {0,1,1,0} };

	for (int i = 0; i < 4; i++) {
		fout << displace + xlist[i] * scale << "," << displace + ylist[i] * scale << " ";
	}
	fout << displace + xlist[0] * scale << "," << displace + ylist[0] * scale;
	fout << "\" fill = \"white\" stroke = \"" + black + "\" stroke-width = \"" + std::to_string(5) + "\" /> " << std::endl;
	/*<svg>
		<rect width = "200" height = "100" fill = "#BBC42A" / >
		< / svg>*/
	fout << "</svg>" << std::endl;
	fout.close();
}
void write_svg_knot_vectors(const std::string& file, const std::vector<double> &U, 
	const std::vector<double>&V) {
	std::ofstream fout;
	fout.open(file);
	double scale = 1000;
	double width1 = 5;
	double width2 = 3;
	double displace = 30;
	std::string color = "#000000";
	std::string color1 = "#7AA20D";
	fout << "<svg>" << std::endl;
	std::vector<double> xlist = { {0,0,1,1} };
	std::vector<double> ylist = { {0,1,1,0} };
	fout << "<polyline points=\"";
	for (int i = 0; i < 4; i++) {
		fout << displace+xlist[i]*scale << "," << displace+ylist[i]*scale << " ";
	}
	fout << displace + xlist[0] * scale << "," << displace + ylist[0] * scale;
	fout << "\" fill = \"white\" stroke = \"" + color + "\" stroke-width = \"" + std::to_string(5) + "\" /> " << std::endl;

	for (int i = 0; i < U.size(); i++) {
		double x1 = displace + U[i]*scale;
		double y1 = displace + 0*scale;
		double x2 = displace + U[i] * scale;
		double y2 = displace + 1 * scale;
		fout << "<line x1=\"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) +
			"\" x2=\"" + std::to_string(x2) + "\" y2=\"" + std::to_string(y2) +
			"\" stroke = \"" + color1 + "\" stroke-width = \"" 
			+std::to_string(width2)+"\" /> "<< std::endl;
	}
	for (int i = 0; i < V.size(); i++) {
		double x1 = displace + 0 * scale;
		double y1 = displace + V[i] * scale;
		double x2 = displace + 1 * scale;
		double y2 = displace + V[i] * scale;
		fout << "<line x1=\"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) +
			"\" x2=\"" + std::to_string(x2) + "\" y2=\"" + std::to_string(y2) +
			"\" stroke = \"" + color1 + "\" stroke-width = \""
			+ std::to_string(width2) + "\" /> " << std::endl;
	}

	fout << "</svg>" << std::endl;
	fout.close();
}
void run_ours(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy=false) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;
	igl::Timer timer;
	double time_knot = 0;
	double time_solve = 0;
	double precision = 0;

	Eigen::MatrixXd ver;
	int nbr = nbr_pts;// nbr of points
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, param_perturbed;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	bool corners = true;
	get_model_sample_points(nbr, ver, F, param, method, corners, path);
	//std::cout << "ver\n" << ver << std::endl;
	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int perturb_itr = 0;
	//double per_ours = 0.9;
	//double per = 0.2;
	int target_steps = 10; 
	bool enable_max_fix_nbr = true;
	igl::opengl::glfw::Viewer viewer;
	if (0) {
		//viewer.data().set_mesh(param_perturbed, F);
		viewer.data().add_points(ver, ecolor);
		viewer.launch();
		exit(0);
	}
	timer.start();
	std::cout << "before generating knot vectors" << std::endl;
	std::cout << "data size " << ver.rows() << std::endl;
	generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param,per_ours,per, target_steps, enable_max_fix_nbr);
	//lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, param_perturbed, F, perturb_itr, per);
	timer.stop();
	time_knot = timer.getElapsedTimeInSec();

	//std::cout << "perturbed max dis " << perturbed_distance(param, param_perturbed) << std::endl;
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "before initialize the basis " << std::endl;
		PartialBasis basis(surface);
		std::cout << "initialize the basis done" << std::endl;
		std::cout<<"before solving control points"<<std::endl;
		timer.start();
		solve_control_points_for_fairing_surface(surface, param, ver,basis);
		timer.stop();
		time_solve = timer.getElapsedTimeInSec();
		surface_visulization(surface, 100, SPs, SFs);
		if (enable_local_energy) {
			double timeitr = 0;
			for (int i = 0; i < 50; i++) {
				timer.start();
				Eigen::MatrixXd energy, euu, evv, euv;
				energy = surface_energy_calculation(surface, basis, 1, euu, evv, euv);
				bool uorv;
				int which;
				double max_energy;
				detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
				std::vector<double> Unew = surface.U;
				std::vector<double> Vnew = surface.V;
				if (!uorv) {// u get updated
					double value = (Unew[which] + Unew[which + 1]) / 2;
					surface.U = knot_vector_insert_one_value(Unew, value);
				}
				else {
					double value = (Vnew[which] + Vnew[which + 1]) / 2;
					surface.V = knot_vector_insert_one_value(Vnew, value);
				}
				std::cout << "knot vector get inserted" << std::endl;
				basis.clear();
				basis.init(surface);
				solve_control_points_for_fairing_surface(surface, param, ver, basis);
				timer.stop();
				timeitr += timer.getElapsedTimeInSec();
				std::cout << " control points solved" << std::endl;
				surface_visulization(surface, 100, SPs, SFs);
				
				igl::write_triangle_mesh(path + "ours_"+"p"+ std::to_string(nbr)+"_refine_"+std::to_string(i)+"_m_"+std::to_string(method)+tail+".obj", SPs, SFs);
				std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps","time_itr","max_energy"} };
				int cps = (surface.nu() + 1)*(surface.nv() + 1);
				precision = max_interpolation_err(ver, param, surface);
				std::vector<double> data = { {time_knot, time_solve,precision, 
					double(surface.nu()),double(surface.nv()), double(cps),timeitr,max_energy} };
				write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv",
					titles, data);
			}
			return;
		}
		
	}
	
	
	std::cout << "final U and V, "<<surface.U.size()<<" "<<surface.V.size() << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	precision = max_interpolation_err(ver, param, surface);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		write_points(path + "pts"+std::to_string(nbr)+"_m_" +std::to_string(method)+".obj", ver);
		igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) +tail+ ".obj", SPs, SFs);
		
		std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps", "time"} };
		int cps = (surface.nu() + 1)*(surface.nv() + 1);
		std::vector<double> data = { {time_knot, time_solve,precision, double(surface.nu()),double(surface.nv()), double(cps),
			time_knot+time_solve} };
		write_csv(path + "ours_" + "p" + std::to_string(nbr)  + "_m_" + std::to_string(method) + tail+".csv",
			titles, data);
	}
	output_timing();
	if (1) {// write svg files
		write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "param.svg", param);
		write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
	}
	std::cout << "total time " << time_knot + time_solve << std::endl;


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
	Eigen::MatrixXd p0(1, 3), p1(1, 3);
	//p0.row(0) = BSplineSurfacePoint(surface, 0, 0);
	//p1.row(0) = BSplineSurfacePoint(surface, 0, 1);
	////std::cout << "p0\n" << p0 << std::endl;
	//viewer.data().add_points(p0, red);
	//viewer.data().add_points(p1, green);
	//viewer.launch();
}

// algorithm = 0, ours;
// algorithm = 1, lofting;
// algorithm = 2, averaging;


void run_other_sampling_strategies(const int model, double &per_ours, const std::string path, const std::string tail,
	const double per, const int algorithm) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;
	igl::Timer timer;
	double time_knot = 0;
	double time_solve = 0;
	double precision = 0;
	int nbr;
	Eigen::MatrixXd ver;
	
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, param_perturbed;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	int sampling_strategy = 3;
	bool enable_local_energy = false;
	///////////////////////////////
	// sampling strategies
	if (sampling_strategy == -1) {// use the random sampling method
		nbr = 100;
		get_model_sample_points(nbr, ver, F, param, method, true, path);
	}
	if (sampling_strategy == 0) {// sample gridded data
		int nbr_row = 9;
		int nbr_col = 30;
		nbr = nbr_row * nbr_col;// nbr of points
		param.resize(nbr, 2);
		int current = 0;
		for (int i = 0; i < nbr_row; i++) {
			for (int j = 0; j < nbr_col; j++) {
				param(current, 0) = double(1) / (nbr_row - 1)*i;
				param(current, 1) = double(1) / (nbr_col - 1)*j;
				current++;
			}
		}
	}
	if (sampling_strategy == 1) {// sample non-uniform density
		srand(1);
		int nbr_other = 150;
		int nbr_local = 100;
		double xmin = 0.3, xmax = 0.5;
		double ymin = 0.3, ymax = 0.5;
		nbr = nbr_other + nbr_local;
		param.resize(nbr, 2);
		int current = 0;
		for (int i = 0; i < INFINITE; i++) {
			Vector3d para3d = Vector3d::Random();
			double u = (1 + para3d[0]) / 2;
			double v = (1 + para3d[1]) / 2;
			if (u>=xmin && u<=xmax && v>=ymin && v<=ymax) {
				continue;
			}
			param.row(current) << u, v;
			current++;
			if (current == nbr_other) {
				break;
			}
		}
		int c1 = 0;
		for (int i = 0; i < INFINITE; i++) {
			Vector3d para3d = Vector3d::Random();
			double u = (1 + para3d[0]) / 2;
			double v = (1 + para3d[1]) / 2;
			if (u<xmin || u>xmax || v<ymin || v>ymax) {
				continue;
			}
			param.row(current+c1) << u, v;
			c1++;
			if (c1 == nbr_local) {
				break;
			}
		}

	}
	if (sampling_strategy == 2) {// sample row data
		int nbr_row = 6;
		int nbr_col = 20;
		nbr = nbr_row * nbr_col;// nbr of points
		param.resize(nbr, 2);
		int current = 0;
		for (int i = 0; i < nbr_row; i++) {
			for (int j = 0; j < nbr_col; j++) {
				param(current, 0) = double(1) / (nbr_row - 1)*i;
				Vector3d spd = Vector3d::Random();
				param(current, 1) = (spd[0]+1)/2;
				current++;
			}
		}
	}
	if (sampling_strategy == 3) {// sample not close data
		double para_noise=0.7;// noise level of the parameters
		int nbr_row = 10;
		int nbr_col = 10;
		double row_dis = double(1) / (nbr_row );
		double col_dis = double(1) / (nbr_col );
		nbr = (nbr_row) * (nbr_col);// nbr of points
		param.resize(nbr, 2);
		int current = 0;
		for (int i = 0; i < nbr_row; i++) {
			for (int j = 0; j < nbr_col; j++) {
				double upo = row_dis * (i + 0.5);
				double vpo = col_dis * (j + 0.5);
				Vector3d noise = Vector3d::Random();
				param(current, 0) = upo + row_dis * 0.5*noise[0] * para_noise;
				param(current, 1) = vpo + col_dis * 0.5*noise[1] * para_noise;
				current++;
			}
		}
		//std::cout << "para\n" << param << std::endl;
	}
	
	//std::cout << "parameters\n" << param << std::endl;
	//////////////////////////////
	if (sampling_strategy >= 0) {
		input_parameters_get_model_sample_points(ver, F, param, method);
	}
	
	//std::cout << "ver\n" << ver << std::endl;
	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int perturb_itr = 0;
	//double per_ours = 0.9;
	//double per = 0.2;
	int target_steps = 10;
	bool enable_max_fix_nbr = true;
	igl::opengl::glfw::Viewer viewer;
	if (0) {
		//viewer.data().set_mesh(param_perturbed, F);
		viewer.data().add_points(ver, ecolor);
		viewer.launch();
		exit(0);
	}
	std::string algo_str;
	std::cout << "before generating knot vectors" << std::endl;
	std::cout << "data size " << ver.rows() << std::endl;
	timer.start();
	
	if (algorithm == 0) {
		algo_str = "ours";
		generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);
	}
	if (algorithm == 1) {
		algo_str = "lofting";
		lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, per);
	}
	if (algorithm == 2) {
		algo_str = "pigel";
		piegl_method_generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per);
	}
	
	
	
	timer.stop();
	time_knot = timer.getElapsedTimeInSec();

	//std::cout << "perturbed max dis " << perturbed_distance(param, param_perturbed) << std::endl;
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "before initialize the basis " << std::endl;
		PartialBasis basis(surface);
		std::cout << "initialize the basis done" << std::endl;
		std::cout << "before solving control points" << std::endl;
		timer.start();
		solve_control_points_for_fairing_surface(surface, param, ver, basis);
		timer.stop();
		time_solve = timer.getElapsedTimeInSec();
		surface_visulization(surface, 100, SPs, SFs);
		if (enable_local_energy) {
			double timeitr = 0;
			for (int i = 0; i < 25; i++) {
				timer.start();
				Eigen::MatrixXd energy, euu, evv, euv;
				energy = surface_energy_calculation(surface, basis, 1, euu, evv, euv);
				bool uorv;
				int which;
				double max_energy;
				detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
				std::vector<double> Unew = surface.U;
				std::vector<double> Vnew = surface.V;
				if (!uorv) {// u get updated
					double value = (Unew[which] + Unew[which + 1]) / 2;
					surface.U = knot_vector_insert_one_value(Unew, value);
				}
				else {
					double value = (Vnew[which] + Vnew[which + 1]) / 2;
					surface.V = knot_vector_insert_one_value(Vnew, value);
				}
				std::cout << "knot vector get inserted" << std::endl;
				basis.clear();
				basis.init(surface);
				solve_control_points_for_fairing_surface(surface, param, ver, basis);
				timer.stop();
				timeitr += timer.getElapsedTimeInSec();
				std::cout << " control points solved" << std::endl;
				surface_visulization(surface, 100, SPs, SFs);

				igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);
				std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps","time_itr","max_energy"} };
				int cps = (surface.nu() + 1)*(surface.nv() + 1);
				precision = max_interpolation_err(ver, param, surface);
				std::vector<double> data = { {time_knot, time_solve,precision,
					double(surface.nu()),double(surface.nv()), double(cps),timeitr,max_energy} };
				write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv",
					titles, data);
			}
			return;
		}

	}


	std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	precision = max_interpolation_err(ver, param, surface);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
		igl::write_triangle_mesh(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);

		std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps", "time"} };
		int cps = (surface.nu() + 1)*(surface.nv() + 1);
		std::vector<double> data = { {time_knot, time_solve,precision, double(surface.nu()),double(surface.nv()), double(cps),
			time_knot + time_solve} };
		write_csv(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".csv",
			titles, data);
	}
	output_timing();
	if (1) {// write svg files
		write_svg_pts(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "param.svg", param);
		write_svg_knot_vectors(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
	}
	std::cout << "total time " << time_knot + time_solve << std::endl;

}
void ours_inserting_redundant(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_insertion = false) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;
	igl::Timer timer;
	double time_knot = 0;
	double time_solve = 0;
	double precision = 0;

	Eigen::MatrixXd ver;
	int nbr = nbr_pts;// nbr of points
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, param_perturbed;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	bool corners = false;
	get_model_sample_points(nbr, ver, F, param, method, corners, path);
	//std::cout << "ver\n" << ver << std::endl;
	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int perturb_itr = 0;
	//double per_ours = 0.9;
	//double per = 0.2;
	int target_steps = 10;
	bool enable_max_fix_nbr = true;
	igl::opengl::glfw::Viewer viewer;
	if (0) {
		//viewer.data().set_mesh(param_perturbed, F);
		viewer.data().add_points(ver, ecolor);
		viewer.launch();
		exit(0);
	}
	timer.start();
	std::cout << "before generating knot vectors" << std::endl;
	std::cout << "data size " << ver.rows() << std::endl;
	generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);
	//lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, param_perturbed, F, perturb_itr, per);
	timer.stop();
	time_knot = timer.getElapsedTimeInSec();

	//std::cout << "perturbed max dis " << perturbed_distance(param, param_perturbed) << std::endl;
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "before initialize the basis " << std::endl;
		PartialBasis basis(surface);
		std::cout << "initialize the basis done" << std::endl;
		std::cout << "before solving control points" << std::endl;
		//timer.start();
		//solve_control_points_for_fairing_surface(surface, param, ver, basis);
		////timer.stop();
		//time_solve = timer.getElapsedTimeInSec();
		//surface_visulization(surface, 100, SPs, SFs);
		bool local_energy = true;
		if (enable_insertion) {
			/*double minvalue = 0.2; or 0.5, 0.015, 20
			double maxvalue = minvalue + 0.015;
			int nbrinsert = 30;*/// best settings

			// we first insert knots to create problems
			double minvalue = 0.4;
			double maxvalue = minvalue + 0.1;
			int nbrinsert = 20;
			double interval = (maxvalue - minvalue) / nbrinsert;
			
			for (int i = 0; i < nbrinsert; i++) {
				std::vector<double> Unew = surface.U;
				std::vector<double> Vnew = surface.V;
				double value = minvalue + interval * i;
				surface.V = knot_vector_insert_one_value(Vnew, value);
				
			}
			timer.start();
			basis.clear();
			basis.init(surface);
			solve_control_points_for_fairing_surface(surface, param, ver, basis);
			timer.stop();
			time_solve = timer.getElapsedTimeInSec();
			surface_visulization(surface, 100, SPs, SFs);
			if (local_energy) {
				double timeitr = 0;
				for (int i = 0; i < 100; i++) {
					timer.start();
					Eigen::MatrixXd energy, euu, evv, euv;
					energy = surface_energy_calculation(surface, basis, 1, euu, evv, euv);
					bool uorv;
					int which;
					double max_energy;
					detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
					std::vector<double> Unew = surface.U;
					std::vector<double> Vnew = surface.V;
					if (!uorv) {// u get updated
						double value = (Unew[which] + Unew[which + 1]) / 2;
						surface.U = knot_vector_insert_one_value(Unew, value);
						std::cout << "##insert u = " << value << std::endl;
					}
					else {
						double value = (Vnew[which] + Vnew[which + 1]) / 2;
						surface.V = knot_vector_insert_one_value(Vnew, value);
						std::cout << "##insert v = " << value << std::endl;
					}
					std::cout << "knot vector get inserted" << std::endl;
					basis.clear();
					basis.init(surface);
					solve_control_points_for_fairing_surface(surface, param, ver, basis);
					timer.stop();
					timeitr += timer.getElapsedTimeInSec();
					std::cout << " control points solved" << std::endl;
					surface_visulization(surface, 100, SPs, SFs);

					igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);
					std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps","time_itr","max_energy"} };
					int cps = (surface.nu() + 1)*(surface.nv() + 1);
					precision = max_interpolation_err(ver, param, surface);
					std::vector<double> data = { {time_knot, time_solve,precision,
						double(surface.nu()),double(surface.nv()), double(cps),timeitr,max_energy} };
					write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv",
						titles, data);
				}
				return;
			}
			
			
		}

	}


	std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	precision = max_interpolation_err(ver, param, surface);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
		igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);

		std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps", "time"} };
		int cps = (surface.nu() + 1)*(surface.nv() + 1);
		std::vector<double> data = { {time_knot, time_solve,precision, double(surface.nu()),double(surface.nv()), double(cps),
			time_knot + time_solve} };
		write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".csv",
			titles, data);
	}
	output_timing();
	if (1) {// write svg files
		write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "param.svg", param);
		write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
	}
	std::cout << "total time " << time_knot + time_solve << std::endl;


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
	Eigen::MatrixXd p0(1, 3), p1(1, 3);
	//p0.row(0) = BSplineSurfacePoint(surface, 0, 0);
	//p1.row(0) = BSplineSurfacePoint(surface, 0, 1);
	////std::cout << "p0\n" << p0 << std::endl;
	//viewer.data().add_points(p0, red);
	//viewer.data().add_points(p1, green);
	//viewer.launch();
}
void run_Seungyong(const int model, const int nbr_pts, const double tolerance, const std::string path) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;
	igl::Timer timer;
	double time_run = 0;
	Eigen::MatrixXd ver;
	int nbr = nbr_pts;// nbr of points
	Eigen::MatrixXi F;
	Eigen::MatrixXd param;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	bool corners = false;
	get_model_sample_points(nbr,  ver, F, param, method, corners, path);

	int degree1 = 3;
	int degree2 = 3;
	double fair_parameter = 1e-7;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	igl::opengl::glfw::Viewer viewer;

	std::vector<Bsurface> surfaces;
	timer.start();
	iteratively_approximate_method(degree1, degree2, Uknot, Vknot, param, ver, tolerance, surfaces, fair_parameter);
	timer.stop();
	time_run = timer.getElapsedTimeInSec();
	std::cout << "run_Seungyong finished, surface levels " << surfaces.size() << std::endl;
	std::cout << "final surface size " << surfaces.back().nu()+1 << " " << surfaces.back().nv()+1
		<< " total control points " << (surfaces.back().nu()+1)*(surfaces.back().nv()+1) << std::endl;
	std::cout << "Seungyong time " << time_run << std::endl;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	
	write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
	for (int i = 0; i < surfaces.size(); i++) {
		surface_visulization(surfaces[i], 100, SPs, SFs);
		igl::write_triangle_mesh(path + "Seungyong_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) +"_level_"+
			std::to_string(i) + ".obj", SPs, SFs);
	}
	
}

void run_piegl(const int model, const int nbr_pts, const double per = 0.2) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;

	igl::Timer timer;
	double time_total = 0;
	Eigen::MatrixXd ver;
	int nbr = nbr_pts;// nbr of points
	Eigen::MatrixXi F;
	Eigen::MatrixXd param;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	bool corners = false;
	get_model_sample_points(nbr, ver, F, param, method, corners, "./");

	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int perturb_itr = 0;


	//mesh_parameter_perturbation(param, F, param_perturbed, perturb_itr);
	//std::cout << "param_perturbed\n" << param_perturbed << std::endl;
	igl::opengl::glfw::Viewer viewer;
	if (0) {
		//viewer.data().set_mesh(param_perturbed, F);
		viewer.data().add_points(ver, ecolor);
		viewer.launch();
		exit(0);
	}
	timer.start();
	piegl_method_generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param,  per);
	//lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, param_perturbed, F, perturb_itr, per);

	Eigen::MatrixXd paramout(ver.rows(), 3), zeros(ver.rows(), 1);
	zeros = Eigen::MatrixXd::Constant(ver.rows(), 1, 0);

	//paramout << param_perturbed, zeros;
	//std::cout << "perturbed max dis " << perturbed_distance(param, param_perturbed) << std::endl;
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "before initialize the basis " << std::endl;
		PartialBasis basis(surface);
		std::cout << "initialize the basis done" << std::endl;
		std::cout << "before solving control points" << std::endl;
		solve_control_points_for_fairing_surface(surface, param, ver, basis);
		surface_visulization(surface, 100, SPs, SFs);

	}
	timer.stop();
	time_total = timer.getElapsedTimeInSec();

	std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	std::cout << "total time " << time_total << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
		write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
		igl::write_triangle_mesh(path + "piegl_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", SPs, SFs);
	}
	output_timing();

	std::cout << "final surface size " << surface.nu()+1 << " " << surface.nv()+1
		<< " total control points " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;



	/*
	viewer.data().set_edges(bdver, edges, fcolor);*/
	//viewer.data().add_points(vector_to_matrix_3d(inter_pts), ecolor);

	// see the linear interpolated surface
	//viewer.data().set_mesh(param, F);
	//viewer.data().add_edges(edge0, edge1, pcolor);

	//viewer.data().set_mesh(param_perturbed, F);
	viewer.data().clear();
	viewer.data().set_mesh(SPs, SFs);
	//viewer.data().add_points(ver, ecolor);
	Eigen::MatrixXd p0(1, 3), p1(1, 3);
	p0.row(0) = BSplineSurfacePoint(surface, 0, 0);
	p1.row(0) = BSplineSurfacePoint(surface, 0, 1);
	//std::cout << "p0\n" << p0 << std::endl;
	viewer.data().add_points(p0, red);
	viewer.data().add_points(p1, green);
	viewer.launch();
}
void run_lofting(const int model, const int nbr_pts, const double per,const std::string path) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;

	igl::Timer timer;
	double time_total = 0;
	Eigen::MatrixXd ver;
	int nbr = nbr_pts;// nbr of points
	Eigen::MatrixXi F;
	Eigen::MatrixXd param;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	bool corners = false;
	get_model_sample_points(nbr, ver, F, param, method, corners, "./");

	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int perturb_itr = 0;


	//mesh_parameter_perturbation(param, F, param_perturbed, perturb_itr);
	//std::cout << "param_perturbed\n" << param_perturbed << std::endl;
	igl::opengl::glfw::Viewer viewer;
	if (0) {
		//viewer.data().set_mesh(param_perturbed, F);
		viewer.data().add_points(ver, ecolor);
		viewer.launch();
		exit(0);
	}
	timer.start();
	//piegl_method_generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per);
	lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, per);

	Eigen::MatrixXd paramout(ver.rows(), 3), zeros(ver.rows(), 1);
	zeros = Eigen::MatrixXd::Constant(ver.rows(), 1, 0);

	//paramout << param_perturbed, zeros;
	//std::cout << "perturbed max dis " << perturbed_distance(param, param_perturbed) << std::endl;
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "before initialize the basis " << std::endl;
		PartialBasis basis(surface);
		std::cout << "initialize the basis done" << std::endl;
		std::cout << "before solving control points" << std::endl;
		solve_control_points_for_fairing_surface(surface, param, ver, basis);
		surface_visulization(surface, 100, SPs, SFs);

	}
	timer.stop();
	time_total = timer.getElapsedTimeInSec();

	std::cout << "nu and nv " << surface.nu() << " " << surface.nv() << std::endl;
	std::cout << "#cp " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	std::cout << "total time " << time_total << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		
		write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
		igl::write_triangle_mesh(path + "lofting_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", SPs, SFs);
	}
	output_timing();

	std::cout << "final surface size " << surface.nu() + 1 << " " << surface.nv() + 1
		<< " total control points " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;



	/*
	viewer.data().set_edges(bdver, edges, fcolor);*/
	//viewer.data().add_points(vector_to_matrix_3d(inter_pts), ecolor);

	// see the linear interpolated surface
	//viewer.data().set_mesh(param, F);
	//viewer.data().add_edges(edge0, edge1, pcolor);

	//viewer.data().set_mesh(param_perturbed, F);
	viewer.data().clear();
	viewer.data().set_mesh(SPs, SFs);
	//viewer.data().add_points(ver, ecolor);
	Eigen::MatrixXd p0(1, 3), p1(1, 3);
	p0.row(0) = BSplineSurfacePoint(surface, 0, 0);
	p1.row(0) = BSplineSurfacePoint(surface, 0, 1);
	//std::cout << "p0\n" << p0 << std::endl;
	viewer.data().add_points(p0, red);
	viewer.data().add_points(p1, green);
	viewer.launch();
}
void read_mesh_series(std::string path, std::string namebase, int end) {
	std::string file0 = path + namebase + std::to_string(0) + ".obj";
	Eigen::MatrixXd ver;
	Eigen::MatrixXi f;
	igl::read_triangle_mesh(file0, ver, f);
	for (int i = 1; i < end + 1; i++) {// 0 has already been included
		std::cout << "reading " << i << std::endl;
		std::string file = path + namebase + std::to_string(i) + ".obj";
		Eigen::MatrixXd ver_t;
		Eigen::MatrixXi f_t;
		igl::read_triangle_mesh(file, ver_t, f_t);
		ver = ver + ver_t;
	}
	file0= path + namebase+ "sum.obj";
	igl::write_triangle_mesh(file0, ver, f);
}
void run_mesh_reconstruction(const std::string inpath, const std::string modelname, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy = false) {
	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.9, 0.9, 0.9;; pcolor << 0, 0.9, 0.5;
	red << 1, 0, 0; green << 0, 1, 0; blue << 0, 0, 1;
	igl::Timer timer;
	double time_knot = 0;
	double time_solve = 0;
	double precision = 0;

	Eigen::MatrixXd ver;
	
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, paramout;
	//std::string modelname = "mask3kf.obj";
	//std::string modelname = "simplified_2500_TigerMask4.obj";
	std::string meshfile = inpath + modelname;
		//"D:\\vs\\sparse_data_interpolation\\meshes\\meshmodels\\" + modelname;
	
	//get_mesh_vertices_and_parametrization(ver, F, param);
	bool corners = false;
	mesh_parameterization(meshfile, ver, param, F);
	paramout.resize(param.rows(), 3);
	Eigen::VectorXd param_zero= Eigen::VectorXd::Zero(param.rows());
	paramout << param, param_zero;
	igl::write_triangle_mesh(path + "param_" + modelname + tail, paramout, F);

	//get_model_sample_points(nbr, ver, F, param, method, corners);
	//std::cout << "ver\n" << ver << std::endl;
	int nbr = param.rows();
	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot = { {0,0,0,0,1,1,1,1} };
	std::vector<double> Vknot = Uknot;
	int target_steps = 10;
	bool enable_max_fix_nbr = true;
	igl::opengl::glfw::Viewer viewer;
	if (0) {
		//viewer.data().set_mesh(param_perturbed, F);
		viewer.data().add_points(ver, ecolor);
		viewer.launch();
		exit(0);
	}
	timer.start();
	std::cout << "before generating knot vectors" << std::endl;
	std::cout << "data size " << ver.rows() << std::endl;
	generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);
	
	timer.stop();
	time_knot = timer.getElapsedTimeInSec();

	//std::cout << "perturbed max dis " << perturbed_distance(param, param_perturbed) << std::endl;
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	int visual_nbr = 200;
	if (1) {
		surface.degree1 = 3;
		surface.degree2 = 3;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "before initialize the basis " << std::endl;
		PartialBasis basis(surface);
		std::cout << "initialize the basis done" << std::endl;
		std::cout << "before solving control points" << std::endl;
		timer.start();
		solve_control_points_for_fairing_surface(surface, param, ver, basis);
		timer.stop();
		time_solve = timer.getElapsedTimeInSec();
		surface_visulization(surface, visual_nbr, SPs, SFs);
		if (enable_local_energy) {
			double timeitr = 0;
			for (int i = 0; i < 20; i++) {
				timer.start();
				Eigen::MatrixXd energy, euu, evv, euv;
				energy = surface_energy_calculation(surface, basis, 1, euu, evv, euv);
				bool uorv;
				int which;
				double max_energy;
				detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
				std::vector<double> Unew = surface.U;
				std::vector<double> Vnew = surface.V;
				if (!uorv) {// u get updated
					double value = (Unew[which] + Unew[which + 1]) / 2;
					surface.U = knot_vector_insert_one_value(Unew, value);
				}
				else {
					double value = (Vnew[which] + Vnew[which + 1]) / 2;
					surface.V = knot_vector_insert_one_value(Vnew, value);
				}
				std::cout << "knot vector get inserted" << std::endl;
				basis.clear();
				basis.init(surface);
				solve_control_points_for_fairing_surface(surface, param, ver, basis);
				timer.stop();
				timeitr += timer.getElapsedTimeInSec();
				std::cout << " control points solved" << std::endl;
				surface_visulization(surface, visual_nbr, SPs, SFs);

				igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + modelname + tail, SPs, SFs);
				std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps","time_itr","max_energy"} };
				int cps = (surface.nu() + 1)*(surface.nv() + 1);
				precision = max_interpolation_err(ver, param, surface);
				std::vector<double> data = { {time_knot, time_solve,precision,
					double(surface.nu()),double(surface.nv()), double(cps),timeitr,max_energy} };
				write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + modelname + tail + ".csv",
					titles, data);
			}
			return;
		}

	}


	std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	precision = max_interpolation_err(ver, param, surface);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		write_points(path + "pts" + std::to_string(nbr) + "_m_" + modelname, ver);
		igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail, SPs, SFs);

		write_control_pts(surface.control_points, path + "cps" + std::to_string(nbr) + "_m_" + modelname);
		std::vector<std::string> titles = { {"time_knot","time_solve","precision","nu","nv", "cps"} };
		int cps = (surface.nu() + 1)*(surface.nv() + 1);
		std::vector<double> data = { {time_knot, time_solve,precision, double(surface.nu()),double(surface.nv()), double(cps)} };
		write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail + ".csv",
			titles, data);
	}
	output_timing();
	if (1) {// write svg files
		write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail + "param.svg", param);
		write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail + "knots.svg", surface.U, surface.V);
	}
	std::cout << "total time " << time_knot + time_solve << std::endl;


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
	Eigen::MatrixXd p0(1, 3), p1(1, 3);
	//p0.row(0) = BSplineSurfacePoint(surface, 0, 0);
	//p1.row(0) = BSplineSurfacePoint(surface, 0, 1);
	////std::cout << "p0\n" << p0 << std::endl;
	//viewer.data().add_points(p0, red);
	//viewer.data().add_points(p1, green);
	//viewer.launch();
}
void run_pia(const int model, const int nbr_pts, const int max_itr,const double threadshold,
	const int cp_nbr_sqrt, const std::string &path) {
	igl::Timer timer;
	double time_total = 0;
	Eigen::MatrixXd ver;
	int nbr = nbr_pts;// nbr of points
	Eigen::MatrixXi F;
	Eigen::MatrixXd param;
	//get_mesh_vertices_and_parametrization(ver, F, param);
	int method = model;
	bool corners = false;
	get_model_sample_points(nbr, ver, F, param, method, corners, "./");

	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot;
	std::vector<double> Vknot;
	int nu = cp_nbr_sqrt-1;
	int nv = cp_nbr_sqrt-1;
	std::cout << "starting pia method" << std::endl;
	timer.start();
	initialize_pia_knot_vectors(degree1, degree2, Uknot, Vknot, nu, nv);
	std::cout << "pia knot vectors get initialized" << std::endl;
	print_vector(Uknot);
	print_vector(Vknot);
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = degree1;
		surface.degree2 = degree2;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "nu " << surface.nu() << std::endl;
		std::cout << "nv " << surface.nv() << std::endl;
		progressive_iterative_approximation(surface, param, ver, max_itr, threadshold);
		surface_visulization(surface, 100, SPs, SFs);

	}
	timer.stop();
	time_total = timer.getElapsedTimeInSec();

	std::cout << "nu and nv " << surface.nu() << " " << surface.nv() << std::endl;
	std::cout << "#cp " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	std::cout << "total time " << time_total << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;
		write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
		igl::write_triangle_mesh(path + "pia_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", SPs, SFs);
	}
	output_timing();

	std::cout << "final surface size " << surface.nu() + 1 << " " << surface.nv() + 1
		<< " total control points " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;

}

void run_pia_mesh_reconstruct(const std::string meshfile, const int max_itr, const double threadshold,
	const int cp_nbr_sqrt, const std::string &path) {
	igl::Timer timer;
	double time_total = 0;
	Eigen::MatrixXd ver;
	
	Eigen::MatrixXi F;
	Eigen::MatrixXd param;
	mesh_parameterization(meshfile, ver, param, F);
	int degree1 = 3;
	int degree2 = 3;
	std::vector<double> Uknot;
	std::vector<double> Vknot;
	int nu = cp_nbr_sqrt - 1;
	int nv = cp_nbr_sqrt - 1;
	std::cout << "starting pia method" << std::endl;
	timer.start();
	initialize_pia_knot_vectors(degree1, degree2, Uknot, Vknot, nu, nv);
	std::cout << "pia knot vectors get initialized" << std::endl;
	print_vector(Uknot);
	print_vector(Vknot);
	Bsurface surface;
	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	if (1) {
		surface.degree1 = degree1;
		surface.degree2 = degree2;
		surface.U = Uknot;
		surface.V = Vknot;
		std::cout << "nu " << surface.nu() << std::endl;
		std::cout << "nv " << surface.nv() << std::endl;
		progressive_iterative_approximation(surface, param, ver, max_itr, threadshold);
		surface_visulization(surface, 100, SPs, SFs);

	}
	timer.stop();
	time_total = timer.getElapsedTimeInSec();

	std::cout << "nu and nv " << surface.nu() << " " << surface.nv() << std::endl;
	std::cout << "#cp " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;
	print_vector(surface.U);
	print_vector(surface.V);
	std::cout << "maximal interpolation error " << max_interpolation_err(ver, param, surface) << std::endl;
	std::cout << "total time " << time_total << std::endl;
	bool write_file = true;
	if (write_file)
	{
		Eigen::MatrixXi Pf;

		igl::write_triangle_mesh(meshfile+"_pia_cp_"+std::to_string(cp_nbr_sqrt*cp_nbr_sqrt)+".obj", SPs, SFs);
	}
	output_timing();

	std::cout << "final surface size " << surface.nu() + 1 << " " << surface.nv() + 1
		<< " total control points " << (surface.nu() + 1)*(surface.nv() + 1) << std::endl;

}