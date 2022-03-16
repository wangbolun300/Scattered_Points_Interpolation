#pragma once
// This file contains the benchmark surfaces.
//
//
//
//
#include <cmath>
#include <sparse_interp/Types.hpp>
#include <limits>
#include <math.h>
#include <sparse_interp/mesh_processing.h>

namespace SIBSplines
{

	namespace examples
	{
		double peak_function(const double x, const double y)
		{
			double r = 3 * pow(1 - x, 2) * exp(-x * x - (y + 1) * (y + 1)) - 10 * (0.2 * x - pow(x, 3) - pow(y, 5)) * exp(-x * x - y * y) - 1 / 3 * exp(-pow(x + 1, 2) - y * y);
			if (fabs(r) < SCALAR_ZERO)
			{
				r = 0;
			}
			r /= 2;
			return r;
		}

		Eigen::MatrixXd get_peak_sample_points(const int nbr, const int skip)
		{
			Eigen::MatrixXd ver;
			ver.resize(nbr - skip, 3);
			for (int i = 0; i < nbr - 4; i++)
			{
				Vector3d para3d = Vector3d::Random();
				if (i < skip)
				{
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

		double contour_function(const double x, const double y)
		{
			return 4 * x * exp(-x * x - y * y);
		}
		Eigen::MatrixXd get_contour_sample_points(const int nbr, const int skip)
		{
			Eigen::MatrixXd ver;
			ver.resize(nbr - skip, 3);
			double s = 2;
			for (int i = 0; i < nbr - 4; i++)
			{
				Vector3d para3d = Vector3d::Random();
				if (i < skip)
				{
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
		double hyperbolic_function(const double x, const double y)
		{
			return x * x - y * y;
		}
		Eigen::MatrixXd get_hyperbolic_sample_points(const int nbr, const int skip)
		{
			Eigen::MatrixXd ver;
			ver.resize(nbr - skip, 3);
			double s = 1;
			for (int i = 0; i < nbr - 4; i++)
			{
				Vector3d para3d = Vector3d::Random();
				if (i < skip)
				{
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
		double sinus_function(const double x, const double y)
		{
			return sin(3 * 3.1415926 * (x * x + y * y)) / 10;
		}
		Eigen::MatrixXd get_sinus_sample_points(const int nbr, const int skip)
		{
			Eigen::MatrixXd ver;
			ver.resize(nbr - skip, 3);
			double s = 1;
			for (int i = 0; i < nbr - 4; i++)
			{
				Vector3d para3d = Vector3d::Random();
				if (i < skip)
				{
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

		Vector3d bilinear_function(const Vector3d &v0s, const Vector3d &v0e, const Vector3d &v1s, const Vector3d &v1e,
								   const double u, const double v)
		{
			Vector3d m0 = (v0e - v0s) * u + v0s;
			Vector3d m1 = (v1e - v1s) * u + v1s;
			Vector3d m = (m1 - m0) * v + m0;
			return m;
		}
		Eigen::MatrixXd get_bilinear_sample_points(const int nbr, const int skip, Eigen::MatrixXd &param)
		{
			Eigen::MatrixXd ver;
			ver.resize(nbr - skip, 3);
			param.resize(nbr - skip, 2);
			double s = 1; // domain scale is [-s, s]x[-s, s]
			Vector3d v0s(0, 0, 1);
			Vector3d v0e(1, 1, 0);
			Vector3d v1s(0, 1, 1);
			Vector3d v1e(0, 0, 0);
			for (int i = 0; i < nbr - 4; i++)
			{
				Vector3d para3d = Vector3d::Random();
				if (i < skip)
				{
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

		// u0 v0 are in [0,1].
		// the function u in [1,2], v in [-pi/2, pi/2]
		Vector3d snail_function(double u0, double v0)
		{
			const double pi = 3.1415926;
			double u = 1 + u0;
			double v = -pi / 2 + pi * v0;
			double x = u * cos(v) * sin(u);
			double y = u * cos(u) * cos(v);
			double z = -u * sin(v);
			return Vector3d(x, y, z);
		}

		// direct project x y as u v, and rescale the scene into [0, 1]x[0, 1]
		void direct_project_x_y_and_parametrization(const Eigen::MatrixXd &ver, Eigen::MatrixXd &param, Eigen::MatrixXi &F)
		{

			Eigen::MatrixXd paras(ver.rows(), 2), para_new; // the 2d projection of the points
			paras << ver.col(0), ver.col(1);
			para_new = paras;
			Eigen::VectorXi loop;
			std::cout << "start find border" << std::endl;
			// find the convex hull on 2d,
			find_border_loop(paras, loop);
			if (0)
			{ // if user is interested in this part
				Eigen::MatrixXd bdver(loop.size(), 3);
				for (int i = 0; i < loop.size(); i++)
				{
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
			for (int i = 0; i < paras.rows(); i++)
			{
				double x = paras(i, 0);
				double y = paras(i, 1);
				if (xmin > x)
				{
					xmin = x;
				}
				if (xmax < x)
				{
					xmax = x;
				}
				if (ymin > y)
				{
					ymin = y;
				}
				if (ymax < y)
				{
					ymax = y;
				}
			}
			std::cout << "bounding box: (" << xmin << ", " << ymin << "), (" << xmax << ", " << ymax << ")" << std::endl;

			double xsize = xmax - xmin;
			double ysize = ymax - ymin;
			for (int i = 0; i < paras.rows(); i++)
			{
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
			for (int i = 0; i < param.rows(); i++)
			{
				assert(param(i, 0) <= 1 && param(i, 0) >= 0 && param(i, 1) <= 1 && param(i, 1) >= 0);
			}
			// triangulation
			constrained_delaunay_triangulation(para_new, loop, F);
			// std::cout << "check directly parametrization\n" << param << std::endl;
		}

		//#define B_SPLINE_WITH_NOISE
		double noise_scale = 0.05;
		// method = 0, peak;
		// method = 1, contour;
		// method = 2, hyperbolic;
		// method = 3, sinus;
		// method = 4, bilinear
		void input_parameters_get_model_sample_points(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
													  const Eigen::MatrixXd &param, const int method)
		{
			const int nbr_models = 6;

			const auto mfunction_value = [](const int method, const double x, const double y)
			{
				Vector3d v0s(0, 0, 1);
				Vector3d v0e(1, 1, 0);
				Vector3d v1s(0, 1, 1);
				Vector3d v1e(0, 0, 0);
				switch (method)
				{
				case 0:
					return Vector3d(x, y, peak_function(x, y));
				case 1:
#ifdef B_SPLINE_WITH_NOISE
				{
					Vector3d norm(4 * exp(-x * x - y * y) - 8 * x * x * exp(-x * x - y * y), -8 * x * y * exp(-x * x - y * y), -1);
					norm = norm.normalized();
					Vector3d noise = norm * noise_scale * Vector3d::Random()[0];
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
					Vector3d noise = norm * noise_scale * Vector3d::Random()[0];
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
			std::vector<double> scale = {{3, 2, 1, 1, 1}};
			std::vector<int> seed = {{3, 7, 33, 10, 5, 0}};
			Eigen::MatrixXd ver;
			ver.resize(param.rows(), 3);

			double s = scale[method]; // domain scale is [-s, s]x[-s, s]
			int end = param.rows();
			srand(seed[method]);
			bool need_param = true;
			if (method == 4 || method == 5)
			{
				for (int i = 0; i < end; i++)
				{
					ver.row(i) = mfunction_value(method, param(i, 0), param(i, 1));
				}
				Eigen::VectorXi loop;

				find_border_loop(param, loop);
				constrained_delaunay_triangulation(param, loop, F);
			}
			else
			{
				for (int i = 0; i < end; i++)
				{
					double x = s * (param(i, 0) * 2 - 1);
					double y = s * (param(i, 1) * 2 - 1);
					ver.row(i) = mfunction_value(method, x, y);
				}
			}
			V = ver;
			if (1)
			{
				double xmin = INFINITY, ymin = INFINITY, zmin = INFINITY;
				double xmax = -INFINITY, ymax = -INFINITY, zmax = -INFINITY;
				int pnbr = 100;
				Eigen::MatrixXi faces;
				std::vector<std::vector<Vector3d>> pts;
				ver.resize(pnbr * pnbr, 3);
				int verline = 0;
				for (int i = 0; i < pnbr; i++)
				{
					for (int j = 0; j < pnbr; j++)
					{
						double upara = double(i) / (pnbr - 1);
						double vpara = double(j) / (pnbr - 1);
						double u, v;
						if (method == 4 || method == 5)
						{
							u = upara;
							v = vpara;
						}
						else
						{
							u = 2 * s * upara - s;
							v = 2 * s * vpara - s;
						}
						ver.row(verline) = mfunction_value(method, u, v);
						if (xmin > ver(verline, 0))
						{
							xmin = ver(verline, 0);
						}
						if (ymin > ver(verline, 1))
						{
							ymin = ver(verline, 1);
						}
						if (zmin > ver(verline, 2))
						{
							zmin = ver(verline, 2);
						}
						if (xmax < ver(verline, 0))
						{
							xmax = ver(verline, 0);
						}
						if (ymax < ver(verline, 1))
						{
							ymax = ver(verline, 1);
						}
						if (zmax < ver(verline, 2))
						{
							zmax = ver(verline, 2);
						}
						verline++;
					}
				}
				std::cout << "model bounding box \n"
						  << xmin << ", " << xmax << "\n"
						  << ymin << ", " << ymax << "\n"
						  << zmin << ", " << zmax << std::endl;
				faces.resize(2 * (pnbr - 1) * (pnbr - 1), 3);
				int fline = 0;
				for (int i = 0; i < pnbr - 1; i++)
				{
					for (int j = 0; j < pnbr - 1; j++)
					{
						faces.row(fline) = Vector3i(i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
						faces.row(fline + 1) = Vector3i(i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
						fline += 2;
					}
				}
				write_triangle_mesh("./model_" + std::to_string(method) + ".obj", ver, faces);
				std::cout << "write original model" << std::endl;
			}
		}

		void get_model_sample_points(const int nbr, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
									 Eigen::MatrixXd &param, const int method, bool corners, std::string path)
		{
			const int nbr_models = 6;

			const auto mfunction_value = [](const int method, const double x, const double y)
			{
				Vector3d v0s(0, 0, 1);
				Vector3d v0e(1, 1, 0);
				Vector3d v1s(0, 1, 1);
				Vector3d v1e(0, 0, 0);
				switch (method)
				{
				case 0:
					return Vector3d(x, y, peak_function(x, y));
				case 1:
#ifdef B_SPLINE_WITH_NOISE
				{
					Vector3d norm(4 * exp(-x * x - y * y) - 8 * x * x * exp(-x * x - y * y), -8 * x * y * exp(-x * x - y * y), -1);
					norm = norm.normalized();
					Vector3d noise = norm * noise_scale * Vector3d::Random()[0];
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
					Vector3d noise = norm * noise_scale * Vector3d::Random()[0];
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
			std::vector<double> scale = {{3, 2, 1, 1, 1}};
			std::vector<int> seed = {{3, 7, 33, 10, 5, 0}};
			Eigen::MatrixXd ver;
			ver.resize(nbr, 3);
			param.resize(nbr, 2);
			double s = scale[method]; // domain scale is [-s, s]x[-s, s]
			int end = corners ? nbr - 4 : nbr;
			srand(seed[method]);
			bool need_param = true;
			if (method == 4 || method == 5)
			{
				for (int i = 0; i < end; i++)
				{
					Vector3d para3d = Vector3d::Random();
					double u = (1 + para3d[0]) / 2;
					double v = (1 + para3d[1]) / 2;
					param.row(i) << u, v;
					ver.row(i) = mfunction_value(method, u, v);
				}
				if (corners)
				{
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
			else
			{
				for (int i = 0; i < end; i++)
				{
					Vector3d para3d = Vector3d::Random();
					double x = s * para3d[0];
					double y = s * para3d[1];
					ver.row(i) = mfunction_value(method, x, y);
				}
				if (corners)
				{
					ver.row(nbr - 4) = mfunction_value(method, -s, -s);
					ver.row(nbr - 3) = mfunction_value(method, -s, s);
					ver.row(nbr - 2) = mfunction_value(method, s, s);
					ver.row(nbr - 1) = mfunction_value(method, s, -s);
				}
			}
			V = ver;
			if (need_param)
			{
				direct_project_x_y_and_parametrization(V, param, F);
			}
			if (1)
			{
				double xmin = INFINITY, ymin = INFINITY, zmin = INFINITY;
				double xmax = -INFINITY, ymax = -INFINITY, zmax = -INFINITY;
				int pnbr = 100;
				Eigen::MatrixXi faces;
				std::vector<std::vector<Vector3d>> pts;
				ver.resize(pnbr * pnbr, 3);
				int verline = 0;
				for (int i = 0; i < pnbr; i++)
				{
					for (int j = 0; j < pnbr; j++)
					{
						double upara = double(i) / (pnbr - 1);
						double vpara = double(j) / (pnbr - 1);
						double u, v;
						if (method == 4 || method == 5)
						{
							u = upara;
							v = vpara;
						}
						else
						{
							u = 2 * s * upara - s;
							v = 2 * s * vpara - s;
						}
						ver.row(verline) = mfunction_value(method, u, v);
						if (xmin > ver(verline, 0))
						{
							xmin = ver(verline, 0);
						}
						if (ymin > ver(verline, 1))
						{
							ymin = ver(verline, 1);
						}
						if (zmin > ver(verline, 2))
						{
							zmin = ver(verline, 2);
						}
						if (xmax < ver(verline, 0))
						{
							xmax = ver(verline, 0);
						}
						if (ymax < ver(verline, 1))
						{
							ymax = ver(verline, 1);
						}
						if (zmax < ver(verline, 2))
						{
							zmax = ver(verline, 2);
						}
						verline++;
					}
				}
				std::cout << "model bounding box \n"
						  << xmin << ", " << xmax << "\n"
						  << ymin << ", " << ymax << "\n"
						  << zmin << ", " << zmax << std::endl;
				faces.resize(2 * (pnbr - 1) * (pnbr - 1), 3);
				int fline = 0;
				for (int i = 0; i < pnbr - 1; i++)
				{
					for (int j = 0; j < pnbr - 1; j++)
					{
						faces.row(fline) = Vector3i(i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
						faces.row(fline + 1) = Vector3i(i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
						fline += 2;
					}
				}
				write_triangle_mesh(path + "model_" + std::to_string(method) + ".obj", ver, faces);
				std::cout << "write original model" << std::endl;
			}
		}
	}
}