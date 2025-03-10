#include "test.h"
#include <sparse_interp/curve.h>
#include "functions.hpp"
#include <iostream>
#include <sparse_interp/surface.h>
#include <cmath>
#include <sparse_interp/mesh_processing.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/harmonic.h>
#include <sparse_interp/energy.h>
#include <igl/write_triangle_mesh.h>
#include <igl/Timer.h>
#include <sparse_interp/comparison.h>

namespace SIBSplines
{
	std::string example_root_path(SI_MESH_DIR);
	using namespace examples;
	using namespace comparison;

	
	void write_control_pts(std::vector<std::vector<Vector3d>> &cps, std::string file)
	{
		std::ofstream fout;
		fout.open(file);
		for (int i = 0; i < cps.size(); i++)
		{
			for (int j = 0; j < cps[i].size(); j++)
			{
				fout << "v " << cps[i][j][0] << " " << cps[i][j][1] << " " << cps[i][j][2] << std::endl;
			}
		}
		fout.close();
	}
	void write_csv(const std::string &file, const std::vector<std::string> titles, const std::vector<double> data)
	{
		std::ofstream fout;
		fout.open(file);
		for (int i = 0; i < titles.size() - 1; i++)
		{
			fout << titles[i] << ",";
		}
		fout << titles.back() << std::endl;
		for (int i = 0; i < data.size() - 1; i++)
		{
			fout << data[i] << ",";
		}
		fout << data.back() << std::endl;
		fout.close();
	}
	void write_svg_pts(const std::string &file, const Eigen::MatrixXd &param)
	{
		std::ofstream fout;
		fout.open(file);
		double scale = 1000;
		double dot = 10;
		double displace = 30;
		std::string color = "#ED6E46";
		std::string black = "#000000";
		fout << "<svg>" << std::endl;
		for (int i = 0; i < param.rows(); i++)
		{
			double x = displace + param(i, 0) * scale;
			double y = displace + param(i, 1) * scale;
			fout << " <circle cx=\"" << x << "\" cy = \"" << y << "\" r = \"" << dot << "\" fill = \"" << color << "\" />" << std::endl;
		}
		fout << "<polyline points=\"";

		std::vector<double> xlist = {{0, 0, 1, 1}};
		std::vector<double> ylist = {{0, 1, 1, 0}};

		for (int i = 0; i < 4; i++)
		{
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
	void write_svg_knot_vectors(const std::string &file, const std::vector<double> &U,
								const std::vector<double> &V)
	{
		std::ofstream fout;
		fout.open(file);
		double scale = 1000;
		double width1 = 5;
		double width2 = 3;
		double displace = 30;
		std::string color = "#000000";
		std::string color1 = "#7AA20D";
		fout << "<svg>" << std::endl;
		std::vector<double> xlist = {{0, 0, 1, 1}};
		std::vector<double> ylist = {{0, 1, 1, 0}};
		fout << "<polyline points=\"";
		for (int i = 0; i < 4; i++)
		{
			fout << displace + xlist[i] * scale << "," << displace + ylist[i] * scale << " ";
		}
		fout << displace + xlist[0] * scale << "," << displace + ylist[0] * scale;
		fout << "\" fill = \"white\" stroke = \"" + color + "\" stroke-width = \"" + std::to_string(5) + "\" /> " << std::endl;

		for (int i = 0; i < U.size(); i++)
		{
			double x1 = displace + U[i] * scale;
			double y1 = displace + 0 * scale;
			double x2 = displace + U[i] * scale;
			double y2 = displace + 1 * scale;
			fout << "<line x1=\"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) +
						"\" x2=\"" + std::to_string(x2) + "\" y2=\"" + std::to_string(y2) +
						"\" stroke = \"" + color1 + "\" stroke-width = \"" + std::to_string(width2) + "\" /> "
				 << std::endl;
		}
		for (int i = 0; i < V.size(); i++)
		{
			double x1 = displace + 0 * scale;
			double y1 = displace + V[i] * scale;
			double x2 = displace + 1 * scale;
			double y2 = displace + V[i] * scale;
			fout << "<line x1=\"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) +
						"\" x2=\"" + std::to_string(x2) + "\" y2=\"" + std::to_string(y2) +
						"\" stroke = \"" + color1 + "\" stroke-width = \"" + std::to_string(width2) + "\" /> "
				 << std::endl;
		}

		fout << "</svg>" << std::endl;
		fout.close();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void run_ours(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
				  const double per, const bool enable_local_energy = false)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;
		igl::Timer timer;
		double time_knot = 0;
		double time_solve = 0;
		double precision = 0;

		Eigen::MatrixXd ver;
		int nbr = nbr_pts; // nbr of points
		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		int method = model;
		bool corners = true;
		get_model_sample_points(nbr, ver, F, param, method, corners, path);
		// std::cout << "ver\n" << ver << std::endl;
		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		int perturb_itr = 0;
		// double per_ours = 0.9;
		// double per = 0.2;
		int target_steps = 10;
		bool enable_max_fix_nbr = true;
		igl::opengl::glfw::Viewer viewer;
		if (0)
		{
			viewer.data().add_points(ver, ecolor);
			viewer.launch();
			exit(0);
		}
		timer.start();
		std::cout << "before generating knot vectors" << std::endl;
		std::cout << "data size " << ver.rows() << std::endl;
		Bsurface surface;
		surface.generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);

		timer.stop();
		time_knot = timer.getElapsedTimeInSec();

		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		if (1)
		{
			surface.degree1 = 3;
			surface.degree2 = 3;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "before initialize the basis " << std::endl;
			PartialBasis basis(surface);
			std::cout << "initialize the basis done" << std::endl;
			std::cout << "before solving control points" << std::endl;
			timer.start();
			surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
			timer.stop();
			time_solve = timer.getElapsedTimeInSec();
			surface.surface_visulization(surface, 100, SPs, SFs);
			if (enable_local_energy)
			{
				double timeitr = 0;
				for (int i = 0; i < 50; i++)
				{
					timer.start();
					Eigen::MatrixXd energy, euu, evv, euv;
					energy = surface.surface_energy_calculation(surface, basis, 1, euu, evv, euv);
					bool uorv;
					int which;
					double max_energy;
					surface.detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
					std::vector<double> Unew = surface.U;
					std::vector<double> Vnew = surface.V;
					if (!uorv)
					{ // u get updated
						double value = (Unew[which] + Unew[which + 1]) / 2;
						surface.U = knot_vector_insert_one_value(Unew, value);
					}
					else
					{
						double value = (Vnew[which] + Vnew[which + 1]) / 2;
						surface.V = knot_vector_insert_one_value(Vnew, value);
					}
					std::cout << "knot vector get inserted" << std::endl;
					basis.clear();
					basis.init(surface);
					surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
					timer.stop();
					timeitr += timer.getElapsedTimeInSec();
					std::cout << " control points solved" << std::endl;
					surface.surface_visulization(surface, 100, SPs, SFs);

					igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);
					std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time_itr", "max_energy"}};
					int cps = (surface.nu() + 1) * (surface.nv() + 1);
					precision = surface.max_interpolation_err(ver, param, surface);
					std::vector<double> data = {{time_knot, time_solve, precision,
												 double(surface.nu()), double(surface.nv()), double(cps), timeitr, max_energy}};
					write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv",
							  titles, data);
				}
				return;
			}
		}

		std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		precision = surface.max_interpolation_err(ver, param, surface);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
		bool write_file = true;
		if (write_file)
		{
			Eigen::MatrixXi Pf;
			write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
			igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);

			std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time"}};
			int cps = (surface.nu() + 1) * (surface.nv() + 1);
			std::vector<double> data = {{time_knot, time_solve, precision, double(surface.nu()), double(surface.nv()), double(cps),
										 time_knot + time_solve}};
			write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".csv",
					  titles, data);
		}
		output_timing();
		if (1)
		{ // write svg files
			write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "param.svg", param);
			write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
		}
		std::cout << "total time " << time_knot + time_solve << std::endl;

		viewer.data().clear();
		viewer.data().set_mesh(SPs, SFs);
		viewer.data().add_points(ver, ecolor);
		Eigen::MatrixXd p0(1, 3), p1(1, 3);
	}

	// algorithm = 0, ours;
	// algorithm = 1, lofting;
	// algorithm = 2, averaging;

	void run_other_sampling_strategies(const int model, double &per_ours, const std::string path, const std::string tail,
									   const double per, const int algorithm)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;
		igl::Timer timer;
		double time_knot = 0;
		double time_solve = 0;
		double precision = 0;
		int nbr;
		Eigen::MatrixXd ver;

		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		// get_mesh_vertices_and_parametrization(ver, F, param);
		int method = model;
		int sampling_strategy = 3;
		bool enable_local_energy = false;
		///////////////////////////////
		// sampling strategies
		if (sampling_strategy == -1)
		{ // use the random sampling method
			nbr = 100;
			get_model_sample_points(nbr, ver, F, param, method, true, path);
		}
		if (sampling_strategy == 0)
		{ // sample gridded data
			int nbr_row = 9;
			int nbr_col = 30;
			nbr = nbr_row * nbr_col; // nbr of points
			param.resize(nbr, 2);
			int current = 0;
			for (int i = 0; i < nbr_row; i++)
			{
				for (int j = 0; j < nbr_col; j++)
				{
					param(current, 0) = double(1) / (nbr_row - 1) * i;
					param(current, 1) = double(1) / (nbr_col - 1) * j;
					current++;
				}
			}
		}
		if (sampling_strategy == 1)
		{ // sample non-uniform density
			srand(1);
			int nbr_other = 150;
			int nbr_local = 100;
			double xmin = 0.3, xmax = 0.5;
			double ymin = 0.3, ymax = 0.5;
			nbr = nbr_other + nbr_local;
			param.resize(nbr, 2);
			int current = 0;
			for (int i = 0; i < INFINITY; i++)
			{
				Vector3d para3d = Vector3d::Random();
				double u = (1 + para3d[0]) / 2;
				double v = (1 + para3d[1]) / 2;
				if (u >= xmin && u <= xmax && v >= ymin && v <= ymax)
				{
					continue;
				}
				param.row(current) << u, v;
				current++;
				if (current == nbr_other)
				{
					break;
				}
			}
			int c1 = 0;
			for (int i = 0; i < INFINITY; i++)
			{
				Vector3d para3d = Vector3d::Random();
				double u = (1 + para3d[0]) / 2;
				double v = (1 + para3d[1]) / 2;
				if (u < xmin || u > xmax || v < ymin || v > ymax)
				{
					continue;
				}
				param.row(current + c1) << u, v;
				c1++;
				if (c1 == nbr_local)
				{
					break;
				}
			}
		}
		if (sampling_strategy == 2)
		{ // sample row data
			int nbr_row = 6;
			int nbr_col = 20;
			nbr = nbr_row * nbr_col; // nbr of points
			param.resize(nbr, 2);
			int current = 0;
			for (int i = 0; i < nbr_row; i++)
			{
				for (int j = 0; j < nbr_col; j++)
				{
					param(current, 0) = double(1) / (nbr_row - 1) * i;
					Vector3d spd = Vector3d::Random();
					param(current, 1) = (spd[0] + 1) / 2;
					current++;
				}
			}
		}
		if (sampling_strategy == 3)
		{							 // sample not close data
			double para_noise = 0.7; // noise level of the parameters
			int nbr_row = 10;
			int nbr_col = 10;
			double row_dis = double(1) / (nbr_row);
			double col_dis = double(1) / (nbr_col);
			nbr = (nbr_row) * (nbr_col); // nbr of points
			param.resize(nbr, 2);
			int current = 0;
			for (int i = 0; i < nbr_row; i++)
			{
				for (int j = 0; j < nbr_col; j++)
				{
					double upo = row_dis * (i + 0.5);
					double vpo = col_dis * (j + 0.5);
					Vector3d noise = Vector3d::Random();
					param(current, 0) = upo + row_dis * 0.5 * noise[0] * para_noise;
					param(current, 1) = vpo + col_dis * 0.5 * noise[1] * para_noise;
					current++;
				}
			}
			// std::cout << "para\n" << param << std::endl;
		}

		// std::cout << "parameters\n" << param << std::endl;
		//////////////////////////////
		if (sampling_strategy >= 0)
		{
			input_parameters_get_model_sample_points(ver, F, param, method);
		}

		// std::cout << "ver\n" << ver << std::endl;
		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		int perturb_itr = 0;
		// double per_ours = 0.9;
		// double per = 0.2;
		int target_steps = 10;
		bool enable_max_fix_nbr = true;
		igl::opengl::glfw::Viewer viewer;
		if (0)
		{

			viewer.data().add_points(ver, ecolor);
			viewer.launch();
			exit(0);
		}
		std::string algo_str;
		std::cout << "before generating knot vectors" << std::endl;
		std::cout << "data size " << ver.rows() << std::endl;
		timer.start();
		Bsurface surface;
		if (algorithm == 0)
		{
			algo_str = "ours";
			surface.generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);
		}
		if (algorithm == 1)
		{
			algo_str = "lofting";
			lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, per);
		}
		if (algorithm == 2)
		{
			algo_str = "pigel";
			piegl_method_generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per);
		}

		timer.stop();
		time_knot = timer.getElapsedTimeInSec();

		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		if (1)
		{
			surface.degree1 = 3;
			surface.degree2 = 3;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "before initialize the basis " << std::endl;
			PartialBasis basis(surface);
			std::cout << "initialize the basis done" << std::endl;
			std::cout << "before solving control points" << std::endl;
			timer.start();
			surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
			timer.stop();
			time_solve = timer.getElapsedTimeInSec();
			surface.surface_visulization(surface, 100, SPs, SFs);
			if (enable_local_energy)
			{
				double timeitr = 0;
				for (int i = 0; i < 25; i++)
				{
					timer.start();
					Eigen::MatrixXd energy, euu, evv, euv;
					energy = surface.surface_energy_calculation(surface, basis, 1, euu, evv, euv);
					bool uorv;
					int which;
					double max_energy;
					surface.detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
					std::vector<double> Unew = surface.U;
					std::vector<double> Vnew = surface.V;
					if (!uorv)
					{ // u get updated
						double value = (Unew[which] + Unew[which + 1]) / 2;
						surface.U = knot_vector_insert_one_value(Unew, value);
					}
					else
					{
						double value = (Vnew[which] + Vnew[which + 1]) / 2;
						surface.V = knot_vector_insert_one_value(Vnew, value);
					}
					std::cout << "knot vector get inserted" << std::endl;
					basis.clear();
					basis.init(surface);
					surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
					timer.stop();
					timeitr += timer.getElapsedTimeInSec();
					std::cout << " control points solved" << std::endl;
					surface.surface_visulization(surface, 100, SPs, SFs);

					igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);
					std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time_itr", "max_energy"}};
					int cps = (surface.nu() + 1) * (surface.nv() + 1);
					precision = surface.max_interpolation_err(ver, param, surface);
					std::vector<double> data = {{time_knot, time_solve, precision,
												 double(surface.nu()), double(surface.nv()), double(cps), timeitr, max_energy}};
					write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv",
							  titles, data);
				}
				return;
			}
		}

		std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		precision = surface.max_interpolation_err(ver, param, surface);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
		bool write_file = true;
		if (write_file)
		{
			Eigen::MatrixXi Pf;
			write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
			igl::write_triangle_mesh(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);

			std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time"}};
			int cps = (surface.nu() + 1) * (surface.nv() + 1);
			std::vector<double> data = {{time_knot, time_solve, precision, double(surface.nu()), double(surface.nv()), double(cps),
										 time_knot + time_solve}};
			write_csv(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".csv",
					  titles, data);
		}
		output_timing();
		if (1)
		{ // write svg files
			write_svg_pts(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "param.svg", param);
			write_svg_knot_vectors(path + algo_str + "_p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
		}
		std::cout << "total time " << time_knot + time_solve << std::endl;
	}
	void ours_inserting_redundant(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
								  const double per, const bool enable_insertion = false)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;
		igl::Timer timer;
		double time_knot = 0;
		double time_solve = 0;
		double precision = 0;

		Eigen::MatrixXd ver;
		int nbr = nbr_pts; // nbr of points
		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		// get_mesh_vertices_and_parametrization(ver, F, param);
		int method = model;
		bool corners = false;
		get_model_sample_points(nbr, ver, F, param, method, corners, path);
		// std::cout << "ver\n" << ver << std::endl;
		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		int perturb_itr = 0;
		// double per_ours = 0.9;
		// double per = 0.2;
		int target_steps = 10;
		bool enable_max_fix_nbr = true;
		igl::opengl::glfw::Viewer viewer;
		if (0)
		{
			viewer.data().add_points(ver, ecolor);
			viewer.launch();
			exit(0);
		}
		timer.start();
		std::cout << "before generating knot vectors" << std::endl;
		std::cout << "data size " << ver.rows() << std::endl;
		Bsurface surface;
		surface.generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);
		timer.stop();
		time_knot = timer.getElapsedTimeInSec();

		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		if (1)
		{
			surface.degree1 = 3;
			surface.degree2 = 3;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "before initialize the basis " << std::endl;
			PartialBasis basis(surface);
			std::cout << "initialize the basis done" << std::endl;
			std::cout << "before solving control points" << std::endl;
			// timer.start();
			// solve_control_points_for_fairing_surface(surface, param, ver, basis);
			////timer.stop();
			// time_solve = timer.getElapsedTimeInSec();
			// surface_visulization(surface, 100, SPs, SFs);
			bool local_energy = true;
			if (enable_insertion)
			{
				/*double minvalue = 0.2; or 0.5, 0.015, 20
				double maxvalue = minvalue + 0.015;
				int nbrinsert = 30;*/
				// best settings

				// we first insert knots to create problems
				double minvalue = 0.4;
				double maxvalue = minvalue + 0.1;
				int nbrinsert = 20;
				double interval = (maxvalue - minvalue) / nbrinsert;

				for (int i = 0; i < nbrinsert; i++)
				{
					std::vector<double> Unew = surface.U;
					std::vector<double> Vnew = surface.V;
					double value = minvalue + interval * i;
					surface.V = knot_vector_insert_one_value(Vnew, value);
				}
				timer.start();
				basis.clear();
				basis.init(surface);
				surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
				timer.stop();
				time_solve = timer.getElapsedTimeInSec();
				surface.surface_visulization(surface, 100, SPs, SFs);
				if (local_energy)
				{
					double timeitr = 0;
					for (int i = 0; i < 100; i++)
					{
						timer.start();
						Eigen::MatrixXd energy, euu, evv, euv;
						energy = surface.surface_energy_calculation(surface, basis, 1, euu, evv, euv);
						bool uorv;
						int which;
						double max_energy;
						surface.detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
						std::vector<double> Unew = surface.U;
						std::vector<double> Vnew = surface.V;
						if (!uorv)
						{ // u get updated
							double value = (Unew[which] + Unew[which + 1]) / 2;
							surface.U = knot_vector_insert_one_value(Unew, value);
							std::cout << "##insert u = " << value << std::endl;
						}
						else
						{
							double value = (Vnew[which] + Vnew[which + 1]) / 2;
							surface.V = knot_vector_insert_one_value(Vnew, value);
							std::cout << "##insert v = " << value << std::endl;
						}
						std::cout << "knot vector get inserted" << std::endl;
						basis.clear();
						basis.init(surface);
						surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
						timer.stop();
						timeitr += timer.getElapsedTimeInSec();
						std::cout << " control points solved" << std::endl;
						surface.surface_visulization(surface, 100, SPs, SFs);

						igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);
						std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time_itr", "max_energy"}};
						int cps = (surface.nu() + 1) * (surface.nv() + 1);
						precision = surface.max_interpolation_err(ver, param, surface);
						std::vector<double> data = {{time_knot, time_solve, precision,
													 double(surface.nu()), double(surface.nv()), double(cps), timeitr, max_energy}};
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
		precision = surface.max_interpolation_err(ver, param, surface);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
		bool write_file = true;
		if (write_file)
		{
			Eigen::MatrixXi Pf;
			write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
			igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);

			std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time"}};
			int cps = (surface.nu() + 1) * (surface.nv() + 1);
			std::vector<double> data = {{time_knot, time_solve, precision, double(surface.nu()), double(surface.nv()), double(cps),
										 time_knot + time_solve}};
			write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + ".csv",
					  titles, data);
		}
		output_timing();
		if (1)
		{ // write svg files
			write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "param.svg", param);
			write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
		}
		std::cout << "total time " << time_knot + time_solve << std::endl;

		viewer.data().clear();
		viewer.data().set_mesh(SPs, SFs);
		viewer.data().add_points(ver, ecolor);
		Eigen::MatrixXd p0(1, 3), p1(1, 3);
	}
	void run_Seungyong(const int model, const int nbr_pts, const double tolerance, const std::string path)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;
		igl::Timer timer;
		double time_run = 0;
		Eigen::MatrixXd ver;
		int nbr = nbr_pts; // nbr of points
		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		// get_mesh_vertices_and_parametrization(ver, F, param);
		int method = model;
		bool corners = false;
		get_model_sample_points(nbr, ver, F, param, method, corners, path);

		int degree1 = 3;
		int degree2 = 3;
		double fair_parameter = 1e-7;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		igl::opengl::glfw::Viewer viewer;

		std::vector<Bsurface> surfaces;
		timer.start();
		iteratively_approximate_method(degree1, degree2, Uknot, Vknot, param, ver, tolerance, surfaces, fair_parameter);
		timer.stop();
		time_run = timer.getElapsedTimeInSec();
		std::cout << "run_Seungyong finished, surface levels " << surfaces.size() << std::endl;
		std::cout << "final surface size " << surfaces.back().nu() + 1 << " " << surfaces.back().nv() + 1
				  << " total control points " << (surfaces.back().nu() + 1) * (surfaces.back().nv() + 1) << std::endl;
		std::cout << "Seungyong time " << time_run << std::endl;
		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;

		write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj", ver);
		for (int i = 0; i < surfaces.size(); i++)
		{
			surfaces[i].surface_visulization(surfaces[i], 100, SPs, SFs);
			igl::write_triangle_mesh(path + "Seungyong_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + "_level_" +
										 std::to_string(i) + ".obj",
									 SPs, SFs);
		}
	}

	void run_piegl(const int model, const int nbr_pts, const double per = 0.2)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;

		igl::Timer timer;
		double time_total = 0;
		Eigen::MatrixXd ver;
		int nbr = nbr_pts; // nbr of points
		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		// get_mesh_vertices_and_parametrization(ver, F, param);
		int method = model;
		bool corners = false;
		get_model_sample_points(nbr, ver, F, param, method, corners, "./");

		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		int perturb_itr = 0;

		igl::opengl::glfw::Viewer viewer;
		if (0)
		{
			viewer.data().add_points(ver, ecolor);
			viewer.launch();
			exit(0);
		}
		timer.start();
		piegl_method_generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per);

		Eigen::MatrixXd paramout(ver.rows(), 3), zeros(ver.rows(), 1);
		zeros = Eigen::MatrixXd::Constant(ver.rows(), 1, 0);

		Bsurface surface;
		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		if (1)
		{
			surface.degree1 = 3;
			surface.degree2 = 3;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "before initialize the basis " << std::endl;
			PartialBasis basis(surface);
			std::cout << "initialize the basis done" << std::endl;
			std::cout << "before solving control points" << std::endl;
			surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
			surface.surface_visulization(surface, 100, SPs, SFs);
		}
		timer.stop();
		time_total = timer.getElapsedTimeInSec();

		std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
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

		std::cout << "final surface size " << surface.nu() + 1 << " " << surface.nv() + 1
				  << " total control points " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;

		viewer.data().clear();
		viewer.data().set_mesh(SPs, SFs);
		// viewer.data().add_points(ver, ecolor);
		Eigen::MatrixXd p0(1, 3), p1(1, 3);
		p0.row(0) = surface.BSplineSurfacePoint(surface, 0, 0);
		p1.row(0) = surface.BSplineSurfacePoint(surface, 0, 1);
		// std::cout << "p0\n" << p0 << std::endl;
		viewer.data().add_points(p0, red);
		viewer.data().add_points(p1, green);
		viewer.launch();
	}
	void run_lofting(const int model, const int nbr_pts, const double per, const std::string path)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;

		igl::Timer timer;
		double time_total = 0;
		Eigen::MatrixXd ver;
		int nbr = nbr_pts; // nbr of points
		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		// get_mesh_vertices_and_parametrization(ver, F, param);
		int method = model;
		bool corners = false;
		get_model_sample_points(nbr, ver, F, param, method, corners, "./");

		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		int perturb_itr = 0;

		igl::opengl::glfw::Viewer viewer;
		if (0)
		{
			viewer.data().add_points(ver, ecolor);
			viewer.launch();
			exit(0);
		}
		timer.start();
		// piegl_method_generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per);
		lofting_method_generate_interpolation_knot_vectors(false, degree1, degree2, Uknot, Vknot, param, per);

		Eigen::MatrixXd paramout(ver.rows(), 3), zeros(ver.rows(), 1);
		zeros = Eigen::MatrixXd::Constant(ver.rows(), 1, 0);

		Bsurface surface;
		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		if (1)
		{
			surface.degree1 = 3;
			surface.degree2 = 3;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "before initialize the basis " << std::endl;
			PartialBasis basis(surface);
			std::cout << "initialize the basis done" << std::endl;
			std::cout << "before solving control points" << std::endl;
			surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
			surface.surface_visulization(surface, 100, SPs, SFs);
		}
		timer.stop();
		time_total = timer.getElapsedTimeInSec();

		std::cout << "nu and nv " << surface.nu() << " " << surface.nv() << std::endl;
		std::cout << "#cp " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
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
				  << " total control points " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;

		viewer.data().clear();
		viewer.data().set_mesh(SPs, SFs);
		// viewer.data().add_points(ver, ecolor);
		Eigen::MatrixXd p0(1, 3), p1(1, 3);
		p0.row(0) = surface.BSplineSurfacePoint(surface, 0, 0);
		p1.row(0) = surface.BSplineSurfacePoint(surface, 0, 1);
		// std::cout << "p0\n" << p0 << std::endl;
		viewer.data().add_points(p0, red);
		viewer.data().add_points(p1, green);
		viewer.launch();
	}
	void read_mesh_series(std::string path, std::string namebase, int end)
	{
		std::string file0 = path + namebase + std::to_string(0) + ".obj";
		Eigen::MatrixXd ver;
		Eigen::MatrixXi f;
		igl::read_triangle_mesh(file0, ver, f);
		for (int i = 1; i < end + 1; i++)
		{ // 0 has already been included
			std::cout << "reading " << i << std::endl;
			std::string file = path + namebase + std::to_string(i) + ".obj";
			Eigen::MatrixXd ver_t;
			Eigen::MatrixXi f_t;
			igl::read_triangle_mesh(file, ver_t, f_t);
			ver = ver + ver_t;
		}
		file0 = path + namebase + "sum.obj";
		igl::write_triangle_mesh(file0, ver, f);
	}
	void run_mesh_reconstruction(const std::string inpath, const std::string modelname, double &per_ours, const std::string path, const std::string tail,
								 const double per, const bool enable_local_energy = false)
	{
		Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3), pcolor(1, 3), red(1, 3), green(1, 3), blue(1, 3);
		fcolor << 1, 0, 0;
		ecolor << 0.9, 0.9, 0.9;
		;
		pcolor << 0, 0.9, 0.5;
		red << 1, 0, 0;
		green << 0, 1, 0;
		blue << 0, 0, 1;
		igl::Timer timer;
		double time_knot = 0;
		double time_solve = 0;
		double precision = 0;

		Eigen::MatrixXd ver;

		Eigen::MatrixXi F;
		Eigen::MatrixXd param, paramout;
		// std::string modelname = "mask3kf.obj";
		// std::string modelname = "simplified_2500_TigerMask4.obj";
		std::string meshfile = inpath + modelname;
		//"D:\\vs\\sparse_data_interpolation\\meshes\\meshmodels\\" + modelname;

		// get_mesh_vertices_and_parametrization(ver, F, param);
		bool corners = false;
		mesh_parameterization(meshfile, ver, param, F);
		paramout.resize(param.rows(), 3);
		Eigen::VectorXd param_zero = Eigen::VectorXd::Zero(param.rows());
		paramout << param, param_zero;
		igl::write_triangle_mesh(path + "param_" + modelname + tail, paramout, F);

		// get_model_sample_points(nbr, ver, F, param, method, corners);
		// std::cout << "ver\n" << ver << std::endl;
		int nbr = param.rows();
		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
		std::vector<double> Vknot = Uknot;
		int target_steps = 10;
		bool enable_max_fix_nbr = true;
		igl::opengl::glfw::Viewer viewer;
		if (0)
		{
			viewer.data().add_points(ver, ecolor);
			viewer.launch();
			exit(0);
		}
		timer.start();
		std::cout << "before generating knot vectors" << std::endl;
		std::cout << "data size " << ver.rows() << std::endl;
		Bsurface surface;
		surface.generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours, per, target_steps, enable_max_fix_nbr);

		timer.stop();
		time_knot = timer.getElapsedTimeInSec();

		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		int visual_nbr = 200;
		if (1)
		{
			surface.degree1 = 3;
			surface.degree2 = 3;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "before initialize the basis " << std::endl;
			PartialBasis basis(surface);
			std::cout << "initialize the basis done" << std::endl;
			std::cout << "before solving control points" << std::endl;
			timer.start();
			surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
			timer.stop();
			time_solve = timer.getElapsedTimeInSec();
			surface.surface_visulization(surface, visual_nbr, SPs, SFs);
			if (enable_local_energy)
			{
				double timeitr = 0;
				for (int i = 0; i < 20; i++)
				{
					timer.start();
					Eigen::MatrixXd energy, euu, evv, euv;
					energy = surface.surface_energy_calculation(surface, basis, 1, euu, evv, euv);
					bool uorv;
					int which;
					double max_energy;
					surface.detect_max_energy_interval(surface, energy, euu, evv, uorv, which, max_energy);
					std::vector<double> Unew = surface.U;
					std::vector<double> Vnew = surface.V;
					if (!uorv)
					{ // u get updated
						double value = (Unew[which] + Unew[which + 1]) / 2;
						surface.U = knot_vector_insert_one_value(Unew, value);
					}
					else
					{
						double value = (Vnew[which] + Vnew[which + 1]) / 2;
						surface.V = knot_vector_insert_one_value(Vnew, value);
					}
					std::cout << "knot vector get inserted" << std::endl;
					basis.clear();
					basis.init(surface);
					surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
					timer.stop();
					timeitr += timer.getElapsedTimeInSec();
					std::cout << " control points solved" << std::endl;
					surface.surface_visulization(surface, visual_nbr, SPs, SFs);

					igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + modelname + tail, SPs, SFs);
					std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time_itr", "max_energy"}};
					int cps = (surface.nu() + 1) * (surface.nv() + 1);
					precision = surface.max_interpolation_err(ver, param, surface);
					std::vector<double> data = {{time_knot, time_solve, precision,
												 double(surface.nu()), double(surface.nv()), double(cps), timeitr, max_energy}};
					write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" + std::to_string(i) + "_m_" + modelname + tail + ".csv",
							  titles, data);
				}
				return;
			}
		}

		std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		precision = surface.max_interpolation_err(ver, param, surface);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
		bool write_file = true;
		if (write_file)
		{
			Eigen::MatrixXi Pf;
			write_points(path + "pts" + std::to_string(nbr) + "_m_" + modelname, ver);
			igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail, SPs, SFs);

			write_control_pts(surface.control_points, path + "cps" + std::to_string(nbr) + "_m_" + modelname);
			std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv", "cps"}};
			int cps = (surface.nu() + 1) * (surface.nv() + 1);
			std::vector<double> data = {{time_knot, time_solve, precision, double(surface.nu()), double(surface.nv()), double(cps)}};
			write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail + ".csv",
					  titles, data);
		}
		output_timing();
		if (1)
		{ // write svg files
			write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail + "param.svg", param);
			write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + modelname + tail + "knots.svg", surface.U, surface.V);
		}
		std::cout << "total time " << time_knot + time_solve << std::endl;

		viewer.data().clear();
		viewer.data().set_mesh(SPs, SFs);
		viewer.data().add_points(ver, ecolor);
		Eigen::MatrixXd p0(1, 3), p1(1, 3);
		// p0.row(0) = BSplineSurfacePoint(surface, 0, 0);
		// p1.row(0) = BSplineSurfacePoint(surface, 0, 1);
		////std::cout << "p0\n" << p0 << std::endl;
		// viewer.data().add_points(p0, red);
		// viewer.data().add_points(p1, green);
		// viewer.launch();
	}
	void run_pia(const int model, const int nbr_pts, const int max_itr, const double threadshold,
				 const int cp_nbr_sqrt, const std::string &path)
	{
		igl::Timer timer;
		double time_total = 0;
		Eigen::MatrixXd ver;
		int nbr = nbr_pts; // nbr of points
		Eigen::MatrixXi F;
		Eigen::MatrixXd param;
		// get_mesh_vertices_and_parametrization(ver, F, param);
		int method = model;
		bool corners = false;
		get_model_sample_points(nbr, ver, F, param, method, corners, "./");

		int degree1 = 3;
		int degree2 = 3;
		std::vector<double> Uknot;
		std::vector<double> Vknot;
		int nu = cp_nbr_sqrt - 1;
		int nv = cp_nbr_sqrt - 1;
		std::cout << "starting pia method" << std::endl;
		timer.start();
		comparison::initialize_pia_knot_vectors(degree1, degree2, Uknot, Vknot, nu, nv);
		std::cout << "pia knot vectors get initialized" << std::endl;
		print_vector(Uknot);
		print_vector(Vknot);
		Bsurface surface;
		Eigen::MatrixXd SPs;
		Eigen::MatrixXi SFs;
		if (1)
		{
			surface.degree1 = degree1;
			surface.degree2 = degree2;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "nu " << surface.nu() << std::endl;
			std::cout << "nv " << surface.nv() << std::endl;
			progressive_iterative_approximation(surface, param, ver, max_itr, threadshold);
			surface.surface_visulization(surface, 100, SPs, SFs);
		}
		timer.stop();
		time_total = timer.getElapsedTimeInSec();

		std::cout << "nu and nv " << surface.nu() << " " << surface.nv() << std::endl;
		std::cout << "#cp " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
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
				  << " total control points " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;
	}

	void run_pia_mesh_reconstruct(const std::string meshfile, const int max_itr, const double threadshold,
								  const int cp_nbr_sqrt, const std::string &path)
	{
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
		if (1)
		{
			surface.degree1 = degree1;
			surface.degree2 = degree2;
			surface.U = Uknot;
			surface.V = Vknot;
			std::cout << "nu " << surface.nu() << std::endl;
			std::cout << "nv " << surface.nv() << std::endl;
			progressive_iterative_approximation(surface, param, ver, max_itr, threadshold);
			surface.surface_visulization(surface, 100, SPs, SFs);
		}
		timer.stop();
		time_total = timer.getElapsedTimeInSec();

		std::cout << "nu and nv " << surface.nu() << " " << surface.nv() << std::endl;
		std::cout << "#cp " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;
		print_vector(surface.U);
		print_vector(surface.V);
		std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;
		std::cout << "total time " << time_total << std::endl;
		bool write_file = true;
		if (write_file)
		{
			Eigen::MatrixXi Pf;

			igl::write_triangle_mesh(meshfile + "_pia_cp_" + std::to_string(cp_nbr_sqrt * cp_nbr_sqrt) + ".obj", SPs, SFs);
		}
		output_timing();

		std::cout << "final surface size " << surface.nu() + 1 << " " << surface.nv() + 1
				  << " total control points " << (surface.nu() + 1) * (surface.nv() + 1) << std::endl;
	}
	
	void run_temperature_fitting()
	{
		std::vector<std::vector<double>> stations, border;
		read_csv_data_lbl(example_root_path + "ShanxiStations.csv", stations);
		read_csv_data_lbl(example_root_path + "ShanxiBorder.csv", border);
		
		std::cout << "there are " << stations.size() << " stations\n";
		std::vector<int> stationIds(stations.size());
		std::vector<std::vector<double>> temps(stations.size()), tempsRaw;
		read_csv_data_lbl(example_root_path + "ShanxiTemperatures.csv", tempsRaw);
		std::cout<<"there are "<<tempsRaw.size()<<" rows of temp data, on avg each station has "<<double(tempsRaw.size()) / stations.size()<<" data\n";
		std::cout<<"the cols "<<tempsRaw[0].size()<<"\n";
		for(int i=0;i<10;i++)
		{
			for(auto ele : tempsRaw[i])
			{
				std::cout<<ele<<",";
			}
			std::cout<<"\n";
		}
		Eigen::MatrixXd param(stations.size(),2), paramb(border.size(), 3);
		Eigen::Vector2d cmin(stations[0][1], stations[0][2]), cmax(stations[0][1], stations[0][2]);
		for (int i=0;i<stations.size();i++)
		{
			stationIds[i] = stations[i][0];
			param(i,0) = stations[i][1];
			param(i,1) = stations[i][2];
			
		}
		for (int i = 0; i < border.size(); i++)
		{
			paramb(i, 0) = border[i][1]; // note that the x and y for border is swapped
			paramb(i, 1) = border[i][0];
			paramb(i, 2) = 0;
			if (paramb(i, 0) < cmin[0])
			{
				cmin[0] = paramb(i, 0);
			}
			if (paramb(i, 1) < cmin[1])
			{
				cmin[1] = paramb(i, 1);
			}
			if (paramb(i, 0) > cmax[0])
			{
				cmax[0] = paramb(i, 0);
			}
			if (paramb(i, 1) > cmax[1])
			{
				cmax[1] = paramb(i, 1);
			}
		}
		std::cout<<"min and max of the scene "<<cmin.transpose()<<" ,,, "<<cmax.transpose()<<"\n";
		// rescale the model into [0,1]x[0,1]
		double length = std::max(cmax[0]-cmin[0], cmax[1]-cmin[1]); // make the longer side fit in [0,1]
		Eigen::Vector2d offset(0.05,0.05); // the model not too close to the boundary
		// rescale the parameters
		for (int i = 0; i < param.rows(); i++)
		{
			param.row(i) = (Eigen::Vector2d(param.row(i)) - cmin) * 0.9 / length + offset;
		}
		// for (int i = 0; i < paramb.rows(); i++)
		// {
		// 	paramb.row(i) = (Eigen::Vector2d(paramb.row(i)) - cmin) * 0.9 / length + offset; 
		// }
		// std::cout<<"param\n"<<param<<"\nparamb\n"<<paramb<<"\n";

		// collect the data list
		for (int i=0;i<tempsRaw.size();i++)
		{
			int stid = tempsRaw[i][0]; // the station id
			int whichId = -1;
			for(int j=0;j<stationIds.size();j++)
			{
				if(stid == stationIds[j])
				{
					whichId = j;
					break;
				}
			}
			if (whichId >= 0)
			{
				temps[whichId].push_back(tempsRaw[i].back());
				if(i<10)
				{
					std::cout<<"print temp, "<<temps[whichId].back()<<"\n";
				}
			}
		}
		for(int i=0;i<temps.size()-1;i++)
		{
			if(temps[i].size()!=temps[i+1].size())
			{
				std::cout<<"error here 1 "<<temps[i].size()<<" "<<temps[i+1].size()<<"\n";
				exit(0);
			}
		}
		// there are 3 knot vectors: U, V, and T.
		// compute the U and V knot vectors:
		int degree = 3;
		Bsurface surface;
		int nbr = param.rows();					// the number of data points
		surface.degree1 = degree;					// degree of u direction
		surface.degree2 = degree;					// degree of v direction
		surface.U = {{0, 0, 0, 0, 1, 1, 1, 1}}; // the initial U knot vector
		surface.V = surface.U;					// the initial V knot vector
		int target_steps = 10;					// the number of iterations for constructing lists $L$.
		bool enable_max_fix_nbr = true;			// progressively update the knot vectors to make the two knot vectors balanced in length.
		double delta = 0.9;						// the parameter to improve the solving stability
		double per = 0.5;						// the parameter inherited from [Wen-Ke Wang et al, 2008, CAD]
		
		// generate knot vectors to make sure the data points can be interpolated
		surface.generate_interpolation_knot_vectors(surface.degree1, surface.degree2, surface.U, surface.V, param, delta, per, target_steps, enable_max_fix_nbr);
		int refinement = std::max(surface.nu() + 1, surface.nv() + 1) * 0.5;
		std::cout<<"Refine the knot vectors by adding "<<refinement<<" knots.\n";
		// refine the shape 
		surface.RefineKnots(refinement);
		std::cout << "knot vectors generated" << std::endl;
		PartialBasis basis(surface);
		int nbrTime = temps[0].size(); // the number of time steps
		std::vector<std::vector<std::vector<Vector3d>>> CPs(nbrTime), CPAll; // the control points each layer, and the final control points
		for (int i = 0; i < nbrTime; i++)
		{
			// clean the data, remove suspecious data
			std::vector<Vector3d> datapoints;
			std::vector<Vector2d> dataparams;
			for(int j=0;j<temps.size();j++)
			{
				double temperature = temps[j][i];
				if (temperature > 100 || temperature < -100) 
				{
					continue;
				}
				datapoints.push_back(Vector3d(param(j,0), param(j,1), temperature));
				dataparams.push_back(Vector2d(param(j,0), param(j,1)));
			}
			Eigen::MatrixXd ver = vec_list_to_matrix(datapoints);
			Eigen::MatrixXd par(dataparams.size(),2);
			for (int j = 0; j < dataparams.size(); j++)
			{
				par.row(j) = dataparams[j];
			}
			// compute the control points of this layer
			surface.solve_control_points_for_fairing_surface(surface, par, ver, basis);
			CPs[i] = surface.control_points;
		}
		int ncpi = CPs[0].size();
		int ncpj = CPs[0][0].size();
		std::vector<double> T, parT;
		// construct the third knot vector T using even parametrization
		for(int i=0;i<nbrTime;i++)
		{
			parT.push_back(i*1.0 / (nbrTime - 1));
		}
		// make sure the last parameter is precisely 1
		parT.back()=1;

		// generate the knot vector T
		T = {{0, 0, 0, 0, 1, 1, 1, 1}};
		bool fullFixed;
		T = fix_knot_vector_to_interpolate_curve_WKW(degree, T, parT, per, fullFixed, -1);
		if (!fullFixed)
		{
			std::cout << "The knot vector is not properly generated!\n";
			exit(0);
		}
		std::cout<<"T knots: ";
		for(auto ele : T)
		{
			std::cout<<ele<<",";
		}
		std::cout<<"\n";
		// solve the final control points
		Bcurve curveTool;
		curveTool.degree = degree;
		curveTool.U = T;
		CPAll.resize(ncpi);
		for (int i = 0; i < ncpi; i++)
		{
			CPAll[i].resize(ncpj);
			for (int j = 0; j < ncpj; j++)
			{
				std::vector<Vector3d> localdata(nbrTime);
				for (int k = 0; k < nbrTime; k++)
				{
					localdata[k] = CPs[k][i][j];
				}
				curveTool.solve_control_points_for_fairing_curve(curveTool, parT, localdata, 0.5, 0.5);
				CPAll[i][j] = curveTool.control_points;
			}
		}
		std::cout << "all the control points are solved\n";
		PolynomialBasis tBasisTool;
		auto tbf = tBasisTool.calculate_single(degree, T); // obtain the basis functions of the knot vector T
		std::cout << "basis functions T is computed\n";
		// check interpolation errors
		double error = 0;
		for (int i = 0; i < temps.size(); i++)
		{
			for (int j = 0; j < temps[i].size(); j++)
			{
				double temperature = temps[i][j];
				if (temperature < -100 || temperature > 100)
				{
					continue;
				}
				double u = param(i, 0);
				double v = param(i, 1);
				double t = parT[j];
				Eigen::Vector3d fitting = surface.BSplineSurfacePoint(basis.Ubasis, basis.Vbasis, tbf, u, v, t, surface.U, surface.V, T, CPAll, degree, 
				degree, degree);
				double e = temperature - fitting[2];
				if (e > error)
				{
					error = e;
				}
			}
		}
		std::cout << "maximal distance error is " << error << "\n";
		// make surface of a given moment, plot the mesh, and the original data
		auto shiftCoordBack=[](const Eigen::Vector2d &par, const Eigen::Vector2d& offset, const double length, const Eigen::Vector2d& cmin)
		{
			Eigen::Vector2d result;
			result = (par-offset)*length/0.9+cmin;
			return result;
		};
		Eigen::MatrixXi quads;
		Eigen::MatrixXd vertices, points;
		std::vector<Vector3d> vers;
		int discretization = 150; // the mesh size
		int moment = 10; // the moment-th step
		double scale = 0.1; // re-scale the height
		double itv = 1.0/discretization;
		double t = parT[moment];
		
		for(int i=0;i<discretization;i++)
		{
			double u = i * itv;
			for (int j = 0; j < discretization; j++)
			{
				double v = j * itv;
				Eigen::Vector3d fitting = surface.BSplineSurfacePoint(basis.Ubasis, basis.Vbasis, tbf, u, v, t, surface.U, surface.V, T, CPAll, degree, 
				degree, degree);
				Vector2d par(fitting(0), fitting(1));
				Vector2d shifted = shiftCoordBack(par, offset, length, cmin);
				fitting[0] = shifted[0];
				fitting[1] = shifted[1];
				fitting[2] *= scale; // rescale for better visulize the temperature
				vers.push_back(fitting);
			}
		}
		vertices = vec_list_to_matrix(vers);
		surface.constructRegularF(vertices.rows(), discretization, quads);
		vers.clear();
		for (int i = 0; i < temps.size(); i++)
		{
			double temperature = temps[i][moment];
			if (temperature < -100 || temperature > 100)
			{
				continue;
			}
			double u = stations[i][1];
			double v = stations[i][2];
			temperature *= scale; // rescale for better visulize the temperature
			vers.push_back(Eigen::Vector3d(u,v,temperature)); // push the original data
		}
		points = vec_list_to_matrix(vers);
		// cut the quad mesh using the boundary info, and convert it into a triangular mesh
		Eigen::MatrixXd Vclean;
		Eigen::MatrixXi Fclean;
		cutBoundaryGenerateTopology(vertices, quads, paramb, Vclean, Fclean);
		std::vector<double> colors;
		for (int i = 0; i < Vclean.rows(); i++)
		{
			colors.push_back(Vclean(i, 2) / scale);
		}

		// // export the ground truth
		quads.resize(0,0);
		// remind that the longitude is in an reversed direction, so we need to invert the tempeture
		points.col(2) = -1 * points.col(2);
		Vclean.col(2) = -1 * Vclean.col(2);
		std::cout<<"Predicting the temperature in Shanxi province at 18:00, 02, March 2025...\n";
		writeMesh(example_root_path+"Y2025M03D02H1800GroundTruth.obj",points, quads);
		// the fitting result
		writeMesh(example_root_path+"Y2025M03D02H1800.obj",Vclean, Fclean);
		// the temperature as the color
		write_csv(example_root_path+"Y2025M03D02H1800Color.csv", colors);
		std::cout << "fitting result min tem " << *std::min_element(colors.begin(), colors.end()) << " max tem " << *std::max_element(colors.begin(), colors.end()) << "\n";
		
		// given a location, check the fitting result
		int location = 10;
		std::cout<<"checking a given location, Ningwu County, #"<<stations[location][0]<<", that locates at ("<<stations[location][1]<<", "<<stations[location][2]<<")\n";
		double u = param(location, 0);
		double v = param(location, 1);
		colors.clear();
		std::vector<double> times;
		for(int i=0;i<discretization;i++)
		{
			double t = itv * i;
			Eigen::Vector3d fitting = surface.BSplineSurfacePoint(basis.Ubasis, basis.Vbasis, tbf, u, v, t, surface.U, surface.V, T, CPAll, degree, 
				degree, degree);
			colors.push_back(fitting[2]);
			times.push_back(t);
		}
		// the temperature at this location
		write_csv(example_root_path+"FittingResultsGivenLocation.csv", colors);
		write_csv(example_root_path+"FittingResultsGivenLocation_times.csv", times);
		colors.clear();
		times.clear();
		for (int i = 0; i < temps[location].size(); i++)
		{
			double temperature = temps[location][i];
			if (temperature < -100 || temperature > 100)
			{
				continue;
			}
			colors.push_back(temperature);
			times.push_back(parT[i]);
		}
		// the ground truth temperature at this location
		write_csv(example_root_path+"groundTruthGivenLocation.csv", colors);
		write_csv(example_root_path+"groundTruthGivenLocation_times.csv", times);
	}
}