#include <sparse_interp/comparison.h>
#include <sparse_interp/curve.h>
#include <sparse_interp/mesh_processing.h>
#include <sparse_interp/energy.h>

namespace SIBSplines
{
	//
	//
	// Lofting
	//
	//
	//
	//
	namespace comparison
	{
		void lofting_method_generate_interpolation_knot_vectors(const bool start_from_v_direction, int degree1, int degree2,
																std::vector<double> &Uknot, std::vector<double> &Vknot,
																const Eigen::MatrixXd &param_original,
																const double per)
		{

			std::vector<double> Ugrid, Vgrid;
			Eigen::MatrixXi grid_map;
			generate_UV_grid(param_original, Ugrid, Vgrid, grid_map);
			std::cout << "** UV grid sizes, " << Ugrid.size() << ", " << Vgrid.size() << std::endl;
			bool fully_fixed;
			for (int i = 0; i < Vgrid.size(); i++)
			{
				std::vector<double> paras = get_iso_line_parameters(degree1, degree2, true, i, Ugrid, Vgrid, grid_map);
				// std::cout << "\nthe " << i << "th iso line parameters " << std::endl;
				// print_vector(paras);
				Uknot = fix_knot_vector_to_interpolate_curve_WKW(degree1, Uknot, paras, per, fully_fixed);
				assert(fully_fixed == true);
			}
			std::cout << "finished initialize Uknot" << std::endl;
			print_vector(Uknot);
			/*std::cout << "\n** the fixed U knot" << std::endl;
			 */

			// fix iso-u lines knot vector
			for (int i = 0; i < Ugrid.size(); i++)
			{ // for each u parameter
				std::vector<double> paras = get_iso_line_parameters(degree1, degree2, false, i, Ugrid, Vgrid, grid_map);
				/*std::cout << "\nthe " << i << "th iso line parameters " << std::endl;
				print_vector(paras);*/
				Vknot = fix_knot_vector_to_interpolate_curve_WKW(degree2, Vknot, paras, per, fully_fixed);
				assert(fully_fixed == true);
			}
			std::cout << "finished initialize Vknot" << std::endl;
			print_vector(Vknot);
			if (start_from_v_direction)
			{
				Vknot = fix_knot_vector_to_interpolate_curve_WKW(degree2, Vknot, Vgrid, per, fully_fixed);
			}
			else
			{
				Uknot = fix_knot_vector_to_interpolate_curve_WKW(degree1, Uknot, Ugrid, per, fully_fixed);
			}
		}

		//
		//
		// Averaging
		//
		//
		//
		//
		void piegl_method_generate_interpolation_knot_vectors(int degree1, int degree2,
															  std::vector<double> &Uknot, std::vector<double> &Vknot,
															  const Eigen::MatrixXd &param_original,
															  const double per)
		{
			std::vector<double> Ugrid, Vgrid;
			Eigen::MatrixXi grid_map;
			generate_UV_grid(param_original, Ugrid, Vgrid, grid_map);
			bool ff;
			Uknot = fix_knot_vector_to_interpolate_curve_WKW(degree1, Uknot, Ugrid, per, ff);
			Vknot = fix_knot_vector_to_interpolate_curve_WKW(degree2, Vknot, Vgrid, per, ff);
		}

		//
		//
		// Multilevel B-spline
		//
		//

		std::vector<double> knot_vector_level(const int degree, const int level)
		{
			std::vector<double> result;
			for (int i = 0; i < degree + 1; i++)
			{
				result.push_back(0);
			}
			double current = 0;
			int nbr_intervals = 1;
			for (int i = 0; i < level; i++)
			{
				nbr_intervals *= 2;
			}
			if (level > 0)
			{
				double length = 1.0 / nbr_intervals;
				for (int i = 0; i < nbr_intervals - 1; i++)
				{
					current += length;
					result.push_back(current);
				}
			}
			for (int j = 0; j < degree + 1; j++)
			{
				result.push_back(1);
			}
			return result;
		}
		Eigen::MatrixXd equality_part_of_Seungyong(Bsurface &surface,
												   const std::vector<std::vector<std::vector<int>>> &overlaps, const int nbr_related)
		{

			int psize = (surface.nu() + 1) * (surface.nv() + 1); // total number of control points.

			Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nbr_related, psize);
			int degree1 = surface.degree1;
			int degree2 = surface.degree2;
			std::vector<double> U = surface.U;
			std::vector<double> V = surface.V;

			for (int i = 0; i < result.rows(); i++)
			{
				int counter = 0;
				for (int j = 0; j < result.cols(); j++)
				{
					// figure out the jth control point corresponding to which Pij
					int coff_i = j / (surface.nv() + 1);
					int coff_j = j - coff_i * (surface.nv() + 1);
					if (overlaps[coff_i][coff_j].size() > 0)
					{
						if (counter >= i)
						{
							result(i, j) = 1;
							break;
						}
						counter++;
					}
				}
			}
			return result;
		}

		void none_zero_basis_ids(const int interval, const std::vector<double> &U, const double u, const int nu,
								 const int degree, int &id0, int &id1)
		{
			if (interval < 0)
			{
				id0 = nu;
				id1 = nu;
				return;
			}
			double left = U[interval];
			if (u > left)
			{
				id0 = interval - degree;
				id1 = interval;
				return;
			}
			// u==left
			if (left == U[0])
			{
				id0 = 0;
				id1 = 0;
				return;
			}
			else
			{
				id0 = interval - degree;
				id1 = interval - 1;
				return;
			}
		}
		Eigen::MatrixXd left_part_of_Seungyong(Bsurface &surface, PartialBasis &basis,
											   const Eigen::MatrixXd &param, std::vector<std::vector<std::vector<int>>> &overlaps,
											   std::vector<int> &uinterval, std::vector<int> &vinterval, int &nbr_related)
		{
			std::vector<double> U = surface.U;
			std::vector<double> V = surface.V;
			uinterval.resize(param.rows());
			vinterval.resize(param.rows());
			int degree1 = surface.degree1;
			int degree2 = surface.degree2;
			// std::cout << "before 01" << std::endl;
			for (int i = 0; i < param.rows(); i++)
			{
				double u = param(i, 0);
				double v = param(i, 1);
				int which = -1;
				for (int j = 0; j < U.size() - 1; j++)
				{
					if (u >= U[j] && u < U[j + 1])
					{
						which = j;
						break;
					}
				}
				uinterval[i] = which;

				which = -1;
				for (int j = 0; j < V.size() - 1; j++)
				{
					if (v >= V[j] && v < V[j + 1])
					{
						which = j;
						break;
					}
				}
				vinterval[i] = which;
			}
			// std::cout << "before 02" << std::endl;
			int nu = surface.nu();
			int nv = surface.nv();
			// (nu+1)x(nv+1) control points
			overlaps.resize(nu + 1);
			nbr_related = 0;
			for (int i = 0; i < overlaps.size(); i++)
			{
				overlaps[i].resize(nv + 1);
			}
			// std::cout << "before 03" << std::endl;
			for (int i = 0; i < param.rows(); i++)
			{
				double u = param(i, 0);
				double v = param(i, 1);
				int k1 = uinterval[i];
				int k2 = vinterval[i];
				int id0, id1, id2, id3; // id0=i-p, id1=i.
				none_zero_basis_ids(k1, U, u, nu, degree1, id0, id1);
				none_zero_basis_ids(k2, V, v, nv, degree2, id2, id3);

				assert(id0 >= 0 && id0 <= nu + 1);
				assert(id1 >= 0 && id1 <= nu + 1);
				assert(id2 >= 0 && id2 <= nv + 1);
				assert(id3 >= 0 && id3 <= nv + 1);
				for (int j = id0; j < id1 + 1; j++)
				{
					for (int k = id2; k < id3 + 1; k++)
					{
						assert(j < overlaps.size());
						assert(k < overlaps[j].size());
						overlaps[j][k].push_back(i);
					}
				}
			}
			// std::cout << "before 04" << std::endl;
			for (int i = 0; i < nu + 1; i++)
			{
				for (int j = 0; j < nv + 1; j++)
				{
					if (overlaps[i][j].size() > 0)
					{
						nbr_related++;
					}
				}
			}
			int cnbr = (nu + 1) * (nv + 1); // nbr of the control points

			// int sz = nbr_related + cnbr; // matrix size
			// Eigen::MatrixXd result(sz, sz);
			int sz = cnbr; // matrix size
			Eigen::MatrixXd result = Eigen::MatrixXd::Identity(1, 1);
			// Eigen::MatrixXd lu = energy_part_of_surface_least_square(surface, basis);
			// Eigen::MatrixXd rd = Eigen::MatrixXd::Zero(nbr_related, nbr_related);// right down corner part
			// Eigen::MatrixXd ld = equality_part_of_Seungyong(surface, overlaps,nbr_related);
			// Eigen::MatrixXd ru = -ld.transpose();
			// result << lu, ru, ld, rd;
			return result;
		}
		double right_part_element_for_Seungyong(const int rid, const int cid,
												const std::vector<std::vector<std::vector<int>>> &overlaps, const std::vector<int> &uinterval,
												const std::vector<int> &vinterval, const std::vector<double> &U,
												const std::vector<double> &V, const int degree1, const int degree2, int nu, int nv,
												const Eigen::MatrixXd &ver, const Eigen::MatrixXd &param, const int dimension)
		{
			std::vector<double> bottoms(ver.rows());
			for (int i = 0; i < param.rows(); i++)
			{
				int k1 = uinterval[i];
				int k2 = vinterval[i];
				double u = param(i, 0);
				double v = param(i, 1);
				double btm = 0;
				int id0, id1, id2, id3; // id0=i-p, id1=i.
				none_zero_basis_ids(k1, U, u, nu, degree1, id0, id1);
				none_zero_basis_ids(k2, V, v, nv, degree2, id2, id3);

				for (int j = id0; j < id1 + 1; j++)
				{
					for (int k = id2; k < id3 + 1; k++)
					{
						double value = Nip(j, degree1, u, U) * Nip(k, degree2, v, V);
						btm += value * value;
					}
				}
				assert(btm > 0);
				bottoms[i] = btm;
			}

			std::vector<int> list = overlaps[rid][cid];
			std::vector<double> tops(list.size());
			std::vector<double> elements(list.size());
			for (int i = 0; i < list.size(); i++)
			{
				int idtmp = list[i];
				double u = param(idtmp, 0);
				double v = param(idtmp, 1);
				double top = Nip(rid, degree1, u, U) * Nip(cid, degree2, v, V);
				tops[i] = top;
				if (top <= 0)
				{
					std::cout << "reason, N " << Nip(rid, degree1, u, U) << " u " << u << std::endl;
					std::cout << "reason, N " << Nip(cid, degree2, v, V) << " v " << v << std::endl;
				}
				elements[i] = top * ver(idtmp, dimension) / bottoms[idtmp];
			}

			double sum = 0;
			for (int i = 0; i < list.size(); i++)
			{
				sum += tops[i] * tops[i];
			}
			double result = 0;
			for (int i = 0; i < list.size(); i++)
			{
				result += tops[i] * tops[i] * elements[i];
			}
			assert(sum > 0);
			result = result / sum;
			return result;
		}

		Eigen::VectorXd right_part_of_Seungyong(Bsurface &surface, const Eigen::MatrixXd &paras,
												const Eigen::MatrixXd &ver, const int dimension,
												const std::vector<std::vector<std::vector<int>>> &overlaps,
												const std::vector<int> &uinterval, const std::vector<int> &vinterval, const int nbr_related)
		{

			int nu = surface.nu();
			int nv = surface.nv();
			int psize = (surface.nu() + 1) * (surface.nv() + 1); // total number of control points.
			Eigen::VectorXd result(psize);
			int dealing = 0;
			for (int i = 0; i < nu + 1; i++)
			{
				for (int j = 0; j < nv + 1; j++)
				{
					if (overlaps[i][j].size() > 0)
					{
						double value = right_part_element_for_Seungyong(i, j, overlaps, uinterval, vinterval, surface.U, surface.V, surface.degree1,
																		surface.degree2, nu, nv, ver, paras, dimension);
						assert(value < std::numeric_limits<double>::infinity() &&
							   value > -std::numeric_limits<double>::infinity());
						result[dealing] = value;
					}
					else
					{
						result[dealing] = 0;
					}
					dealing++;
				}
			}
			return result;

			for (int i = 0; i < psize; i++)
			{
				result[i] = 0;
			}
			for (int i = psize; i < result.size(); i++)
			{
				int checking = i - psize;
				int counter = 0;
				bool jump = false;
				double value = 0;
				for (int j = 0; j < nu + 1; j++)
				{
					for (int k = 0; k < nv + 1; k++)
					{
						if (overlaps[j][k].size() > 0)
						{
							if (counter >= checking)
							{
								value = right_part_element_for_Seungyong(j, k, overlaps, uinterval, vinterval, surface.U, surface.V, surface.degree1,
																		 surface.degree2, nu, nv, ver, paras, dimension);
								assert(value < std::numeric_limits<double>::infinity() &&
									   value > -std::numeric_limits<double>::infinity());
								jump = true;
								break;
							}
							counter++;
						}
					}
					if (jump)
					{
						break;
					}
				}
				result[i] = value;
			}
			return result;
		}
		void iteratively_approximate_method(int degree1, int degree2,
											std::vector<double> &Uknot, std::vector<double> &Vknot,
											const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver,
											const double tolerance,
											std::vector<Bsurface> &surfaces, const double per)
		{
			typedef Eigen::SparseMatrix<double> SparseMatrixXd;
			int level = 0;
			Eigen::MatrixXd err = ver;
			while (1)
			{
				Bsurface slevel;
				slevel.degree1 = degree1;
				slevel.degree2 = degree2;
				slevel.U = knot_vector_level(degree1, level);
				slevel.V = knot_vector_level(degree2, level);
				PartialBasis basis(slevel);
				std::vector<std::vector<std::vector<int>>> overlaps;
				std::vector<int> uinterval, vinterval;
				int nbr_related;
				std::cout << "before get left" << std::endl;
				Eigen::MatrixXd left = left_part_of_Seungyong(slevel, basis, param, overlaps, uinterval, vinterval, nbr_related);
				std::cout << "get left" << std::endl;
				/////////////////

				// SparseMatrixXd matB;
				// Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

				// matB = left.sparseView();

				// solver.compute(matB);
				// if (solver.info() != Eigen::Success) {
				//	// decomposition failed
				//	std::cout << "solving failed" << std::endl;
				//	return;
				// }

				/////////////////
				int sz = (slevel.nu() + 1) * (slevel.nv() + 1);
				std::vector<Vector3d> cps(sz);
				for (int i = 0; i < 3; i++)
				{
					Eigen::VectorXd right = right_part_of_Seungyong(slevel, param, err, i, overlaps, uinterval, vinterval, nbr_related);
					Eigen::MatrixXd p_lambda;

					p_lambda = right;
					/*	if (solver.info() != Eigen::Success) {
							std::cout << "solving failed" << std::endl;
							return;
						}*/
					for (int j = 0; j < cps.size(); j++)
					{
						cps[j][i] = p_lambda(j, 0);
					}
					// std::cout << "right part\n" << right << std::endl;
				}
				push_control_point_list_into_surface(slevel, cps);
				surfaces.push_back(slevel);
				double max_error;
				int max_overlap = 0;
				err = slevel.interpolation_err_for_apprximation(err, param, slevel, max_error);
				std::cout << "the ith iteration " << level << std::endl;
				std::cout << "sizes " << slevel.U.size() << " " << slevel.V.size() << std::endl;
				std::cout << "max error is " << max_error << std::endl;
				for (int i = 0; i < overlaps.size(); i++)
				{
					for (int j = 0; j < overlaps[i].size(); j++)
					{
						if (overlaps[i][j].size() > max_overlap)
						{
							max_overlap = overlaps[i][j].size();
						}
					}
				}
				std::cout << "max overlap is " << max_overlap << std::endl;

				// std::cout << "error matrix\n" << err << std::endl;
				if (max_error <= tolerance)
				{
					break;
				}
				level++;
			}
			return;
		}

		//
		//
		//
		// IGA
		//
		//
		//
		// this is an implementation of "B-spline surface fitting by iterative geometric interpolation/approximation
		// algorithms, 2012, CAD"
		double para_distance(const double uvalue, const double vvalue, const Eigen::MatrixXd &param, const int row)
		{
			Eigen::MatrixXd cornor(1, 2);
			cornor << uvalue, vvalue;
			return (cornor.row(0) - param.row(row)).norm();
		}

		// set up the initial surface as a bilinear surface interpolating the four conors of the surface
		void set_up_initial_surface(Bsurface &surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver)
		{
			int nbr = ver.rows();
			double dis00 = 2, dis01 = 2, dis10 = 2, dis11 = 2;
			int p00, p01, p10, p11;
			for (int i = 0; i < nbr; i++)
			{
				double dis = para_distance(0, 0, param, i);
				if (dis < dis00)
				{
					dis00 = dis;
					p00 = i;
				}

				dis = para_distance(0, 1, param, i);
				if (dis < dis01)
				{
					dis01 = dis;
					p01 = i;
				}

				dis = para_distance(1, 0, param, i);
				if (dis < dis10)
				{
					dis10 = dis;
					p10 = i;
				}

				dis = para_distance(1, 1, param, i);
				if (dis < dis11)
				{
					dis11 = dis;
					p11 = i;
				}
			}
			int c1 = surface.nu() + 1;
			int c2 = surface.nv() + 1;
			surface.control_points.resize(surface.nu() + 1);
			for (int i = 0; i < surface.nu() + 1; i++)
			{
				surface.control_points[i].resize(surface.nv() + 1);
			}
			surface.control_points[0][0] = ver.row(p00);
			surface.control_points[0][surface.nv()] = ver.row(p01);
			surface.control_points[surface.nu()][surface.nv()] = ver.row(p11);
			surface.control_points[surface.nu()][0] = ver.row(p10);

			for (int i = 0; i < surface.nu() + 1; i++)
			{
				for (int j = 0; j < surface.nv() + 1; j++)
				{
					double u = double(i) / surface.nu();
					double v = double(j) / surface.nv();
					Vector3d p0 = surface.control_points[0][0] +
								  (surface.control_points[surface.nu()][0] - surface.control_points[0][0]) * u;
					Vector3d p1 = surface.control_points[0][surface.nv()] +
								  (surface.control_points[surface.nu()][surface.nv()] - surface.control_points[0][surface.nv()]) * u;
					surface.control_points[i][j] = p0 + (p1 - p0) * v;
				}
			}
		}
		std::vector<Vector3d> get_the_delta_k(Bsurface &surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver)
		{
			std::vector<Vector3d> result;
			int size = param.rows();
			result.resize(size);
			for (int i = 0; i < size; i++)
			{
				Vector3d vs = surface.BSplineSurfacePoint(surface, param(i, 0), param(i, 1));
				Vector3d data = ver.row(i);
				result[i] = data - vs;
			}
			return result;
		}

		// returned ids means the parameter effects P_{id0}, ...
		std::vector<int> get_parameter_corresponded_control_points(const std::vector<double> &U, const int degree,
																   const double parameter)
		{
			int n = U.size() - degree - 2; // n+1 is the nbr of control points.
			std::vector<int> ids;
			if (parameter == U.back())
			{ // if parameter==1, then
				ids.push_back(n);
				return ids;
			}
			int interval;
			for (int i = 0; i < U.size() - 1; i++)
			{
				if (parameter >= U[i] && parameter < U[i + 1])
				{
					interval = i;
					break;
				}
			}
			for (int i = interval - degree; i < interval + 1; i++)
			{
				ids.push_back(i);
			}
			return ids;
		}

		void get_the_data_ids_and_weights_for_each_control_point(
			std::vector<std::vector<std::vector<int>>> &ids, std::vector<std::vector<std::vector<double>>> &weights,
			Bsurface &surface, const Eigen::MatrixXd &param)
		{
			int cnbr_0 = surface.nu() + 1; // nbr of control points in u direction
			int cnbr_1 = surface.nv() + 1;
			int pnbr = param.rows(); // nbr of data points
			ids.resize(cnbr_0);
			weights.resize(cnbr_0);
			for (int i = 0; i < cnbr_0; i++)
			{
				ids[i].resize(cnbr_1);
				weights[i].resize(cnbr_1);
			}
			for (int i = 0; i < pnbr; i++)
			{
				double u = param(i, 0);
				double v = param(i, 1);
				std::vector<int> u_ids = get_parameter_corresponded_control_points(surface.U, surface.degree1, u);
				std::vector<int> v_ids = get_parameter_corresponded_control_points(surface.V, surface.degree2, v);
				for (int r1 = 0; r1 < u_ids.size(); r1++)
				{
					for (int r2 = 0; r2 < v_ids.size(); r2++)
					{
						ids[u_ids[r1]][v_ids[r2]].push_back(i);
						double N1 = Nip(u_ids[r1], surface.degree1, u, surface.U);
						double N2 = Nip(v_ids[r2], surface.degree2, v, surface.V);
						weights[u_ids[r1]][v_ids[r2]].push_back(N1 * N2);
					}
				}
			}
		}

		void get_Delta(const std::vector<Vector3d> &delta, const std::vector<std::vector<std::vector<int>>> &ids,
					   const std::vector<std::vector<std::vector<double>>> &weights, std::vector<std::vector<Vector3d>> &Delta)
		{
			int pnbr = delta.size();
			int cnbr_u = ids.size();
			int cnbr_v = ids[0].size();
			Delta.resize(cnbr_u);
			for (int i = 0; i < cnbr_u; i++)
			{
				Delta[i].resize(cnbr_v);
			}
			for (int i = 0; i < cnbr_u; i++)
			{
				for (int j = 0; j < cnbr_v; j++)
				{
					Delta[i][j] = Vector3d(0, 0, 0);
					int nbr_vecs = ids[i][j].size();
					double weight_sum = 0;
					for (int k = 0; k < nbr_vecs; k++)
					{
						int delta_id = ids[i][j][k];
						Vector3d vec = delta[delta_id];
						double dweight = weights[i][j][k];
						weight_sum += dweight;
						Delta[i][j] += vec * dweight;
					}
					if (weight_sum == 0)
					{
						// std::cout << "WARNING: IPA method need fairing" << std::endl;
					}
					else
					{
						Delta[i][j] /= weight_sum;
					}
				}
			}
		}

		void control_points_add_Delta(std::vector<std::vector<Vector3d>> &cps, const std::vector<std::vector<Vector3d>> &Delta)
		{
			for (int i = 0; i < cps.size(); i++)
			{
				for (int j = 0; j < cps[i].size(); j++)
				{
					cps[i][j] += Delta[i][j];
				}
			}
		}
		double get_maximal_step_length(const std::vector<std::vector<Vector3d>> &Delta)
		{
			double result = 0;
			for (int i = 0; i < Delta.size(); i++)
			{
				for (int j = 0; j < Delta[i].size(); j++)
				{
					if (Delta[i][j].norm() > result)
					{
						result = Delta[i][j].norm();
					}
				}
			}
			return result;
		}
		bool Delta_step_length_larger_than_threadshold(const std::vector<std::vector<Vector3d>> &Delta, const double th)
		{
			double l = get_maximal_step_length(Delta);
			std::cout << "step length " << l << std::endl;
			if (l > th)
			{
				return true;
			}
			return false;
		}

		// gamma is an input parameter for the weight.
		void get_matrix_laplacian_fairing_of_control_points(double gamma, const std::vector<std::vector<Vector3d>> &cp,
															Eigen::MatrixXd &matrix)
		{

			auto get_mid_and_corners = [](std::array<std::array<int, 2>, 4> &cid1, std::array<std::array<int, 2>, 4> &cid2,
										  int uid, int vid)
			{
				cid1[0] = {uid - 1, vid};
				cid1[1] = {uid + 1, vid};
				cid1[2] = {uid, vid - 1};
				cid1[3] = {uid, vid + 1};

				cid2[0] = {uid - 1, vid - 1};
				cid2[1] = {uid + 1, vid - 1};
				cid2[2] = {uid + 1, vid + 1};
				cid2[3] = {uid - 1, vid + 1};
			};
			int unbr = cp.size();
			int vnbr = cp[0].size();
			auto to_matrix_id = [](const int uid, const int vid, const int unbr, const int vnbr, int &id)
			{
				id = uid * vnbr + vid;
			};
			std::array<std::array<int, 2>, 4> cid1; // four mid-control points
			std::array<std::array<int, 2>, 4> cid2; // four corner-control points

			matrix = Eigen::MatrixXd::Zero(unbr * vnbr, unbr * vnbr);
			for (int i = 0; i < unbr; i++)
			{
				for (int j = 0; j < vnbr; j++)
				{
					get_mid_and_corners(cid1, cid2, i, j); // get the ids of mid points and corner points
					int id_row;
					to_matrix_id(i, j, unbr, vnbr, id_row); // we are dealing with the ith row of the matrix
					double weight = 1;
					std::vector<int> mid_valid_col;
					std::vector<int> cor_valid_col;

					for (int k = 0; k < 4; k++)
					{ // get the denominator of the weithts, and the ids of the midpoints and corpoints.
						if (cid1[k][0] >= 0 && cid1[k][0] < unbr && cid1[k][1] >= 0 && cid1[k][1] < vnbr)
						{
							weight += exp(-gamma);
							int id_col;
							to_matrix_id(cid1[k][0], cid1[k][1], unbr, vnbr, id_col);
							mid_valid_col.push_back(id_col);
						}
						if (cid2[k][0] >= 0 && cid2[k][0] < unbr && cid2[k][1] >= 0 && cid2[k][1] < vnbr)
						{
							weight += exp(-2 * gamma);
							int id_col;
							to_matrix_id(cid2[k][0], cid2[k][1], unbr, vnbr, id_col);
							cor_valid_col.push_back(id_col);
						}
					}

					for (int k = 0; k < mid_valid_col.size(); k++)
					{
						matrix(id_row, mid_valid_col[k]) = exp(-gamma) / weight;
					}
					for (int k = 0; k < cor_valid_col.size(); k++)
					{
						matrix(id_row, cor_valid_col[k]) = exp(-2 * gamma) / weight;
					}
					matrix(id_row, id_row) = 1 / weight;
				}
			}
		}
		void laplacian_fairing_of_control_points(const Eigen::MatrixXd &matrix, std::vector<std::vector<Vector3d>> &cp)
		{
			int unbr = cp.size();
			int vnbr = cp[0].size();
			int current = 0;
			Eigen::MatrixXd right(unbr * vnbr, 3), result(unbr * vnbr, 3);

			current = 0;
			for (int i = 0; i < unbr; i++)
			{
				for (int j = 0; j < vnbr; j++)
				{
					right.row(current) = cp[i][j];
					current++;
				}
			}
			// get right part.

			current = 0;
			result = matrix * (right);
			for (int i = 0; i < unbr; i++)
			{
				for (int j = 0; j < vnbr; j++)
				{
					cp[i][j] = result.row(current);
					current++;
				}
			}
		}

		// threads hold is the maximal step length
		void progressive_iterative_approximation(Bsurface &surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver,
												 const int max_itr, const double threadshold)
		{
			int nu = surface.nu();
			int nv = surface.nv();
			set_up_initial_surface(surface, param, ver); // get the surface of iteration 0
			std::cout << "initial surface set up " << std::endl;
			std::vector<std::vector<std::vector<int>>> ids; // map from control points to data points ids.
			std::vector<std::vector<std::vector<double>>> weights;
			get_the_data_ids_and_weights_for_each_control_point(ids, weights, surface, param);
			std::cout << "weights calculated" << std::endl;
			double gamma = 2;
			Eigen::MatrixXd matrix;
			get_matrix_laplacian_fairing_of_control_points(gamma, surface.control_points, matrix); // prepare fairing
			// std::cout << "matrix\n" << matrix << std::endl; exit(0);
			for (int i = 0; i < max_itr; i++)
			{
				std::cout << "iteration " << i << std::endl;
				std::vector<Vector3d> delta_k = get_the_delta_k(surface, param, ver); // get the error vector for each data point
				std::vector<std::vector<Vector3d>> Delta;
				get_Delta(delta_k, ids, weights, Delta);
				if (!Delta_step_length_larger_than_threadshold(Delta, threadshold))
				{
					std::cout << "step length less than given threadshold" << std::endl;
					return; // if the step length is too short, do not update the surface.
				}
				control_points_add_Delta(surface.control_points, Delta); // update surface

				if (i % 19 == 0)
				{
					laplacian_fairing_of_control_points(matrix, surface.control_points);
				}
			}
		}

		void initialize_pia_knot_vector_single(const int degree1, const int nu, std::vector<double> &U)
		{
			U.resize(nu + degree1 + 2);
			for (int i = 0; i < degree1 + 1; i++)
			{
				U[i] = 0;
			}
			for (int i = nu + 1; i < nu + degree1 + 2; i++)
			{
				U[i] = 1;
			}
			int nbr_intervals = nu + 1 - degree1;
			double length_interval = double(1) / nbr_intervals;
			for (int i = degree1 + 1; i < nu + 1; i++)
			{
				U[i] = U[i - 1] + length_interval;
			}
		}

		void initialize_pia_knot_vectors(int degree1, int degree2,
										 std::vector<double> &U, std::vector<double> &V, int nu, int nv)
		{

			initialize_pia_knot_vector_single(degree1, nu, U);
			initialize_pia_knot_vector_single(degree2, nv, V);
			return;
		}

	}
}