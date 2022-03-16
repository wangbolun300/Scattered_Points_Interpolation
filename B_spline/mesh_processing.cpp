#include <sparse_interp/mesh_processing.h>
#include <array>
#include <sparse_interp/curve.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/write_triangle_mesh.h>
#include <queue>
#include <igl/predicates/predicates.h>
#include <igl/triangle/cdt.h>
#include <igl/triangle/triangulate.h>
#include <igl/cotmatrix.h>
#include <igl/point_mesh_squared_distance.h>

namespace SIBSplines
{
	Eigen::MatrixXd vector_to_matrix_3d(const std::vector<Vector3d> &v)
	{
		Eigen::MatrixXd result(v.size(), 3);
		for (int i = 0; i < v.size(); i++)
		{
			result.row(i) = v[i];
		}
		return result;
	}
	std::vector<Vector3d> matrix3d_to_vector(const Eigen::MatrixXd &v)
	{
		std::vector<Vector3d> result(v.rows());
		for (int i = 0; i < v.rows(); i++)
		{
			result[i] = v.row(i);
		}
		return result;
	}

	void vertices_to_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXi &edges)
	{
		edges.resize(pts.rows() - 1, 2);
		for (int i = 0; i < edges.rows(); i++)
		{
			edges(i, 0) = i;
			edges(i, 1) = i + 1;
		}
	}

	void read_and_visual_mesh(const std::string &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
	{
		bool read = igl::read_triangle_mesh(filename, V, F);
	}

	// map the vertices of the mesh boundary to a [0, 1] x [0, 1] domain
	void map_vertices_to_square(const Eigen::MatrixXd &V, const Eigen::VectorXi &bnd, Eigen::MatrixXd &bnd_uv)
	{
		bnd_uv.resize(bnd.size(), 2);
		double total_dist = 0;
		std::vector<double> pdis(bnd.size());
		pdis[0] = (V.row(bnd[1]) - V.row(bnd[0])).norm();
		for (int i = 1; i < pdis.size() - 1; i++)
		{
			pdis[i] = (V.row(bnd[i + 1]) - V.row(bnd[i])).norm() + pdis[i - 1];
		}
		pdis[pdis.size() - 1] = pdis[pdis.size() - 2] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();
		total_dist = pdis[pdis.size() - 1];

		int k = 1;
		double quarter = total_dist / 4;
		std::array<double, 4> corners;
		corners[0] = 0;
		for (int i = 0; i < pdis.size() - 1; i++)
		{
			if (pdis[i] <= quarter * k && pdis[i + 1] > quarter * k)
			{
				if (i != corners[k - 1])
				{
					corners[k] = i;
				}
				else
				{
					corners[k] = i + 1;
				}
				k += 1;
				if (k == 4)
					break;
			}
		}
		// u=0, v increase
		bnd_uv.row(0) << 0, 0;
		for (int i = 0; i < corners[1]; i++)
		{
			bnd_uv.row(i + 1) << 0, pdis[i] / pdis[corners[1] - 1];
		}

		// v=1, u increase
		for (int i = corners[1]; i < corners[2]; i++)
		{
			bnd_uv.row(i + 1) << (pdis[i] - pdis[corners[1] - 1]) / (pdis[corners[2] - 1] - pdis[corners[1] - 1]), 1;
		}

		// u=1, v decrease
		for (int i = corners[2]; i < corners[3]; i++)
		{
			bnd_uv.row(i + 1) << 1, 1 - (pdis[i] - pdis[corners[2] - 1]) / (pdis[corners[3] - 1] - pdis[corners[2] - 1]);
		}

		// v=0, u decrease
		for (int i = corners[3]; i < pdis.size() - 1; i++)
		{
			bnd_uv.row(i + 1) << 1 - (pdis[i] - pdis[corners[3] - 1]) / (pdis[pdis.size() - 1] - pdis[corners[3] - 1]), 0;
		}
		// std::cout << "the parameters\n"<<bnd_uv << std::endl;
	}
	// remove the vertices not used in contructing the surface
	// this function does not consider if there are duplicated vertices (two vertices in
	// exact the same position)
	void remove_redundent_mesh_vertices(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
										Eigen::MatrixXd &Vout, Eigen::MatrixXi &Fout)
	{

		std::vector<bool> located(V.rows(), false);
		std::vector<int> map(V.rows());
		int nbr = 0;
		for (int i = 0; i < F.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (located[F(i, j)] == false)
				{
					located[F(i, j)] = true;
					map[F(i, j)] = nbr;
					nbr++;
				}
			}
		}
		Fout.resize(F.rows(), 3);
		for (int i = 0; i < F.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				Fout(i, j) = map[F(i, j)];
			}
		}
		Vout.resize(nbr, 3);
		int k = 0;
		for (int i = 0; i < V.rows(); i++)
		{
			if (located[i])
				Vout.row(map[i]) = V.row(i);
		}
	}

	// detect the duplicated vertices, and re-orginaze the faces to avoid open
	// boundary caused by vertices duplication
	void remove_duplicated_vertices(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
									Eigen::MatrixXd &Vout, Eigen::MatrixXi &Fout)
	{
		std::vector<int> map(V.rows(), -1);
		std::vector<bool> duplicated(V.rows(), false);
		int to_delete = 0;
		for (int i = 0; i < V.rows(); i++)
		{
			for (int j = i; j < V.rows(); j++)
			{
				if (i < j)
				{
					if (V.row(i) == V.row(j))
					{
						duplicated[j] = true;
						if (duplicated[i] == true)
						{
							map[j] = map[i];
						}
						else
						{
							map[j] = i;
							to_delete++;
						}
					}
				}
			}
		}
		Vout = V;
		Fout.resize(F.rows(), 3);
		for (int i = 0; i < Fout.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (duplicated[F(i, j)] == true)
				{
					Fout(i, j) = map[F(i, j)];
				}
				else
				{
					Fout(i, j) = F(i, j);
				}
			}
		}
	}

	void generate_clean_mesh_data_for_parametrization(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
													  Eigen::MatrixXd &Vout, Eigen::MatrixXi &Fout)
	{

		// the duplicated vertices will be located, and faces will not use all the vertices
		remove_duplicated_vertices(V, F, Vout, Fout);

		//  remove the unrelated vertices
		remove_redundent_mesh_vertices(Vout, Fout, Vout, Fout);
	}
	// this is to remove some faces whose x, y or z coordinates (the parameter 'axis') of vertices larger than 'value'
	// axis select from 0, 1, 2.
	void remove_some_faces(const int axis, const double value, const bool remove_larger,
						   const Eigen::MatrixXd &V, const Eigen::MatrixXi F, Eigen::MatrixXi &newF)
	{
		int rows = F.rows();
		newF.resize(1, 3);
		int already_updated = 0;
		for (int i = 0; i < rows; i++)
		{
			if (remove_larger)
			{
				if (V(F(i, 0), axis) > value || V(F(i, 1), axis) > value || V(F(i, 2), axis) > value)
				{
					continue;
				}
			}
			else
			{
				if (V(F(i, 0), axis) < value || V(F(i, 1), axis) < value || V(F(i, 2), axis) < value)
				{
					continue;
				}
			}
			newF.conservativeResize(already_updated + 1, 3);
			newF.row(already_updated) = F.row(i);
			already_updated++;
		}
	}

	// use igl::harmonic to map the mesh (0 genus) to a u-v domain.
	void mesh_parameterization(
		const std::string &meshfile, Eigen::MatrixXd &V, Eigen::MatrixXd &param, Eigen::MatrixXi &F)
	{
		/*const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
		const std::string filename = path + "camel_smallest.obj";*/
		Eigen::MatrixXd Vori;
		Eigen::MatrixXi Fori;
		read_and_visual_mesh(meshfile, Vori, Fori);

		V = Vori;
		F = Fori;
		/////////////////////////////////////////
		// parameterization part
		Eigen::VectorXi bnd;
		igl::boundary_loop(F, bnd); // boundary vertices detection
		Eigen::MatrixXd bnd_uv;
		/*igl::map_vertices_to_circle(V, bnd, bnd_uv);*/
		map_vertices_to_square(V, bnd, bnd_uv);

		igl::harmonic(V, F, bnd, bnd_uv, 1, param);
	}

	// find the shortest edge length of the mesh.
	double shortest_edge_length_of_parametric_domain(const Eigen::MatrixXd &paras, const Eigen::MatrixXi &F)
	{
		double rst = 1;
		for (int i = 0; i < F.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				int id0 = F(i, j);
				int id1 = F(i, (j + 1) % 3);
				/*if (paras(id0, 0) == paras(id1, 0)) {
					if (paras(id0, 0) == 0 || paras(id0, 0) == 1) {
						continue;
					}
				}
				*/
				double dist = (paras.row(id0) - paras.row(id1)).norm();
				if (dist < rst)
				{
					rst = dist;
				}
			}
		}
		return rst;
	}
	int orient_2d(const Vector2d &a, const Vector2d &b, const Vector2d &c)
	{
		auto result = igl::predicates::orient2d(a, b, c);
		if (result == igl::predicates::Orientation::POSITIVE)
		{
			return 1;
		}
		if (result == igl::predicates::Orientation::NEGATIVE)
		{
			return -1;
		}
		return 0;
	}

	Eigen::MatrixXi update_info_list(const int info1, const int info2, const int ref, const int k,
									 const int j, const Eigen::MatrixXi &list)
	{
		Eigen::MatrixXi result = list;
		int rows = result.rows();
		result(k, j) = ref;
		result(k + 1, j) = ref;

		if (info1 != -1 && info2 != -1)
		{
			for (int i = 0; i < rows; i++)
			{
				if (result(i, j) == info2 || result(i, j) == info1)
				{
					result(i, j) = ref;
				}
			}
		}
		return result;
	}


	int find_a_border_point(const Eigen::MatrixXd &V, const bool is_v, const bool smallest)
	{
		double value = V(0, is_v);
		int id = 0;
		if (smallest)
		{
			for (int i = 0; i < V.rows(); i++)
			{
				if (value > V(i, is_v))
				{
					value = V(i, is_v);
					id = i;
				}
			}
		}
		else
		{
			for (int i = 0; i < V.rows(); i++)
			{
				if (value < V(i, is_v))
				{
					value = V(i, is_v);
					id = i;
				}
			}
		}
		return id;
	}

	// when co-linear, check if p is on closed s0-s1
	bool point_on_segment_2d(const Vector2d &s0, const Vector2d &s1, const Vector2d &p)
	{
		double umin = std::min(s0[0], s1[0]);
		double umax = std::max(s0[0], s1[0]);
		double vmin = std::min(s0[1], s1[1]);
		double vmax = std::max(s0[1], s1[1]);
		double u = p[0];
		double v = p[1];
		if (u <= umax && u >= umin)
		{
			if (v <= vmax && v >= vmin)
			{
				return true;
			}
		}
		return false;
	}

	// given a 2d point set V, detect its borders loop
	void find_border_loop(const Eigen::MatrixXd &V, Eigen::VectorXi &loop)
	{
		// first find the point with smallest u
		assert(V.cols() == 2);
		int id0 = find_a_border_point(V, false, true); // umin
		int id1 = find_a_border_point(V, true, true);  // vmin
		if (id0 == id1)
		{
			id1 = find_a_border_point(V, true, false);
		}
		// use id0,id1 as a skeleton
		// now uminid is the first point of loop

		bool id1_used = false;					  // when find id1. we use this flag to switch to another temporary edge
		std::vector<bool> border_flags(V.rows()); // if a point is a border point, set it as true;
		for (int i = 0; i < border_flags.size(); i++)
		{
			border_flags[i] = false;
		}
		// std::cout << "start push border" << std::endl;
		std::vector<int> vec;
		vec.push_back(id0);
		border_flags[id0] = true;
		for (int i = 0; i < vec.size(); i++)
		{
			// std::cout << "ith border point "<<i <<" "<<vec[i]<< std::endl;
			int current = vec[i];
			if (vec.size() >= 2 && vec.front() == vec.back())
			{ // the first element is the last element, searching finished
				break;
			}
			int tmpid = id1_used ? id0 : id1;

			// now check all the orientations against segment [vec[i],tmpid], and update tmpid
			for (int j = 0; j < V.rows(); j++)
			{
				if (border_flags[j])
				{
					continue;
				}
				// TODO
				Vector2d s0 = V.row(vec[i]), s1 = V.row(tmpid), query = V.row(j);
				int ori = orient_2d(s0, s1, query);
				if (ori == -1)
				{
					tmpid = j;
				}
				if (ori == 0)
				{
					if (point_on_segment_2d(s0, s1, query))
					{
						tmpid = j;
					}
				}
			}
			vec.push_back(tmpid); //
			border_flags[tmpid] = true;
			if (tmpid == id1)
			{
				id1_used = true;
			}
			if (tmpid == id0)
			{
				break;
			}
		}
		if (vec.front() == vec.back())
		{
			vec.pop_back(); // avoid duplication
		}
		loop.resize(vec.size());
		for (int i = 0; i < vec.size(); i++)
		{
			loop[i] = vec[i];
		}
	}

	// the inputs are cleaned(no duplicated) data points.
	// if V has 3 cols, then truncate it to 2d.
	void constrained_delaunay_triangulation(
		const Eigen::MatrixXd &V, const Eigen::VectorXi &Edge_ids,
		Eigen::MatrixXi &F)
	{

		Eigen::MatrixXd Vin(V.rows(), 2);
		if (V.cols() == 2)
		{
			Vin = V;
		}
		else
		{
			Vin << V.col(0), V.col(1);
		}
		Eigen::MatrixXi Ein(Edge_ids.size(), 2);
		for (int i = 0; i < Edge_ids.size() - 1; i++)
		{
			Ein(i, 0) = Edge_ids[i];
			Ein(i, 1) = Edge_ids[i + 1];
		}
		Ein(Edge_ids.size() - 1, 0) = Edge_ids[Edge_ids.size() - 1];
		Ein(Edge_ids.size() - 1, 1) = Edge_ids[0];

		Eigen::MatrixXd WV;
		Eigen::MatrixXi WF, WE;
		Eigen::VectorXi J;
		Eigen::MatrixXd H;
		H.resize(0, 0); // there is no hole
		// std::string flags = "-c";
		// igl::triangle::cdt(Vin, Ein, flags, WV, WF, WE, J);

		igl::triangle::triangulate(Vin, Ein, H, "", WV, WF);

		F = WF;
		assert(WV.rows() == Vin.rows());
	}

	// given parameters and the connectivity, get the U and V parameters, and a map showing the positions of the points
	// in U and V
	void generate_UV_grid(const Eigen::MatrixXd &param,
						  std::vector<double> &U, std::vector<double> &V, Eigen::MatrixXi &map)
	{
		std::vector<int> Umap, Vmap;
		int pnbr = param.rows();
		Umap.resize(pnbr); // show the u parameter correspond to which U
		Vmap.resize(pnbr);
		U.clear();
		V.clear();
		U.reserve(pnbr);
		V.reserve(pnbr);
		auto cmp = [](std::pair<double, int> i1, std::pair<double, int> i2)
		{
			return i1.first > i2.first;
		};

		std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(cmp)> queue_u(cmp);
		std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(cmp)> queue_v(cmp);
		for (int i = 0; i < pnbr; i++)
		{
			std::pair<double, int> upair;
			std::pair<double, int> vpair;
			upair.first = param(i, 0);
			upair.second = i;
			vpair.first = param(i, 1);
			vpair.second = i;
			queue_u.push(upair);
			queue_v.push(vpair);
		}
		double ut = -1, vt = -1;
		int uorder = -1, vorder = -1;
		for (int i = 0; i < pnbr; i++)
		{
			bool udiff = false;
			bool vdiff = false;
			if (queue_u.top().first != ut)
			{
				udiff = true;
				ut = queue_u.top().first;
				U.push_back(ut);
			}
			if (queue_v.top().first != vt)
			{
				vdiff = true;
				vt = queue_v.top().first;
				V.push_back(vt);
			}

			Umap[queue_u.top().second] = U.size() - 1;
			Vmap[queue_v.top().second] = V.size() - 1;
			queue_u.pop();
			queue_v.pop();
		}
		map = Eigen::MatrixXi::Constant(U.size(), V.size(), -1);
		for (int i = 0; i < pnbr; i++)
		{
			int urefer = Umap[i];
			int vrefer = Vmap[i];
			if (map(urefer, vrefer) != -1)
			{
				std::cout << "PARAMETRIZATION ERROR, DUPLICATION" << std::endl;
			}
			map(urefer, vrefer) = i;
		}
	}


	// 0, in the middle of the triangle
	// 1, on edge F0, F1
	// 2, on edge F1, F2
	// 3, on edge F2, F0
	// since in our algorithm, the point is not possible on the vertices of the mesh
	bool is_point_in_triangle_2d(const Vector2d &v0, const Vector2d &v1, const Vector2d &v2, const Vector2d &p,
								 int &orientation, int &statement)
	{

		int o0 = orient_2d(v0, v1, p);
		int o1 = orient_2d(v1, v2, p);
		int o2 = orient_2d(v2, v0, p);
		if (o0 == o1 && o1 == o2)
		{
			orientation = o0;
			statement = 0;
			return true;
		}
		if (o0 == o1 && o2 == 0)
		{
			orientation = o0;
			statement = 3;
			return true;
		}
		if (o1 == o2 && o0 == 0)
		{
			orientation = o1;
			statement = 1;
			return true;
		}
		if (o2 == o0 && o1 == 0)
		{
			orientation = o0;
			statement = 2;
			return true;
		}
		return false;
	}
	// given a parameter u and v, we find one ring triangle in F, and returns the orientation of F
	// the status is:
	// 0, in the middle of the triangle
	// 1, on edge F0, F1
	// 2, on edge F1, F2
	// 3, on edge F2, F0
	// since in our algorithm, the point is not possible on the vertices of the mesh
	int find_one_ring_for_param(const double u, const double v, const Eigen::MatrixXd &param, const Eigen::MatrixXi &F,
								int &orientation, int &statement)
	{
		Vector2d p(u, v);
		for (int i = 0; i < F.rows(); i++)
		{
			Vector2d v0 = param.row(F(i, 0));
			Vector2d v1 = param.row(F(i, 1));
			Vector2d v2 = param.row(F(i, 2));

			if (is_point_in_triangle_2d(v0, v1, v2, p, orientation, statement))
			{
				return i;
			}
		}
		assert(1);
		std::cout << " The Code Should Not Go Here" << std::endl;
		return -1;
	}
	double area2d(const Vector2d &p0, const Vector2d &p1, const Vector2d &p2)
	{
		Eigen::MatrixXd m(2, 2);
		m << p0 - p1, p0 - p2;
		double v = m.determinant();
		// Vector2d vec = (p0 - p1).cross(p0 - p2);
		// double v = vec.norm();
		return fabs(v / 2);
		return 0;
	}
	bool check_UV_validation(const std::vector<double> &U)
	{
		for (int i = 0; i < U.size() - 1; i++)
		{
			if (U[i] >= U[i + 1])
			{
				return false;
			}
		}
		return true;
	}

	

	double point_mesh_distance(const Eigen::MatrixXd ver, const Eigen::MatrixXd &Vmesh,
							   const Eigen::MatrixXi &Fmesh, Eigen::VectorXd &sqrD)
	{
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		igl::point_mesh_squared_distance(ver, Vmesh, Fmesh, sqrD, I, C);

		return sqrt(sqrD.maxCoeff());
	}

	void write_points(const std::string &file, const Eigen::MatrixXd &ver)
	{
		std::ofstream fout;
		fout.open(file);
		for (int i = 0; i < ver.rows(); i++)
		{
			fout
				//<< std::setprecision(17)
				<< "v " << ver(i, 0) << " " << ver(i, 1) << " " << ver(i, 2) << std::endl;
		}
		fout.close();
	}
	bool write_triangle_mesh(const std::string filename, const Eigen::MatrixXd &ver, Eigen::MatrixXi &f)
	{
		return igl::write_triangle_mesh(filename, ver, f);
	}
}