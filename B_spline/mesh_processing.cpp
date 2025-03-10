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
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <igl/segment_segment_intersect.h>

namespace SIBSplines
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<> CGMesh;
	std::vector<Eigen::Vector3d> mat_to_vec_list(const Eigen::MatrixXd &m)
	{
		std::vector<Eigen::Vector3d> result(m.rows());
		for (int i = 0; i < m.rows(); i++)
		{
			result[i] = m.row(i);
		}
		return result;
	}
	Eigen::MatrixXd vec_list_to_matrix(const std::vector<Eigen::Vector3d> &vec)
	{
		Eigen::MatrixXd mat;
		if (vec.size() > 0)
		{
			mat.resize(vec.size(), 3);
			for (int i = 0; i < mat.rows(); i++)
			{
				mat.row(i) = vec[i];
			}
		}
		else
		{
			std::cout << "The vertex vector is empty" << std::endl;
		}
		return mat;
	}
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
	void writeMesh(const std::string filename, const Eigen::MatrixXd &ver, Eigen::MatrixXi &F)
	{
		std::ofstream fout;
		fout.open(filename);
		for (int i = 0; i < ver.rows(); i++)
		{
			fout << "v " << ver(i, 0) << " " << ver(i, 1) << " " << ver(i, 2) << "\n";
		}
		for (int i = 0; i < F.rows(); i++)
		{
			fout << "f ";
			for (int j = 0; j < F.cols(); j++)
			{
				fout << F(i, j) + 1 << " ";
			}
			fout << "\n";
		}
		fout.close();
	}

	inline int orient2d(
		const Eigen::Vector2d &p, const Eigen::Vector2d &q, const Eigen::Vector2d &r)
	{
		igl::predicates::exactinit();
		return (int)igl::predicates::orient2d(p, q, r);
	}

	inline bool point_on_segment2d(const Eigen::Vector2d &p, const Eigen::Vector2d &a, const Eigen::Vector2d &b, bool know_colinear)
	{
		if (!know_colinear)
		{
			if (orient2d(p, a, b) != 0)
			{
				return false;
			}
		}
		// colinear

		// on one of the point.
		if (p == a || p == b)
		{
			return true;
		}
		// degeneration
		if (a == b)
		{
			return false;
		}
		// find a point not on the line
		Eigen::Vector2d r = Eigen::Vector2d::Random();
		if (orient2d(r, a, b) == 0)
		{
			while (1)
			{
				r = Eigen::Vector2d::Random();
				if (orient2d(r, a, b) != 0)
				{
					break;
				}
			}
		}
		int o1 = orient2d(p, r, a);
		int o2 = orient2d(p, b, r);
		if (o1 == o2) // p is on the same side of the two edges
		{
			return true;
		}
		return false;
	}

	inline bool segment_segment_intersection2d(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
											   const Eigen::Vector2d &t0, const Eigen::Vector2d &t1)
	{
		if (p0 == t0 || p0 == t1 || p1 == t0 || p1 == t1)
		{
			return true;
		}
		std::array<int, 4> os;
		os[0] = orient2d(p0, t0, t1);
		os[1] = orient2d(p1, t0, t1);
		os[2] = orient2d(t0, p0, p1);
		os[3] = orient2d(t1, p0, p1);

		if (os[0] == 0)
		{
			if (point_on_segment2d(p0, t0, t1, true))
			{
				return true;
			}
		}
		if (os[1] == 0)
		{
			if (point_on_segment2d(p1, t0, t1, true))
			{
				return true;
			}
		}
		if (os[2] == 0)
		{
			if (point_on_segment2d(t0, p0, p1, true))
			{
				return true;
			}
		}
		if (os[3] == 0)
		{
			if (point_on_segment2d(t1, p0, p1, true))
			{
				return true;
			}
		}
		if (os[0] * os[1] * os[2] * os[3] == 0)
		{
			return false;
		}
		// check if the two segments crosses
		if (os[0] * os[1] == -1 && os[2] * os[3] == -1)
		{
			return true;
		}
		return false;
	}
	// 0: no intersectin, 1: intersect interior. 2. intersect q0 or q1. 3. the point p0 is on the segment
	// the code is only for the purpose of our code, not exactly solving the generic problems
	// the input p0 p1 defines the ray, p1 is far away enough
	int ray_segment_intersection2d(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
								   const Eigen::Vector2d &q0, const Eigen::Vector2d &q1)
	{
		if (p0 == q0 || p0 == q1)
		{
			return 3;
		}
		int o1 = orient2d(p0, p1, q0);
		int o2 = orient2d(p0, p1, q1);
		if (o1 * o2 > 0) // q0, q2 on the same side of the ray. no intersection
		{
			return 0;
		}
		if (o1 * o2 == 0)
		{
			return 2;
		}
		if (point_on_segment2d(p0, q0, q1, false))
		{
			return 3;
		}
		return segment_segment_intersection2d(p0, p1, q0, q1);
	}
	// p0, p1 defines a ray (actually it is a segment, p1 is far enough from p0),
	// -1: need to shoot another ray.
	// 0: not intersect
	// 1: intersect
	int ray_source_in_loop(const Eigen::Vector3d &p0in, const Eigen::Vector3d &p1in, const Eigen::MatrixXd &Vl)
	{
		Eigen::Vector2d p0(p0in[0], p0in[1]);
		Eigen::Vector2d p1(p1in[0], p1in[1]);
		int ni = 0;
		int vnbr = Vl.rows();
		for (int i = 0; i < vnbr; i++)
		{
			int id0 = i;
			int id1 = (i + 1) % vnbr;
			Eigen::Vector2d q0(Vl(id0, 0), Vl(id0, 1));
			Eigen::Vector2d q1(Vl(id1, 0), Vl(id1, 1));
			int check = ray_segment_intersection2d(p0, p1, q0, q1);
			if (!check)
			{
				continue;
			}
			if (check == 1)
			{
				ni++;
			}
			if (check == 2)
			{
				return -1; // shoot on the vertices of the segment. need another shoot.
			}
			if (check == 3) // the sorce point p0 is on the segment, the point is inside the loop
			{
				return 1;
			}
		}
		return ni % 2;
	}
	Eigen::Vector2d findRandomPointOutCircle(const Eigen::Vector3d &bcIn, const double r)
	{
		double radius = (Eigen::Vector3d::Random()[0] + 1) / 2 * 3.1415926;
		Eigen::Vector2d result;
		Eigen::Vector2d bc(bcIn[0], bcIn[1]);
		return Eigen::Vector2d(r * (1.001) * cos(radius), r * (1.001) * sin(radius)) + bc;
	}
	std::vector<int> pointsInLoop(const Eigen::MatrixXd &V, const Eigen::MatrixXd &Vl)
	{
		Eigen::Vector3d vmin(Vl.col(0).minCoeff(), Vl.col(1).minCoeff(), Vl.col(2).minCoeff());
		Eigen::Vector3d vmax(Vl.col(0).maxCoeff(), Vl.col(1).maxCoeff(), Vl.col(2).maxCoeff());
		std::cout << "min and max of the scene: \n"
				  << vmin.transpose() << ", \n"
				  << vmax.transpose() << "\n";
		Eigen::Vector3d centroid = (vmin + vmax) / 2;
		double r = (vmax - vmin).norm() / 2 * (1.001); // radius
		std::vector<int> inLoop;
		for (int i = 0; i < V.rows(); i++)
		{
			Eigen::Vector3d p = V.row(i);
			if (p[0] < vmin[0] || p[1] < vmin[1] || p[0] > vmax[0] || p[1] > vmax[1])
			{
				continue;
			}
			int check = 0;
			while (1)
			{
				Eigen::Vector2d p12d = findRandomPointOutCircle(centroid, r);
				check = ray_source_in_loop(p, Eigen::Vector3d(p12d[0], p12d[1], 0), Vl);
				if (check == -1)
				{
					continue;
				}
				break;
			}
			if (check)
			{
				inLoop.push_back(i);
			}
		}
		return inLoop;
	}
	// classify the original quads into preserved full quads or quads need to be cutted
	// fvin records if the four vertices of a quad are inside or outside.
	// flags records the inner vertices as true.
	void classifyOriginalQuads(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<int> &pIn,
							   std::vector<int> &quadFull, std::vector<int> &quadCut, std::vector<std::array<bool, 4>> &fvin,
							   std::vector<bool> &flags)
	{
		int fnbr = F.rows();
		int vnbr = V.rows();
		flags = std::vector<bool>(vnbr, false);
		for (int i = 0; i < pIn.size(); i++)
		{
			flags[pIn[i]] = true;
		}
		// std::cout<<"cutted fids\n";
		for (int i = 0; i < fnbr; i++)
		{
			int counter = 0;
			std::array<bool, 4> inside;
			for (int j = 0; j < 4; j++)
			{
				if (flags[F(i, j)])
				{
					counter++;
					inside[j] = true;
				}
				else
				{
					inside[j] = false;
				}
			}
			if (counter == 4)
			{
				quadFull.push_back(i);
				continue;
			}
			if (counter > 0 && counter < 4)
			{
				quadCut.push_back(i);
				fvin.push_back(inside);
				// std::cout<<i<<", ";
				// if(i==304)
				// {
				//     std::cout << "checking " << i << ", the counter, " << counter << " verin, " << inside[0] << ", " << inside[1] << ", " << inside[2]
				//               << ", " << inside[3] << "\n";
				//     std::cout<<"the face, "<< F.row(i)<<"\n";
				// }
			}
		}
		std::cout << "\n";
	}
	void matrix2Mesh(CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
	{
		mesh.clear();
		int nv = V.rows();
		std::vector<CGMesh::VertexHandle> vhs;
		for (int i = 0; i < nv; i++)
			vhs.push_back(mesh.add_vertex(CGMesh::Point(V(i, 0), V(i, 1), V(i, 2))));

		int nf = F.rows();
		int fv = F.cols();
		for (int i = 0; i < nf; i++)
		{
			std::vector<CGMesh::VertexHandle> face_vhandles;
			face_vhandles.push_back(vhs[F(i, 0)]);
			face_vhandles.push_back(vhs[F(i, 1)]);
			face_vhandles.push_back(vhs[F(i, 2)]);
			if (fv == 4)
			{
				face_vhandles.push_back(vhs[F(i, 3)]);
			}
			mesh.add_face(face_vhandles);
		}
	}

	void remapVertices(const int vnbr, const int vnbrAll, const Eigen::MatrixXd &Vall,
					   const std::vector<bool> &flags, std::vector<int> &map, Eigen::MatrixXd &Vclean)
	{
		int counter = 0;
		map = std::vector<int>(vnbrAll, -1);
		for (int i = 0; i < flags.size(); i++)
		{
			if (flags[i])
			{
				map[i] = counter;
				counter++;
			}
		}
		for (int i = flags.size(); i < vnbrAll; i++)
		{
			map[i] = counter;
			counter++;
		}
		Vclean.resize(counter, 3);
		for (int i = 0; i < map.size(); i++)
		{
			if (map[i] >= 0)
			{
				Vclean.row(map[i]) = Vall.row(i);
			}
		}
	}

	void cleanFaceList(const std::vector<int> &map, const std::vector<int> &quadFull, const Eigen::MatrixXi &F,
					   std::vector<std::array<int, 4>> &Fclean)
	{
		Fclean.resize(quadFull.size());
		for (int i = 0; i < quadFull.size(); i++)
		{
			Eigen::Vector4i qd = F.row(quadFull[i]);
			Fclean[i] = {map[qd[0]], map[qd[1]],
						 map[qd[2]], map[qd[3]]};
			if (map[qd[0]] < 0 || map[qd[1]] < 0 ||
				map[qd[2]] < 0 || map[qd[3]] < 0)
			{
				std::cout << "wrong in map!\n";
				exit(0);
			}
		}
	}
	template <typename T>
	void cleanFaceList(const std::vector<int> &map, std::vector<T> &faces, const int order)
	{
		std::vector<T> fcp = faces;
		for (int i = 0; i < faces.size(); i++)
		{
			for (int j = 0; j < order; j++)
			{
				fcp[i][j] = map[faces[i][j]];
				if (map[faces[i][j]] < 0)
				{
					std::cout << "error in templated face cleaner!\n";
					exit(0);
				}
			}
		}
		faces = fcp;
	}
	template <typename T>
	void cleanFaceList(const std::vector<int> &map, std::vector<T> &faces, const std::vector<int> &order)
	{
		std::vector<T> fcp = faces;
		for (int i = 0; i < faces.size(); i++)
		{
			for (int j = 0; j < order[i]; j++)
			{
				fcp[i][j] = map[faces[i][j]];
				if (map[faces[i][j]] < 0)
				{
					std::cout << "error in templated face cleaner!\n";
					exit(0);
				}
			}
		}
		faces = fcp;
	}
	void cleanLoopList(const std::vector<int> &map, std::vector<int> &loop)
	{
		std::vector<int> fcp = loop;
		for (int i = 0; i < loop.size(); i++)
		{
			fcp[i] = map[loop[i]];
		}
		loop = fcp;
	}
	Eigen::MatrixXd cleanVWithMap(const Eigen::MatrixXd &V, const std::vector<int> &vmap)
	{
		std::vector<Eigen::Vector3d> result;
		for (int i = 0; i < V.rows(); i++)
		{
			if (vmap[i] >= 0)
				result.push_back(Eigen::Vector3d(V.row(i)));
		}
		return vec_list_to_matrix(result);
	}
	template <typename T>
	void writeObjPly(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
					 const std::vector<T> &F, const int order)
	{
		std::ofstream fout;
		fout.open(filename);
		for (int i = 0; i < V.rows(); i++)
		{
			if (vmap[i] >= 0)
				fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
		}
		for (int i = 0; i < F.size(); i++)
		{
			fout << "f";

			for (int j = 0; j < F[i].size(); j++)
			{
				fout << " " << F[i][j] + 1;
			}
			fout << "\n";
		}
		fout.close();
	}
	void writeObjPly(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
					 const std::vector<int> &F)
	{
		std::ofstream fout;
		fout.open(filename);
		for (int i = 0; i < V.rows(); i++)
		{
			if (vmap[i] >= 0)
				fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
		}
		for (int i = 0; i < F.size(); i++)
		{
			fout << "f";

			fout << "," << F[i] + 1;

			fout << "\n";
		}
		fout.close();
	}

	void writeObjLoop(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
					  const std::vector<int> &Loop)
	{
		std::ofstream fout;
		fout.open(filename);
		for (int i = 0; i < V.rows(); i++)
		{
			if (vmap[i] >= 0)
				fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
		}
		for (int i = 0; i < Loop.size() - 1; i++)
		{
			fout << "l " << Loop[i] + 1 << " " << Loop[i + 1] + 1;

			fout << "\n";
		}
		fout << "l " << Loop.back() + 1 << " " << Loop.front() + 1;
		fout.close();
	}
	Eigen::Vector3d edge_loop_intersectionPoint(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
												const Eigen::Vector3d &pp0, const Eigen::Vector3d &pp1,
												const Eigen::MatrixXd &Vl, bool &found, int &round)
	{
		int vnbr = Vl.rows();
		int eid = -1;
		found = false;
		Eigen::Vector3d result(-1, -1, -1);
		round = -1;
		for (int i = 0; i < vnbr; i++)
		{
			int id0 = i;
			int id1 = (i + 1) % vnbr;
			Eigen::Vector2d q0(Vl(id0, 0), Vl(id0, 1));
			Eigen::Vector2d q1(Vl(id1, 0), Vl(id1, 1));
			bool check = segment_segment_intersection2d(p0, p1, q0, q1);
			if (check)
			{
				if (eid >= 0)
				{
					std::cout << "this edge is intersecting the loop with more than 1 point!\n";
					exit(0);
				}
				eid = i;
				double t0, t1;
				Eigen::Vector3d p03(p0[0], p0[1], 0);
				Eigen::Vector3d q03(q0[0], q0[1], 0);
				Eigen::Vector3d p13(p1[0], p1[1], 0);
				Eigen::Vector3d q13(q1[0], q1[1], 0);
				bool intersect = igl::segment_segment_intersect(p03, p13 - p03, q03, q13 - q03, t0, t1);
				if (!intersect)
				{
					std::cout << "numerical problem happens, we better implement this function by ourselves\n";
					exit(0);
				}
				if (t0 < 1e-6)
				{
					t0 = 0;
					round = 0;
				}
				if (t0 > 1 - 1e-6)
				{
					t0 = 1;
					round = 1;
				}
				Eigen::Vector3d i3 = pp0 + t0 * (pp1 - pp0);
				result = i3;
				found = true;
				////////////// TODO the following "break" may cause problems, please check the logic here if problem happens
				break;
			}
		}
		return result;
	}
	void findALlIntersections(const CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
							  const Eigen::MatrixXd &Vl, std::vector<CGMesh::EdgeHandle> &ehs,
							  std::vector<Eigen::Vector3d> &intersections)
	{
		for (CGMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
		{
			OpenMesh::EdgeHandle eh = e_it.handle();
			OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(eh, 0);
			int vid0 = mesh.from_vertex_handle(hh).idx();
			int vid1 = mesh.to_vertex_handle(hh).idx();
			Eigen::Vector2d p0(V(vid0, 0), V(vid0, 1)), p1(V(vid1, 0), V(vid1, 1));
			bool found = false;
			int round;
			Eigen::Vector3d intersection = edge_loop_intersectionPoint(p0, p1, V.row(vid0), V.row(vid1), Vl, found, round);
			if (found)
			{
				ehs.push_back(eh);
				intersections.push_back(intersection);
			}
		}
	}

	bool verInQuad(const int vid, const Eigen::Vector4i &face)
	{
		if (vid == face[0] || vid == face[1] || vid == face[2] || vid == face[3])
		{
			return true;
		}
		return false;
	}

	// -1: the edge is not cutted.
	// other: the intersection point id (in the intersection list) on the edge
	int edgeCutted(const int vf, const int vt, const Eigen::Vector4i &quad, CGMesh &mesh, const std::vector<CGMesh::EdgeHandle> &ehs,
				   const std::array<bool, 4> &fvin)
	{
		int result = -1;
		bool in0 = false, in1 = false;
		bool f0 = false, f1 = false;
		for (int i = 0; i < 4; i++)
		{
			if (vf == quad[i])
			{
				in0 = fvin[i];
				f0 = true;
			}
			if (vt == quad[i])
			{
				in1 = fvin[i];
				f1 = true;
			}
		}
		if (!f0 || !f1)
		{
			std::cout << "the edge is not in the quad!!!\n";
			std::cout << "the edge: " << vf << ", " << vt << ", the quad, " << quad.transpose() << "\n";
			exit(0);
		}
		if (in0 == in1)
		{
			return -1;
		}
		// this edge is cutted, find the corresponding halfedge.
		bool found = false;
		int counter = 0;
		for (auto eh : ehs)
		{
			CGMesh::HalfedgeHandle hh = mesh.halfedge_handle(eh, 0);
			int v0 = mesh.from_vertex_handle(hh).idx();
			int v1 = mesh.to_vertex_handle(hh).idx();
			bool c1 = vf == v0 && vt == v1;
			bool c2 = vf == v1 && vt == v0;
			if (c1 || c2)
			{
				found = true;
				return counter;
			}
			counter++;
		}
		if (!found)
		{
			std::cout << "ERROR: the cutted edge is not found in the list\n";
			exit(0);
		}
	}

	// from a cut point to another cut point
	std::vector<int> fromCutPtToCutPt(const std::vector<int> &vloop, const int start, const int vnbr)
	{
		std::vector<int> result;
		if (vloop[start] < vnbr)
		{
			std::cout << "please take a proper start point!\n";
			exit(0);
		}
		result.push_back(vloop[start]);
		bool found = false;
		for (int i = 1; i < vloop.size(); i++)
		{
			int id = (i + start) % vloop.size();
			result.push_back(vloop[id]);
			if (vloop[id] >= vnbr)
			{
				found = true;
				break;
			}
		}
		if (!found)
		{
			std::cout << "the loop end point is not found!\n";
			exit(0);
		}
		return result;
	}

	// pick a inner polygon from the two polygons of the cutted quad
	// Vflags: size is vnbr, showing which is inside.
	// be: boundary edge
	void classifyVloop(const std::vector<bool> &Vflags, const std::vector<int> &vloop, const int vnbr,
					   std::vector<int> &poly, std::array<int, 2> &be)
	{
		std::array<std::vector<int>, 2> ps;
		int counter = 0;
		for (int i = 0; i < vloop.size(); i++)
		{
			if (vloop[i] >= vnbr) // this is a cut point
			{
				ps[counter] = fromCutPtToCutPt(vloop, i, vnbr);
				be[counter] = vloop[i];
				counter++;
			}
		}
		if (Vflags[ps[0][1]])
		{
			poly = ps[0];
		}
		else
		{
			poly = ps[1];
			if (!Vflags[ps[1][1]])
			{
				std::cout << "Both the two polygons are outside!\n";
				std::cout << "the vnbr " << vnbr << "\n";
				for (auto r : ps[0])
				{
					std::cout << r << ", (" << Vflags[r] << "), ";
				}
				std::cout << "\n";
				for (auto r : ps[1])
				{
					std::cout << r << ", (" << Vflags[r] << "), ";
				}
				std::cout << "\n";
				exit(0);
			}
		}
		return;
	}

	void sortBoundaryEdgesIntoLoop(const std::vector<std::array<int, 2>> &bes, std::vector<int> &loop)
	{
		int enbr = bes.size();
		std::vector<bool> echecked(enbr, false);
		loop.push_back(bes[0][0]);
		loop.push_back(bes[0][1]);
		echecked[0] = true;
		int nbrChecked = 1;
		while (1)
		{
			if (nbrChecked == enbr - 1)
			{
				break;
			}
			int id = loop.back();
			for (int i = 0; i < enbr; i++)
			{
				if (echecked[i])
				{
					continue;
				}
				if (id == bes[i][0] || id == bes[i][1])
				{
					if (id == bes[i][0])
					{
						loop.push_back(bes[i][1]);
					}
					else
					{
						loop.push_back(bes[i][0]);
					}
					echecked[i] = true;
					nbrChecked++;
				}
			}
		}
	}
	// bLoop is the boundary loop
	// split the quads into
	void splitIntersectedQuads(CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
							   const Eigen::MatrixXd &Vl, const std::vector<CGMesh::EdgeHandle> &ehs,
							   const std::vector<int> &quadCut, const std::vector<std::array<bool, 4>> &fvin,
							   const std::vector<bool> &Vflags, std::vector<std::array<int, 4>> &rQuads,
							   std::vector<std::array<int, 3>> &rTris, std::vector<std::array<int, 5>> &rPenta,
							   std::vector<int> &bLoop)
	{
		int vnbr = V.rows();
		// the generated quads, triangles and pentagons.

		std::vector<std::array<int, 2>> bEdges; // boundary edges

		for (int i = 0; i < quadCut.size(); i++)
		{
			CGMesh::FaceHandle fh = mesh.face_handle(quadCut[i]);
			std::vector<int> vloop; // the loop for the face
			// iterate over the 4 edges, if there is a cut, add the vertex into the loop.
			CGMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(fh);
			OpenMesh::HalfedgeHandle hh = fh_it.handle();

			for (int j = 0; j < 4; j++)
			{
				hh = mesh.next_halfedge_handle(hh);
				int vfrom = mesh.from_vertex_handle(hh).idx();
				int vto = mesh.to_vertex_handle(hh).idx();
				for (int c : vloop)
				{
					if (c == vfrom)
					{
						std::cout << "the vertex is already in \nwhich edge " << j << "\n";
						exit(0);
					}
				}

				vloop.push_back(vfrom);
				bool notin = false;
				if (!verInQuad(vfrom, F.row(quadCut[i])))
				{
					std::cout << "ver is not in!\n"
							  << vfrom << "\n";
					notin = true;
				}
				if (!verInQuad(vto, F.row(quadCut[i])))
				{
					std::cout << "ver is not in!\n"
							  << vfrom << "\n";
					notin = true;
				}
				if (notin)
				{
					std::cout << "face , " << F.row(quadCut[i]) << "\nwhich edge, " << j << "\n";
					std::cout << "the face vertices itr:\n";
					for (CGMesh::FaceVertexIter fvitr = mesh.fv_begin(fh); fvitr != mesh.fv_end(fh); ++fvitr)
					{
						std::cout << fvitr.handle().idx() << ", ";
					}
					std::cout << "\n";
				}

				int iid = edgeCutted(vfrom, vto, F.row(quadCut[i]), mesh, ehs, fvin[i]);
				if (iid < 0)
				{ // not cutted, continue;
					continue;
				}
				else
				{
					vloop.push_back(iid + vnbr);
				}
			}

			// for (CGMesh::FaceEdgeIter fh_it = mesh.fe_begin(fh); fh_it != mesh.fe_end(fh); ++fh_it)
			// {
			//     OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(fh_it, 0);
			//     int vfrom = mesh.from_vertex_handle(hh).idx();
			//     int vto = mesh.to_vertex_handle(hh).idx();
			//     for(int c : vloop)
			//     {
			//         if(c==vfrom)
			//         {
			//             std::cout<<"the vertex is already in \n";
			//             exit(0);
			//         }
			//     }

			//     vloop.push_back(vfrom);

			//     if(!verInQuad(vfrom, F.row(quadCut[i])))
			//     {
			//         std::cout<<"ver is not in!\n";
			//         exit(0);
			//     }
			//     if(!verInQuad(vto, F.row(quadCut[i])))
			//     {
			//         std::cout<<"ver is not in!\n";
			//         exit(0);
			//     }

			//     int iid = edgeCutted(vfrom, vto, F.row(quadCut[i]), mesh, ehs, fvin[i]);
			//     if (iid < 0)
			//     { // not cutted, continue;
			//         continue;
			//     }
			//     else
			//     {
			//         vloop.push_back(iid + vnbr);
			//     }
			// }
			if (vloop.size() != 6)
			{
				std::cout << "The vloop has " << vloop.size() << " vertices!\n";
				exit(0);
				continue; // TODO This may skip some intersected shapes on the boundary, but resolution wise it is fine.
			}
			// classify the cutted quad into two polygons. keep the inner one.
			std::vector<int> poly;
			std::array<int, 2> be;
			bool foundin = false;
			for (int vid : vloop)
			{
				if (vid < vnbr)
				{
					if (Vflags[vid])
					{
						foundin = true;
					}
				}
			}
			if (!foundin)
			{
				std::cout << "error: the cutted quad doesn't have any inside ver, fid, " << quadCut[i] << "\n";
				std::cout << "face, " << F.row(quadCut[i]) << "\n";
				std::cout << "loop: \n";
				for (int vid : vloop)
				{
					std::cout << vid << ", ";
				}
				std::cout << "\n";

				// exit(0);
			}
			classifyVloop(Vflags, vloop, vnbr, poly, be);
			bEdges.push_back(be);
			int psize = poly.size();
			if (psize < 3 || psize > 5)
			{
				std::cout << "psize wrong!!!\n";
				exit(0);
			}
			if (psize == 3)
			{
				rTris.push_back({poly[0], poly[1], poly[2]});
			}
			if (psize == 4)
			{
				rQuads.push_back({poly[0], poly[1], poly[2], poly[3]});
			}
			if (psize == 5)
			{
				rPenta.push_back({poly[0], poly[1], poly[2], poly[3], poly[4]});
			}
		}
		// sort the boundary edges into a loop
		sortBoundaryEdgesIntoLoop(bEdges, bLoop);
	}
	// returns the halfedge id of the edge v0-v1
	int whichHalfedge(const int v0, const int v1, CGMesh &mesh, const std::vector<CGMesh::EdgeHandle> &ehs)
	{

		for (int i = 0; i < ehs.size(); i++)
		{
			CGMesh::HalfedgeHandle hh = mesh.halfedge_handle(ehs[i], 0);
			int vf = mesh.from_vertex_handle(hh).idx();
			int vt = mesh.to_vertex_handle(hh).idx();
			if ((v0 == vf && v1 == vt) || (v0 == vt && v1 == vf))
			{
				return i;
			}
		}
		std::cout << "We cannot find the halfedge!\n";
		exit(0);
	}
	void findConnectivityForInnerVertices(const std::vector<bool> &Vflags, const std::vector<int> &pin, const int vnbr,
										  CGMesh &mesh, const std::vector<CGMesh::EdgeHandle> &ehs,
										  std::vector<std::array<int, 5>> &connect)
	{
		for (int i = 0; i < pin.size(); i++)
		{
			CGMesh::VertexHandle vh = mesh.vertex_handle(pin[i]);
			std::vector<int> nbs; // neighbours
			for (CGMesh::VertexOHalfedgeIter voh_it = mesh.voh_begin(vh); voh_it != mesh.voh_end(vh); ++voh_it)
			{
				CGMesh::HalfedgeHandle hh = voh_it.handle();
				int idout = mesh.to_vertex_handle(hh).idx();
				if (idout == pin[i])
				{
					std::cout << "from and to are the same!\n";
					exit(0);
				}
				if (Vflags[idout]) // this vertex is inside the boundary loop, add it into the neighbours
				{
					nbs.push_back(idout);
				}
				else // this vertex is outside, then check which halfedge it is.
				{
					int vid = whichHalfedge(pin[i], idout, mesh, ehs) + vnbr;
					nbs.push_back(vid);
				}
			}
			if (nbs.size() == 4) // record only the valence 4 vertices
				connect.push_back({pin[i], nbs[0], nbs[1], nbs[2], nbs[3]});
		}
	}

	Eigen::MatrixXd connectVertexList(const Eigen::MatrixXd &V0, const std::vector<Eigen::Vector3d> &V1)
	{
		int vnbr = V0.rows() + V1.size();
		Eigen::MatrixXd V(vnbr, 3);
		V.topRows(V0.rows()) = V0;
		V.bottomRows(V1.size()) = vec_list_to_matrix(V1);
		return V;
	}
	Eigen::MatrixXi polygons2Triangles(const std::vector<std::vector<int>> &polys)
	{
		auto splitPly = [](const std::vector<int> &p)
		{
			std::vector<std::vector<int>> result;
			if(p.size()==3)
			{
				result.push_back(p);
				return result;
			}
			for (int i = 1; i < p.size() - 1; i++)
			{
				std::vector<int> tmp;
				tmp.push_back(p[0]);
				tmp.push_back(p[i]);
				tmp.push_back(p[i + 1]);
				result.push_back(tmp);
			}
			return result;
		};
		std::vector<Eigen::Vector3i> faces;
		faces.reserve(polys.size() * 3);
		for (int i = 0; i < polys.size(); i++)
		{
			std::vector<std::vector<int>> splitted = splitPly(polys[i]);
			for(auto tri : splitted)
			{
				faces.push_back(Eigen::Vector3i(tri[0], tri[1],tri[2]));
			}
		}
		Eigen::MatrixXi F(faces.size(), 3);
		for(int i=0;i<faces.size();i++)
		{
			F.row(i) = faces[i];
		}
		return F;
	}
	void cutBoundaryGenerateTopology(const Eigen::MatrixXd &Vquad, const Eigen::MatrixXi &Fquad, const Eigen::MatrixXd &Vcurve, Eigen::MatrixXd &Vout,
									 Eigen::MatrixXi &Fout)
	{
		int vnbr = Vquad.rows();
		int fnbr = Fquad.rows();
		Eigen::MatrixXd Vproj = Vquad; // prject the vertices onto 2d
		Vproj.col(2) = Eigen::VectorXd::Zero(vnbr);
		std::cout << "The quad has " << Vquad.rows() << " vertices\n";
		// mark the points that are inside the loop.
		std::vector<int> pin = pointsInLoop(Vproj, Vcurve);
		std::cout << "there are " << pin.size() << " vertices in the circle\n";

		Eigen::MatrixXd Vin(pin.size(), 3);
		for (int i = 0; i < pin.size(); i++)
		{
			Vin.row(i) = Vquad.row(pin[i]);
		}
		std::vector<int> quadFull;
		std::vector<int> quadCut;
		std::vector<std::array<bool, 4>> fvin;
		std::vector<bool> vFlags;
		// record the quads that are cutted.
		classifyOriginalQuads(Vquad, Fquad, pin, quadFull, quadCut, fvin, vFlags);
		std::cout << "there are " << quadFull.size() << " quads are full, and " << quadCut.size() << " are cutted\n";
		CGMesh mesh;
		matrix2Mesh(mesh, Vquad, Fquad);
		std::vector<CGMesh::EdgeHandle> ehs;
		std::vector<Eigen::Vector3d> intersections;
		// find all the intersections on the edges.
		findALlIntersections(mesh, Vquad, Fquad, Vcurve, ehs, intersections);
		std::cout << "find intersection point nbr, " << ehs.size() << "\n";
		// get the cutted polygons and the boundary loop
		// the triangles, the quads and the pentagons
		std::vector<std::array<int, 3>> rT;
		std::vector<std::array<int, 4>> rQ;
		std::vector<std::array<int, 5>> rP;
		std::vector<std::vector<int>> polys;
		std::vector<int> bLoop;
		splitIntersectedQuads(mesh, Vquad, Fquad, Vcurve, ehs, quadCut, fvin, vFlags, rQ, rT, rP, bLoop);
		std::cout << "new generated quads, " << rQ.size() << ", trians, " << rT.size() << ", pentas, " << rP.size() << "\n";
		std::vector<std::array<int, 5>> connect;
		// get the connection for all the inner vertices.
		findConnectivityForInnerVertices(vFlags, pin, vnbr, mesh, ehs, connect);
		Eigen::MatrixXd Vall = connectVertexList(Vquad, intersections);

		// clean the data
		std::vector<int> map;
		Eigen::MatrixXd Vclean;
		remapVertices(vnbr, Vall.rows(), Vall, vFlags, map, Vclean);
		std::vector<std::array<int, 4>> quadsClean;
		// clean the face lists
		cleanFaceList(map, quadFull, Fquad, quadsClean);
		cleanFaceList(map, rQ, 4);
		cleanFaceList(map, rT, 3);
		cleanFaceList(map, rP, 5);
		// clean the connectivity
		cleanFaceList(map, connect, 5);
		// clean the boundary loop
		cleanLoopList(map, bLoop);

		auto extendPolygon = []<typename T>(std::vector<T> &inshape, std::vector<std::vector<int>> &extended)
		{
			for (auto ply : inshape)
			{
				std::vector<int> plg;
				for (auto ver : ply)
				{
					plg.push_back(ver);
				}
				extended.push_back(plg);
			}
		};
		// obtain all the polygons
		extendPolygon(quadsClean, polys);
		extendPolygon(rQ, polys);
		extendPolygon(rT, polys);
		extendPolygon(rP, polys);
		Eigen::MatrixXi Tri = polygons2Triangles(polys);
		Vout = cleanVWithMap(Vall, map);
		Fout = Tri;
		std::cout<<"After cutting, there are "<<Vout.rows()<<" vertices\n";
		// write them

		// writeObjPly(crpcFolder + "rQF.obj", Vall, map, quadsClean, 4);
		// writeObjPly(crpcFolder + "rQ.obj", Vall, map, rQ, 4);
		// writeObjPly(crpcFolder + "rT.obj", Vall, map, rT, 3);
		// writeObjPly(crpcFolder + "rP.obj", Vall, map, rP, 5);
		// writeObjPly(crpcFolder + "Poly.obj", Vall, map, polys, -1);
		// writeObjPly(crpcFolder + "connect.txt", Vall, map, connect, 5);
		// writeObjLoop(crpcFolder + "loop.obj", Vall, map, bLoop);
		// writeObjPly(crpcFolder + "loop.txt", Vall, map, bLoop);
		// std::cout << "connection size " << connect.size() << "\n";
		// writeSomePoints(foldername + "test.obj", Vall, map, connect, 766);

		// test and plot the connectivity

		// write the points inside
		// igl::writeOBJ(foldername + "pIn.obj", Vin, Eigen::MatrixXi());
	}
	void write_csv(const std::string &file, const std::vector<double> data)
	{
		std::ofstream fout;
		fout.open(file);
		for (int i = 0; i < data.size() - 1; i++)
		{
			fout << data[i] << ",";
		}
		fout << data.back() << std::endl;
		fout.close();
	}
}