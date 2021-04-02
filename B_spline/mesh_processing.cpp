#include <mesh_processing.h>
#include<array>
#include<curve.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/write_triangle_mesh.h>
#include <queue>
#include <igl/predicates/predicates.h>
//#include <igl/triangle/cdt.h>
//#include<igl/triangle/triangulate.h>
Eigen::MatrixXd list_to_matrix_3d(const std::vector<std::vector<double>>& v) {
	Eigen::MatrixXd result(v.size(), 3);
	for (int i = 0; i < v.size(); i++) {
		result(i, 0) = v[i][0];
		result(i, 1) = v[i][1];
		result(i, 2) = v[i].size() >= 3 ? v[i][2] : 0;
	}
	return result;
}
Eigen::MatrixXd list_to_matrix_3d(const std::vector<Vector3d>& v, const std::vector<int>& selected) {
	Eigen::MatrixXd result(selected.size(), 3);
	for (int i = 0; i < selected.size(); i++) {
		result(i, 0) = v[selected[i]][0];
		result(i, 1) = v[selected[i]][1];
		result(i, 2) = v[selected[i]].size() >= 3 ? v[selected[i]][2] : 0;
	}
	return result;

}

Eigen::MatrixXd vector_to_matrix_3d(const std::vector<Vector3d>& v) {
	Eigen::MatrixXd result(v.size(), 3);
	for (int i = 0; i < v.size(); i++) {
		result.row(i) = v[i];
	}
	return result;
}
std::vector<Vector3d> matrix3d_to_vector(const Eigen::MatrixXd& v) {
	std::vector<Vector3d> result(v.rows());
	for (int i = 0; i < v.rows(); i++) {
		result[i] = v.row(i);
	}
	return result;
}


void vertices_to_edges(const Eigen::MatrixXd& pts, Eigen::MatrixXi &edges) {
	edges.resize(pts.rows() - 1, 2);
	for (int i = 0; i < edges.rows(); i++) {
		edges(i, 0) = i; edges(i, 1) = i + 1;
	}
}
void test_read_mesh(const std::string &filename) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	bool read = igl::read_triangle_mesh(filename, V, F);

}

void read_and_visual_mesh(const std::string &filename,Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	bool read = igl::read_triangle_mesh(filename, V, F);
}

void map_vertices_to_square(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd &bnd_uv){
	bnd_uv.resize(bnd.size(), 2);
	double total_dist = 0;
	std::vector<double> pdis(bnd.size());
	pdis[0] = (V.row(bnd[1]) - V.row(bnd[0])).norm();
	for (int i = 1; i < pdis.size() - 1; i++) {
		pdis[i] = (V.row(bnd[i + 1]) - V.row(bnd[i])).norm() + pdis[i - 1];
		
	}
	pdis[pdis.size() - 1] = pdis[pdis.size() - 2] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();
	total_dist = pdis[pdis.size() - 1];

	int k = 1;
	double quarter = total_dist / 4;
	std::array<double, 4> corners;
	corners[0] = 0;
	for (int i = 0; i < pdis.size()-1; i++) {
		if (pdis[i] <= quarter * k&&pdis[i + 1] > quarter*k) {
			if (i != corners[k - 1]) {
				corners[k] = i;
			}
			else {
				corners[k] = i + 1;
			}
			k += 1;
			if (k == 4) break;
		}
	}
	// u=0, v increase
	bnd_uv.row(0) << 0, 0;
	for (int i = 0; i < corners[1]; i++) {
		bnd_uv.row(i + 1) << 0, pdis[i] / pdis[corners[1]-1];
	}

	// v=1, u increase
	for (int i =corners[1]; i < corners[2]; i++) {
		bnd_uv.row(i + 1) << (pdis[i] - pdis[corners[1] - 1])/( pdis[corners[2] - 1]- pdis[corners[1] - 1]), 1;
	}

	// u=1, v decrease
	for (int i = corners[2]; i < corners[3]; i++) {
		bnd_uv.row(i + 1) << 1, 1 - (pdis[i] - pdis[corners[2] - 1]) / (pdis[corners[3] - 1] - pdis[corners[2] - 1]);
	}

	// v=0, u decrease
	for (int i = corners[3]; i < pdis.size()-1; i++) {
		bnd_uv.row(i + 1) << 1 - (pdis[i] - pdis[corners[3] - 1]) / (pdis[pdis.size()-1] - pdis[corners[3] - 1]), 0;
	}
	//std::cout << "the parameters\n"<<bnd_uv << std::endl;
	

}
// remove the vertices not used in contructing the surface
// this function does not consider if there are duplicated vertices (two vertices in 
// exact the same position)
void remove_redundent_mesh_vertices(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
	Eigen::MatrixXd& Vout, Eigen::MatrixXi& Fout) {


	std::vector<bool> located(V.rows(), false);
	std::vector<int> map(V.rows());
	int nbr = 0;
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			if (located[F(i, j)] == false) {
				located[F(i, j)] = true;
				map[F(i, j)] = nbr;
				nbr++;
			}

		}
	}
	Fout.resize(F.rows(), 3);
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			Fout(i, j) = map[F(i, j)];
		}
	}
	Vout.resize(nbr, 3);
	int k = 0;
	for (int i = 0; i < V.rows(); i++) {
		if (located[i])
			Vout.row(map[i]) = V.row(i);
	}
}

// detect the duplicated vertices, and re-orginaze the faces to avoid open
// boundary caused by vertices duplication
void remove_duplicated_vertices(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
	Eigen::MatrixXd& Vout, Eigen::MatrixXi& Fout) {
	std::vector<int> map(V.rows(), -1);
	std::vector<bool> duplicated(V.rows(), false);
	int to_delete = 0;
	for (int i = 0; i < V.rows(); i++) {
		for (int j = i; j < V.rows(); j++) {
			if (i < j) {
				if (V.row(i) == V.row(j)) {
					duplicated[j] = true;
					if (duplicated[i] == true) {
						map[j] = map[i];
					}
					else {
						map[j] = i;
						to_delete++;
					}

				}
			}
		}
	}
	Vout = V;
	Fout.resize(F.rows(), 3);
	for (int i = 0; i < Fout.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			if (duplicated[F(i, j)] == true) {
				Fout(i, j) = map[F(i, j)];
			}
			else {
				Fout(i, j) = F(i, j);
			}
		}
	}
}

void generate_clean_mesh_data_for_parametrization(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
	Eigen::MatrixXd& Vout, Eigen::MatrixXi& Fout) {

	// the duplicated vertices will be located, and faces will not use all the vertices
	remove_duplicated_vertices(V, F, Vout, Fout);

	//  remove the unrelated vertices
	remove_redundent_mesh_vertices(Vout, Fout, Vout, Fout);

}
// this is to remove some faces of a mesh accroding to the vertices coordinates.
// axis select from 0, 1, 2. 
void remove_some_faces(const int axis, const double value, const bool remove_larger,
	const Eigen::MatrixXd& V, const Eigen::MatrixXi F, Eigen::MatrixXi& newF) {
	int rows = F.rows();
	newF.resize(1, 3);
	int already_updated = 0;
	for (int i = 0; i < rows; i++) {
		if (remove_larger) {
			if (V(F(i, 0), axis) > value || V(F(i, 1), axis) > value || V(F(i, 2), axis) > value) {
				continue;
			}
		}
		else {
			if (V(F(i, 0), axis) < value || V(F(i, 1), axis) < value || V(F(i, 2), axis) < value) {
				continue;
			}
		}
		newF.conservativeResize(already_updated + 1, 3);
		newF.row(already_updated) = F.row(i);
		already_updated++;

	}

}
void mesh_parameterization(
	const std::string &meshfile, Eigen::MatrixXd& V, Eigen::MatrixXd &param, Eigen::MatrixXi &F) {
	/*const std::string path = "D:\\vs\\sparse_data_interpolation\\meshes\\";
	const std::string filename = path + "camel_smallest.obj";*/
	Eigen::MatrixXd Vori; Eigen::MatrixXi Fori;
	read_and_visual_mesh(meshfile, Vori, Fori);
	
	
	V = Vori; F = Fori;
	/////////////////////////////////////////
	// parameterization part
	Eigen::VectorXi bnd;
	igl::boundary_loop(F, bnd);// boundary vertices detection
	Eigen::MatrixXd bnd_uv;
	/*igl::map_vertices_to_circle(V, bnd, bnd_uv);*/
	map_vertices_to_square(V, bnd, bnd_uv);

	igl::harmonic(V, F, bnd, bnd_uv, 1, param);

}

double shortest_edge_length_of_parametric_domain(const Eigen::MatrixXd& paras, const Eigen::MatrixXi &F) {
	double rst = 1;
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int id0 = F(i, j);
			int id1 = F(i, (j + 1) % 3);
			/*if (paras(id0, 0) == paras(id1, 0)) {
				if (paras(id0, 0) == 0 || paras(id0, 0) == 1) {
					continue;
				}
			}
			*/
			double dist = (paras.row(id0) - paras.row(id1)).norm();
			if (dist < rst) {
				rst = dist;
			}
		}
	}
	return rst;
}
int orient_2d(const Vector2d& a, const Vector2d &b, const Vector2d &c) {
	auto result= igl::predicates::orient2d(a, b, c);
	if (result == igl::predicates::Orientation::POSITIVE) {
		return 1;
	}
	if (result == igl::predicates::Orientation::NEGATIVE) {
		return -1;
	}
	return 0;
}
// check if the perturbed triangles are still hold the orientation
// if the orientations are the same, then this parameter is safe to be perturbed.
bool check_the_perturbed_one_parameter(const Eigen::MatrixXd& para_in, const Eigen::MatrixXd& perturbed, const Eigen::MatrixXi &F,
	const int parid) {
	for (int i = 0; i < F.rows(); i++) {
		int f0 = F(i, 0), f1 = F(i, 1), f2 = F(i, 2);
		if (f0 == parid || f1 == parid || f2 == parid) {
			// processing is below
		}
		else {
			continue;
		}
		
		int original_orientation = orient_2d(Vector2d(para_in.row(f0)), Vector2d(para_in.row(f1)), Vector2d(para_in.row(f2)));
		int orientation = orient_2d(Vector2d(perturbed.row(f0)), Vector2d(perturbed.row(f1)), Vector2d(perturbed.row(f2)));
		return original_orientation == orientation;
	}
}
bool perturb_paras_and_check_validity(Eigen::MatrixXd &paras, const Eigen::MatrixXi& F, const int u_or_v,
	const Eigen::VectorXi order, const Eigen::MatrixXi &info, const int info_ref, const double target_value) {
	Eigen::MatrixXd para_tmp = paras;
	for (int i = 0; i < info.rows(); i++) {
		if (info(i, u_or_v) == info_ref) {
			para_tmp(order[i], u_or_v) = target_value;
		}
	}
	for (int i = 0; i < info.rows(); i++) {
		if (info(i, u_or_v) == info_ref) {
			int para_id = order[i];
			bool able = check_the_perturbed_one_parameter(paras, para_tmp, F, para_id);
			if (!able) {
				return false;
			}
		}
	}
	paras = para_tmp;
	return true;
}


// itr is the number of iterations;
void mesh_parameter_perturbation(const Eigen::MatrixXd &para_in, const Eigen::MatrixXi &F, 
	Eigen::MatrixXd &para_out, int itr) {
	std::cout << "inside perturbation" << std::endl;
	// the the first and second elements are u and v, the third is the reference
	auto cmp_u = [](Vector3d i1, Vector3d i2) {
		return i1[0] > i2[0];
	};
	auto cmp_v = [](Vector3d i1, Vector3d i2) {
		return i1[1] > i2[1];
	};
	std::cout << "check 1" << std::endl;
	std::priority_queue<Vector3d, std::vector<Vector3d>, decltype(cmp_u)> queue_u(cmp_u);
	std::priority_queue<Vector3d, std::vector<Vector3d>, decltype(cmp_v)> queue_v(cmp_v);
	std::cout << "check 1.5" << std::endl;
	for (int i = 0; i < para_in.rows(); i++) {
		Vector3d vecu = Vector3d(para_in(i, 0), para_in(i, 1), i);
		Vector3d vecv = Vector3d(para_in(i, 0), para_in(i, 1), i);
		queue_u.push(vecu);
		queue_v.push(vecv);
	}
	std::cout << "check 2" << std::endl;
	
	Eigen::VectorXi u_order(para_in.rows()), v_order(para_in.rows());
	std::cout << "check 2.5,parasize "<<para_in.rows() << std::endl;
	for (int i = 0; i < para_in.rows(); i++) {
		Vector3d ut = queue_u.top();
		u_order(i) = ut[2];
		Vector3d vt = queue_v.top();
		v_order(i) = vt[2];
		queue_u.pop();
		queue_v.pop();
		std::cout << "i " << i << std::endl;
	}
	std::cout << "check 3" << std::endl;
	//// next get the shortest edge of parameters
	Eigen::MatrixXd para_modified = para_in;
	int rows = para_in.rows();
	 Eigen::MatrixXi pinfo= Eigen::MatrixXi::Constant(rows, 2, -1);// perturbed info
	 std::cout << "before perturbation" << std::endl;
	for (int i = 0; i < itr; i++) {
		double perturb_dist = shortest_edge_length_of_parametric_domain(para_modified, F);
		perturb_dist /= 2;// the perturbation threshold
		
		
		for (int j = 0; j < 2; j++) {// perturb u or v
			for (int k = 1; k < rows - 1; k++) {
				Eigen::VectorXi order;
				if (j == 0) {
					order = u_order;
				}
				else {
					order = v_order;
				}
				int paraid1 = order[k];
				int paraid2 = order[k + 1];
				double value1 = para_modified(paraid1, j);
				double value2 = para_modified(paraid2, j);
				if (value1 == 0 || value1 == 1 || value2 == 0 || value2 == 1) {
					continue;// border points don't perturb
				}
				if (fabs(value1 - value2) < perturb_dist) {
					Eigen::MatrixXd para_tmp = para_modified;
					Eigen::MatrixXi update_info = pinfo;// when trying to clap, we modify update_info. if ok, we update pinfo
					int info_ref;
					int info_small = std::min(update_info(k, j), update_info(k + 1, j));
					int info_big= std::max(update_info(k, j), update_info(k + 1, j));
					if (info_small == -1) {
						if (info_big != -1) {
							info_ref = info_big;
						}
						else {
							info_ref = k;
						}
					}
					else {
						info_ref = info_small;
					}
					double target_value = (value1 + value2) / 2;
					update_info(k, j) = info_ref;
					update_info(k + 1, j) = info_ref;
					bool check = perturb_paras_and_check_validity(para_tmp, F, j, order, update_info, info_ref, target_value);
					if (check) {// if perturbable, reset parameters and info
						para_modified = para_tmp;
						pinfo = update_info;
						std::cout << "perturbed" << std::endl;
					}
				}
			}
		}
		std::cout << "itr finished" << std::endl;

	}
	std::cout << "finished perturbation" << std::endl;
	std::cout << "p info\n" << pinfo << std::endl;
	para_out = para_modified;

}
// the inputs are cleaned(no duplicated) data points.
// if V has 3 cols, then truncate it to 2d.
//void constrained_delaunay_triangulation(
//	const Eigen::MatrixXd& V, const Eigen::VectorXi& Edge_ids,
//	Eigen::MatrixXi& F) {
//	
//	Eigen::MatrixXd Vin(V.rows(), 2);
//	if (V.cols() == 2) {
//		Vin = V;
//	}
//	else {
//		Vin << V.col(0), V.col(1);
//	}
//	Eigen::MatrixXi Ein(Edge_ids.size(), 2);
//	for (int i = 0; i < Edge_ids.size() - 1; i++) {
//		Ein(i, 0) = Edge_ids[i];
//		Ein(i, 1) = Edge_ids[i + 1];
//	}
//	Ein(Edge_ids.size() - 1, 0) = Edge_ids[Edge_ids.size() - 1];
//	Ein(Edge_ids.size() - 1, 1) = Edge_ids[0];
//
//	Eigen::MatrixXd WV;
//	Eigen::MatrixXi WF,WE;
//	Eigen::VectorXi J;
//	std::string flags = "-c";
//	igl::triangle::cdt(Vin, Ein, flags, WV, WF, WE, J);
//	F = WF;
//	assert(WV.rows() == Vin.rows());
//}
