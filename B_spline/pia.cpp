// this is an implementation of "B-spline surface fitting by iterative geometric interpolation/approximation
// algorithms, 2012, CAD"

#include<basis.h>
#include<surface.h>

double para_distance(const double uvalue, const double vvalue, const Eigen::MatrixXd &param, const int row) {
	Eigen::MatrixXd cornor(1, 2);
	cornor << uvalue, vvalue;
	return (cornor.row(0) - param.row(row)).norm();
}
void set_up_initial_surface(Bsurface &surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver) {
	int nbr = ver.rows();
	double dis00 = 2, dis01 = 2, dis10 = 2, dis11 = 2;
	int p00, p01, p10, p11;
	for (int i = 0; i < nbr; i++) {
		double dis = para_distance(0, 0, param, i);
		if (dis < dis00) {
			dis00=dis;
			p00 = i;
		}

		dis = para_distance(0, 1, param, i);
		if (dis < dis01) {
			dis01 = dis;
			p01 = i;
		}
		
		dis = para_distance(1, 0, param, i);
		if (dis < dis10) {
			dis10 = dis;
			p10 = i;
		}

		dis = para_distance(1, 1, param, i);
		if (dis < dis11) {
			dis11 = dis;
			p11 = i;
		}
	}
	int c1 = surface.nu() + 1;
	int c2 = surface.nv() + 1;
	surface.control_points.resize(surface.nu() + 1);
	for (int i = 0; i < surface.nu() + 1; i++) {
		surface.control_points[i].resize(surface.nv() + 1);
	}
	surface.control_points[0][0] = ver.row(p00);
	surface.control_points[0][surface.nv()] = ver.row(p01);
	surface.control_points[surface.nu()][surface.nv()] = ver.row(p11);
	surface.control_points[surface.nu()][0] = ver.row(p10);

	for (int i = 0; i < surface.nu()+1; i++) {
		for (int j = 0; j < surface.nv()+1; j++) {
			double u = double(i) / surface.nu();
			double v = double(j) / surface.nv();
			Vector3d p0 = surface.control_points[0][0] +
				(surface.control_points[surface.nu()][0] - surface.control_points[0][0])*u;
			Vector3d p1= surface.control_points[0][surface.nv()] +
				(surface.control_points[surface.nu()][surface.nv()] - surface.control_points[0][surface.nv()])*u;
			surface.control_points[i][j] = p0 + (p1 - p0)*v;
		}
	}
}
std::vector<Vector3d> get_the_delta_k(Bsurface& surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver) {
	std::vector<Vector3d> result;
	int size = param.rows();
	result.resize(size);
	for (int i = 0; i < size; i++) {
		Vector3d vs = BSplineSurfacePoint(surface, param(i,0), param(i, 1));
		Vector3d data = ver.row(i);
		result[i] = data - vs;
	}
	return result;
}

// returned ids means the parameter effects P_{id0}, ...
std::vector<int> get_parameter_corresponded_control_points(const std::vector<double>& U, const int degree,
	const double parameter) {
	int n = U.size() - degree - 2; // n+1 is the nbr of control points.
	std::vector<int> ids; 
	if (parameter == U.back()) { // if parameter==1, then 
		ids.push_back(n);
		return ids;
	}
	int interval;
	for (int i = 0; i < U.size() - 1; i++) {
		if (parameter >= U[i] && parameter < U[i + 1]) {
			interval = i;
			break;
		}
	}
	for (int i = interval - degree; i < interval + 1; i++) {
		ids.push_back(i);
	}
	return ids;
}


void get_the_data_ids_and_weights_for_each_control_point(
	std::vector<std::vector<std::vector<int>>> &ids, std::vector<std::vector<std::vector<double>>> &weights,
	Bsurface &surface, const Eigen::MatrixXd &param
) {
	int cnbr_0 = surface.nu() + 1;// nbr of control points in u direction
	int cnbr_1 = surface.nv() + 1;
	int pnbr = param.rows();// nbr of data points
	ids.resize(cnbr_0); 
	weights.resize(cnbr_0);
	for (int i = 0; i < cnbr_0; i++) {
		ids[i].resize(cnbr_1);
		weights[i].resize(cnbr_1);
	}
	for (int i = 0; i < pnbr; i++) {
		double u = param(i, 0);
		double v = param(i, 1);
		std::vector<int> u_ids = get_parameter_corresponded_control_points(surface.U, surface.degree1, u);
		std::vector<int> v_ids = get_parameter_corresponded_control_points(surface.V, surface.degree2, v);
		for (int r1 = 0; r1 < u_ids.size(); r1++) {
			for (int r2 = 0; r2 < v_ids.size(); r2++){
				ids[u_ids[r1]][v_ids[r2]].push_back(i);
				double N1 = Nip(u_ids[r1], surface.degree1, u, surface.U);
				double N2 = Nip(v_ids[r2], surface.degree2, v, surface.V);
				weights[u_ids[r1]][v_ids[r2]].push_back(N1*N2);

			} 
		}
	}
}

void get_Delta(const std::vector<Vector3d> &delta, const std::vector<std::vector<std::vector<int>>> &ids,
	const std::vector<std::vector<std::vector<double>>> &weights, std::vector<std::vector<Vector3d>>& Delta) {
	int pnbr = delta.size();
	int cnbr_u = ids.size();
	int cnbr_v = ids[0].size();
	Delta.resize(cnbr_u);
	for (int i = 0; i < cnbr_u; i++) {
		Delta[i].resize(cnbr_v);
	}
	for (int i = 0; i < cnbr_u; i++) {
		for (int j = 0; j < cnbr_v; j++) {
			Delta[i][j] = Vector3d(0, 0, 0);
			int nbr_vecs = ids[i][j].size();
			double weight_sum = 0;
			for (int k = 0; k < nbr_vecs; k++) {
				int delta_id = ids[i][j][k];
				Vector3d vec = delta[delta_id];
				double dweight = weights[i][j][k];
				weight_sum+=dweight;
				Delta[i][j] += vec * dweight;
			}
			if (weight_sum == 0) {
				//std::cout << "WARNING: IPA method need fairing" << std::endl;
			}
			else {
				Delta[i][j] /= weight_sum;
			}
				
		}
	}
}

void control_points_add_Delta(std::vector<std::vector<Vector3d>>& cps, const std::vector<std::vector<Vector3d>>&Delta) {
	for (int i = 0; i < cps.size(); i++) {
		for (int j = 0; j < cps[i].size(); j++) {
			cps[i][j] += Delta[i][j];
		}
	}
}
double get_maximal_step_length(const std::vector<std::vector<Vector3d>> &Delta) {
	double result = 0;
	for (int i = 0; i < Delta.size(); i++) {
		for (int j = 0; j < Delta[i].size(); j++) {
			if (Delta[i][j].norm() > result) {
				result = Delta[i][j].norm();
			}
		}
	}
	return result;
}
bool Delta_step_length_larger_than_threadshold(const std::vector<std::vector<Vector3d>> &Delta, const double th) {
	double l = get_maximal_step_length(Delta);
	std::cout << "step length " << l << std::endl;
	if (l > th) {
		return true;
	}
	return false;
}

// gamma is an input parameter for the weight.
void get_matrix_laplacian_fairing_of_control_points(double gamma, const std::vector<std::vector<Vector3d>>&cp, 
	Eigen::MatrixXd &matrix) {
	
	auto get_mid_and_corners=[](std::array<std::array<int, 2>, 4> &cid1, std::array<std::array<int, 2>, 4> &cid2,
		int uid, int vid) {
		cid1[0] = { uid - 1,vid };
		cid1[1] = { uid + 1,vid };
		cid1[2] = { uid,vid - 1 };
		cid1[3] = { uid,vid + 1 };

		cid2[0] = { uid - 1,vid - 1 };
		cid2[1] = { uid + 1,vid - 1 };
		cid2[2] = { uid + 1,vid + 1 };
		cid2[3] = { uid - 1,vid + 1 };

	};
	int unbr = cp.size();
	int vnbr = cp[0].size();
	auto to_matrix_id = [](const int uid, const int vid, const int unbr, const int vnbr, int& id) {
		id = uid * vnbr + vid;
	};
	std::array<std::array<int, 2>, 4> cid1;// four mid-control points
	std::array<std::array<int, 2>, 4> cid2;// four corner-control points
	
	matrix=Eigen::MatrixXd::Zero(unbr*vnbr, unbr*vnbr);
	for (int i = 0; i < unbr; i++) {
		for (int j = 0; j < vnbr; j++) {
			get_mid_and_corners(cid1, cid2, i, j);// get the ids of mid points and corner points
			int id_row;
			to_matrix_id(i, j, unbr, vnbr, id_row);// we are dealing with the ith row of the matrix
			double weight = 1;
			std::vector<int> mid_valid_col;
			std::vector<int> cor_valid_col;

			for (int k = 0; k < 4; k++) {// get the denominator of the weithts, and the ids of the midpoints and corpoints.
				if (cid1[k][0] >= 0 && cid1[k][0] < unbr&&cid1[k][1] >= 0 && cid1[k][1] < vnbr) {
					weight += exp(-gamma);
					int id_col;
					to_matrix_id(cid1[k][0], cid1[k][1], unbr, vnbr, id_col);
					mid_valid_col.push_back(id_col);
				}
				if (cid2[k][0] >= 0 && cid2[k][0] < unbr&&cid2[k][1] >= 0 && cid2[k][1] < vnbr) {
					weight += exp(-2*gamma);
					int id_col;
					to_matrix_id(cid2[k][0], cid2[k][1], unbr, vnbr, id_col);
					cor_valid_col.push_back(id_col);
				}
			}
			
			for (int k = 0; k < mid_valid_col.size(); k++) {
				matrix(id_row, mid_valid_col[k]) = exp(-gamma) / weight;
			}
			for (int k = 0; k < cor_valid_col.size(); k++) {
				matrix(id_row, cor_valid_col[k]) = exp(-2*gamma) / weight;
			}
			matrix(id_row, id_row) = 1 / weight;

		}
	}
}
void laplacian_fairing_of_control_points(const Eigen::MatrixXd& matrix, std::vector<std::vector<Vector3d>>&cp) {
	int unbr = cp.size();
	int vnbr = cp[0].size();
	int current = 0;
	Eigen::MatrixXd right(unbr*vnbr, 3), result(unbr*vnbr, 3);

	
	current = 0;
	for (int i = 0; i < unbr; i++) {
		for (int j = 0; j < vnbr; j++) {
			right.row(current) = cp[i][j];
			current++;
		}
	}
	// get right part.

	current = 0;
	result = matrix*(right);
	for (int i = 0; i < unbr; i++) {
		for (int j = 0; j < vnbr; j++) {
			cp[i][j] = result.row(current);
			current++;
		}
	}
}

// threads hold is the maximal step length
void progressive_iterative_approximation(Bsurface &surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver,
	const int max_itr, const double threadshold) {
	int nu = surface.nu();
	int nv = surface.nv();
	set_up_initial_surface(surface, param, ver);// get the surface of iteration 0
	std::cout << "initial surface set up " << std::endl;
	std::vector<std::vector<std::vector<int>>> ids;// map from control points to data points ids.
	std::vector<std::vector<std::vector<double>>> weights;
	get_the_data_ids_and_weights_for_each_control_point(ids, weights, surface, param);
	std::cout << "weights calculated" << std::endl;
	double gamma = 2;
	Eigen::MatrixXd matrix;
	get_matrix_laplacian_fairing_of_control_points(gamma, surface.control_points, matrix);// prepare fairing
	//std::cout << "matrix\n" << matrix << std::endl; exit(0);
	for (int i = 0; i < max_itr; i++) {
		std::cout << "iteration " << i << std::endl;
		std::vector<Vector3d> delta_k = get_the_delta_k(surface, param, ver); // get the error vector for each data point
		std::vector<std::vector<Vector3d>> Delta;
		get_Delta(delta_k, ids, weights, Delta);
		if (!Delta_step_length_larger_than_threadshold(Delta, threadshold)) {
			std::cout << "step length less than given threadshold" << std::endl;
			return; // if the step length is too short, do not update the surface.
		}
		control_points_add_Delta(surface.control_points, Delta);// update surface

		
		if (i % 19 == 0) {
			laplacian_fairing_of_control_points(matrix, surface.control_points);
		}
	}


	
}
