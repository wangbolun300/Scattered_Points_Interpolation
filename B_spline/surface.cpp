#include<surface.h>
#include<curve.h>
#include <iomanip> 
#include<mesh_processing.h>

namespace SIBSplines{
int Bsurface::nu() {
	return U.size() - 2 - degree1;
}
int Bsurface::nv() {
	return V.size() - 2 - degree2;
}
// TODO make this function more efficient by removing 0 valued elements
Vector3d BSplineSurfacePoint_(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V, const double upara,
	const double vpara, const std::vector<std::vector<Vector3d>>& control) {
	Eigen::Vector3d result = Eigen::Vector3d(0, 0, 0);
	int nu = U.size() - 2 - degree1;// n + 1 = number of control points
	assert(nu + 1 == control.size());
	int nv = V.size() - 2 - degree2;// n + 1 = number of control points
	assert(nv + 1 == control[0].size());

	for (int i = 0; i < nu + 1; i++) {
		double base1 = Nip(i, degree1, upara, U);
		for (int j = 0; j < nv + 1; j++) {
			double base2 = Nip(j, degree2, vpara, V);
			result += base1 * base2 * control[i][j];
		}
		//std::cout << "base " << base << std::endl;

	}
	return result;
}
Vector3d Bsurface::BSplineSurfacePoint(const Bsurface& surface, const double upara, const double vpara) {
	return BSplineSurfacePoint_(surface.degree1, surface.degree2, surface.U, surface.V, upara, vpara, surface.control_points);
}
// id_list may contain id, if contains, return the one right after id; if the one is the id is the last of id_list,
	// still return id. this happens when:
/*
t0: oxoo
t1: ooxx
t2:     xooo
t3:     oxoo
*/
// where the xs are the choosen ones. x in the first row is the choosen control point of t0, and then t1 can choose two 
// control points

// if not contain id, return the first element of id_list(t2 row);

int find_the_id_after_this_one(const int id, const std::vector<int>& id_list) {
	int which = -1;
	for (int i = 0; i < id_list.size(); i++) {
		if (id == id_list[i]) {
			which = i;
			break;
		}
	}
	if (which == -1) {
		return id_list[0];
	}
	else {
		if (which == id_list.size() - 1) {
			return id;
		}
		else {
			return id_list[which + 1];
		}
	}

	std::cout << "ERROR: SHOUD NOT GO HERE IN int find_the_id_after_this_one()!!!" << std::endl;
	return -1;
}
// get feasible control point matrix. if checking v direction control points (v_direction=true), make sure that the Uin knot vector
// is already fixed. for the explanation of Ugrid, Vgrid and UVmap, see function generate_UV_grid() in 'mesh_processing.h'
Eigen::MatrixXi get_feasible_control_point_matrix(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin, const bool v_direction,
	 const std::vector<double>& Ugrid, const std::vector<double>&Vgrid, const Eigen::MatrixXi& UVmap,
	const int nbr_para, std::vector<std::vector<std::array<int, 2>>>&para_to_feasible, const double per, 
	per_too_large& per_flag, std::vector<int> &feasible_order) {
	// if checking v_direction, the control points distribution depends on Uin
	para_to_feasible.resize(nbr_para);
	feasible_order.resize(nbr_para);
	assert(UVmap.rows() == Ugrid.size() && UVmap.cols() == Vgrid.size());
	std::vector<double> kv;// knot vector
	int degree;
	std::vector<double> grid, other_grid;
	if (v_direction) {
		kv = Uin;
		degree = degree1;
		grid = Vgrid;//  the v parameter of each iso-v line
		other_grid = Ugrid;
	}
	else {
		kv = Vin;
		degree = degree2;
		grid = Ugrid;
		other_grid = Vgrid;
	}
	int csize = kv.size() - degree - 1;// Nip is from 0 to csize-1. 
	int rows = grid.size();
	Eigen::MatrixXi result= Eigen::MatrixXi::Constant(rows, csize, -1); // initialize the matrix with -1
	std::vector<std::vector<std::vector<int>>> feasible;
	feasible.resize(rows);
	for (int i = 0; i < rows; i++) {
		feasible[i].resize(csize);
	}

	// get basic feasible points, there are overlaps meaning this control points affect more than one parameter
	for (int i = 0; i < rows; i++) {// for each iso-line
		for (int j = 0; j < other_grid.size(); j++) {
			int id = v_direction ? UVmap(j, i) : UVmap(i, j);
			if (id < 0) {
				continue;
			}
			double para = other_grid[j];
			// get the bacis feasible points of parameter whose index is id
			
			std::vector<int> basic_fp = feasible_control_point_of_given_parameter(para, kv, degree,per);
			
			for (int r = 0; r < basic_fp.size(); r++) {
				feasible[i][basic_fp[r]].push_back(id);
			}
		}
	}
	// next deal with duplication to generate a reduced_feasible matrix, for each ith row and jth column, there
	// are at most 1 element
	std::vector<std::vector<std::vector<int>>> reduced_feasible;
	reduced_feasible.resize(rows);
	for (int i = 0; i < rows; i++) {
		reduced_feasible[i].resize(csize);
	}
	for (int i = 0; i < rows; i++) {
		int now_checking = -1;
		for (int j = 0; j < csize; j++) {
			int ij_size = feasible[i][j].size();
			if (ij_size == 0) {
				now_checking = -1;
				continue;
			}
			if (ij_size == 1) {// if no overlape, directly select this one
				now_checking = feasible[i][j][0];// the current id
				reduced_feasible[i][j].push_back(now_checking);
				para_to_feasible[now_checking].push_back({ {i,j} });
			}
			if (ij_size > 1) {
				if (now_checking == -1) {// should take the first element 
					now_checking= feasible[i][j][0];// the current id
					reduced_feasible[i][j].push_back(now_checking);
					para_to_feasible[now_checking].push_back({ {i,j} });
				}
				else {
					now_checking = find_the_id_after_this_one(now_checking, feasible[i][j]);
					reduced_feasible[i][j].push_back(now_checking);
					para_to_feasible[now_checking].push_back({ {i,j} });
				}
			}
		}
	}
	// now there are at most one feasible point in each grid 
	
	for (int i = 0; i < result.rows(); i++) {
		for (int j = 0; j < result.cols(); j++) {
			int rf_size = reduced_feasible[i][j].size();
			if (rf_size == 0) { 
				continue; 
			}
			assert(rf_size == 1);
			result(i,j) = reduced_feasible[i][j][0];
		}
	}
	per_flag.flag = true;
	for (int i = 0; i < para_to_feasible.size(); i++) {
		if (para_to_feasible[i].empty()) {
			per_flag.flag = false;
			std::cout << "AN ERROR OCCURED: NO FEASIBLE POINT IS CHOOSEN FOR VERTEX #" << i << "\nPLEASE MAKE per OR per_ours SMALLER" << std::endl;
			return result;
		}
		assert(!para_to_feasible[i].empty());
	}

	int tt = 0;
	std::vector<bool> oflag(nbr_para, false);
	for (int i = 0; i < result.rows(); i++) {
		for (int j = 0; j < result.cols(); j++) {
			int tid = result(i, j);
			if (tid < 0) {
				continue;
			}
			if (oflag[tid] == false) {
				oflag[tid] = true;
				feasible_order[tt] = tid;
				tt++;
			}
		}
	}
	if (tt != nbr_para) {
		std::cout << "ERROR IN feasible_order" << std::endl;
		exit(0);
	}
	return result;
}
// bool updated shows if there is difference between the input and output
// maximal_selection is the maximal number of fcps we select in this step
Eigen::MatrixXi select_FCP_based_on_weight_naive(const Eigen::MatrixXi& fcp, std::vector<std::vector<std::array<int, 2>>> &para_to_feasible,
	const Eigen::MatrixXd &weight_matrix//, bool& updated, const int maximal_selection
) {
	//updated = false;

	// reorder the weight
	std::vector<std::vector<std::vector<std::array<int, 2>>>> selected(para_to_feasible.size());
	int selected_nbr = 0;
	for (int i = 0; i < para_to_feasible.size(); i++) {
		double weight = 0;
		for (int j = 0; j < para_to_feasible[i].size(); j++) {
			int id0 = para_to_feasible[i][j][0];
			int id1 = para_to_feasible[i][j][1];
			double w = weight_matrix(id0, id1);
			if (w > weight) {
				std::vector<std::array<int, 2>> newlist;
				newlist.push_back({ {id0,id1} });
				selected[i].push_back(newlist);
				weight = w;
			}
			else {
				if (w == weight) {
					selected[i].back().push_back({ {id0,id1} });
				}
			}
		}
	}

	// select highest weighted fcp
	bool total_copy = false;
	Eigen::MatrixXi result = Eigen::MatrixXi::Constant(fcp.rows(), fcp.cols(), -1);
	for (int i = 0; i < para_to_feasible.size(); i++) {
		assert(!selected[i].empty());

		//if (selected[i].size() > 1 && !total_copy) {

		//	selected_nbr += 1;
		//	if (selected_nbr > maximal_selection) {
		//		total_copy = true;
		//	}
		//	//updated = true;
		//}
		if (total_copy) {
			for (int j = 0; j < para_to_feasible[i].size(); j++) {
				int id0 = para_to_feasible[i][j][0];
				int id1 = para_to_feasible[i][j][1];
				result(id0, id1) = i;
			}
		}
		else {
			std::vector<std::array<int, 2>> highest = selected[i].back();
			para_to_feasible[i] = highest;
			for (int j = 0; j < highest.size(); j++) {
				int id0 = highest[j][0];
				int id1 = highest[j][1];
				result(id0, id1) = i;
				if (weight_matrix(id0, id1) == 1) {// if no one shares fcp, meaning the parameter is 0 or 1, choose one fcp is enough
					break;
				}
			}
		}

	}
	return result;

}
bool para_to_feasible_validation(std::vector<std::vector<std::array<int, 2>>> &para_to_feasible) {
	for (int i = 0; i < para_to_feasible.size(); i++) {
		if (para_to_feasible[i].size() == 0) {
			return false;
		}
	}
	return true;
}

bool fcp_weight_validation(const Eigen::MatrixXi& fcp,
	const Eigen::MatrixXd &weight_matrix) {
	int row = fcp.rows();
	int col = weight_matrix.cols();
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			int v1 = fcp(i, j);
			double v2 = weight_matrix(i, j);
			if (v1 >= 0) {
				if (v2 < 0) {
					std::cout << "<,<, fcp "<<v1 << std::endl;
					std::cout << "value " << v2 << std::endl;
					return false;
				}
				if (v2 == 0) {
					std::cout << "==" << std::endl;
					return false;
				}
			}
		}
	}
	return true;
}

bool paras_to_feasible_and_fcp_validation(std::vector<std::vector<std::array<int, 2>>> &ptf,
	const Eigen::MatrixXi& fcp) {

	for (int i = 0; i < ptf.size(); i++) {
		bool found = false;
		for (int j = 0; j < ptf[i].size(); j++) {
			int id0 = ptf[i][j][0];
			int id1 = ptf[i][j][1];
			if (fcp(id0, id1) >= 0) {
				found = true;
			}
		}
		if (!found) {
			std::cout << "the i= " << i << "th element not find corresponding in fcp" << std::endl;
			return false;
		}
	}
	return true;
}
// bool updated shows if there is difference between the input and output
// maximal_selection is the maximal number of fcps we select in this step
Eigen::MatrixXi select_FCP_based_on_weight(const Eigen::MatrixXi& fcp, 
	std::vector<std::vector<std::array<int, 2>>> &para_to_feasible,
	const Eigen::MatrixXd &weight_matrix, bool& updated, const int maximal_selection,
	const std::vector<int> &feasible_order) {
	assert(fcp_weight_validation(fcp, weight_matrix));
	updated = false;
	assert(para_to_feasible_validation(para_to_feasible));
	assert(paras_to_feasible_and_fcp_validation(para_to_feasible, fcp));
	// reorder the weight
	std::vector<std::vector<std::vector<std::array<int, 2>>>> selected(para_to_feasible.size());
	int selected_nbr = 0;
	for (int i = 0; i < para_to_feasible.size(); i++) {
		double weight = 0;
		for (int j = 0; j < para_to_feasible[i].size(); j++) {
			int id0 = para_to_feasible[i][j][0];
			int id1 = para_to_feasible[i][j][1];
			double w = weight_matrix(id0, id1);
			
			if (w > weight) {
				std::vector<std::array<int,2>> newlist;
				newlist.push_back({ {id0,id1} });
				selected[i].push_back(newlist);
				weight = w;
			}
			else {
				if (w == weight) {
					assert(selected[i].size() > 0);
					selected[i].back().push_back({ {id0,id1} });
				}
			}
		}
		if (selected[i].empty()) {
			std::cout << "problematic i " << i << std::endl;
			for (int j = 0; j < para_to_feasible[i].size(); j++) {
				int id0 = para_to_feasible[i][j][0];
				int id1 = para_to_feasible[i][j][1];
				std::cout << para_to_feasible[i][j][0] << " "<<para_to_feasible[i][j][1] << std::endl;
				double w = weight_matrix(id0, id1);
				std::cout << "w " << w << std::endl;
			}
		}
		assert(selected[i].size() > 0);
	}

	// select highest weighted fcp
	bool total_copy = false;
	Eigen::MatrixXi result = Eigen::MatrixXi::Constant(fcp.rows(), fcp.cols(), -1);
	for (int i = 0; i < para_to_feasible.size(); i++) {
		if (selected[i].empty()) {
			std::cout << "problematic i " << i << std::endl;
		}
		assert(!selected[i].empty());
		
		if (selected[i].size() > 1&&!total_copy) {
			
			selected_nbr += 1;
			if (selected_nbr > maximal_selection) {
				total_copy = true;
			}
			updated = true;
		}
		if (total_copy) {
			for (int j = 0; j < para_to_feasible[i].size(); j++) {
				int id0 = para_to_feasible[i][j][0];
				int id1 = para_to_feasible[i][j][1];
				result(id0, id1) = i;
			}
		}
		else {
			std::vector<std::array<int, 2>> highest = selected[i].back();
			std::vector<std::array<int, 2>> sout;
			for (int j = 0; j < highest.size(); j++) {
				int id0 = highest[j][0];
				int id1 = highest[j][1];
				result(id0, id1) = i;
				sout.push_back({ {id0,id1} });
				if (weight_matrix(id0, id1) == 1) {// if no one shares fcp, meaning the parameter is 0 or 1, choose one fcp is enough
					break;
				}
			}
			para_to_feasible[i] = sout;
		}
		
	} 
	if (updated == false) {// the weight matrix will not renew, now force to remove some redundant fcps, 
		for (int i = 0; i < feasible_order.size(); i++) {
			int which = feasible_order[i];
			if (para_to_feasible[which].size() > 1) {// this is a redundant fcp
				selected_nbr += 1;
				
				if (selected_nbr > maximal_selection) {
					break;
				}
				for (int j = 0; j < para_to_feasible[which].size() - 1; j++) {
					int id0 = para_to_feasible[which][j][0];
					int id1 = para_to_feasible[which][j][1];
					result(id0, id1) = -1;
					//std::cout << "\nForcing deleting " << id0 << " " << id1 << std::endl;
				}
				std::array<int, 2> last_element = para_to_feasible[which].back();
				para_to_feasible[which].resize(1);// delete other redundant fcps
				para_to_feasible[which][0] = last_element;
				updated = true;
			}
		}
	}
	return result;

}

// pick one fcp for each parameter. this function is used only when the selecting based on weight is finished but
// still have redundant fcps
Eigen::MatrixXi remove_redundant_FCP(const Eigen::MatrixXi& fcp, std::vector<std::vector<std::array<int, 2>>> &para_to_feasible) {
	int nbr_re = 0;
	Eigen::MatrixXi result = Eigen::MatrixXi::Constant(fcp.rows(), fcp.cols(), -1);
	for (int i = 0; i < para_to_feasible.size(); i++) {
		if (para_to_feasible[i].size() > 1) {
			nbr_re++;
		}
		int id0 = para_to_feasible[i].front()[0];
		int id1 = para_to_feasible[i].front()[1];
		result(id0, id1) = i;
	}
	//std::cout << "** redundant fcp nbr " << nbr_re << std::endl;
	return result;

}
Eigen::MatrixXd weight_matrix_calculation(const Eigen::MatrixXi& fcp,const Eigen::MatrixXi& interval_matrix,
	const int degree) {
	// we define weight here. when more fcp share control points, the weight is smaller;
	// when there are more points in this row, the weight is smaller

	int rows = fcp.rows();
	int cols = fcp.cols();

	const auto weight_strategy = [](const std::vector<int>& vec, const int nbr_points) {
		int sum = 1;
		// weight is the 
		for (int i = 1; i < vec.size(); i++) {
			sum += vec[i] * vec[i] * i;
		}
		double result = 1 / double(sum);// *nbr_points);
		assert(result > 0);
		return result;
	};
	// initialize weight matrix
	Eigen::MatrixXd weight_matrix = Eigen::MatrixXd::Constant(rows, cols, -1);
	
	

	for (int j = 0; j < cols; j++) {// for each column
		for (int i = 0; i < rows; i++) {
			
			if (fcp(i, j) < 0) {
				continue;
			}
			int this_id = fcp(i, j);
			int interval = interval_matrix(i, j);
			std::vector<int> interval_info(degree + 2, 0);// h = interval_info[r] means there are h fcps sharing r control points 
			if (interval < 0) {// interval value is -2: no other fcp share control points with it 
				weight_matrix(i, j) = 1;
			}
			else {// search this column if there is overlap
				int col_nbr_pts = 0;// how many points there are in this column (except for whose parameter is 0 or 1)
				int col_before_nbr = 1;// how many points related but before this point
				int row_nbr_pts = 1;
				for (int k = 0; k < j; k++) {
					if (this_id == fcp(k, j)) {
						row_nbr_pts += 1;
					}
				}
				
				for (int k = 0; k < rows; k++) {
					if (fcp(k, j) < 0) {// if this point does not exist
						continue;
					}
					int other_interval = interval_matrix(k, j);
					if (other_interval < 0) {// if -1 or -2, no need to count 
						continue;
					}
					col_nbr_pts += 1;
					int diff = abs(other_interval - interval);
					if (diff > degree) {
						continue;
					}
					int nbr_share = degree + 1 - diff;// meaning that fcp(i,j) share nbr_share points with fcp(k,j).
					if (k < i) {
						col_before_nbr += 1;
					}
					assert(nbr_share > 0);
					interval_info[nbr_share] += 1;
				}
				double weight = 
#ifndef WEIGHT_NAIVE			
					double(row_nbr_pts)
					*
					1/double(col_before_nbr)
					*
#endif
					weight_strategy(interval_info, col_nbr_pts);
				if (weight >= 1) {
					std::cout << "large weight " << weight << std::endl;
				}
				assert(weight > 0);
				
				weight_matrix(i, j) = weight;
			}
		}
	}
	assert(fcp_weight_validation(fcp, weight_matrix));
	return weight_matrix;
}
bool para_to_feasible_is_clean(std::vector<std::vector<std::array<int, 2>>> &para_to_feasible) {
	for (int i = 0; i < para_to_feasible.size(); i++) {
		if (para_to_feasible[i].size() > 1) {
			return false;
		}
	}
	return true;
}
// calculate weights and select ACP according to the weight
Eigen::MatrixXi calculate_active_control_points_from_feasible_control_points(const Eigen::MatrixXi& fcp,const bool v_direction,
	const std::vector<double> &Uknot, const std::vector<double> &Vknot, 
	const Eigen::MatrixXd& paras, const int degree1, const int degree2, 
	std::vector<std::vector<std::array<int, 2>>> &para_to_feasible, const int target_steps, 
	const std::vector<int> &feasible_order) {
	assert(paras_to_feasible_and_fcp_validation(para_to_feasible, fcp));
	assert(para_to_feasible.size() == paras.rows());
	int maximal_processing = 
		std::max(1, int(para_to_feasible.size() / target_steps));

	// if return -2, the parameter is 0 or 1, then no one share control points with it
	const auto para_in_which_interval = [](const double para, const std::vector<double> &knots) {
		int result = -2;
		for (int i = 0; i < knots.size() - 1; i++) {
			if (para == knots.front()) {
				break;
			}
			if (para >= knots[i] && para < knots[i + 1]) {
				result = i;
				break;
			}
		}
		return result;
	};
	

	int rows = fcp.rows();
	int cols = fcp.cols();

	

	// initialize weight matrix
	Eigen::MatrixXd weight_matrix = Eigen::MatrixXd::Constant(rows, cols, -1);

	
	

	std::vector<double> kv;
	int degree;
	if (v_direction) {// if control points belongs to iso-v lines, then we are considering v knot vector
		kv = Vknot;
		degree = degree2;
	}
	else {
		kv = Uknot;
		degree = degree1;
	}
	

	// now in interval_matrix, value = -1 means no fcp here; value = -2 means the parameter is 0 or 1; otherwise, 
	// the value implies the parameter is in [u_value, u_(value+1)]
	// interval_matrix no need updating, fcp need to be updated

	bool updated = true;
	Eigen::MatrixXi selected_fcp = fcp;
	int nbr_rounds = 0;
	while (updated) {
		//std::cout << "before round " << nbr_rounds << std::endl;
		// initialize para_matrix. it will record the u or v parameters of each fcp
		Eigen::MatrixXd para_matrix = Eigen::MatrixXd::Constant(rows, cols, -1);
		int uv = v_direction ? 1 : 0;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (selected_fcp(i, j) < 0) {
					continue;
				}
				para_matrix(i, j) = paras(selected_fcp(i, j), uv);
			}
		}
		// initialize interval matrix. it will record which interval [u_i, u_(i+1)] does the parameters in
		Eigen::MatrixXi interval_matrix = Eigen::MatrixXi::Constant(rows, cols, -1);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (selected_fcp(i, j) < 0) {
					continue;
				}
				double para = para_matrix(i, j);
				int which = para_in_which_interval(para, kv);
				interval_matrix(i, j) = which;
			}
		}
		// calculate weight
		weight_matrix = weight_matrix_calculation(selected_fcp, interval_matrix, degree);
		// calculate weight finished
		assert(fcp_weight_validation(selected_fcp, weight_matrix));
		selected_fcp = 
			select_FCP_based_on_weight(selected_fcp, para_to_feasible, weight_matrix,updated, maximal_processing,
				feasible_order);
			
		nbr_rounds += 1;
		//std::cout << "after round " << nbr_rounds << std::endl;
		
	}
	//std::cout << "select fcp nbr of rounds " << nbr_rounds << std::endl;
	//exit(0);
	//selected_fcp = remove_redundant_FCP(selected_fcp, para_to_feasible);
	if (!para_to_feasible_is_clean(para_to_feasible)) {
		std::cout << "FCP still have rendundant elements" << std::endl;
		exit(0);
	}
	return selected_fcp;
}
std::vector<double> get_iso_line_parameters_from_ACP(const Eigen::MatrixXi&ACP, const int id, const Eigen::MatrixXd& paras, const bool v_direction) {
	int uv = v_direction ? 1 : 0;// if checking iso-v lines, then we are dealing with V parameters
	std::vector<double> result;
	Eigen::VectorXi indices = ACP.col(id);
	for (int i = 0; i < indices.size(); i++) {
		if (indices[i] < 0) {
			continue;
		}
		result.push_back(paras(indices[i], uv));
	}
	return result;
}
bool check_ACP_calidation(const Eigen::MatrixXi& ACP, const int nbr_para) {
	std::vector<bool> check(nbr_para, false);
	for (int i = 0; i < ACP.rows(); i++) {
		for (int j = 0; j < ACP.cols(); j++) {
			int value = ACP(i, j);
			if (value >= 0) {
				check[value] = true;
			}
		}
	}
	for (int i = 0; i < check.size(); i++) {
		if (check[i] == false) {
			std::cout << "ACP ERROR: MISSING #" << i << ", ACP:\n" << ACP << std::endl;
			return false;
		}
	}
	return true;
}

// return true, then the surface knot fixing for both u and v is finished
bool progressively_generate_interpolation_knot_vectors(const bool v_direction, int degree1, int degree2,
	std::vector<double>& Uknot, std::vector<double>& Vknot, const std::vector<double> Ugrid, const std::vector<double>& Vgrid,
	const Eigen::MatrixXi& grid_map, const Eigen::MatrixXd& param, const double per_ours,const double per, 
	const int target_steps, const bool enable_max_fix_nbr, per_too_large& per_flag) {
	const int nbr_para = param.rows();
	std::vector<std::vector<std::array<int, 2>>> para_to_feasible;
	std::vector<int> feasible_order;// the order of the parameters in FCP
	Eigen::MatrixXi FCP = get_feasible_control_point_matrix(3, 3, Uknot, Vknot, v_direction, Ugrid, Vgrid, grid_map, nbr_para,
		para_to_feasible,per_ours,per_flag, feasible_order);
	if (per_flag.flag == false) {
		// this means the per_ours should be smaller
		return false;
	}
	Eigen::MatrixXi ACP =
#ifdef NO_SELECTING_ACP
		FCP;
#else
#ifdef NAIVE_SELECTING_ACP
		remove_redundant_FCP(FCP, para_to_feasible);
#else
		calculate_active_control_points_from_feasible_control_points(FCP, v_direction, Uknot, Vknot, param,
		degree1, degree2, para_to_feasible, target_steps, feasible_order);
#endif
#endif
	assert(check_ACP_calidation(ACP, nbr_para));
	std::vector<double> kv, kv_other;
	if (v_direction) {// if checking iso-v lines, then we are fixing V knot vector
		kv = Vknot;
		kv_other = Uknot;
	}
	else {
		kv = Uknot;
		kv_other = Vknot;
	}
	bool finished = false;
	int nbr_fixed_in_each_itr;
	if (enable_max_fix_nbr) {
		int estimate_nbr_cp = param.rows()*param.rows() / 2;
		nbr_fixed_in_each_itr = std::max(1,std::abs(int(kv.size()-kv_other.size())));// std::max((int)sqrt(estimate_nbr_cp) / target_steps, 2);
		//std::cout << "max fix nbr, " << nbr_fixed_in_each_itr << std::endl;
	}
	else {
		nbr_fixed_in_each_itr = -1;
	}

	for (int i = 0; i < ACP.cols(); i++) {
		bool fully_fixed = false;
		std::vector<double> parameters = get_iso_line_parameters_from_ACP(ACP, i, param, v_direction);
		std::vector<double> kv_new = fix_knot_vector_to_interpolate_curve_WKW(3, kv, parameters,per, 
			fully_fixed, nbr_fixed_in_each_itr);
		int real_inserted = kv_new.size() - kv.size();
		kv = kv_new;
		if (i == ACP.cols() - 1 && fully_fixed) {
			finished = true;
			break;
		}
		
		if (enable_max_fix_nbr){
			if (kv.size() > kv_other.size()) {// the size difference between Uknot and Vknot decides which knot vector to refine
				break;
			}
			if (nbr_fixed_in_each_itr - real_inserted >= 0) {
				nbr_fixed_in_each_itr = nbr_fixed_in_each_itr - real_inserted;
			}
			else {
				break;
			}
		}
	}
	if (v_direction) {// if checking iso-v lines, then we are fixing V knot vector
		Vknot = kv;
		
	}
	else {
		Uknot = kv;
	}
	//std::cout << "** sizes " << Uknot.size() << " " << Vknot.size() << std::endl;
	return finished;
}

// if v_direction is true, then checking iso-v lines to fix U knot vector
std::vector<double> get_iso_line_parameters(const int degree1, const int degree2, const bool v_direction, const int line_id,
	const std::vector<double>& Ugrid, const std::vector<double>& Vgrid, const Eigen::MatrixXi& grid_map) {
	std::vector<double> result;
	std::vector<double> grid;
	Eigen::VectorXi line_map;
	
	if (v_direction) {
		grid = Ugrid;
		line_map = grid_map.col(line_id);
	}
	else {
		grid = Vgrid;
		line_map = grid_map.row(line_id);
		
	}

	assert(grid.size() == line_map.size());
	for (int i = 0; i < grid.size(); i++) {
		if (line_map[i] >= 0) {
			result.push_back(grid[i]);
			
		}
		
	}
	return result;

}

std::vector<double> update_knot_vector_based_on_grid(const int degree1, const int degree2, 
	const bool v_direction, 
	const std::vector<double>& Ugrid, const std::vector<double>& Vgrid, const Eigen::MatrixXi& grid_map,
	const double per, const std::vector<double>& Uin) {
	int size;
	if (v_direction) {
		size = Vgrid.size();
	}
	else {
		size = Ugrid.size();
	}
	std::vector<double> Utemp = Uin;
	for (int i = 0; i < size; i++) {
		std::vector<double> paras = get_iso_line_parameters(degree1, degree2, v_direction, i, Ugrid, Vgrid, grid_map);
		bool fully_fixed=false;
		Utemp = fix_knot_vector_to_interpolate_curve_WKW(degree1, Utemp, paras, per, fully_fixed);
		assert(fully_fixed == true);
	}
	return Utemp;
}
void Bsurface::generate_interpolation_knot_vectors( int degree1, int degree2,
	std::vector<double>& Uknot, std::vector<double>& Vknot,
	const Eigen::MatrixXd& param_original, 
	const double per_ours,const double per, const int target_steps, const bool enable_max_fix_nbr, per_too_large &per_flag) {
	
	
	std::vector<double> Ugrid, Vgrid;
	Eigen::MatrixXi grid_map;
	generate_UV_grid(param_original, Ugrid, Vgrid, grid_map);
	std::cout << "** UV grid sizes, " << Ugrid.size() << ", " << Vgrid.size() << std::endl;
	bool fully_fixed;
	std::vector<double> Utemp = Uknot, Vtemp = Vknot;
	// initialize U
	Utemp=update_knot_vector_based_on_grid(degree1, degree2, true, Ugrid, Vgrid, grid_map, per, Utemp);

	std::cout << "finished initialize Uknot" << std::endl;
	print_vector(Utemp);

	// initialize V
	Vtemp = update_knot_vector_based_on_grid(degree1, degree2, false, Ugrid, Vgrid, grid_map, per, Vtemp);
	std::cout << "finished initialize Vknot" << std::endl;
	print_vector(Vtemp);
	bool finished = false;
	bool v_direction;// = start_from_v_direction;
	if (Utemp.size() > Vtemp.size()) {// iso-v line, update V
		v_direction = true;
	}
	else {
		v_direction = false;
	}
	
	while (!finished) {
		finished = progressively_generate_interpolation_knot_vectors(v_direction, degree1, degree2,
			Utemp, Vtemp, Ugrid, Vgrid, grid_map, param_original, per_ours, per, target_steps, enable_max_fix_nbr,
			per_flag);
		if (per_flag.flag == false) {
			// per_ours need to be reduced
			return;
		}
		v_direction = !v_direction;
		/*if (!finished) {
			std::cout << "switched, UV size "<<Utemp.size()<<" "<<Vtemp.size() << std::endl;

		}*/
	}
	print_vector(Utemp);
	print_vector(Vtemp);

	//std::cout << "Uknot Vknot size "<< Utemp.size()<<" "<<Vtemp.size() << std::endl;

	
	
	std::cout << "knot fixing finished, sizes " << Utemp.size() << " " << Vtemp.size() << std::endl;
	// post processing
	bool post_processing = false;
	if (post_processing) {

		// V direction, fix V
		Utemp = Uknot;
		finished = progressively_generate_interpolation_knot_vectors(true, degree1, degree2,
			Utemp, Vtemp, Ugrid, Vgrid, grid_map, param_original, per_ours, per, target_steps, false,// not enable_max_fix_nbr
			per_flag);
		if (per_flag.flag == false) {
			// per_ours need to be reduced
			return;
		}
		assert(finished);

		// U direction, fix U
		Vtemp = Vknot;
		finished = progressively_generate_interpolation_knot_vectors(false, degree1, degree2,
			Utemp, Vtemp, Ugrid, Vgrid, grid_map, param_original, per_ours, per, target_steps, false,// not enable_max_fix_nbr
			per_flag);
		if (per_flag.flag == false) {
			// per_ours need to be reduced
			return;
		}
		assert(finished);

		std::cout << "post processing results, " << Utemp.size() << " " << Vtemp.size() << std::endl;
		print_vector(Utemp);
		print_vector(Vtemp);
	}

	
	Uknot = Utemp;
	Vknot = Vtemp;
	return;
}


void Bsurface::generate_interpolation_knot_vectors(int degree1, int degree2,
	std::vector<double>& Uknot, std::vector<double>& Vknot,
	const Eigen::MatrixXd& param_original, 
	double &per_ours, const double per, const int target_steps, const bool enable_max_fix_nbr) {
	per_too_large per_flag;
	per_flag.flag = false;
	double per_ours_tmp = per_ours;
	while (1) {
		generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param_original, per_ours_tmp,
			per, target_steps, enable_max_fix_nbr, per_flag);
		if (per_flag.flag) {
			break;
		}
		else {
			per_ours_tmp *= 0.9;
			std::cout<<"delta is reduced to "<<per_ours_tmp<<std::endl;
		}
	}
	per_ours = per_ours_tmp;
	return;	
}
double Bsurface::max_interpolation_err(const Eigen::MatrixXd&ver, const Eigen::MatrixXd& param, Bsurface& surface) {
	double err = 0;
	for (int i = 0; i < ver.rows(); i++) {
		Vector3d v = ver.row(i);
		double parau = param(i, 0);
		double parav = param(i, 1);
		Vector3d vs = BSplineSurfacePoint(surface, parau, parav);
		double newerr = (v - vs).norm();
		if (newerr > err) {
			err = newerr;
		}
	}
	return err;
}
// calculate interpolation error for approximation/interpolation method
Eigen::MatrixXd Bsurface::interpolation_err_for_apprximation(const Eigen::MatrixXd&ver, 
	const Eigen::MatrixXd& param, Bsurface& surface,double &max_err) {
	Eigen::MatrixXd result(ver.rows(), ver.cols());
	double err = 0;
	for (int i = 0; i < ver.rows(); i++) {
		Vector3d v = ver.row(i);
		double parau = param(i, 0);
		double parav = param(i, 1);
		Vector3d vs = BSplineSurfacePoint(surface, parau, parav);
		Vector3d diff = v - vs;
		result.row(i) << diff[0], diff[1], diff[2];
		double newerr = diff.norm();
		if (newerr > err) {
			err = newerr;
		}
	}
	max_err = err;
	return result;
}

void B_spline_surface_to_mesh(Bsurface &surface, const int pnbr, Eigen::MatrixXd &ver, Eigen::MatrixXi& faces) {
	std::vector<std::vector<Vector3d>> pts;
	ver.resize(pnbr*pnbr, 3);
	int verline = 0;
	for (int i = 0; i < pnbr; i++) {
		for (int j = 0; j < pnbr; j++) {
			double upara = double(i) / (pnbr - 1);
			double vpara = double(j) / (pnbr - 1);
			ver.row(verline) = surface.BSplineSurfacePoint(surface, upara, vpara);
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




void Bsurface::surface_visulization(Bsurface& surface, const int nbr, Eigen::MatrixXd & v, Eigen::MatrixXi &f) {
	B_spline_surface_to_mesh(surface, nbr, v, f);
	return;
}
// void surface_visulization(std::vector<Bsurface>& surfaces, const int pnbr, Eigen::MatrixXd & ver, 
// 	Eigen::MatrixXi &faces, const int without) {
// 	std::vector<std::vector<Vector3d>> pts;
// 	ver.resize(pnbr*pnbr, 3);
// 	int verline = 0;
// 	for (int i = 0; i < pnbr; i++) {
// 		for (int j = 0; j < pnbr; j++) {
// 			double upara = double(i) / (pnbr - 1);
// 			double vpara = double(j) / (pnbr - 1);
// 			ver.row(verline) = Vector3d(0, 0, 0);
// 			for (int si = 0; si < surfaces.size()- without; si++) {
// 				ver.row(verline) += BSplineSurfacePoint(surfaces[si], upara, vpara);
// 			}
			
// 			verline++;
// 		}
// 	}
// 	faces.resize(2 * (pnbr - 1)*(pnbr - 1), 3);
// 	int fline = 0;
// 	for (int i = 0; i < pnbr - 1; i++) {
// 		for (int j = 0; j < pnbr - 1; j++) {
// 			faces.row(fline) = Vector3i(i + pnbr * (j + 1), i + pnbr * j, i + pnbr * (1 + j) + 1);
// 			faces.row(fline + 1) = Vector3i(i + pnbr * (1 + j) + 1, i + pnbr * j, i + pnbr * j + 1);
// 			fline += 2;
// 		}
// 	}
// 	return;
// }

}