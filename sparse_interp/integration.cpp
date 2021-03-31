#include<integration.h>
#include<iostream>
#include<fstream>
#include <igl/opengl/glfw/Viewer.h>
#include<test.h>
#include<mesh_processing.h>

igl::opengl::glfw::Viewer global_viewer_1;
void test_intergration() {
	std::cout << "integrated" << std::endl;
}
using namespace std;
std::vector<std::vector<double>> read_CSV_data(const string inputFileName) {

	std::vector<std::vector<double>> result;
	ifstream infile;
	infile.open(inputFileName);
	if (!infile.is_open())
	{
		cout << "Path Wrong!!!!" << endl;
		return result;
	}

	int l = 0;
	while (infile) // there is input overload classfile
	{
		l++;
		string s;
		if (!getline(infile, s)) break;
		if (s[0] != '#') {
			istringstream ss(s);
			vector<double> record;
			record.clear();
			while (ss) {
				string line;
				if (!getline(ss, line, ','))
					break;
				try {
					record.push_back(stod(line));
					
				}
				catch (const std::invalid_argument e) {
					cout << "NaN found in file " << inputFileName << " line " << l
						<< endl;
					e.what();
				}
			}

			result.push_back(record);
		}
	}
	if (!infile.eof()) {
		cerr << "Could not read file " << inputFileName << "\n";
	}
	return result;
}


void write_edges_obj(const std::string filename, Eigen::MatrixXd& points, const Eigen::MatrixXi& edges) {
	std::ofstream fout;
	fout.open(filename);
	for (int i = 0; i < points.rows(); i++) {
		fout << "v " << points(i, 0) << " " << points(i, 1) << " " << points(i, 2) << std::endl;
	}
	for (int i = 0; i < edges.rows(); i++) {
		fout << "l " << edges(i, 0) + 1 << " " << edges(i, 1) + 1 << std::endl;

	}
	fout.close();
}

void go_through_temperature_interpolation() {
	const auto location_list_to_lists = [](
		const std::vector<std::vector<double>>& lct,std::vector<int>& names, std::vector<std::vector<double>>& locations) {
		names.resize(lct.size());
		locations.resize(lct.size());
		for (int i = 0; i < lct.size(); i++) {
			names[i] = lct[i][0];
			locations[i].push_back(lct[i][1]);
			locations[i].push_back(lct[i][2]);
		}
		return 1;
		
	};
	
	
	
	
	std::string borderfile = "D:\\vs\\sparse_data_interpolation\\sparse_data\\data\\Shanxi_border.csv";
	std::vector<std::vector<double>> border = read_CSV_data(borderfile);
	
	// all the border points
	Eigen::MatrixXd points = list_to_matrix_3d(border);
	// all the border points
	Eigen::MatrixXd knotP = points;
	Eigen::MatrixXi knotE;
	vertices_to_edges(points, knotE);
	
	//write_edges_obj("D:\\vs\\sparse_data_interpolation\\sparse_data\\data\\border_edges.obj", points, knotE);
	int expect_p_nbr = 50;
	double tolerance = 0.02;
	std::vector<int> pids;
	
	std::vector<Vector3d> border_list = matrix3d_to_vector(points);
	std::cout << "before remove redundant, size, " << border_list.size() << std::endl;
	border_list = border::border_remove_redundant_points(border_list);// remove redundant points
	std::cout << "after remove redundant, size, " << border_list.size() << std::endl;
	border::get_simp_points(border_list, expect_p_nbr, tolerance, pids);
	std::cout << "final feature points " <<pids.size()<< std::endl;

	// get the simplified points
	Eigen::MatrixXd knotP_simp = list_to_matrix_3d(border_list, pids);
	Eigen::MatrixXi knotE_simp;
	vertices_to_edges(knotP_simp, knotE_simp);
		//write_edges_obj("D:\\vs\\sparse_data_interpolation\\sparse_data\\data\\border_edges_simp.obj", knotP_simp, knotE_simp);
	
	std::vector<std::vector<double>> location_list = 
		read_CSV_data("D:\\vs\\sparse_data_interpolation\\sparse_data\\data\\station.csv");

	std::vector<int> station_ids;
	std::vector<std::vector<double>> location_list_clean;
	location_list_to_lists(location_list, station_ids, location_list_clean);
	int bordersize = knotP_simp.rows();// the size of the border points
	Eigen::MatrixXd location_matrix = list_to_matrix_3d(location_list_clean);
	Eigen::MatrixXd all_points(bordersize + location_matrix.rows(),3);
	
	// get all the points for parameterization
	all_points << knotP_simp, location_matrix;
	Eigen::VectorXi bnd(bordersize);
	for (int i = 0; i < bordersize; i++) {
		bnd[i] = i;
	}



	// get the boundary parameters
	Eigen::MatrixXd boundary_parameters;
	map_vertices_to_square(all_points, bnd, boundary_parameters);

	//Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	//fcolor << 1, 0, 0; ecolor << 0.5, 0.5, 0.5;
	//global_viewer_1.data().set_edges(knotP, knotE, fcolor);
	////global_viewer_1.data().add_points(points, ecolor);
	//global_viewer_1.launch();
}
