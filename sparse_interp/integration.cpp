#include<integration.h>
#include<iostream>
#include<fstream>
#include <igl/opengl/glfw/Viewer.h>
#include<test.h>
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

Eigen::MatrixXd list_to_matrix_3d(const std::vector<std::vector<double>>& v) {
	Eigen::MatrixXd result(v.size(), 3);
	for (int i = 0; i < v.size(); i++) {
		result(i, 0) = v[i][0];
		result(i, 1) = v[i][1];
		result(i, 2) = v[i].size() >= 3 ? v[i][2] : 0;
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

void visual_border() {
	std::string borderfile = "D:\\vs\\sparse_data_interpolation\\sparse_data\\data\\Shanxi_border.csv";
	std::vector<std::vector<double>> border = read_CSV_data(borderfile);
	
	Eigen::MatrixXd points = list_to_matrix_3d(border);
	Eigen::MatrixXd knotP = points;
	Eigen::MatrixXi knotE;
	vertices_to_edges(points, knotE);
	std::cout << "front," << knotP.row(0) << std::endl;
	std::cout << "back," << knotP.row(knotP.rows()-1) << std::endl;
	//write_edges_obj("D:\\vs\\sparse_data_interpolation\\sparse_data\\data\\border_edges.obj", points, knotE);
	int expect_p_nbr = 50;
	double tolerance = 0.1;
	std::vector<int> pids;
	
	std::vector<Vector3d> border_list = matrix3d_to_vector(points);
	std::cout << "before remove redundant, size, " << border_list.size() << std::endl;
	border_list = border::border_remove_redundant_points(border_list);
	std::cout << "after remove redundant, size, " << border_list.size() << std::endl;
	border::get_simp_points(border_list, expect_p_nbr, tolerance, pids);
	std::cout << "final feature points " <<pids.size()<< std::endl;

	Eigen::MatrixXd fcolor(1, 3), ecolor(1, 3);
	fcolor << 1, 0, 0; ecolor << 0.5, 0.5, 0.5;
	//global_viewer_1.data().set_edges(knotP, knotE, fcolor);
	////global_viewer_1.data().add_points(points, ecolor);
	//global_viewer_1.launch();
}
