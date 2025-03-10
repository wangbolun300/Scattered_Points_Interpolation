#include"basis.h"
#include<iostream>
#include<assert.h>
#include<fstream>
namespace SIBSplines{
// the B-spline basis function N_{i,0}(u). U is the knot vector
double Ni0(const int i, const double u, const std::vector<double> &U) {
	if (u >= U[i] && u < U[i + 1]) return 1.0;
	
	return 0.0;
}

double handle_division(const double a, const double b) {
	
	if (b == 0) {
		// if the denominator is 0, then this term is 0
		return 0;
	}
	else return a / b;
}
// the B-spline basis function N_{i,p}(u)
double Nip(const int i, const int p, const double u, const std::vector<double> &U) {
	if (p == 0) {
		return Ni0(i, u, U);
	}
	if (u == U.back()) {
		if (i == U.size() - 2 - p) {
			return 1;
		}
		else {
			return 0;
		}
	}
	double result1 = handle_division(u - U[i], U[i + p] - U[i])*Nip(i, p - 1, u, U);
	double result2 = handle_division(U[i + p + 1] - u, U[i + p + 1] - U[i + 1])*Nip(i + 1, p - 1, u, U);
	return result1 + result2;
}

double Degree_1_Derivative(const int i, const int p, const double u, const std::vector<double> &U) {
	double result1 = handle_division(p, U[i + p] - U[i])*Nip(i, p - 1, u, U);
	double result2 = handle_division(p, U[i + p + 1] - U[i + 1])*Nip(i + 1, p - 1, u, U);
	return result1 - result2;
}

// the kth derivate of basis function
double kth_derivate_Nip(const int k, const int i, 
	const int p, const double u, const std::vector<double> &U) {
	assert(k <= p);
	if (k == 0) {// if k=0, return itself
		return Nip(i, p, u, U);
	}
	double result1 = handle_division(p, U[i + p] - U[i])*kth_derivate_Nip(k-1, i, p - 1, u, U);
	double result2 = handle_division(p, U[i + p + 1] - U[i + 1])*kth_derivate_Nip(k-1, i + 1, p - 1, u, U);
	return result1 - result2;
}
std::vector<double> Centripetal_parameterization(std::vector<Vector3d>& pts) {
	int nbr = pts.size();
	std::vector<double> paras(nbr);
	double total = 0;
	std::vector<double> distances(nbr - 1);
	for (int i = 0; i < nbr - 1; i++) {
		distances[i] = sqrt((pts[i + 1] - pts[i]).norm());
		total += distances[i];
	}
	paras[0] = 0;
	for (int i = 0; i < nbr - 1; i++) {
		paras[i + 1] = paras[i] + distances[i] / total;
	}
	return paras;
}

std::vector<double> Chord_parameterization(std::vector<Vector3d>& pts) {
	int nbr = pts.size();
	std::vector<double> paras(nbr);
	double total = 0;
	std::vector<double> distances(nbr - 1);
	for (int i = 0; i < nbr - 1; i++) {
		// the only difference with Centripetal_parameterization is here
		distances[i] = (pts[i + 1] - pts[i]).norm();
		total += distances[i];
	}
	paras[0] = 0;
	for (int i = 0; i < nbr - 1; i++) {
		paras[i + 1] = paras[i] + distances[i] / total;
	}
	return paras;
}
void print_vector(const std::vector<double>& input) {
	for (int i = 0; i < input.size(); i++) {
		std::cout << input[i] << ", ";
	}
	std::cout << std::endl;
}
void print_vector(const std::vector<int>& input) {
	for (int i = 0; i < input.size(); i++) {
		std::cout << input[i] << ", ";
	}
	std::cout << std::endl;
}
bool read_csv_data_lbl(const std::string fname, std::vector<std::vector<double>> &data){
    std::cout<<"reading "<<fname<<std::endl;
    if (fname.length() == 0)
        return false;
	int count = 0;
    std::ifstream infile;
    std::vector<std::vector<double>> results;

    infile.open(fname);
    if (!infile.is_open())
    {
        std::cout << "Path Wrong!!!!" << std::endl;
        std::cout << "path, " << fname << std::endl;
        return false;
    }

    int l = 0;
    while (infile) // there is input overload classfile
    {
        std::string s;
        if (!getline(infile, s))
            break;

        if (s[0] != '#')
        {
            std::istringstream ss(s);
            std::vector<double> record;
            int c = 0;
            while (ss)
            {
                std::string line;
                if (!getline(ss, line, ','))
                    break;
                try
                {

                    record.push_back(std::stod(line));
                    c++;

                }
                catch (const std::invalid_argument e)
                {
                    // std::cout << "NaN found in file " << fname
                    //           <<  std::endl;
					// 		  std::cout<<"line "<<line<<"\n";
					// 		  std::cout<<"count "<<count<<"\n";
					// 		  count++;
                    e.what();
                }
            }

            
            results.push_back(record);
        }
    }
    data = results;
    if (!infile.eof())
    {
        std::cerr << "Could not read file " << fname << "\n";
    }
    std::cout<<fname<<" get readed"<<std::endl;
    return true;
}
}