#pragma once
#include<Eigen/Core>
#include<array>
#include<iostream>
typedef  Eigen::Vector3d Vector3d;
static const int STAIR_FORWARD = 0;
static const int STAIR_BACKWARD = 1;
static const int STAIR_WHOLE = 2;
void print_vector(const std::vector<double>& input);
void print_vector(const std::vector<int>& input);