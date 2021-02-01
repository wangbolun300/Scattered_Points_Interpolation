#pragma once
#include<Eigen/Core>
#include<array>
#include<iostream>
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector2i Vector2i;
static const int STAIR_FORWARD = 0;
static const int STAIR_BACKWARD = 1;
static const int STAIR_WHOLE = 2;
static const int STAIR_HIGHEST = 3;
void print_vector(const std::vector<double>& input);
void print_vector(const std::vector<int>& input);