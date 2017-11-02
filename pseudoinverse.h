//Copyright (c) 2016-2017 Mr.Zhang
//All rights reserved.
//This code is based on Eigen library

#ifndef PSEUDOINVERSE_H
#define PSEUDOINVERSE_H
#define EPS 0.1e-15
#include "Eigen/Dense"
//#include "eigen-eigen-b9cd8366d4e8\Eigen\Dense"

using namespace Eigen;

/*Function*/
MatrixXd pinv(MatrixXd&,double tol=-1.0);
#endif
