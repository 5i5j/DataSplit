//pseuodinverse.cpp
//Copyright (c) 2016-2017 Dr. Tim Yi
//All rights reserved.
#include "pseudoinverse.h"
/*************************************************************************
 *
 * pinv()
 *
 * Parameters:
 *
 * MatrixXd& A          - input matrix(dynamic double)
 *
 * double tol           - The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS(class(A)).  
 *
 * Return Value:
 *
 * MatrixXd             - return pseudoinverse matrix
 *
 * Description:
 *   1.[U,S,V] = svd(A) % A = U*S*V'
 *   2.Any singular values less than a tolerance are treated as zero.
 *   3. invA=V./S*U'
 ************************************************************************/
MatrixXd pinv(MatrixXd& A, double tol)
{
	int rows,cols;
	int r = 0;
	rows = A.rows();
	cols = A.cols();
	//Singular value decomposition
	JacobiSVD<MatrixXd>svd(A, ComputeThinU | ComputeThinV);
	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();
	MatrixXd S = svd.singularValues();
	//pseudoinverse matrix inv_A
	MatrixXd inv_A(cols,rows);
	//The default tolerance is MAX(SIZE(A))*NORM(S)*EPS
	if(tol <= 0 )
		tol = std::max(rows,cols) * S(0) * EPS;
	for( int i=0; i<S.size(); ++i )
		if(S(i) > tol )
			r++;
	//Any singular values less than a tolerance are treated as zero
	for( int i=0; i<cols; ++i )
	{
		for( int k=0; k<r; ++k )
			V(i,k)/= S(k);
	}
	//Calculation pseudoinverse matrix inv_A
	for( int i=0; i<cols; ++i )
		for( int j=0; j<rows; ++j )
		{
			double sum = 0;
			for( int k=0; k<r; ++k )
				sum += V(i,k)*U(j,k);
			inv_A(i,j)= sum;
		}
	return inv_A;
}
