#ifndef GMNR_COMMON_H
#define GMNR_COMMON_H

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace gmnr{
	typedef float Scalar;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

	typedef Eigen::Matrix<Scalar, 2, 2> Matrix2D;
	typedef Eigen::Matrix<Scalar, 2, 1> Vector2D;

	typedef Vector2D Point2D;
	typedef Vector2D Normal2D;

	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 2> PointSet2D;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 2> NormalSet2D;

	typedef Eigen::Matrix<Scalar, 3, 3> Matrix3D;
	typedef Eigen::Matrix<Scalar, 3, 1> Vector3D;

	typedef Vector3D Point3D;
	typedef Vector3D Normal3D;

	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3> PointSet3D;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3> NormalSet3D;

	typedef Eigen::Transform<Scalar, 3, Eigen::Affine> Affine3D;
	typedef Eigen::Translation<Scalar, 3> Translation3D;
	typedef Eigen::AngleAxis<Scalar> AngleAxis;
	typedef Eigen::Quaternion<Scalar> Quaternion;

	typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
	typedef Eigen::Triplet<Scalar> Triplet;
	typedef std::vector<Triplet> TripletList;
}

#endif