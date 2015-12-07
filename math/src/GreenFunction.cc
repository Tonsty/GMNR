#include <GMNR/common.h>
#include <GMNR/math/GreenFunction.h>

namespace gmnr{
	Scalar greenFunc( Scalar x, int d){
		switch (d){
		case 1: return x*x*x;
		case 2: {
					if (x < Eigen::NumTraits<Scalar>::epsilon()) return 0;
					
					return x*x*log((double)x);
				}
		case 3: return x;
		default: return 0;
		}
	};

	Matrix greenFunc(const Matrix &_Xi, const Matrix &_Xj){

		Matrix M(_Xi.rows(), _Xj.rows());
		for (int i = 0; i < _Xi.rows(); i++ )
		{
			for (int j = 0; j < _Xj.rows(); j++)
			{
				M(i, j) = greenFunc( (_Xi.row(i) - _Xj.row(j)).norm(), _Xi.cols());
				//std::cout << "Xi = " << _Xi.row(i) << std::endl; 
				//std::cout << "Xj = " << _Xj.row(j) << std::endl; 
				//std::cout << "Mij = " << M(i, j) << std::endl; 
			}
		}
		return M;
	}
};