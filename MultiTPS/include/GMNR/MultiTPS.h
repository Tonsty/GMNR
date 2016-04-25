#ifndef GMNR_MULTITPS_MULTITPS_H
#define GMNR_MULTITPS_MULTITPS_H

#ifdef _WIN32
	#ifndef LIBMULTITPS_API
		#ifdef LIBMULTITPS_DYNAMIC
			#if LIBMULTITPS_BUILD
				#define LIBMULTITPS_API __declspec(dllexport)
			#else 
				#define LIBMULTITPS_API __declspec(dllimport)
			#endif
		#else 
			#define LIBMULTITPS_API
		#endif
	#endif
#else
	#define LIBMULTITPS_API
#endif

#ifdef _DEBUG
#include <iostream>
#endif
#include <vector>

#include <GMNR/ThinPlateSplines.h>

namespace gmnr{

	class MultiTPS{

#ifdef _DEBUG
		// used for debug
		friend LIBMULTITPS_API std::ostream& operator<<(std::ostream& os, const MultiTPS& tps);
#endif

	public:
		LIBMULTITPS_API MultiTPS();
		LIBMULTITPS_API MultiTPS(const Matrix &_X, 
			const Matrix &_Y,
			const std::vector<int> &_m,
			const std::vector<int> &_alpha,
			const std::vector<int> &_beta,
			const Vector &_kappa,
			const Vector &_lambda,
			const int &_max_iter_num = 0,
			const Scalar &_iter_rate = -1.0);
		const std::vector<TPSFunction>& getfs() const { return fs_;}
	protected:
		std::vector<TPSFunction> fs_;
	private:
		LIBMULTITPS_API void solve_new_(const Matrix &_X, 
			const Matrix &_Y,
			const std::vector<int> &_m,
			const std::vector<int> &_alpha,
			const std::vector<int> &_beta,
			const Vector &_kappa,
			const Vector &_lambda);
		LIBMULTITPS_API void solve_old_(const Matrix &_X, 
			const Matrix &_Y,
			const std::vector<int> &_m,
			const std::vector<int> &_alpha,
			const std::vector<int> &_beta,
			const Vector &_kappa,
			const Vector &_lambda);

		LIBMULTITPS_API void solve_new_sparse_(const Matrix &_X, 
			const Matrix &_Y,
			const std::vector<int> &_m,
			const std::vector<int> &_alpha,
			const std::vector<int> &_beta,
			const Vector &_kappa,
			const Vector &_lambda);

		LIBMULTITPS_API void solve_iterative_(const Matrix &_X, 
			const Matrix &_Y,
			const std::vector<int> &_m,
			const std::vector<int> &_alpha,
			const std::vector<int> &_beta,
			const Vector &_kappa,
			const Vector &_lambda,
			const int &_max_iter_num,
			const Scalar &_iter_rate);
	};

	class ApproxiMultiTPS : public MultiTPS {

	public:
		LIBMULTITPS_API ApproxiMultiTPS();
		LIBMULTITPS_API ApproxiMultiTPS(const Matrix &_X, 
			const Matrix &_Y,
			const std::vector<int> &_m,
			const std::vector<int> &_alpha,
			const std::vector<int> &_beta,
			const Vector &_kappa,
			const Vector &_lambda,
			const std::vector<int> &_na,
			const std::vector<int> &_nb,
			const int &_max_iter_num = 0,
			const Scalar &_iter_rate = -1.0);
	};
};

#endif