#include <iostream>
#include <fstream>
#include <unordered_map>
#include <math.h>
#include <Eigen/Core>
#include <GMNR/MultiTPS.h>

namespace gmnr{

	void setSmallValueToZero(Matrix &A, Scalar threshold = 5e-6) {
		for (int i = 0; i < A.size(); i++) {
			A(i) = abs(A(i)) > threshold ? A(i) : 0; 
		}
	}

	MultiTPS::MultiTPS(){}

	MultiTPS::MultiTPS(const Matrix &_X, 
		const Matrix &_Y,
		const std::vector<int> &_m,
		const std::vector<int> &_alpha,
		const std::vector<int> &_beta,
		const Vector &_kappa,
		const Vector &_lambda,
		const int &_max_iter_num,
		const Scalar &_iter_rate) {
			//solve_old_(_X, _Y, _m, _alpha, _beta, _kappa, _lambda);

			//std::vector<TPSFunction> old_fs = fs_; 
			//for (int i = 0; i < old_fs.size(); i++){
			//	std::cout << "old_fs[" << i << "]=\n" << old_fs[i] << std::endl;
			//}

			if(_max_iter_num > 0) solve_iterative_(_X, _Y, _m, _alpha, _beta, _kappa, _lambda, _max_iter_num, _iter_rate);
			else solve_new_sparse_(_X, _Y, _m, _alpha, _beta, _kappa, _lambda);
			
			//std::vector<TPSFunction> new_fs = fs_; 
			//for (int i = 0; i < new_fs.size(); i++){
			//	std::cout << "new_fs[" << i << "]=\n" << new_fs[i] << std::endl;
			//}
	}

	template<typename T>
	std::ostream& operator << (std::ostream &os, const std::vector<T> &vec){
		for (int i = 0; i < vec.size(); i++){
			os << vec[i] << std::endl;
		}
		return os;
	}

	template<typename K, typename V>
	std::ostream& operator << (std::ostream &os, const std::unordered_map<K,V> &map){
		for (std::unordered_map<K,V>::const_iterator iter = map.begin(); iter != map.end(); iter++){
			os << (*iter).first << "\t" << (*iter).second << std::endl;
		}
		return os;
	}

	void MultiTPS::solve_new_(const Matrix &_X, 
		const Matrix &_Y,
		const std::vector<int> &_m,
		const std::vector<int> &_alpha,
		const std::vector<int> &_beta,
		const Vector &_kappa,
		const Vector &_lambda) {
			int P = _m.size(), N = _lambda.size();
			int d = _X.cols();
			int m = 0; // or m = _X.rows();
			for (int i = 0; i < P; i++) {
				m += _m[i];
			}

			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
			Matrix Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			std::vector<int> Xpos(P), Ypos(P), Apos(P);

			Xpos[0] = 0;
			for (int mu = 1; mu < P; mu++) {
				Xpos[mu] = Xpos[mu-1] + _m[mu-1];
			}
			Ypos[0] = 0;
			for (int mu = 1; mu < P; mu++) {
				Ypos[mu] = Ypos[mu-1] + _m[mu-1];
			}
			Apos[0] = 0;
			for (int mu = 1; mu < P; mu++) {
				Apos[mu] = Apos[mu-1] + _m[mu-1];
			}

			//std::cout << "Xpos = \n" << Xpos << std::endl;
			//std::cout << "Ypos = \n" << Ypos << std::endl;

			std::vector< std::vector<int> > mu_alpha(N), mu_beta(N);
			std::vector< std::unordered_map<int, int> > bar_t(N), bar_s(N);
			for (int mu = 0; mu < P; mu++) {
				int k = _alpha[mu];
				bar_t[k][mu] = mu_alpha[k].size();
				mu_alpha[k].push_back(mu);
			}
			for (int mu = 0; mu < P; mu++) {
				int k = _beta[mu];
				bar_s[k][mu] = mu_beta[k].size();
				mu_beta[k].push_back(mu);
			}

			//for (int k = 0; k < N; k++){
			//	std::cout << "mu_alpha[" << k << "]= \n" << mu_alpha[k] << std::endl;
			//	std::cout << "mu_beta[" << k << "]= \n" << mu_beta[k] << std::endl;
			//}
			//for (int k = 0; k < N; k++){
			//	std::cout << "bar_t[" << k << "]= \n" << bar_t[k] << std::endl;
			//	std::cout << "bar_s[" << k << "]= \n" << bar_s[k] << std::endl;
			//}

			std::vector<int> sum_m_mu_a(N), sum_m_mu_b(N);
			std::vector< std::vector<int> > M_a_pos(N), M_b_pos(N);

			for (int k = 0; k < N; k++) {
				int tk = mu_alpha[k].size(), sk = mu_beta[k].size();

				M_a_pos[k].resize(tk);
				int sum_m_mu_a_k = 0;
				for (int t = 0; t < tk; t++) {
					M_a_pos[k][t] = sum_m_mu_a_k;
					int mu_a_k_t = mu_alpha[k][t];
					sum_m_mu_a_k += _m[mu_a_k_t];
				}
				sum_m_mu_a[k] = sum_m_mu_a_k;

				M_b_pos[k].resize(sk);
				int sum_m_mu_b_k = 0;
				for (int s = 0; s < sk; s++) {
					M_b_pos[k][s] = sum_m_mu_b_k;
					int mu_b_k_s = mu_beta[k][s];
					sum_m_mu_b_k += _m[mu_b_k_s];
				}
				sum_m_mu_b[k] = sum_m_mu_b_k;
			}

			//std::cout << "sum_m_mu_a \n" << sum_m_mu_a << std::endl;
			//std::cout << "sum_m_mu_b \n" << sum_m_mu_b << std::endl;
			//for (int k = 0; k < N; k++){
			//	std::cout << "M_a_pos[" << k << "]= \n" << M_a_pos[k] << std::endl;
			//	std::cout << "M_b_pos[" << k << "]= \n" << M_b_pos[k] << std::endl;
			//}

			std::vector< Matrix > M(N);
			for (int k = 0; k < N; k++) {
				int tk = mu_alpha[k].size(), sk = mu_beta[k].size();
				int sum_m_mu_a_k = sum_m_mu_a[k], sum_m_mu_b_k = sum_m_mu_b[k];

				Matrix M_k_a_a(sum_m_mu_a_k, sum_m_mu_a_k);
				Matrix M_k_a_b(sum_m_mu_a_k, sum_m_mu_b_k);
				Matrix M_k_b_a(sum_m_mu_b_k, sum_m_mu_a_k);
				Matrix M_k_b_b(sum_m_mu_b_k, sum_m_mu_b_k);

				for (int ti = 0, ipos = 0; ti < tk; ti++) {
					int mu_a_k_ti = mu_alpha[k][ti];
					for (int tj = 0, jpos = 0; tj < tk; tj++) {
						int mu_a_k_tj = mu_alpha[k][tj];
						M_k_a_a.block(ipos, jpos, _m[mu_a_k_ti], _m[mu_a_k_tj]) = 
							greenFunc(_X.block(Xpos[mu_a_k_ti], 0, _m[mu_a_k_ti], d), _X.block(Xpos[mu_a_k_tj], 0, _m[mu_a_k_tj], d));
						jpos += _m[mu_a_k_tj];
					}
					ipos += _m[mu_a_k_ti];
				}
				for (int ti = 0, ipos = 0; ti < tk; ti++) {
					int mu_a_k_ti = mu_alpha[k][ti];
					for (int sj = 0, jpos = 0; sj < sk; sj++) {
						int mu_b_k_sj = mu_beta[k][sj];
						M_k_a_b.block(ipos, jpos, _m[mu_a_k_ti], _m[mu_b_k_sj]) = 
							greenFunc(_X.block(Xpos[mu_a_k_ti], 0, _m[mu_a_k_ti], d), _Y.block(Ypos[mu_b_k_sj], 0, _m[mu_b_k_sj], d));
						jpos += _m[mu_b_k_sj];
					}
					ipos += _m[mu_a_k_ti];
				}
				for (int si = 0, ipos = 0; si < sk; si++) {
					int mu_b_k_si = mu_beta[k][si];
					for (int tj = 0, jpos = 0; tj < tk; tj++) {
						int mu_a_k_tj = mu_alpha[k][tj];
						M_k_b_a.block(ipos, jpos, _m[mu_b_k_si], _m[mu_a_k_tj]) = 
							greenFunc(_Y.block(Ypos[mu_b_k_si], 0, _m[mu_b_k_si], d), _X.block(Xpos[mu_a_k_tj], 0, _m[mu_a_k_tj], d));
						jpos += _m[mu_a_k_tj];
					}
					ipos += _m[mu_b_k_si];
				}
				for (int si = 0, ipos = 0; si < sk; si++) {
					int mu_b_k_si = mu_beta[k][si];
					for (int sj = 0, jpos = 0; sj < sk; sj++) {
						int mu_b_k_sj = mu_beta[k][sj];
						M_k_b_b.block(ipos, jpos, _m[mu_b_k_si], _m[mu_b_k_sj]) = 
							greenFunc(_Y.block(Ypos[mu_b_k_si], 0, _m[mu_b_k_si], d), _Y.block(Ypos[mu_b_k_sj], 0, _m[mu_b_k_sj], d));
						jpos += _m[mu_b_k_sj];
					}
					ipos += _m[mu_b_k_si];
				}
				M[k].resize(sum_m_mu_a_k + sum_m_mu_b_k, sum_m_mu_a_k + sum_m_mu_b_k);
				if(sum_m_mu_a_k == 0){
					M[k] << M_k_b_b;
				}else if(sum_m_mu_b_k == 0){
					M[k] << M_k_a_a;
				}else{
					M[k] << M_k_a_a, M_k_a_b,
						M_k_b_a, M_k_b_b;
				}
			}

			Matrix Ca = Matrix::Zero(P*(d+1), N*(d+1)), Cb = Matrix::Zero(P*(d+1), N*(d+1));
			for (int mu = 0, pos = 0; mu < P; mu++) {
				int ai = _alpha[mu], bi = _beta[mu];
				Ca.block(pos, ai * (d+1), d+1, d+1) = Matrix::Identity(d+1, d+1);
				Cb.block(pos, bi * (d+1), d+1, d+1) = Matrix::Identity(d+1, d+1);
				pos += d+1;
			}

			std::vector<Matrix> D(N);
			for (int k = 0; k < N; k++) {
				int tk = mu_alpha[k].size(), sk = mu_beta[k].size();
				int sum_m_mu_a_k = sum_m_mu_a[k], sum_m_mu_b_k = sum_m_mu_b[k];
				Matrix D_k_a = Matrix::Zero(sum_m_mu_a_k, m), D_k_b = Matrix::Zero(sum_m_mu_b_k, m);
				for (int t = 0, ipos = 0; t < tk; t++) {
					int mu = mu_alpha[k][t];
					D_k_a.block(ipos, Xpos[mu], _m[mu], _m[mu]) = Matrix::Identity(_m[mu], _m[mu]);
					ipos += _m[mu];
				}
				for (int s = 0, ipos = 0; s < sk; s++) {
					int mu = mu_beta[k][s];
					D_k_b.block(ipos, Ypos[mu], _m[mu], _m[mu]) = Matrix::Identity(_m[mu], _m[mu]);
					ipos += _m[mu_beta[k][s]];
				}
				D[k].resize(sum_m_mu_a_k + sum_m_mu_b_k, 2*m);
				if(sum_m_mu_a_k == 0){
					D[k] << Matrix::Zero(sum_m_mu_b_k, m), D_k_b;
				}else if (sum_m_mu_b_k == 0){
					D[k] << D_k_a, Matrix::Zero(sum_m_mu_a_k, m);
				}else {
					D[k] << D_k_a, Matrix::Zero(sum_m_mu_a_k, m),
						Matrix::Zero(sum_m_mu_b_k, m), D_k_b;
				}
			}

			//for (int k = 0; k < N; k++){
			//	std::cout << "D[" << k << "]= \n" << D[k] << std::endl;
			//}

			Matrix Ha(m, 2*m), Hb(m, 2*m), S(m, 2*m);
			for (int mu = 0, pos = 0; mu < P; mu++) {
				int i = _alpha[mu], j = _beta[mu];

				Matrix bar_M_i_mu_t(_m[mu], 2*m), bar_M_j_mu_s(_m[mu], 2*m);

				Matrix D_i;
				int sum_m_mu_a_i = sum_m_mu_a[i], sum_m_mu_b_i = sum_m_mu_b[i];
				D_i = D[i];
				Matrix D_j;
				int sum_m_mu_a_j = sum_m_mu_a[j], sum_m_mu_b_j = sum_m_mu_b[j];
				D_j = D[j];

				Matrix M_i_a_bar_t_i_mu;
				int bar_t_i_mu = bar_t[i][mu];
				M_i_a_bar_t_i_mu = M[i].block( M_a_pos[i][bar_t_i_mu], 0, _m[mu], sum_m_mu_a_i + sum_m_mu_b_i);

				Matrix M_j_b_bar_s_j_mu;
				int bar_s_j_mu = bar_s[j][mu];
				M_j_b_bar_s_j_mu = M[j].block(sum_m_mu_a_j + M_b_pos[j][bar_s_j_mu], 0, _m[mu], sum_m_mu_a_j + sum_m_mu_b_j);

				bar_M_i_mu_t = M_i_a_bar_t_i_mu * D_i;
				bar_M_j_mu_s = M_j_b_bar_s_j_mu * D_j;
				
				Ha.block(pos, 0, _m[mu], 2*m) = bar_M_i_mu_t;
				Hb.block(pos, 0, _m[mu], 2*m) = bar_M_j_mu_s;
				pos += _m[mu];
			}
			S = Ha - Hb;

			//std::cout << "Ha = \n" << Ha << std::endl;
			//std::cout << "Hb = \n" << Hb << std::endl;
			//std::cout << "S = \n" << S << std::endl;

			Matrix Ja(m, N * (d+1)), Jb(m, N * (d+1)), T(m, N * (d+1));
			Matrix Ka = Matrix::Zero(m, m), Kb = Matrix::Zero(m, m);
			for (int mu = 0, pos = 0; mu < P; mu++) {
				Ja.block(pos, 0, _m[mu], N * (d+1)) = X.block(pos, 0, _m[mu], d+1) * Ca.block(mu * (d+1), 0, d+1, N * (d+1));
				Jb.block(pos, 0, _m[mu], N * (d+1)) = Y.block(pos, 0, _m[mu], d+1) * Cb.block(mu * (d+1), 0, d+1, N * (d+1));
				Ka.block(pos, pos, _m[mu], _m[mu]) = Matrix::Identity(_m[mu], _m[mu]) * _kappa(_alpha[mu]);
				Kb.block(pos, pos, _m[mu], _m[mu]) = Matrix::Identity(_m[mu], _m[mu]) * _kappa(_beta[mu]);
				pos += _m[mu];
			}
			T = Ja - Jb;                                                                                                                                                                                                                

			//Matrix diag_lambda_M = Matrix::Zero(2*m, 2*m);
			//for (int k = 0, pos = 0; k < N; k++) {
			//	int sum_m_mu_a_k = sum_m_mu_a[k], sum_m_mu_b_k = sum_m_mu_b[k];
			//	diag_lambda_M.block(pos, pos, sum_m_mu_a_k + sum_m_mu_b_k, sum_m_mu_a_k + sum_m_mu_b_k) = M[k] * _lambda(k);
			//	pos += sum_m_mu_a_k + sum_m_mu_b_k;
			//}
			//Matrix D_matrice(2*m, 2*m);
			//for (int k = 0, pos = 0; k < N; k++){
			//	D_matrice.block(pos, 0, sum_m_mu_a[k] + sum_m_mu_b[k], 2*m) = D[k];
			//	pos += sum_m_mu_a[k] + sum_m_mu_b[k];
			//}
			//Matrix D_T_diag_lambda_M_D;
			//D_T_diag_lambda_M_D = D_matrice.transpose() * diag_lambda_M * D_matrice;
			Matrix D_T_diag_lambda_M_D = Matrix::Zero(2*m, 2*m);
			for (int k = 0; k < N; k++) {
				D_T_diag_lambda_M_D += _lambda[k] * D[k].transpose() * M[k] * D[k];
			}

			Matrix P11(2*m, 2*m), P12(2*m, N * (d+1));
			P11 = S.transpose() * S + Ha.transpose() * Ka * Ha + Hb.transpose() * Kb * Hb;
			P12 = S.transpose() * T + Ha.transpose() * Ka * Ja + Hb.transpose() * Kb * Jb;
			Matrix Z1(2*m, d);
			Z1 = Ha.transpose() * Ka * _X + Hb.transpose() * Kb * _Y;
			Matrix P21(N * (d+1), 2*m), P22(N * (d+1), N * (d+1));
			P21 = T.transpose() * S + Ja.transpose() * Ka * Ha + Jb.transpose() * Kb * Hb;
			P22 = T.transpose() * T + Ja.transpose() * Ka * Ja + Jb.transpose() * Kb * Jb;
			Matrix Z2(N * (d+1), d);
			Z2 = Ja.transpose() * Ka * _X + Jb.transpose() * Kb * _Y;

			Matrix equationL = Matrix::Zero(2*m + N * (d+1), 2*m + N * (d+1));
			equationL << P11 + D_T_diag_lambda_M_D, P12,
						P21, P22;

			Matrix equationR = Matrix::Zero(2*m + N * (d+1), d);
			equationR << Z1, Z2;

			Matrix solution(2*m + N * (d+1), d);

			//solution = (equationL.transpose() * equationL).inverse() * (equationL.transpose() * equationR);

			Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
			std::cout << "Condition Number : " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			std::cout << "Singular Vector : \n" << svd.singularValues() << std::endl;
			std::cout << "Rows = " << equationL.rows() << " Columns = " << equationL.cols() << std::endl;
			std::cout << "before reset threshold, rank = " << svd.rank() << std::endl;
			svd.setThreshold(Eigen::NumTraits<Scalar>::epsilon());
			std::cout << "after reset threshold, rank = " << svd.rank() << std::endl;
			solution = svd.solve(equationR);

			//std::cout << "Ax - b = \n" << equationL * solution - equationR << std::endl; 

			Matrix A(2*m, d), B(N * (d+1), d);
			A = solution.block(0, 0, 2*m, d);
			B = solution.block(2*m, 0, N * (d+1), d);

			Scalar tps_energy1 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy1 += ( T.block(pos, 0, _m[mu], N * (d+1)) * B + S.block(pos, 0, _m[mu], 2*m) * A ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy1 = " << tps_energy1 << std::endl;
			Scalar tps_energy2 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy2 += _kappa(_alpha[mu]) * 
					( Ja.block(pos, 0, _m[mu], N * (d+1) ) * B + Ha.block(pos, 0, _m[mu], 2*m) * A - _X.block(pos, 0, _m[mu], d) ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy2 = " << tps_energy2 << std::endl;
			Scalar tps_energy3 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy3 += _kappa(_beta[mu]) *
					( Jb.block(pos, 0, _m[mu], N * (d+1) ) * B + Hb.block(pos, 0, _m[mu], 2*m) * A - _Y.block(pos, 0, _m[mu], d) ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy3 = " << tps_energy3 << std::endl;
			//Scalar tps_energy4 = (A.transpose() * D_matrice.transpose() * diag_lambda_M * D_matrice * A).trace();
			Scalar tps_energy4 = (A.transpose() * D_T_diag_lambda_M_D * A).trace();
			std::cout << "tps_energy4 = " << tps_energy4 << std::endl;

			fs_.resize(N);
			for (int k = 0; k < N; k++){
				Matrix Ak(sum_m_mu_a[k] + sum_m_mu_b[k], d), Xk(sum_m_mu_a[k] + sum_m_mu_b[k], d);
				for (int ti = 0, pos = 0; ti < mu_alpha[k].size(); ti++){
					Ak.block(pos, 0, _m[mu_alpha[k][ti]], d) = 
						A.block(Apos[mu_alpha[k][ti]], 0, _m[mu_alpha[k][ti]], d);
					Xk.block(pos, 0, _m[mu_alpha[k][ti]], d) = 
						_X.block(Xpos[mu_alpha[k][ti]], 0, _m[mu_alpha[k][ti]], d);
					pos += _m[mu_alpha[k][ti]];
				}
				for (int sj = 0, pos = 0; sj < mu_beta[k].size(); sj++) {
					Ak.block(sum_m_mu_a[k] + pos, 0, _m[mu_beta[k][sj]], d) = 
						A.block(m + Apos[mu_beta[k][sj]], 0, _m[mu_beta[k][sj]], d);
					Xk.block(sum_m_mu_a[k] + pos, 0, _m[mu_beta[k][sj]], d) =
						_Y.block(Ypos[mu_beta[k][sj]], 0, _m[mu_beta[k][sj]], d);
					pos += _m[mu_beta[k][sj]];
				}
				fs_[k].setX(Xk);
				fs_[k].setA(Ak);
			}
			for (int k = 0, kpos = 0; k < N; k++){
				Matrix Bk = B.block(kpos, 0, d+1, d);
				fs_[k].setB(Bk);
				kpos += d+1;
			}
	}

	void MultiTPS::solve_new_sparse_(const Matrix &_X, 
		const Matrix &_Y,
		const std::vector<int> &_m,
		const std::vector<int> &_alpha,
		const std::vector<int> &_beta,
		const Vector &_kappa,
		const Vector &_lambda) {
			int P = _m.size(), N = _lambda.size();
			int d = _X.cols();
			int m = 0; // or m = _X.rows();
			for (int i = 0; i < P; i++) {
				m += _m[i];
			}

			std::cout << "N = " << N << std::endl;
			std::cout << "P = " << P << std::endl;
			std::cout << "m = " << m << std::endl;
			for (int i = 0; i < P; i++) {
				std::cout << "( " << _alpha[i] << ", " << _beta[i] << ") " << _m[i] << std::endl;
			}
			std::cout << "_X : rows = " << _X.rows() << " cols = " << _X.cols() << std::endl;
			std::cout << "_Y : rows = " << _Y.rows() << " cols = " << _Y.cols() << std::endl;

			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
			Matrix Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			std::vector<int> Xpos(P), Ypos(P), Apos(P);

			Xpos[0] = 0;
			for (int mu = 1; mu < P; mu++) {
				Xpos[mu] = Xpos[mu-1] + _m[mu-1];
			}
			Ypos[0] = 0;
			for (int mu = 1; mu < P; mu++) {
				Ypos[mu] = Ypos[mu-1] + _m[mu-1];
			}
			Apos[0] = 0;
			for (int mu = 1; mu < P; mu++) {
				Apos[mu] = Apos[mu-1] + _m[mu-1];
			}

			//std::cout << "Xpos = \n" << Xpos << std::endl;
			//std::cout << "Ypos = \n" << Ypos << std::endl;
			//std::cout << "Apos = \n" << Apos << std::endl;

			std::vector< std::vector<int> > mu_alpha(N), mu_beta(N);
			std::vector< std::unordered_map<int, int> > bar_t(N), bar_s(N);
			for (int mu = 0; mu < P; mu++) {
				int k = _alpha[mu];
				bar_t[k][mu] = mu_alpha[k].size();
				mu_alpha[k].push_back(mu);
			}
			for (int mu = 0; mu < P; mu++) {
				int k = _beta[mu];
				bar_s[k][mu] = mu_beta[k].size();
				mu_beta[k].push_back(mu);
			}

			//for (int k = 0; k < N; k++){
			//	std::cout << "mu_alpha[" << k << "]= \n" << mu_alpha[k] << std::endl;
			//	std::cout << "mu_beta[" << k << "]= \n" << mu_beta[k] << std::endl;
			//}
			//for (int k = 0; k < N; k++){
			//	std::cout << "bar_t[" << k << "]= \n" << bar_t[k] << std::endl;
			//	std::cout << "bar_s[" << k << "]= \n" << bar_s[k] << std::endl;
			//}

			std::vector<int> sum_m_mu_a(N), sum_m_mu_b(N);
			std::vector< std::vector<int> > M_a_pos(N), M_b_pos(N);

			for (int k = 0; k < N; k++) {
				int tk = mu_alpha[k].size(), sk = mu_beta[k].size();

				M_a_pos[k].resize(tk);
				int sum_m_mu_a_k = 0;
				for (int t = 0; t < tk; t++) {
					M_a_pos[k][t] = sum_m_mu_a_k;
					int mu_a_k_t = mu_alpha[k][t];
					sum_m_mu_a_k += _m[mu_a_k_t];
				}
				sum_m_mu_a[k] = sum_m_mu_a_k;

				M_b_pos[k].resize(sk);
				int sum_m_mu_b_k = 0;
				for (int s = 0; s < sk; s++) {
					M_b_pos[k][s] = sum_m_mu_b_k;
					int mu_b_k_s = mu_beta[k][s];
					sum_m_mu_b_k += _m[mu_b_k_s];
				}
				sum_m_mu_b[k] = sum_m_mu_b_k;
			}

			//std::cout << "sum_m_mu_a \n" << sum_m_mu_a << std::endl;
			//std::cout << "sum_m_mu_b \n" << sum_m_mu_b << std::endl;
			//for (int k = 0; k < N; k++){
			//	std::cout << "M_a_pos[" << k << "]= \n" << M_a_pos[k] << std::endl;
			//	std::cout << "M_b_pos[" << k << "]= \n" << M_b_pos[k] << std::endl;
			//}

			std::vector< Matrix > M(N);
			for (int k = 0; k < N; k++) {
				int tk = mu_alpha[k].size(), sk = mu_beta[k].size();
				int sum_m_mu_a_k = sum_m_mu_a[k], sum_m_mu_b_k = sum_m_mu_b[k];

				Matrix M_k_a_a(sum_m_mu_a_k, sum_m_mu_a_k);
				Matrix M_k_a_b(sum_m_mu_a_k, sum_m_mu_b_k);
				Matrix M_k_b_a(sum_m_mu_b_k, sum_m_mu_a_k);
				Matrix M_k_b_b(sum_m_mu_b_k, sum_m_mu_b_k);

				for (int ti = 0, ipos = 0; ti < tk; ti++) {
					int mu_a_k_ti = mu_alpha[k][ti];
					for (int tj = 0, jpos = 0; tj < tk; tj++) {
						int mu_a_k_tj = mu_alpha[k][tj];
						M_k_a_a.block(ipos, jpos, _m[mu_a_k_ti], _m[mu_a_k_tj]) = 
							greenFunc(_X.block(Xpos[mu_a_k_ti], 0, _m[mu_a_k_ti], d), _X.block(Xpos[mu_a_k_tj], 0, _m[mu_a_k_tj], d));
						jpos += _m[mu_a_k_tj];
					}
					ipos += _m[mu_a_k_ti];
				}
				for (int ti = 0, ipos = 0; ti < tk; ti++) {
					int mu_a_k_ti = mu_alpha[k][ti];
					for (int sj = 0, jpos = 0; sj < sk; sj++) {
						int mu_b_k_sj = mu_beta[k][sj];
						M_k_a_b.block(ipos, jpos, _m[mu_a_k_ti], _m[mu_b_k_sj]) = 
							greenFunc(_X.block(Xpos[mu_a_k_ti], 0, _m[mu_a_k_ti], d), _Y.block(Ypos[mu_b_k_sj], 0, _m[mu_b_k_sj], d));
						jpos += _m[mu_b_k_sj];
					}
					ipos += _m[mu_a_k_ti];
				}
				for (int si = 0, ipos = 0; si < sk; si++) {
					int mu_b_k_si = mu_beta[k][si];
					for (int tj = 0, jpos = 0; tj < tk; tj++) {
						int mu_a_k_tj = mu_alpha[k][tj];
						M_k_b_a.block(ipos, jpos, _m[mu_b_k_si], _m[mu_a_k_tj]) = 
							greenFunc(_Y.block(Ypos[mu_b_k_si], 0, _m[mu_b_k_si], d), _X.block(Xpos[mu_a_k_tj], 0, _m[mu_a_k_tj], d));
						jpos += _m[mu_a_k_tj];
					}
					ipos += _m[mu_b_k_si];
				}
				for (int si = 0, ipos = 0; si < sk; si++) {
					int mu_b_k_si = mu_beta[k][si];
					for (int sj = 0, jpos = 0; sj < sk; sj++) {
						int mu_b_k_sj = mu_beta[k][sj];
						M_k_b_b.block(ipos, jpos, _m[mu_b_k_si], _m[mu_b_k_sj]) = 
							greenFunc(_Y.block(Ypos[mu_b_k_si], 0, _m[mu_b_k_si], d), _Y.block(Ypos[mu_b_k_sj], 0, _m[mu_b_k_sj], d));
						jpos += _m[mu_b_k_sj];
					}
					ipos += _m[mu_b_k_si];
				}
				M[k].resize(sum_m_mu_a_k + sum_m_mu_b_k, sum_m_mu_a_k + sum_m_mu_b_k);
				if(sum_m_mu_a_k == 0){
					M[k] << M_k_b_b;
				}else if(sum_m_mu_b_k == 0){
					M[k] << M_k_a_a;
				}else{
					M[k] << M_k_a_a, M_k_a_b,
						M_k_b_a, M_k_b_b;
				}

				std::cout << "M[" << k << "] : rows = " << M[k].rows() << " cols = " << M[k].cols() << std::endl;
			}

			std::vector<SparseMatrix> Ca(P), Cb(P);
			for (int mu = 0; mu < P; mu++) {
				TripletList Ca_coefficients, Cb_coefficients;
				int ai = _alpha[mu], bi = _beta[mu];
				for (int di = 0; di < d+1; di++ ){
					Ca_coefficients.push_back(Triplet(di, ai * (d+1) + di, 1));
					Cb_coefficients.push_back(Triplet(di, bi * (d+1) + di, 1));
				}
				Ca[mu].resize(d+1, N*(d+1));
				Ca[mu].setFromTriplets(Ca_coefficients.begin(), Ca_coefficients.end());
				Cb[mu].resize(d+1, N*(d+1));
				Cb[mu].setFromTriplets(Cb_coefficients.begin(), Cb_coefficients.end());

				std::cout << "Ca[" << mu << "] : rows = " << Ca[mu].rows() << " cols = " << Ca[mu].cols() << " nonzeros = " << Ca[mu].nonZeros() << std::endl;
				std::cout << "Cb[" << mu << "] : rows = " << Cb[mu].rows() << " cols = " << Cb[mu].cols() << " nonzeros = " << Cb[mu].nonZeros() << std::endl;
			}

			//for (int mu = 0; mu < P; mu++){
			//	std::cout << "Ca[" << mu << "]= \n" << Ca[mu] << std::endl;
			//	std::cout << "Cb[" << mu << "]= \n" << Cb[mu]<< std::endl;
			//}

			std::vector<SparseMatrix> D(N);
			for (int k = 0; k < N; k++) {
				int tk = mu_alpha[k].size(), sk = mu_beta[k].size();
				int sum_m_mu_a_k = sum_m_mu_a[k], sum_m_mu_b_k = sum_m_mu_b[k];
				SparseMatrix D_k(sum_m_mu_a_k + sum_m_mu_b_k, 2*m);
				TripletList D_k_coefficients;
				for (int t = 0, ipos = 0; t < tk; t++) {
					int mu = mu_alpha[k][t];
					for (int ti = 0; ti < _m[mu]; ti++) {
						D_k_coefficients.push_back(Triplet(ipos + ti, Xpos[mu] + ti, 1));
					}
					ipos += _m[mu];
				}
				for (int s = 0, ipos = 0; s < sk; s++) {
					int mu = mu_beta[k][s];
					for (int si = 0; si < _m[mu]; si++) {
						D_k_coefficients.push_back(Triplet(sum_m_mu_a_k + ipos + si, m + Ypos[mu] + si , 1));
					}
					ipos += _m[mu_beta[k][s]];
				}

				D_k.setFromTriplets(D_k_coefficients.begin(), D_k_coefficients.end());

				D[k] = D_k;

				std::cout << "D[" << k << "] : rows = " << D[k].rows() << " cols = " << D[k].cols() << " nonzeros = " << D[k].nonZeros() << std::endl;
			}

			//for (int k = 0; k < N; k++){
			//	std::cout << "D[" << k << "]= \n" << D[k] << std::endl;
			//}

			SparseMatrix Ha(m, 2*m), Hb(m, 2*m), S(m, 2*m);
			TripletList Ha_coefficients, Hb_coefficients;
			for (int mu = 0, pos = 0; mu < P; mu++) {
				int i = _alpha[mu], j = _beta[mu];

				SparseMatrix bar_M_i_mu_t(_m[mu], 2*m), bar_M_j_mu_s(_m[mu], 2*m);

				SparseMatrix D_i, D_j;
				int sum_m_mu_a_i = sum_m_mu_a[i], sum_m_mu_b_i = sum_m_mu_b[i];
				int sum_m_mu_a_j = sum_m_mu_a[j], sum_m_mu_b_j = sum_m_mu_b[j];				
				D_i = D[i];
				D_j = D[j];

				Matrix M_i_a_bar_t_i_mu;
				int bar_t_i_mu = bar_t[i][mu];
				M_i_a_bar_t_i_mu = M[i].block(M_a_pos[i][bar_t_i_mu], 0, _m[mu], sum_m_mu_a_i + sum_m_mu_b_i);

				Matrix M_j_b_bar_s_j_mu;
				int bar_s_j_mu = bar_s[j][mu];
				M_j_b_bar_s_j_mu = M[j].block(sum_m_mu_a_j + M_b_pos[j][bar_s_j_mu], 0, _m[mu], sum_m_mu_a_j + sum_m_mu_b_j);

				bar_M_i_mu_t = SparseMatrix(M_i_a_bar_t_i_mu.sparseView()) * D_i;
				bar_M_j_mu_s = SparseMatrix(M_j_b_bar_s_j_mu.sparseView()) * D_j;

				//std::cout << "bar_M_i_mu_t = \n" << bar_M_i_mu_t << std::endl;
				//std::cout << "bar_M_j_mu_s = \n" << bar_M_j_mu_s << std::endl;

				for (int outer_index = 0; outer_index < bar_M_i_mu_t.outerSize(); outer_index++) {
					for (SparseMatrix::InnerIterator it(bar_M_i_mu_t, outer_index); it; ++it) {
						Ha_coefficients.push_back(Triplet(pos + it.row(), it.col(), it.value()));
					}
				}
				for (int outer_index = 0; outer_index < bar_M_j_mu_s.outerSize(); outer_index++) {
					for (SparseMatrix::InnerIterator it(bar_M_j_mu_s, outer_index); it; ++it) {
						Hb_coefficients.push_back(Triplet(pos + it.row(), it.col(), it.value()));
					}
				}

				pos += _m[mu];
			}
			Ha.setFromTriplets(Ha_coefficients.begin(), Ha_coefficients.end());
			Hb.setFromTriplets(Hb_coefficients.begin(), Hb_coefficients.end());
			S = Ha - Hb;

			std::cout << "Ha : rows = " << Ha.rows() << " cols = " << Ha.cols() << " nonzeros = " << Ha.nonZeros() << std::endl;
			std::cout << "Hb : rows = " << Hb.rows() << " cols = " << Hb.cols() << " nonzeros = " << Hb.nonZeros() << std::endl;
			std::cout << "S : rows = " << S.rows() << " cols = " << S.cols() << " nonzeros = " << S.nonZeros() << std::endl;

			//std::cout << "Ha = \n" << Ha << std::endl;
			//std::cout << "Hb = \n" << Hb << std::endl;
			//std::cout << "S = \n" << S << std::endl;

			SparseMatrix Ja(m, N * (d+1)), Jb(m, N * (d+1)), T(m, N * (d+1));
			TripletList Ja_coefficients, Jb_coefficients;
			SparseMatrix Ka(m, m), Kb(m, m);
			TripletList Ka_coefficients, Kb_coefficients;
			for (int mu = 0, pos = 0; mu < P; mu++) {
				SparseMatrix Ja_mu(_m[mu], N * (d+1)), Jb_mu(_m[mu], N * (d+1));
				Ja_mu = SparseMatrix(X.block(pos, 0, _m[mu], d+1).sparseView()) * Ca[mu];
				Jb_mu = SparseMatrix(Y.block(pos, 0, _m[mu], d+1).sparseView()) * Cb[mu];

				for (int outer_index = 0; outer_index < Ja_mu.outerSize(); outer_index++) {
					for (SparseMatrix::InnerIterator it(Ja_mu, outer_index); it; ++it) {
						Ja_coefficients.push_back(Triplet(pos + it.row(),it.col(),it.value()));
					}
				}
				for (int outer_index = 0; outer_index < Jb_mu.outerSize(); outer_index++) {
					for (SparseMatrix::InnerIterator it(Jb_mu, outer_index); it; ++it) {
						Jb_coefficients.push_back(Triplet(pos + it.row(),it.col(),it.value()));
					}
				}

				for (int i = 0; i < _m[mu]; i++) {
					Ka_coefficients.push_back(Triplet(pos + i, pos + i, _kappa(_alpha[mu])));
					Kb_coefficients.push_back(Triplet(pos + i, pos + i, _kappa(_beta[mu])));
				}
				pos += _m[mu];
			}
			Ja.setFromTriplets(Ja_coefficients.begin(), Ja_coefficients.end());
			Jb.setFromTriplets(Jb_coefficients.begin(), Jb_coefficients.end());
			T = Ja - Jb;
			Ka.setFromTriplets(Ka_coefficients.begin(), Ka_coefficients.end());
			Kb.setFromTriplets(Kb_coefficients.begin(), Kb_coefficients.end());

			std::cout << "Ja : rows = " << Ja.rows() << " cols = " << Ja.cols() << " nonzeros = " << Ja.nonZeros() << std::endl;
			std::cout << "Jb : rows = " << Jb.rows() << " cols = " << Jb.cols() << " nonzeros = " << Jb.nonZeros() << std::endl;
			std::cout << "T : rows = " << T.rows() << " cols = " << T.cols() << " nonzeros = " << T.nonZeros() << std::endl;
			std::cout << "Ka : rows = " << Ka.rows() << " cols = " << Ka.cols() << " nonzeros = " << Ka.nonZeros() << std::endl;
			std::cout << "Kb : rows = " << Kb.rows() << " cols = " << Kb.cols() << " nonzeros = " << Kb.nonZeros() << std::endl;

			Matrix D_T_diag_lambda_M_D = Matrix::Zero(2*m, 2*m);
			for (int k = 0; k < N; k++) {
				SparseMatrix temp = _lambda[k] * D[k].transpose() * SparseMatrix(M[k].sparseView()) * D[k];
				D_T_diag_lambda_M_D += Matrix(temp);
			}
			std::cout << "D_T_diag_lambda_M_D : rows = " << D_T_diag_lambda_M_D.rows() << " cols = " << D_T_diag_lambda_M_D.cols() << " nonzeros = " << D_T_diag_lambda_M_D.nonZeros() << std::endl;

			SparseMatrix temp;
			Matrix P11(2*m, 2*m), P12(2*m, N * (d+1));
			temp = S.transpose() * S;
			P11 = Matrix(temp);
			std::cout << "P11 : rows = " << P11.rows() << " cols = " << P11.cols() << " nonzeros = " << P11.nonZeros() << std::endl;
			temp = Ha.transpose() * Ka * Ha;
			P11 += Matrix(temp);
			std::cout << "P11 : rows = " << P11.rows() << " cols = " << P11.cols() << " nonzeros = " << P11.nonZeros() << std::endl;
			temp = Hb.transpose() * Kb * Hb;
			P11 += Matrix(temp);
			std::cout << "P11 : rows = " << P11.rows() << " cols = " << P11.cols() << " nonzeros = " << P11.nonZeros() << std::endl;
			temp = S.transpose() * T;
			P12 = Matrix(temp);
			std::cout << "P12 : rows = " << P12.rows() << " cols = " << P12.cols() << " nonzeros = " << P12.nonZeros() << std::endl;
			temp = Ha.transpose() * Ka * Ja;
			P12 += Matrix(temp);
			std::cout << "P12 : rows = " << P12.rows() << " cols = " << P12.cols() << " nonzeros = " << P12.nonZeros() << std::endl;
			temp = Hb.transpose() * Kb * Jb;
			P12 += Matrix(temp);
			std::cout << "P12 : rows = " << P12.rows() << " cols = " << P12.cols() << " nonzeros = " << P12.nonZeros() << std::endl;

			Matrix Z1(2*m, d);
			temp = Ha.transpose() * Ka * SparseMatrix(_X.sparseView());
			Z1 = Matrix(temp);
			temp = Hb.transpose() * Kb * SparseMatrix(_Y.sparseView());
			Z1 += Matrix(temp);
			std::cout << "Z1 : rows = " << Z1.rows() << " cols = " << Z1.cols() << " nonzeros = " << Z1.nonZeros() << std::endl;

			Matrix P21(N * (d+1), 2*m), P22(N * (d+1), N * (d+1));
			temp = T.transpose() * S;
			P21 = Matrix(temp);
			temp = Ja.transpose() * Ka * Ha;
			P21 += Matrix(temp);
			temp = Jb.transpose() * Kb * Hb;
			P21 += Matrix(temp);
			std::cout << "P21 : rows = " << P21.rows() << " cols = " << P21.cols() << " nonzeros = " << P21.nonZeros() << std::endl;
			temp = T.transpose() * T;
			P22 = Matrix(temp);
			temp = Ja.transpose() * Ka * Ja;
			P22 += Matrix(temp);
			temp = Jb.transpose() * Kb * Jb;
			P22 += Matrix(temp);
			std::cout << "P22 : rows = " << P22.rows() << " cols = " << P22.cols() << " nonzeros = " << P22.nonZeros() << std::endl;

			Matrix Z2(N * (d+1), d);
			temp = Ja.transpose() * Ka * SparseMatrix(_X.sparseView());
			Z2 = Matrix(temp);
			temp = Jb.transpose() * Kb * SparseMatrix(_Y.sparseView());
			Z2 += Matrix(temp);
			std::cout << "Z2 : rows = " << Z2.rows() << " cols = " << Z2.cols() << " nonzeros = " << Z2.nonZeros() << std::endl;

			Matrix equationL = Matrix::Zero(2*m + N * (d+1), 2*m + N * (d+1));
			equationL << P11 + D_T_diag_lambda_M_D, P12,
						P21, P22;
			std::cout << "equationL : rows = " << equationL.rows() << " cols = " << equationL.cols() << " nonzeros = " << equationL.nonZeros() << std::endl;

			Matrix equationR = Matrix::Zero(2*m + N * (d+1), d);
			equationR << Z1, Z2;
			std::cout << "equationR : rows = " << equationR.rows() << " cols = " << equationR.cols() << " nonzeros = " << equationR.nonZeros() << std::endl;

			//std::fstream output_file;
			//output_file.open("equationL.txt", std::ios::out);
			//if (output_file) {
			//	output_file << equationL;
			//}
			//output_file.close();
			//output_file.open("equationR.txt", std::ios::out);
			//if (output_file) {
			//	output_file << equationR;
			//}
			//output_file.close();

			Matrix solution(2*m + N * (d+1), d);

			//solution = (equationL.transpose() * equationL).inverse() * (equationL.transpose() * equationR);

			//Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
			//std::cout << "Condition Number : " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			//std::cout << "Singular Vector : \n" << svd.singularValues() << std::endl;
			//std::cout << "Rows = " << equationL.rows() << " Columns = " << equationL.cols() << std::endl;
			//std::cout << "before reset threshold, rank = " << svd.rank() << std::endl;
			//svd.setThreshold(Eigen::NumTraits<Scalar>::epsilon());
			//std::cout << "after reset threshold, rank = " << svd.rank() << std::endl;
			//solution = svd.solve(equationR);

			Eigen::LDLT<Matrix> ldlt(equationL);
			solution = ldlt.solve(equationR);
			Scalar ldlt_residual = (equationL * solution - equationR).squaredNorm();
			std::cout << "ldlt_residual = "<< ldlt_residual << std::endl;

			if (ldlt_residual > 1.0) {
				Eigen::FullPivLU<Matrix> lu(equationL);
				solution = lu.solve(equationR);
				Scalar lu_residual = (equationL * solution - equationR).squaredNorm();
				if (lu_residual > 1.0) {
					std::cout << "lu_residual =" << lu_residual << std::endl;
					return;
				}
			}

			//std::cout << "Ax - b = \n" << equationL * solution - equationR << std::endl; 

			Matrix A(2*m, d), B(N * (d+1), d);
			A = solution.block(0, 0, 2*m, d);
			B = solution.block(2*m, 0, N * (d+1), d);

			Scalar tps_energy1 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy1 += ( T.block(pos, 0, _m[mu], N * (d+1)) * B + S.block(pos, 0, _m[mu], 2*m) * A ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy1 = " << tps_energy1 << std::endl;
			Scalar tps_energy2 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy2 += _kappa(_alpha[mu]) * 
					( Ja.block(pos, 0, _m[mu], N * (d+1) ) * B + Ha.block(pos, 0, _m[mu], 2*m) * A - _X.block(pos, 0, _m[mu], d) ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy2 = " << tps_energy2 << std::endl;
			Scalar tps_energy3 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy3 += _kappa(_beta[mu]) *
					( Jb.block(pos, 0, _m[mu], N * (d+1) ) * B + Hb.block(pos, 0, _m[mu], 2*m) * A - _Y.block(pos, 0, _m[mu], d) ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy3 = " << tps_energy3 << std::endl;
			Scalar tps_energy4 = (A.transpose() * D_T_diag_lambda_M_D * A).trace();
			std::cout << "tps_energy4 = " << tps_energy4 << std::endl;

			fs_.resize(N);
			for (int k = 0; k < N; k++){
				Matrix Ak(sum_m_mu_a[k] + sum_m_mu_b[k], d), Xk(sum_m_mu_a[k] + sum_m_mu_b[k], d);
				for (int ti = 0, pos = 0; ti < mu_alpha[k].size(); ti++){
					Ak.block(pos, 0, _m[mu_alpha[k][ti]], d) = 
						A.block(Apos[mu_alpha[k][ti]], 0, _m[mu_alpha[k][ti]], d);
					Xk.block(pos, 0, _m[mu_alpha[k][ti]], d) = 
						_X.block(Xpos[mu_alpha[k][ti]], 0, _m[mu_alpha[k][ti]], d);
					pos += _m[mu_alpha[k][ti]];
				}
				for (int sj = 0, pos = 0; sj < mu_beta[k].size(); sj++) {
					Ak.block(sum_m_mu_a[k] + pos, 0, _m[mu_beta[k][sj]], d) = 
						A.block(m + Apos[mu_beta[k][sj]], 0, _m[mu_beta[k][sj]], d);
					Xk.block(sum_m_mu_a[k] + pos, 0, _m[mu_beta[k][sj]], d) =
						_Y.block(Ypos[mu_beta[k][sj]], 0, _m[mu_beta[k][sj]], d);
					pos += _m[mu_beta[k][sj]];
				}
				fs_[k].setX(Xk);
				fs_[k].setA(Ak);
			}
			for (int k = 0, kpos = 0; k < N; k++){
				Matrix Bk = B.block(kpos, 0, d+1, d);
				fs_[k].setB(Bk);
				kpos += d+1;
			}
	}

	void MultiTPS::solve_old_(const Matrix &_X, 
		const Matrix &_Y,
		const std::vector<int> &_m,
		const std::vector<int> &_alpha,
		const std::vector<int> &_beta,
		const Vector &_kappa,
		const Vector &_lambda)	{
			int P = _m.size(), N = _lambda.size();
			int d = _X.cols();
			int m = 0; // or m = _X.rows();
			for (int i = 0; i < P; i++) {
				m += _m[i];
			}

			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
			Matrix Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			Matrix Maa(m, m), Mab(m, m), Mba(m, m), Mbb(m, m);
			for (int i = 0, ipos = 0; i < P; i++) {
				for (int j = 0, jpos = 0; j < P; j++) {
					Maa.block(ipos, jpos, _m[i], _m[j]) = greenFunc(_X.block(ipos, 0, _m[i], d), _X.block(jpos, 0, _m[j], d));
					Mab.block(ipos, jpos, _m[i], _m[j]) = greenFunc(_X.block(ipos, 0, _m[i], d), _Y.block(jpos, 0, _m[j], d));
					Mba.block(ipos, jpos, _m[i], _m[j]) = greenFunc(_Y.block(ipos, 0, _m[i], d), _X.block(jpos, 0, _m[j], d));
					Mbb.block(ipos, jpos, _m[i], _m[j]) = greenFunc(_Y.block(ipos, 0, _m[i], d), _Y.block(jpos, 0, _m[j], d));
					jpos += _m[j];
				}
				ipos += _m[i];
			}
			Matrix Ma(m, 2*m), Mb(m, 2*m);
			Ma << Maa, Mab;
			Mb << Mba, Mbb;

			Matrix Ca = Matrix::Zero(P*(d+1), N*(d+1)), Cb = Matrix::Zero(P*(d+1), N*(d+1));
			for (int i = 0, pos = 0; i < P; i++) {
				int ai = _alpha[i], bi = _beta[i];
				Ca.block(pos, ai * (d+1), d+1, d+1) = Matrix::Identity(d+1, d+1);
				Cb.block(pos, bi * (d+1), d+1, d+1) = Matrix::Identity(d+1, d+1);
				pos += d+1;
			}
			//std::cout << "Ca = \n" << Ca << std::endl;
			//std::cout << "Cb = \n" << Cb << std::endl;

			Matrix Da = Matrix::Zero(P * 2*m, N * 2*m), Db = Matrix::Zero(P * 2*m, N * 2*m);
			for (int i = 0, pos = 0; i < P; i++) {
				int ai = _alpha[i], bi = _beta[i];
				Da.block(pos, ai * 2*m, 2*m, 2*m) = Matrix::Identity(2*m, 2*m);
				Db.block(pos, bi * 2*m, 2*m, 2*m) = Matrix::Identity(2*m, 2*m);
				pos += 2*m;
			}
			//std::cout << "Da = \n" << Da << std::endl;
			//std::cout << "Db = \n" << Db << std::endl;

			Matrix Ha(m, N * 2*m), Hb(m, N * 2*m), S(m, N * 2*m);
			Matrix Ja(m, N * (d+1)), Jb(m, N * (d+1)), T(m, N * (d+1));
			Matrix Ka = Matrix::Zero(m, m), Kb = Matrix::Zero(m, m);
			for (int i = 0, pos = 0; i < P; i++) {
				Ha.block(pos, 0, _m[i], N * 2*m) = Ma.block(pos, 0, _m[i], 2*m) * Da.block(i * 2*m, 0, 2*m, N * 2*m);
				Hb.block(pos, 0, _m[i], N * 2*m) = Mb.block(pos, 0, _m[i], 2*m) * Db.block(i * 2*m, 0, 2*m, N * 2*m);
				Ja.block(pos, 0, _m[i], N * (d+1)) = X.block(pos, 0, _m[i], d+1) * Ca.block(i * (d+1), 0, d+1, N * (d+1));
				Jb.block(pos, 0, _m[i], N * (d+1)) = Y.block(pos, 0, _m[i], d+1) * Cb.block(i * (d+1), 0, d+1, N * (d+1));
				Ka.block(pos, pos, _m[i], _m[i]) = Matrix::Identity(_m[i], _m[i]) * _kappa(_alpha[i]);
				Kb.block(pos, pos, _m[i], _m[i]) = Matrix::Identity(_m[i], _m[i]) * _kappa(_beta[i]);
				pos += _m[i];
			}
			S = Ha - Hb;
			T = Ja - Jb;

			Matrix M(2*m, 2*m);
			M << Ma, Mb;
			Matrix diag_lambda_M = Matrix::Zero(N * 2*m, N * 2*m);
			for (int k = 0, pos = 0; k < N; k++) {
				diag_lambda_M.block(pos, pos, 2*m, 2*m) = M * _lambda(k);
				pos += 2*m;
			}

			Matrix P11(N * 2*m, N * 2*m), P12(N * 2*m, N * (d+1));
			P11 = S.transpose() * S + Ha.transpose() * Ka * Ha + Hb.transpose() * Kb * Hb;
			P12 = S.transpose() * T + Ha.transpose() * Ka * Ja + Hb.transpose() * Kb * Jb;
			Matrix Z1(N * 2*m, d);
			Z1 = Ha.transpose() * Ka * _X + Hb.transpose() * Kb * _Y;
			Matrix P21(N * (d+1), N * 2*m), P22(N * (d+1), N * (d+1));
			P21 = T.transpose() * S + Ja.transpose() * Ka * Ha + Jb.transpose() * Kb * Hb;
			P22 = T.transpose() * T + Ja.transpose() * Ka * Ja + Jb.transpose() * Kb * Jb;
			Matrix Z2(N * (d+1), d);
			Z2 = Ja.transpose() * Ka * _X + Jb.transpose() * Kb * _Y;

			Matrix U = Matrix::Zero(N * 2*m, N * 2*m);
			for (int k = 0, kpos = 0; k < N; k++){
				Matrix Uka = Matrix::Zero(m, m), Ukb = Matrix::Zero(m, m);
				for (int i = 0, ipos = 0; i < P; i++){
					if(_alpha[i] != k) Uka.block(ipos, ipos, _m[i], _m[i]) = Matrix::Identity(_m[i], _m[i]);
					if(_beta[i] != k) Ukb.block(ipos, ipos, _m[i], _m[i]) = Matrix::Identity(_m[i], _m[i]);
					ipos += _m[i];
				}
				Matrix Uk = Matrix::Zero(2*m, 2*m);
				Uk.block(0, 0, m, m) = Uka;
				Uk.block(m, m, m, m) = Ukb;
				U.block(kpos, kpos, 2*m, 2*m) = Uk;
				kpos += 2*m;
			}
			//std::cout << "U = \n" << U << std::endl;

			Matrix equationL = Matrix::Zero(N * 2*m + N * (d+1) + N * 2*m, N * 2*m + N * (d+1));
			equationL << P11 + diag_lambda_M, P12,
				P21, P22,
				U, Matrix::Zero(N * 2*m, N *(d+1));
			Matrix equationR = Matrix::Zero(N * 2*m + N * (d+1) + N * 2*m, d);
			equationR << Z1, Z2, Matrix::Zero(N * 2*m, d);

			Matrix solution(N * 2*m + N * (d+1), d);

			//solution = (equationL.transpose() * equationL).inverse() * (equationL.transpose() * equationR);

			Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
			std::cout << "Condition Number : " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			std::cout << "Singular Vector : \n" << svd.singularValues() << std::endl;
			std::cout << "Rows = " << equationL.rows() << " Columns = " << equationL.cols() << std::endl;
			std::cout << "before reset threshold, rank = " << svd.rank() << std::endl;
			svd.setThreshold(Eigen::NumTraits<Scalar>::epsilon());
			std::cout << "after reset threshold, rank = " << svd.rank() << std::endl;
			solution = svd.solve(equationR);

			//std::cout << "Ax - b = \n" << equationL * solution - equationR << std::endl; 

			Matrix A(N * 2*m, d), B(N * (d+1), d);
			A = solution.block(0, 0, N * 2*m, d);
			B = solution.block(N * 2*m, 0, N * (d+1), d);

			Scalar tps_energy1 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy1 += ( T.block(pos, 0, _m[mu], N * (d+1)) * B + S.block(pos, 0, _m[mu], 2*m*N) * A ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy1 = " << tps_energy1 << std::endl;
			Scalar tps_energy2 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy2 += _kappa(_alpha[mu]) * 
					( Ja.block(pos, 0, _m[mu], N * (d+1) ) * B + Ha.block(pos, 0, _m[mu], 2*m*N) * A - _X.block(pos, 0, _m[mu], d) ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy2 = " << tps_energy2 << std::endl;
			Scalar tps_energy3 = 0;
			for (int mu = 0, pos = 0; mu < P; mu++){
				tps_energy3 += _kappa(_beta[mu]) *
					( Jb.block(pos, 0, _m[mu], N * (d+1) ) * B + Hb.block(pos, 0, _m[mu], 2*m*N) * A - _Y.block(pos, 0, _m[mu], d) ).squaredNorm();
				pos += _m[mu];
			}
			std::cout << "tps_energy3 = " << tps_energy3 << std::endl;
			Scalar tps_energy4 = (A.transpose() * diag_lambda_M * A).trace();
			std::cout << "tps_energy4 = " << tps_energy4 << std::endl;

			fs_.resize(N);
			for (int k = 0, kpos = 0; k < N; k++){
				Matrix Ak = A.block(kpos, 0, 2*m, d);
				std::vector<Matrix> Xks, Aks;
				int num_pt = 0;
				//for (int i = 0, pos = 0; i < P; i++) {
				//	if (_alpha[i] == k) {
				//		Xks.push_back(_X.block(pos, 0, _m[i], d));
				//		Aks.push_back(Ak.block(pos, 0, _m[i], d));
				//		num_pt += _m[i];
				//	}else if (_beta[i] == k) {
				//		Xks.push_back(_Y.block(pos, 0, _m[i], d));
				//		Aks.push_back(Ak.block(m + pos, 0, _m[i], d));
				//		num_pt += _m[i];
				//	}
				//	pos += _m[i];
				//}
				for (int i = 0, pos = 0; i < P; i++) {
					if (_alpha[i] == k) {
						Xks.push_back(_X.block(pos, 0, _m[i], d));
						Aks.push_back(Ak.block(pos, 0, _m[i], d));
						num_pt += _m[i];
					}
					pos += _m[i];
				}
				for (int i = 0, pos = 0; i < P; i++){
					if (_beta[i] == k) {
						Xks.push_back(_Y.block(pos, 0, _m[i], d));
						Aks.push_back(Ak.block(m + pos, 0, _m[i], d));
						num_pt += _m[i];
					}
					pos += _m[i];
				}
				Matrix Xk(num_pt, d), Ak_new(num_pt, d);
				for (int i = 0, pos = 0; i < Xks.size(); i++){
					Xk.block(pos, 0, Xks[i].rows(), d) = Xks[i];
					Ak_new.block(pos, 0, Xks[i].rows(), d) = Aks[i];
					pos += Xks[i].rows();
				}
				fs_[k].setX(Xk);
				fs_[k].setA(Ak_new);
				kpos += 2*m;
			}
			for (int k = 0, kpos = 0; k < N; k++){
				Matrix Bk = B.block(kpos, 0, d+1, d);
				fs_[k].setB(Bk);
				kpos += d+1;
			}

			//for (int k = 0, kpos = 0; k < N; k++) {
			//	Matrix Xk(2*m, d);
			//	Xk << _X, _Y;
			//	fs_[k].setX(Xk);
			//	Matrix Ak = A.block(kpos, 0, 2*m, d);
			//	fs_[k].setA(Ak);
			//	kpos += 2*m;
			//}
	}

	void MultiTPS::solve_iterative_(const Matrix &_X, 
		const Matrix &_Y,
		const std::vector<int> &_m,
		const std::vector<int> &_alpha,
		const std::vector<int> &_beta,
		const Vector &_kappa,
		const Vector &_lambda,
		const int &_max_iter_num,
		const Scalar &_iter_rate) {
			int P = _m.size(), N = _lambda.size();
			int d = _X.cols();
			int m = 0; // or m = _X.rows();
			for (int i = 0; i < P; i++) {
				m += _m[i];
			}

			std::cout << "N = " << N << std::endl;
			std::cout << "P = " << P << std::endl;
			std::cout << "m = " << m << std::endl;
			for (int i = 0; i < P; i++) {
				std::cout << "( " << _alpha[i] << ", " << _beta[i] << ") " << _m[i] << std::endl;
			}
			std::cout << "_X : rows = " << _X.rows() << " cols = " << _X.cols() << std::endl;
			std::cout << "_Y : rows = " << _Y.rows() << " cols = " << _Y.cols() << std::endl;

			//Initialization
			std::vector<int> An(N, 0);
			for (int i = 0, j = 0; i < P; j+= _m[i], i++) {
				An[_alpha[i]] += _m[i];
				An[_beta[i]] += _m[i];
			}
			fs_.resize(N);
			Matrix current_X(_X), current_Y(_Y);

			Scalar min_lambda = _lambda[0], min_kappa = _kappa[0];
			for (int i = 1; i < N; i++) {
				min_lambda = min_lambda < _lambda[i] ? min_lambda : _lambda[i];
				min_kappa = min_kappa < _kappa[i] ? min_kappa : _kappa[i];
			}
			Scalar iter_rate = _iter_rate > 0 ? _iter_rate : pow(min_kappa, 1.0/_max_iter_num);
			Scalar current_rate = 1.0;
			for(int iter_num = 0;  iter_num < _max_iter_num; iter_num++) {
				for (int i = 0; i < N; i++) {
					Matrix X(An[i], d), Y(An[i], d);
					for (int j = 0, k1 = 0, k2 = 0; j < P; k1 += _m[j], j++) {
						if (_alpha[j] == i) {
							X.block(k2, 0, _m[j], d) = _X.block(k1, 0, _m[j], d);
							if(iter_num == 0) Y.block(k2, 0, _m[j], d) = _Y.block(k1, 0, _m[j], d);
							else Y.block(k2, 0, _m[j], d) = current_Y.block(k1, 0, _m[j], d);
							k2 += _m[j];
						} else if (_beta[j] == i) {
							X.block(k2, 0, _m[j], d) = _Y.block(k1, 0, _m[j], d);
							if(iter_num == 0) Y.block(k2, 0, _m[j], d) = _X.block(k1, 0, _m[j], d);
							else Y.block(k2, 0, _m[j], d) = current_X.block(k1, 0, _m[j], d);
							k2 += _m[j];
						}
					}
					Scalar lambda_i;
					if(An[i] > 100) lambda_i = _lambda[i];
					else lambda_i = _lambda[i] * ((1.0 / min_kappa * current_rate > 1.0) ? (1.0 / min_kappa * current_rate) : 1.0); // _lambda[i] * (min_lambda / min_kappa) * (1.0 / min_lambda * current_rate)
					Scalar kappa_i =  _kappa[i] * ((1.0 / min_kappa * current_rate > 1.0) ? (1.0 / min_kappa * current_rate) : 1.0);
					std::cout << "view " << i << " : kappa = " << kappa_i << ", lambda = " << lambda_i << std::endl;
					if (iter_num == 0 && d == 3) {
						std::stringstream ss_viewX;
						ss_viewX << "viewX_" << i << ".xyz";
						std::fstream fs_viewX(ss_viewX.str().c_str(), std::ios::out);
						if(fs_viewX) {
							for(int p = 0; p < X.rows(); p++) {
								for (int q = 0; q < X.cols(); q++) {
									fs_viewX << X(p, q) << " ";
								}
								fs_viewX << std::endl;
							}
						}
						fs_viewX.close();
						std::stringstream ss_viewY;
						ss_viewY << "viewY_" << i << ".xyz";
						std::fstream fs_viewY(ss_viewY.str().c_str(), std::ios::out);
						if(fs_viewY) {
							for(int p = 0; p < Y.rows(); p++) {
								for (int q = 0; q < Y.cols(); q++) {
									fs_viewY << Y(p, q) << " ";
								}
								fs_viewY << std::endl;
							}
						}
						fs_viewY.close();
					}
					fs_[i] = TPSFunction(X, Y, lambda_i, kappa_i);
					for (int j = 0, k1 = 0; j < P; k1 += _m[j], j++) {
						if (_alpha[j] == i) {
							current_X.block(k1, 0, _m[j], d) = fs_[i].evaluate(_X.block(k1, 0, _m[j], d));
						} else if (_beta[j] == i) {
							current_Y.block(k1, 0, _m[j], d) = fs_[i].evaluate(_Y.block(k1, 0, _m[j], d));
						}
					}
					if (iter_num == _max_iter_num - 1 && d == 3) {

						Matrix current_viewX(An[i], d), current_viewY(An[i], d);
						for (int j = 0, k1 = 0, k2 = 0; j < P; k1 += _m[j], j++) {
							if (_alpha[j] == i) {
								current_viewX.block(k2, 0, _m[j], d) = current_X.block(k1, 0, _m[j], d);
								current_viewY.block(k2, 0, _m[j], d) = current_Y.block(k1, 0, _m[j], d);
								k2 += _m[j];
							} else if (_beta[j] == i) {
								current_viewX.block(k2, 0, _m[j], d) = current_Y.block(k1, 0, _m[j], d);
								current_viewY.block(k2, 0, _m[j], d) = current_X.block(k1, 0, _m[j], d);
								k2 += _m[j];
							}
						}
						std::stringstream ss_current_viewX;
						ss_current_viewX << "current_viewX_" << i << ".xyz";
						std::fstream fs_current_viewX(ss_current_viewX.str().c_str(), std::ios::out);
						if(fs_current_viewX) {
							for(int p = 0; p < current_viewX.rows(); p++) {
								for (int q = 0; q < current_viewX.cols(); q++) {
									fs_current_viewX << current_viewX(p, q) << " ";
								}
								fs_current_viewX << std::endl;
							}
						}
						fs_current_viewX.close();
						std::stringstream ss_current_viewY;
						ss_current_viewY << "current_viewY_" << i << ".xyz";
						std::fstream fs_current_viewY(ss_current_viewY.str().c_str(), std::ios::out);
						if(fs_current_viewY) {
							for(int p = 0; p < current_viewY.rows(); p++) {
								for (int q = 0; q < current_viewY.cols(); q++) {
									fs_current_viewY << current_viewY(p, q) << " ";
								}
								fs_current_viewY << std::endl;
							}
						}
						fs_current_viewY.close();
					}
				}
				current_rate *= iter_rate;
			}
	}

	ApproxiMultiTPS::ApproxiMultiTPS() : MultiTPS() {}

	ApproxiMultiTPS::ApproxiMultiTPS(const Matrix &_X, 
		const Matrix &_Y,
		const std::vector<int> &_m,
		const std::vector<int> &_alpha,
		const std::vector<int> &_beta,
		const Vector &_kappa,
		const Vector &_lambda,
		const std::vector<int> &_na,
		const std::vector<int> &_nb,
		const int &_max_iter_num,
		const Scalar &_iter_rate) {
			int P = _m.size(), N = _lambda.size();
			int d = _X.cols();
			int m = 0; // or m = _X.rows();
			for (int i = 0; i < P; i++) {
				m += _m[i];
			}

			std::cout << "N = " << N << std::endl;
			std::cout << "P = " << P << std::endl;
			std::cout << "m = " << m << std::endl;
			for (int i = 0; i < P; i++) {
				std::cout << "( " << _alpha[i] << ", " << _beta[i] << ") " << _m[i] << std::endl;
			}
			std::cout << "_X : rows = " << _X.rows() << " cols = " << _X.cols() << std::endl;
			std::cout << "_Y : rows = " << _Y.rows() << " cols = " << _Y.cols() << std::endl;

			//Initialization
			std::vector<int> An(N, 0), Xn(N, 0);
			for (int i = 0; i < P; i++) {
				An[_alpha[i]] += _na[i];
				An[_beta[i]] += _nb[i];
				Xn[_alpha[i]] += _m[i];
				Xn[_beta[i]] += _m[i];
			}
			fs_.resize(N);
			Matrix current_X(_X), current_Y(_Y);

			Scalar min_lambda = _lambda[0], min_kappa = _kappa[0];
			for (int i = 1; i < N; i++) {
				min_lambda = min_lambda < _lambda[i] ? min_lambda : _lambda[i];
				min_kappa = min_kappa < _kappa[i] ? min_kappa : _kappa[i];
			}
			Scalar iter_rate = _iter_rate > 0 ? _iter_rate : pow(min_kappa, 1.0/_max_iter_num);
			Scalar current_rate = 1.0;
			for(int iter_num = 0;  iter_num < _max_iter_num; iter_num++) {
				for (int i = 0; i < N; i++) {
					Matrix X(Xn[i], d), Y(Xn[i], d);
					for (int j = 0, k1 = 0, k2 = 0; j < P; k1 += _m[j], j++) {
						if (_alpha[j] == i) {
							X.block(k2, 0, _na[j], d) = _X.block(k1, 0, _na[j], d);
							if(iter_num == 0) Y.block(k2, 0, _na[j], d) = _Y.block(k1, 0, _na[j], d);
							else Y.block(k2, 0, _na[j], d) = current_Y.block(k1, 0, _na[j], d);
							k2 += _na[j];
						} else if (_beta[j] == i) {
							X.block(k2, 0, _nb[j], d) = _Y.block(k1 + _m[j] - _nb[j], 0, _nb[j], d);
							if(iter_num == 0) Y.block(k2, 0, _nb[j], d) = _X.block(k1 + _m[j] - _nb[j], 0, _nb[j], d);
							else Y.block(k2, 0, _nb[j], d) = current_X.block(k1 + _m[j] - _nb[j], 0, _nb[j], d);
							k2 += _nb[j];
						}
					}
					for (int j = 0, k1 = 0, k2 = An[i]; j < P; k1 += _m[j], j++) {
						if (_alpha[j] == i) {
							X.block(k2, 0, _m[j] - _na[j], d) = _X.block(k1 + _na[j], 0, _m[j] - _na[j], d);
							if(iter_num == 0) Y.block(k2, 0, _m[j] - _na[j], d) = _Y.block(k1 + _na[j], 0, _m[j] - _na[j], d);
							else Y.block(k2, 0, _m[j] - _na[j], d) = current_Y.block(k1 + _na[j], 0, _m[j] - _na[j], d);
							k2 += _m[j] - _na[j];
						} else if (_beta[j] == i) {
							X.block(k2, 0, _m[j] - _nb[j], d) = _Y.block(k1, 0, _m[j] - _nb[j], d);
							if(iter_num == 0) Y.block(k2, 0, _m[j] - _nb[j], d) = _X.block(k1, 0, _m[j] - _nb[j], d);
							else Y.block(k2, 0, _m[j] - _nb[j], d) = current_X.block(k1, 0, _m[j] - _nb[j], d);
							k2 += _m[j] - _nb[j];
						}
					}
					Scalar lambda_i;
					if(An[i] > 100) lambda_i = _lambda[i];
					else lambda_i = _lambda[i] * ((1.0 / min_kappa * current_rate > 1.0) ? (1.0 / min_kappa * current_rate) : 1.0); // _lambda[i] * (min_lambda / min_kappa) * (1.0 / min_lambda * current_rate)
					Scalar kappa_i =  _kappa[i] * ((1.0 / min_kappa * current_rate > 1.0) ? (1.0 / min_kappa * current_rate) : 1.0);
					std::cout << "view " << i << " : kappa = " << kappa_i << ", lambda = " << lambda_i << std::endl;
					if (iter_num == 0 && d == 3) {
						std::stringstream ss_viewX;
						ss_viewX << "viewX_" << i << ".xyz";
						std::fstream fs_viewX(ss_viewX.str().c_str(), std::ios::out);
						if(fs_viewX) {
							for(int p = 0; p < X.rows(); p++) {
								for (int q = 0; q < X.cols(); q++) {
									fs_viewX << X(p, q) << " ";
								}
								fs_viewX << std::endl;
							}
						}
						fs_viewX.close();
						std::stringstream ss_viewY;
						ss_viewY << "viewY_" << i << ".xyz";
						std::fstream fs_viewY(ss_viewY.str().c_str(), std::ios::out);
						if(fs_viewY) {
							for(int p = 0; p < Y.rows(); p++) {
								for (int q = 0; q < Y.cols(); q++) {
									fs_viewY << Y(p, q) << " ";
								}
								fs_viewY << std::endl;
							}
						}
						fs_viewY.close();
					}
					//fs_[i] = TPSFunction(X, Y, lambda_i, kappa_i);
					fs_[i] = ApproxiTPSFunction(X, Y, lambda_i, kappa_i, An[i]);
					for (int j = 0, k1 = 0; j < P; k1 += _m[j], j++) {
						if (_alpha[j] == i) {
							current_X.block(k1, 0, _m[j], d) = fs_[i].evaluate(_X.block(k1, 0, _m[j], d));
						} else if (_beta[j] == i) {
							current_Y.block(k1, 0, _m[j], d) = fs_[i].evaluate(_Y.block(k1, 0, _m[j], d));
						}
					}
					if (iter_num == _max_iter_num - 1 && d == 3) {

						Matrix current_viewX(Xn[i], d), current_viewY(Xn[i], d);
						for (int j = 0, k1 = 0, k2 = 0; j < P; k1 += _m[j], j++) {
							if (_alpha[j] == i) {
								current_viewX.block(k2, 0, _m[j], d) = current_X.block(k1, 0, _m[j], d);
								current_viewY.block(k2, 0, _m[j], d) = current_Y.block(k1, 0, _m[j], d);
								k2 += _m[j];
							} else if (_beta[j] == i) {
								current_viewX.block(k2, 0, _m[j], d) = current_Y.block(k1, 0, _m[j], d);
								current_viewY.block(k2, 0, _m[j], d) = current_X.block(k1, 0, _m[j], d);
								k2 += _m[j];
							}
						}
						std::stringstream ss_current_viewX;
						ss_current_viewX << "current_viewX_" << i << ".xyz";
						std::fstream fs_current_viewX(ss_current_viewX.str().c_str(), std::ios::out);
						if(fs_current_viewX) {
							for(int p = 0; p < current_viewX.rows(); p++) {
								for (int q = 0; q < current_viewX.cols(); q++) {
									fs_current_viewX << current_viewX(p, q) << " ";
								}
								fs_current_viewX << std::endl;
							}
						}
						fs_current_viewX.close();
						std::stringstream ss_current_viewY;
						ss_current_viewY << "current_viewY_" << i << ".xyz";
						std::fstream fs_current_viewY(ss_current_viewY.str().c_str(), std::ios::out);
						if(fs_current_viewY) {
							for(int p = 0; p < current_viewY.rows(); p++) {
								for (int q = 0; q < current_viewY.cols(); q++) {
									fs_current_viewY << current_viewY(p, q) << " ";
								}
								fs_current_viewY << std::endl;
							}
						}
						fs_current_viewY.close();
					}
				}
				current_rate *= iter_rate;
			}
	}
};