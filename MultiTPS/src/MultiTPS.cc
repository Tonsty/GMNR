#include <iostream>
#include <unordered_map>
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
		const Vector &_lambda) {
			//solve_old_(_X, _Y, _m, _alpha, _beta, _kappa, _lambda);

			//std::vector<TPSFunction> old_fs = fs_; 
			//for (int i = 0; i < old_fs.size(); i++){
			//	std::cout << "old_fs[" << i << "]=\n" << old_fs[i] << std::endl;
			//}

			solve_new_(_X, _Y, _m, _alpha, _beta, _kappa, _lambda);

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

	template<typename T>
	std::ostream& operator << (std::ostream &os, const std::unordered_map<T,T> &map){
		for (std::unordered_map<T,T>::const_iterator iter = map.begin(); iter != map.end(); iter++){
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

				Matrix D_i_a, D_i_b;
				int sum_m_mu_a_i = sum_m_mu_a[i], sum_m_mu_b_i = sum_m_mu_b[i];
				D_i_a = D[i].block(0, 0, sum_m_mu_a_i, m);
				D_i_b = D[i].block(sum_m_mu_a_i, m, sum_m_mu_b_i, m);
				Matrix D_j_a, D_j_b;
				int sum_m_mu_a_j = sum_m_mu_a[j], sum_m_mu_b_j = sum_m_mu_b[j];
				D_j_a = D[j].block(0, 0, sum_m_mu_a_j, m);
				D_j_b = D[j].block(sum_m_mu_a_j, m, sum_m_mu_b_j, m);

				Matrix M_i_aa_bar_t_i_mu, M_i_ab_bar_t_i_mu;
				int bar_t_i_mu = bar_t[i][mu];
				M_i_aa_bar_t_i_mu = M[i].block( M_a_pos[i][bar_t_i_mu], 0, _m[mu], sum_m_mu_a_i);
				M_i_ab_bar_t_i_mu = M[i].block( M_a_pos[i][bar_t_i_mu], sum_m_mu_a_i, _m[mu], sum_m_mu_b_i);

				Matrix M_j_ba_bar_s_j_mu, M_j_bb_bar_s_j_mu;
				int bar_s_j_mu = bar_s[j][mu];
				M_j_ba_bar_s_j_mu = M[j].block(sum_m_mu_a_j + M_b_pos[j][bar_s_j_mu], 0, _m[mu], sum_m_mu_a_j);
				M_j_bb_bar_s_j_mu = M[j].block(sum_m_mu_a_j + M_b_pos[j][bar_s_j_mu], sum_m_mu_a_j, _m[mu], sum_m_mu_b_j);

				if (sum_m_mu_b_i == 0) {
					bar_M_i_mu_t << M_i_aa_bar_t_i_mu * D_i_a, Matrix::Zero(_m[mu], m);
				}else{
					bar_M_i_mu_t << M_i_aa_bar_t_i_mu * D_i_a, M_i_ab_bar_t_i_mu * D_i_b;
				}
				if (sum_m_mu_a_j == 0){
					bar_M_j_mu_s << Matrix::Zero(_m[mu], m), M_j_bb_bar_s_j_mu * D_j_b;
				}else{
					bar_M_j_mu_s << M_j_ba_bar_s_j_mu * D_j_a, M_j_bb_bar_s_j_mu * D_j_b;
				}
				
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

			Matrix diag_lambda_M = Matrix::Zero(2*m, 2*m);
			for (int k = 0, pos = 0; k < N; k++) {
				int sum_m_mu_a_k = sum_m_mu_a[k], sum_m_mu_b_k = sum_m_mu_b[k];
				diag_lambda_M.block(pos, pos, sum_m_mu_a_k + sum_m_mu_b_k, sum_m_mu_a_k + sum_m_mu_b_k) = M[k] * _lambda(k);
				pos += sum_m_mu_a_k + sum_m_mu_b_k;
			}
			Matrix D_matrice(2*m, 2*m);
			for (int k = 0, pos = 0; k < N; k++){
				D_matrice.block(pos, 0, sum_m_mu_a[k] + sum_m_mu_b[k], 2*m) = D[k];
				pos += sum_m_mu_a[k] + sum_m_mu_b[k];
			}
			Matrix D_T_diag_lambda_M_D;
			D_T_diag_lambda_M_D = D_matrice.transpose() * diag_lambda_M * D_matrice;
			//Matrix D_T_diag_lambda_M_D = Matrix::Zero(2*m, 2*m);
			//for (int k = 0; k < N; k++) {
			//	D_T_diag_lambda_M_D += _lambda[k] * D[k].transpose() * M[k] * D[k];
			//}

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
			Scalar tps_energy4 = (A.transpose() * D_matrice.transpose() * diag_lambda_M * D_matrice * A).trace();
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

};