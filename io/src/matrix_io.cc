#include <fstream>
#include <string>

#include <GMNR/io/matrix_io.h>

namespace gmnr {
	void write_matrix(const std::string &file_name, const Matrix& mat) {
		std::fstream fs_mat(file_name.c_str(), std::ios::out);
		if(fs_mat) {
			for(int p = 0; p < mat.rows(); p++) {
				for (int q = 0; q < mat.cols(); q++) {
					fs_mat << mat(p, q) << " ";
				}
				fs_mat << std::endl;
			}
		}
		fs_mat.close();
	}
}
