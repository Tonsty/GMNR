#ifndef TPS_ICP_H
#define TPS_ICP_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/search/kdtree.h>

#include <GMNR/ThinPlateSplines.h>

namespace gmnr {

	struct Data_Pack {

		typedef boost::shared_ptr<Data_Pack> Ptr;

		Data_Pack(pcl::PointCloud<pcl::PointXYZ>::Ptr pts,
			pcl::search::KdTree<pcl::PointXYZ>::Ptr kd,
			pcl::PointCloud<pcl::Normal>::Ptr norms,
			pcl::PointCloud<pcl::Boundary>::Ptr bd) : pts_(pts), kd_(kd), norms_(norms), bd_(bd) {}

		pcl::PointCloud<pcl::PointXYZ>::Ptr pts_;
		pcl::search::KdTree<pcl::PointXYZ>::Ptr kd_;
		pcl::PointCloud<pcl::Normal>::Ptr norms_;
		pcl::PointCloud<pcl::Boundary>::Ptr bd_;

	};

	struct Transform : boost::enable_shared_from_this<Transform> {

		typedef boost::shared_ptr<Transform> Ptr;

		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0) = 0;

		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0) = 0;
	};

	struct ICP {

		typedef boost::shared_ptr<ICP> Ptr;

		ICP(Data_Pack::Ptr target_data) : target_data_(target_data) {}

		virtual bool is(std::string name) = 0;
		virtual std::string get_classname() = 0;
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf = 0) = 0;
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf = 0) = 0; 

		Data_Pack::Ptr target_data_;
	};

	struct Correspondences {

		typedef boost::shared_ptr<Correspondences> Ptr;

		struct Parameters {
			float dist_thresh_;
			float angle_thresh_;
			bool accept_nonoverlap_corrs_;
			int sub_level_;
			Parameters(float dist_thresh = 0.02f, 
				float angle_thresh = 90.0f, 
				bool accept_nonoverlap_corrs = false, 
				int sub_level = 0) : 
				dist_thresh_(dist_thresh), 
				angle_thresh_(angle_thresh), 
				accept_nonoverlap_corrs_(accept_nonoverlap_corrs), 
				sub_level_(sub_level) {}
		};

		pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts_;
		pcl::PointCloud<pcl::Normal>::Ptr src_norms_;
		pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_pts_;
		pcl::PointCloud<pcl::Normal>::Ptr tgt_norms_;
		pcl::PointIndices::Ptr src_ids_;
		pcl::PointIndices::Ptr tgt_ids_;
		std::vector<float> distance2s_;
		std::vector<float> stabiblites_;
		std::vector<int> sub_levels_;

		Correspondences() : src_pts_(new pcl::PointCloud<pcl::PointXYZ>), tgt_pts_(new pcl::PointCloud<pcl::PointXYZ>),
		src_norms_(new pcl::PointCloud<pcl::Normal>), tgt_norms_(new pcl::PointCloud<pcl::Normal>),
		src_ids_(new pcl::PointIndices), tgt_ids_(new pcl::PointIndices) {}

		void clear();
		void compute_correspondences(Data_Pack::Ptr source_data, Data_Pack::Ptr target_data, Transform::Ptr tf, Parameters params);
		void compute_correspondences(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Data_Pack::Ptr target_data, Transform::Ptr tf, Parameters params);
		void extract_larger_than_deviation(pcl::PointIndices::Ptr &extracted_src_ids, float deviation_scale = 1.0f);
	};

	struct LinearTransform : virtual Transform {

		typedef boost::shared_ptr<LinearTransform> Ptr;

		LinearTransform(Eigen::Matrix4f mat = Eigen::Matrix4f::Identity()) : mat_(mat) {}

		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0);

		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr src_indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0);

		Eigen::Matrix4f mat_;
	};

	struct LinearICP : ICP {
		
		typedef boost::shared_ptr<LinearICP> Ptr;

		enum LinearType {RIGID, SCALE_RIGID, AFFINE};

		struct Parameters {
			LinearType linear_type_;
			int max_iter_;
			Parameters(LinearType linear_type = RIGID, int max_iter = 10) : 
				linear_type_(linear_type), max_iter_(max_iter) {};
		} linear_params_;

		Correspondences::Parameters corrs_params_;

		LinearICP(Data_Pack::Ptr target_data) : ICP(target_data) {}

		static LinearTransform::Ptr compute_linear(Correspondences::Ptr corrs, LinearType linear_type);

		virtual bool is(std::string name) {return (name == "LinearICP");}
		virtual std::string get_classname() {return "LinearICP";}
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf = 0);
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf = 0); 
	}; 

	struct TPSTransform : virtual Transform{

		typedef boost::shared_ptr<TPSTransform> Ptr;

		TPSTransform(TPSFunction tps = TPSFunction::Identity(3)) : tps_(tps) {}

		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0);

		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0);

		TPSFunction tps_;
	};

	struct ThinPlateSplinesICP : ICP {

		typedef boost::shared_ptr<ThinPlateSplinesICP> Ptr;

		struct Parameters {
			float lambda_;
			float kappa_;
			int max_iter_;
			float iter_rate_;
			Parameters(float lambda = 1e-5f, float kappa = 10.0f, int max_iter = 10, float iter_rate = 0.5f) : 
				lambda_(lambda), kappa_(kappa), max_iter_(max_iter), iter_rate_(iter_rate) {};
		} tps_params_;

		Correspondences::Parameters corrs_params_;

		ThinPlateSplinesICP(Data_Pack::Ptr target_data) : ICP(target_data) {}

		static TPSTransform::Ptr compute_tps(Correspondences::Ptr corrs, float lambda = 1e-5f, float kappa = 0.0f, int n = -1);

		virtual bool is(std::string name) {return (name == "ThinPlateSplinesICP");}
		virtual std::string get_classname() {return "ThinPlateSplinesICP";}
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf = 0);
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf = 0); 
	};

	struct LinearTPSTransform : LinearTransform, TPSTransform {

		typedef boost::shared_ptr<LinearTPSTransform> Ptr;

		LinearTPSTransform(Eigen::Matrix4f mat = Eigen::Matrix4f::Identity(), TPSFunction tps = TPSFunction::Identity(3)) : LinearTransform(mat), TPSTransform(tps) {}
		
		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0);
		
		virtual void transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
			pcl::PointCloud<pcl::Normal>::Ptr in_norms = 0, pcl::PointCloud<pcl::Normal>::Ptr out_norms = 0);
	};

	struct LinearThinPlateSplinesICP : ICP  {

		typedef boost::shared_ptr<LinearThinPlateSplinesICP> Ptr;

		struct Parameters {
			LinearICP::Parameters linear_params_;
			ThinPlateSplinesICP::Parameters tps_params_;

			float ec_tolerance_;

			Parameters(LinearICP::LinearType linear_type = LinearICP::RIGID, int linear_max_iter = 3, 
				float lambda = 1e-5f, float kappa = 10.0f, int tps_max_iter = 10, float iter_rate = 0.5f, float ec_tolerance = 0.08f) : 
			linear_params_(linear_type, linear_max_iter), tps_params_(lambda, kappa, tps_max_iter, iter_rate), ec_tolerance_(ec_tolerance) {};
		} lineartps_params_;

		Correspondences::Parameters corrs_params_;

		LinearThinPlateSplinesICP(Data_Pack::Ptr target_data) : ICP(target_data) {}

		virtual bool is(std::string name) {return (name == "LinearThinPlateSplinesICP");}
		virtual std::string get_classname() {return "LinearThinPlateSplinesICP";}
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf = 0);
		virtual Transform::Ptr do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf = 0); 
	};

}

#endif