#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <stdlib.h>

#include <pcl/common/transforms.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/filters/extract_indices.h>

#include <GMNR/LeastSquares/PointToPoint.h>
#include <GMNR/LeastSquares/PointToPlane.h>
#include "tps_icp.h"

using namespace gmnr;

std::string src_fmesh, tgt_fmesh;
int src_mesh_id, tgt_mesh_id;

std::string flinks;
std::set<int> linked_tgt_ids;

std::vector<std::string> tgt_fmeshs;
pcl::PolygonMesh::Ptr src_mesh, tgt_mesh;

std::string src_path, src_name, src_extension;
std::string tgt_path, tgt_name, tgt_extension;

std::string tgt_fbd;
pcl::PointCloud<pcl::Boundary>::Ptr tgt_bd;

pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr src_original, tgt_original;
pcl::PointCloud<pcl::Normal>::Ptr src_norms, tgt_norms;
pcl::PointCloud<pcl::PointXYZ>::Ptr src_cloud, tgt_cloud;

pcl::search::KdTree<pcl::PointXYZ>::Ptr tgt_kdtree;

pcl::PointCloud<pcl::PointXYZ>::Ptr src_sub_cloud;
pcl::PointCloud<pcl::Normal>::Ptr src_sub_norms;

pcl::PLYReader reader_ply;
pcl::PLYWriter writer_ply;
pcl::PCDReader reader_pcd;

pcl::visualization::PCLVisualizer::Ptr visualizer;

float sample_rate = 0.02f;
float dist_thresh = 0.01f;
float angle_thresh = 90;
LinearICP::LinearType linear_type = LinearICP::RIGID;
int linear_max_iter = 3;
float lambda = 1e-5f;
float kappa = 10.0f;
int tps_max_iter = 10;
float iter_rate = 0.5f;
float ec_tolerance = 0.08f;

std::string output_directory = ".";
bool write_sub_src = false;
bool write_warped_src = false;

int spin_time = 100;

struct MultiCorrespondence{
	pcl::PointXYZ src_pt_;
	int src_mesh_id_;
	std::vector<pcl::PointXYZ> tgt_pts_;
	std::vector<int> tgt_mesh_ids_;
};
std::vector<MultiCorrespondence> multicorrs;

std::string replace_file_extension(std::string file_name, std::string new_extension) {
	int len = file_name.length();
	int pos = file_name.find_last_of(".");
	std::string rt;
	rt.resize(pos + new_extension.length());
	int i;
	for(i = 0; i < pos; i++) rt[i] = file_name[i];
	for(; i < rt.length(); i++) rt[i] = new_extension[i-pos];
	return rt;
}

void get_file_path_name_extension(std::string file_name, std::string &path, std::string &name, std::string &extension) {
	int len = file_name.length();
	int begin_pos = file_name.find_last_of("/\\");
	if(begin_pos < 0) path = ".";
	else path = file_name.substr(0, begin_pos);
	int end_pos = file_name.find_last_of('.');
	if(end_pos < 0) extension = "";
	else extension = file_name.substr(end_pos+1, len - end_pos - 1);
	end_pos = end_pos >= 0 ? end_pos : len;
	name = file_name.substr(begin_pos+1, end_pos - begin_pos - 1);
}

void load_target_list() {

	if (!tgt_fmesh.empty()) {
		src_mesh_id = 0;
		tgt_fmeshs.push_back(src_fmesh);
		tgt_fmeshs.push_back(tgt_fmesh);
	}else {
		std::string file_name;
		while(std::getline(std::cin, file_name)) {
			tgt_fmeshs.push_back(file_name);
		}
		if(tgt_fmeshs.empty()) {
			PCL_ERROR ("No input tgt PLY file\n");
			exit (-1);
		}
		src_mesh_id = -1;
		for (int i = 0; i < tgt_fmeshs.size(); i++) {
			if(src_fmesh == tgt_fmeshs[i]) {
				src_mesh_id = i;
				break;
			}
		}
		if(src_mesh_id < 0) {
			PCL_ERROR ("src PLY file NOT IN tgt list\n");
			exit (-1);
		}
	}
}

void loadlinks() {

	std::ifstream links_ifs(flinks, std::ios::in);
	if(links_ifs) {
		std::string line;
		while(std::getline(links_ifs, line)) {
			std::stringstream ss(line);
			int a, b;
			ss >> a >> b;
			if(a == src_mesh_id) linked_tgt_ids.insert(b);
			else if(b == src_mesh_id) linked_tgt_ids.insert(a);
		}
		links_ifs.close();
	}else {
		for (int i = 0; i < tgt_fmeshs.size(); i++) {
			linked_tgt_ids.insert(i);
		}
	}
}

void normalize_norms(pcl::PointCloud<pcl::Normal>::Ptr norms) {
	for (int i = 0; i < norms->size(); i++) {
		if (norms->points[i].getNormalVector3fMap().squaredNorm() > 0.0f) {
			norms->points[i].getNormalVector3fMap().normalize();
		}
	}
}

void load_src_data() {
	src_mesh.reset(new pcl::PolygonMesh);
	if (reader_ply.read(src_fmesh, *src_mesh) == 0) {
		std::cerr << "height: " << src_mesh->cloud.height << " , " << "width: " <<  src_mesh->cloud.width << std::endl;
		std::cerr << "polygons: " << src_mesh->polygons.size() << std::endl;
	} else {
		PCL_ERROR ("Couldn't read file %s \n", src_fmesh);
		exit(-1);
	}
	if (visualizer) {
		visualizer->addPolygonMesh(*src_mesh, "source_mesh");
		visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.3, "source_mesh");
		visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.8, 0.0, 0.0, "source_mesh");
		visualizer->spinOnce(spin_time);
	}

	src_original.reset(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	pcl::fromPCLPointCloud2(src_mesh->cloud, *src_original);
	src_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::copyPointCloud(*src_original, *src_cloud);
	src_norms.reset(new pcl::PointCloud<pcl::Normal>);
	pcl::copyPointCloud(*src_original, *src_norms);
	normalize_norms(src_norms);
}

void load_tgt_data() {
	tgt_mesh.reset(new pcl::PolygonMesh);
	if (reader_ply.read(tgt_fmesh, *tgt_mesh) == 0) {
		std::cerr << "height: " << tgt_mesh->cloud.height << " , " << "width: " <<  tgt_mesh->cloud.width << std::endl;
		std::cerr << "polygons: " << tgt_mesh->polygons.size() << std::endl;
	} else {
		PCL_ERROR ("Couldn't read file %s \n", tgt_fmesh);
		exit (-1);
	}
	if (visualizer) {
		visualizer->addPolygonMesh(*tgt_mesh, "target_mesh");
		visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.3, "target_mesh");
		visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.8, "target_mesh");
		visualizer->spinOnce(spin_time);
	}

	tgt_fbd = replace_file_extension(tgt_fmesh, ".bd");
	tgt_bd.reset(new pcl::PointCloud<pcl::Boundary>);
    if (reader_pcd.read(tgt_fbd, *tgt_bd) != 0) { //* load the file
       PCL_ERROR ("Couldn't read file %s \n", tgt_fbd);
       exit (-1);
    }

	tgt_original.reset(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	pcl::fromPCLPointCloud2(tgt_mesh->cloud, *tgt_original);
	tgt_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::copyPointCloud(*tgt_original, *tgt_cloud);
	tgt_norms.reset(new pcl::PointCloud<pcl::Normal>);
	pcl::copyPointCloud(*tgt_original, *tgt_norms);
	normalize_norms(tgt_norms);
	tgt_kdtree.reset(new pcl::search::KdTree<pcl::PointXYZ>);
	tgt_kdtree->setInputCloud(tgt_cloud);
}

void process_args(int argc, char**argv) {

	std::vector<int> argv_ply_indices = pcl::console::parse_file_extension_argument(argc, argv, ".ply");
	if(argv_ply_indices.size() > 0) {
		src_fmesh = argv[argv_ply_indices[0]];
		get_file_path_name_extension(src_fmesh, src_path, src_name, src_extension);
	}
	else {
		PCL_ERROR ("No input src PLY file\n");
		exit (-1);
	}
	if(argv_ply_indices.size() > 1) {
		tgt_fmesh = argv[argv_ply_indices[1]];
		get_file_path_name_extension(tgt_fmesh, tgt_path, tgt_name, tgt_extension);
	}

	std::vector<int> argv_links_indices = pcl::console::parse_file_extension_argument(argc, argv, ".links");
	if(argv_links_indices.size() > 0) flinks = argv[argv_links_indices[0]];

	pcl::console::parse_argument(argc, argv, "--sample_rate", sample_rate);
	pcl::console::parse_argument(argc, argv, "--dist_thresh", dist_thresh);
	pcl::console::parse_argument(argc, argv, "--angle_thresh", angle_thresh);
	std::string linear_type_str;
	pcl::console::parse_argument(argc, argv, "--linear_type", linear_type_str);
	if(linear_type_str == "rigid") linear_type = LinearICP::RIGID;
	if(linear_type_str == "scale_rigid") linear_type = LinearICP::SCALE_RIGID;
	if(linear_type_str == "affine") linear_type = LinearICP::AFFINE;
	pcl::console::parse_argument(argc, argv, "--linear_max_iter", linear_max_iter);
	pcl::console::parse_argument(argc, argv, "--lambda", lambda);
	pcl::console::parse_argument(argc, argv, "--kappa", kappa);
	pcl::console::parse_argument(argc, argv, "--tps_max_iter", tps_max_iter);
	pcl::console::parse_argument(argc, argv, "--iter_rate", iter_rate);
	pcl::console::parse_argument(argc, argv, "--ec_tolerance", ec_tolerance);

	pcl::console::parse_argument(argc, argv, "--output_directory", output_directory);

	write_sub_src = pcl::console::find_switch(argc, argv, "--write_sub_src");
	write_warped_src = pcl::console::find_switch(argc, argv, "--write_warped_src");

	if(pcl::console::find_switch(argc, argv, "--visualize")) visualizer.reset(new pcl::visualization::PCLVisualizer());
	if(visualizer) visualizer->setBackgroundColor(0.5, 0.5, 0.5);
	pcl::console::parse_argument(argc, argv, "--spin_time", spin_time);

}

void create_sub_samples() {
	src_sub_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	src_sub_cloud->resize(sample_rate * src_cloud->size());
	src_sub_norms.reset(new pcl::PointCloud<pcl::Normal>);
	src_sub_norms->resize(src_sub_cloud->size());

	//srand((int)time(0));
	srand(0);
	//random sampling
	std::unordered_set<int> sampled_indices;
	int sub_count = 0;
	while (sub_count < src_sub_cloud->size()) {
		int id = rand() % src_cloud->size();
		if(sampled_indices.find(id) == sampled_indices.end()) {
			sampled_indices.insert(id);
			src_sub_cloud->points[sub_count] = src_cloud->points[id];
			src_sub_norms->points[sub_count] = src_norms->points[id];
			sub_count++;
		}
	}
	if (visualizer) {
		visualizer->addPointCloud(src_sub_cloud, "source_sub");
		visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "source_sub");
		visualizer->spinOnce(spin_time);
	}
	if(write_sub_src) writer_ply.write(output_directory + "/" + src_name + "_sub.ply", *src_sub_cloud);

	multicorrs.resize(src_sub_cloud->size());
	for (int i = 0; i < src_sub_cloud->size(); i++){
		multicorrs[i].src_mesh_id_ = src_mesh_id;
		multicorrs[i].src_pt_ = src_sub_cloud->points[i];
	}
	//return 0;
}

namespace gmnr {

	void Correspondences::clear() {
		src_pts_->clear();
		tgt_pts_->clear();

		src_ids_->indices.clear();
		tgt_ids_->indices.clear();

		src_norms_->clear();
		tgt_norms_->clear();

		distance2s_.clear();
	}

	void Correspondences::compute_correspondences(Data_Pack::Ptr source_data, Data_Pack::Ptr target_data, Transform::Ptr tf, Parameters params) {
		pcl::PointIndices::Ptr src_indices(new pcl::PointIndices);
		for (int i = 0; i < source_data->pts_->size(); i++) {
			src_indices->indices.push_back(i);
		}
		compute_correspondences(source_data, src_indices, target_data, tf, params);
	}

	void Correspondences::compute_correspondences(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Data_Pack::Ptr target_data, Transform::Ptr tf, Parameters params) {

		clear();

		float dist_thresh = params.dist_thresh_;
		float dist_thresh2 = dist_thresh * dist_thresh;
		float angle_thresh = params.angle_thresh_;
		float cos_angle_thresh = cosf(angle_thresh/180 * M_PI);
		int sub_level = params.sub_level_;

		pcl::PointCloud<pcl::PointXYZ>::Ptr source_pts(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::Normal>::Ptr source_norms(new pcl::PointCloud<pcl::Normal>);

		tf->transform(source_data->pts_, src_indices, source_pts, source_data->norms_, source_norms);

		pcl::PointCloud<pcl::PointXYZ>::Ptr target_pts = target_data->pts_;
		pcl::PointCloud<pcl::Normal>::Ptr target_norms = target_data->norms_;

		pcl::search::KdTree<pcl::PointXYZ>::Ptr target_kd = target_data->kd_;
		pcl::PointCloud<pcl::Boundary>::Ptr target_bd = target_data->bd_;

		int real_size = 0;
		for (int i = 0; i < source_pts->size(); i++) {
			std::vector<int> nearset_indices;
			std::vector<float> nearest_distance2s;

			if( target_kd->nearestKSearch(source_pts->points[i], 1, nearset_indices, nearest_distance2s) > 0 ) {

				if(nearest_distance2s[0] > dist_thresh2) continue;
				if(target_bd && target_bd->points[nearset_indices[0]].boundary_point) continue;
				if(target_norms) {
					Eigen::Vector3f tn = target_norms->points[nearset_indices[0]].getNormalVector3fMap();
					if(tn.squaredNorm() == 0.0f) continue;
					if(source_norms) {
						Eigen::Vector3f sn = source_norms->points[i].getNormalVector3fMap();
						if(sn.squaredNorm() == 0.0f) continue;
						if(sn.dot(tn) < cos_angle_thresh) continue;
					}

					Eigen::Vector3f sp = source_pts->points[i].getVector3fMap();
					Eigen::Vector3f tp = target_pts->points[ nearset_indices[0] ].getVector3fMap();
					float iterpola = acosf( (sp - tp).normalized().dot(tn) )  / (M_PI * 0.5f);
					float new_dist_thresh = (1.0f - 0.5f * iterpola) * dist_thresh;
					float new_dist_thresh2 = new_dist_thresh * new_dist_thresh;
					if(nearest_distance2s[0] > dist_thresh2) continue;

					Eigen::Vector3f project_sp = sp - (sp-tp).dot(tn) * tn;
					pcl::PointXYZ tgt_pt;
					tgt_pt.getVector3fMap() = project_sp;
					tgt_pts_->push_back(tgt_pt);

				} else tgt_pts_->push_back(target_pts->points[ nearset_indices[0] ]);
				//std::cerr << "tgt : " << tgt_pts_->points[real_size].x << " " << tgt_pts_->points[real_size].y << " " << tgt_pts_->points[real_size].z << std::endl;
				tgt_ids_->indices.push_back(nearset_indices[0]);

				src_pts_->push_back(source_pts->points[i]);
				//std::cerr << "src : " << src_pts_->points[real_size].x << " " << src_pts_->points[real_size].y << " " << src_pts_->points[real_size].z << std::endl;
				src_ids_->indices.push_back(src_indices->indices[i]);

				distance2s_.push_back( (src_pts_->back().getVector3fMap() - tgt_pts_->back().getVector3fMap()).squaredNorm() ) ;
				sub_levels_.push_back(sub_level);

				if(source_norms) {
					pcl::Normal sn = source_norms->points[i];
					src_norms_->push_back(sn);
				} 
				if(target_norms) {
					pcl::Normal tn = target_norms->points[nearset_indices[0]];
					tgt_norms_->push_back(tn);
				} 
 				real_size++;
			}
		}

		if(visualizer) {

			pcl::Correspondences s_t_corrs(real_size);
			for (int i = 0; i < real_size; i++) {
				s_t_corrs[i].index_query = i;
				s_t_corrs[i].index_match = i;
			}
			pcl::Correspondences s_s_corrs(real_size);
			for (int i = 0; i < real_size; i++) {
				s_s_corrs[i].index_query = i;
				s_s_corrs[i].index_match = src_ids_->indices[i];
			}
			visualizer->addPointCloud(src_pts_, "src_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "src_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "src_pts_");

			visualizer->addPointCloud(tgt_pts_, "tgt_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0, "tgt_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "tgt_pts_");

			visualizer->addCorrespondences<pcl::PointXYZ>(src_pts_, tgt_pts_, s_t_corrs, "s_t_correspondences");
			visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_t_correspondences");
			visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.0, "s_t_correspondences");

			visualizer->addCorrespondences<pcl::PointXYZ>(src_pts_, source_data->pts_, s_s_corrs, "s_s_correspondences");
			visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_s_correspondences");
			visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 1.0, "s_s_correspondences");

			visualizer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(src_pts_, src_norms_, 1, 0.003f, "src_norms_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "src_norms_");

			visualizer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(tgt_pts_, tgt_norms_, 1, 0.003f, "tgt_norms_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "tgt_norms_");

			visualizer->spinOnce(spin_time, true);

			visualizer->removePointCloud("src_pts_");
			visualizer->removePointCloud("tgt_pts_");
			visualizer->removeCorrespondences("s_t_correspondences");
			visualizer->removeCorrespondences("s_s_correspondences");
			visualizer->removePointCloud("src_norms_");
			visualizer->removePointCloud("tgt_norms_");
		}
	}

	void Correspondences::extract_larger_than_deviation(pcl::PointIndices::Ptr &extracted_src_ids, float deviation_scale) {

		extracted_src_ids->indices.clear();

		float sum_distance2 = 0.0f;
		for (int i = 0; i < distance2s_.size(); i++) {
			sum_distance2 += distance2s_[i];
		}
		float deviation_distance2 = sum_distance2 / distance2s_.size();
		float scaled_deviation_distance2 = deviation_distance2 * deviation_scale;

		int extracted_size = 0;

		for (int i = 0; i < distance2s_.size(); i++) {
			if (distance2s_[i] > scaled_deviation_distance2) {
				extracted_src_ids->indices.push_back(src_ids_->indices[i]);
				extracted_size++;
			}
		}
	}

	void extract_euclidean_cluster(pcl::PointCloud<pcl::PointXYZ>::Ptr pts, pcl::PointIndices::Ptr indices, std::vector<pcl::PointIndices::Ptr> &extracted_ids_cluster, float tolerance, int min_size, int max_size) {
		
		pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
		ec.setClusterTolerance(tolerance);
		ec.setMinClusterSize(min_size);
		ec.setMaxClusterSize(max_size);
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
		ec.setSearchMethod(tree);
		ec.setIndices(indices);
		ec.setInputCloud(pts);

		std::vector<pcl::PointIndices> ids_cluster;
		ec.extract (ids_cluster);
		for (int i = 0; i < ids_cluster.size(); i++) {
			extracted_ids_cluster.push_back(pcl::PointIndices::Ptr(new pcl::PointIndices));
			extracted_ids_cluster.back()->indices.swap(ids_cluster[i].indices);
		}

		for (int i = 0; i < extracted_ids_cluster.size(); i++) {
			std::cerr << "cluster " << i << " : " << extracted_ids_cluster[i]->indices.size() << std::endl;
		}
	}

	void extend_pointcloud(pcl::PointCloud<pcl::PointXYZ>::Ptr pts, pcl::search::KdTree<pcl::PointXYZ>::Ptr kd, pcl::PointIndices::Ptr seed_indices, float extend_width, pcl::PointIndices::Ptr &extend_indices) {

		std::unordered_set<int> extend_ids;
		std::vector<int> &indices = seed_indices->indices;
		for (int i = 0; i < indices.size(); i++) {
			pcl::PointXYZ pt = pts->points[indices[i]];
			std::vector<int> nearest_indices;
			std::vector<float> nearest_distance2s;
			if(kd->radiusSearch(pt, extend_width, nearest_indices, nearest_distance2s) > 0) {
				for (int j = 0; j < nearest_indices.size(); j++) extend_ids.insert(nearest_indices[j]);
			}
		}
		for(std::unordered_set<int>::iterator iter = extend_ids.begin(); iter != extend_ids.end(); iter++) {
			extend_indices->indices.push_back(*iter);
		}

		if(visualizer) {
			pcl::PointCloud<pcl::PointXYZ>::Ptr extracted_pts(new pcl::PointCloud<pcl::PointXYZ>);
			pcl::ExtractIndices<pcl::PointXYZ> pts_extractor;
			pts_extractor.setIndices(seed_indices);
			pts_extractor.setInputCloud(pts);
			pts_extractor.filter(*extracted_pts);

			visualizer->addPointCloud(extracted_pts, "large_deviation_src_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "large_deviation_src_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "large_deviation_src_pts_");	
			visualizer->spinOnce(spin_time);	
			visualizer->removePointCloud("large_deviation_src_pts_");

			pts_extractor.setIndices(extend_indices);
			pts_extractor.filter(*extracted_pts);

			visualizer->addPointCloud(extracted_pts, "extended_large_deviation_src_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "extended_large_deviation_src_pts_");
			visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "extended_large_deviation_src_pts_");	
			visualizer->spinOnce(spin_time);	
			visualizer->removePointCloud("extended_large_deviation_src_pts_");
		}

	}

	void combine_corrs_cluster(std::vector<Correspondences::Ptr> corrs_cluster, Correspondences::Ptr combined_corrs) {

		std::unordered_map<int, std::pair<int,int> > record;

		for (int i = 0; i < corrs_cluster.size(); i++) {
			std::vector<int> &src_ids = corrs_cluster[i]->src_ids_->indices;
			for (int j = 0; j < src_ids.size(); j++) {
				int src_id = src_ids[j];
				if(record.find(src_id) != record.end()) {
					std::pair<int,int> pair = record[src_id];
					int orig_i = pair.first, orig_j = pair.second;
					if(corrs_cluster[orig_i]->distance2s_[orig_j] > corrs_cluster[i]->distance2s_[j]) record[src_id] = std::pair<int,int>(i, j);
				}else record[src_id] = std::pair<int,int>(i, j);
			}
		}

		for(std::unordered_map<int, std::pair<int,int> >::iterator iter = record.begin(); iter != record.end(); iter++) {
			std::pair<int,int> pair = (*iter).second;
			int i = pair.first, j = pair.second;

			combined_corrs->src_ids_->indices.push_back(corrs_cluster[i]->src_ids_->indices[j]);
			combined_corrs->tgt_ids_->indices.push_back(corrs_cluster[i]->tgt_ids_->indices[j]);

			combined_corrs->src_pts_->push_back(corrs_cluster[i]->src_pts_->points[j]);
			combined_corrs->tgt_pts_->push_back(corrs_cluster[i]->tgt_pts_->points[j]);

			combined_corrs->src_norms_->push_back(corrs_cluster[i]->src_norms_->points[j]);
			combined_corrs->tgt_norms_->push_back(corrs_cluster[i]->tgt_norms_->points[j]);

			combined_corrs->distance2s_.push_back(corrs_cluster[i]->distance2s_[j]);
		}
	}

	void LinearTransform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
		pcl::PointCloud<pcl::Normal>::Ptr in_norms, pcl::PointCloud<pcl::Normal>::Ptr out_norms) {
		
		pcl::transformPointCloud(*in, *out, mat_);

		if(in_norms && out_norms) {
			out_norms->resize(in_norms->size());
			for (int i = 0; i < in_norms->size(); i++) {
				pcl::Normal &in_norm = in_norms->points[i];
				pcl::Normal &out_norm = out_norms->points[i];
				out_norm.getNormalVector3fMap() = (mat_.block(0,0,3,3) * in_norm.getNormalVector3fMap()).normalized(); 
			}
		}
	}

	void LinearTransform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
		pcl::PointCloud<pcl::Normal>::Ptr in_norms, pcl::PointCloud<pcl::Normal>::Ptr out_norms) {

		pcl::ExtractIndices<pcl::PointXYZ> pts_extractor;
		pts_extractor.setIndices(indices);
		pts_extractor.setInputCloud(in);
		pts_extractor.filter(*out);

		if (in_norms && out_norms) {
			pcl::ExtractIndices<pcl::Normal> norms_extractor;
			norms_extractor.setIndices(indices);
			norms_extractor.setInputCloud(in_norms);
			norms_extractor.filter(*out_norms);
		}

		transform(out, out, out_norms, out_norms);
	}

	LinearTransform::Ptr LinearICP::compute_linear(Correspondences::Ptr corrs, LinearType linear_type) {

		pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = corrs->src_pts_;
		pcl::PointCloud<pcl::Normal>::Ptr src_norms = corrs->src_norms_;
		pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_pts = corrs->tgt_pts_;
		pcl::PointCloud<pcl::Normal>::Ptr tgt_norms = corrs->tgt_norms_;

		PointSet3D src(src_pts->size(), 3), tgt(tgt_pts->size(), 3);
		for (int i = 0; i < src_pts->size(); i++) src.row(i) = src_pts->points[i].getVector3fMap().cast<gmnr::Scalar>();
		for (int i = 0; i < tgt_pts->size(); i++) tgt.row(i) = tgt_pts->points[i].getVector3fMap().cast<gmnr::Scalar>();

		switch (linear_type) {
		case RIGID: {
				if (tgt_norms) {
					NormalSet3D tnorms(tgt_norms->size(), 3);
					for (int i = 0; i < tgt_norms->size(); i++) tnorms.row(i) = tgt_norms->points[i].getNormalVector3fMap().cast<gmnr::Scalar>();
					PointToPlaneLinear point_to_plane(src, tgt, tnorms, true);
					return LinearTransform::Ptr(new LinearTransform(point_to_plane.transfomation().cast<float>()));
				}
				PointToPointUmeyama point_to_point(src, tgt, false, true);
				return LinearTransform::Ptr(new LinearTransform(point_to_point.transfomation().cast<float>()));
			}
		case SCALE_RIGID: {
				PointToPointUmeyama point_to_point(src, tgt, true, true);
				return LinearTransform::Ptr(new LinearTransform(point_to_point.transfomation().cast<float>()));
				//Eigen::Matrix<float, 3, Eigen::Dynamic> src(3, src_pts->size()), tgt(3, tgt_pts->size());
				//for (int i = 0; i < src_pts->size(); i++) src.col(i) = src_pts->points[i].getVector3fMap();
				//for (int i = 0; i < tgt_pts->size(); i++) tgt.col(i) = tgt_pts->points[i].getVector3fMap();
				//return LinearTransform::Ptr(new LinearTransform(pcl::umeyama(src, tgt, linear_type == SCALE_RIGID)));
			}
		case AFFINE: 
		default: return LinearTransform::Ptr(new LinearTransform());
		}
	}

	Transform::Ptr LinearICP::do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf) {
		pcl::PointIndices::Ptr src_indices(new pcl::PointIndices);
		for (int i = 0; i < source_data->pts_->size(); i++) {
			src_indices->indices.push_back(i);
		}
		return do_icp(source_data, src_indices, init_tf);
	} 

	Transform::Ptr LinearICP::do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf) {

		Correspondences::Ptr corrs(new Correspondences);

		LinearTransform::Ptr linear_tf(new LinearTransform);
		if (init_tf) linear_tf->mat_ = boost::dynamic_pointer_cast<LinearTransform>(init_tf)->mat_;
		else init_tf.reset(new LinearTransform);

		int max_iter = linear_params_.max_iter_;
		LinearType linear_type = linear_params_.linear_type_;

		int iter = 0;
		while(iter < max_iter) {
			corrs->compute_correspondences(source_data, src_indices, target_data_, linear_tf, corrs_params_);
			if(corrs->src_pts_->size() < 10) return init_tf;
			LinearTransform::Ptr temp_linear_tf = compute_linear(corrs, linear_type);
			linear_tf->mat_ = temp_linear_tf->mat_ * linear_tf->mat_;
			iter++;
		}
		return linear_tf;
	} 

	void TPSTransform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out,
		pcl::PointCloud<pcl::Normal>::Ptr in_norms, pcl::PointCloud<pcl::Normal>::Ptr out_norms) {

		out->resize(in->size());
		for (int i = 0; i < in->size(); i++) {
			Eigen::Vector3f pt = in->points[i].getVector3fMap();
			Eigen::Vector3f pt_tps = (tps_.evaluate(pt.transpose().cast<gmnr::Scalar>())).transpose().cast<float>();
			out->points[i].getVector3fMap() = pt_tps;
		}

		if (in_norms && out_norms) {
			static float delta = 1e-5f;
			out_norms->resize(in_norms->size());
			for (int i = 0; i < in_norms->size(); i++) {
				pcl::Normal &in_norm = in_norms->points[i];
				pcl::Normal &out_norm = out_norms->points[i];

				Eigen::Vector3f pt_delta = in->points[i].getVector3fMap() + in_norm.getNormalVector3fMap() * delta;
				Eigen::Vector3f pt_delta_tps = (tps_.evaluate(pt_delta.transpose().cast<gmnr::Scalar>())).transpose().cast<float>();

				out_norm.getNormalVector3fMap() = (pt_delta_tps - out->points[i].getVector3fMap()).normalized(); 
			}		
		}
	}

	void TPSTransform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
		pcl::PointCloud<pcl::Normal>::Ptr in_norms, pcl::PointCloud<pcl::Normal>::Ptr out_norms) {

			pcl::ExtractIndices<pcl::PointXYZ> pts_extractor;
			pts_extractor.setIndices(indices);
			pts_extractor.setInputCloud(in);
			pts_extractor.filter(*out);

			if (in_norms && out_norms) {
				pcl::ExtractIndices<pcl::Normal> norms_extractor;
				norms_extractor.setIndices(indices);
				norms_extractor.setInputCloud(in_norms);
				norms_extractor.filter(*out_norms);
			}

			transform(out, out, out_norms, out_norms);
	}

	TPSTransform::Ptr  ThinPlateSplinesICP::compute_tps(Correspondences::Ptr corrs, float lambda, float kappa, int n) {

		pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = corrs->src_pts_;
		pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_pts = corrs->tgt_pts_;

		PointSet3D src(src_pts->size(), 3), tgt(tgt_pts->size(), 3);
		for (int i = 0; i < src_pts->size(); i++) src.row(i) = src_pts->points[i].getVector3fMap().cast<gmnr::Scalar>();
		for (int i = 0; i < tgt_pts->size(); i++) tgt.row(i) = tgt_pts->points[i].getVector3fMap().cast<gmnr::Scalar>();

		std::cerr << "lambda = " << lambda << " kappa = " << kappa << std::endl;

		if(n > 0) return TPSTransform::Ptr(new TPSTransform(ApproxiTPSFunction(src, tgt, lambda, kappa, n)));
		return TPSTransform::Ptr(new TPSTransform(TPSFunction(src, tgt, lambda, kappa)));
	}

	Transform::Ptr ThinPlateSplinesICP::do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf) {
		pcl::PointIndices::Ptr src_indices(new pcl::PointIndices);
		for (int i = 0; i < source_data->pts_->size(); i++) {
			src_indices->indices.push_back(i);
		}
		return do_icp(source_data, src_indices, init_tf);		
	}

	Transform::Ptr ThinPlateSplinesICP::do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf) {

		Correspondences::Ptr corrs(new Correspondences);

		TPSTransform::Ptr tps_tf;
		if (init_tf) tps_tf->tps_ = boost::dynamic_pointer_cast<TPSTransform>(init_tf)->tps_;
		else init_tf.reset(new TPSTransform);

		float lambda = tps_params_.lambda_;
		float kappa = tps_params_.kappa_;
		int max_iter = tps_params_.max_iter_;
		float iter_rate = tps_params_.iter_rate_;
		float frac = 1.0f / max_iter;
		int iter = 0;
		while(iter < max_iter) {

			corrs->compute_correspondences(source_data, src_indices, target_data_, tps_tf, corrs_params_);

			std::vector<int> &src_ids = corrs->src_ids_->indices;
			pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = corrs->src_pts_;
			for (int i = 0; i < src_ids.size(); i++) src_pts->points[i] = source_data->pts_->points[ src_ids[i] ];

			tps_tf = compute_tps(corrs, lambda, kappa, (int)(corrs->src_pts_->size() * frac * (iter+1)));
			kappa *= iter_rate;
			iter++;
		}
		return tps_tf;
	}

	void LinearTPSTransform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
		pcl::PointCloud<pcl::Normal>::Ptr in_norms, pcl::PointCloud<pcl::Normal>::Ptr out_norms) {

		//Eigen::Matrix4f mat = LinearTransform::mat_;
		//TPSFunction tps = TPSTransform::tps_;

		pcl::PointCloud<pcl::PointXYZ>::Ptr out_tmp(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::Normal>::Ptr out_norms_tmp(new pcl::PointCloud<pcl::Normal>);
		LinearTransform::transform(in, out_tmp, in_norms, out_norms_tmp);
		TPSTransform::transform(out_tmp, out, out_norms_tmp, out_norms);
	}

	void LinearTPSTransform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
		pcl::PointCloud<pcl::Normal>::Ptr in_norms, pcl::PointCloud<pcl::Normal>::Ptr out_norms) {

			pcl::ExtractIndices<pcl::PointXYZ> pts_extractor;
			pts_extractor.setIndices(indices);
			pts_extractor.setInputCloud(in);
			pts_extractor.filter(*out);

			if (in_norms && out_norms) {
				pcl::ExtractIndices<pcl::Normal> norms_extractor;
				norms_extractor.setIndices(indices);
				norms_extractor.setInputCloud(in_norms);
				norms_extractor.filter(*out_norms);
			}

			transform(out, out, out_norms, out_norms);
	}

	Transform::Ptr LinearThinPlateSplinesICP::do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf) {
		pcl::PointIndices::Ptr src_indices(new pcl::PointIndices);
		for (int i = 0; i < source_data->pts_->size(); i++) {
			src_indices->indices.push_back(i);
		}
		return do_icp(source_data, src_indices, init_tf);
	}

	Transform::Ptr LinearThinPlateSplinesICP::do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf) {

		Correspondences::Ptr corrs(new Correspondences);
		pcl::PointCloud<pcl::PointXYZ>::Ptr source_pts(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::Normal>::Ptr source_norms(new pcl::PointCloud<pcl::Normal>);

		LinearTPSTransform::Ptr lineartps_tf(new LinearTPSTransform);
		if (init_tf) {
			lineartps_tf->mat_ = boost::dynamic_pointer_cast<LinearTPSTransform>(init_tf)->mat_;
			lineartps_tf->tps_ = boost::dynamic_pointer_cast<LinearTPSTransform>(init_tf)->tps_;
		} else init_tf.reset(new LinearTPSTransform);

		LinearICP::Parameters linear_params = lineartps_params_.linear_params_;
		ThinPlateSplinesICP::Parameters tps_params = lineartps_params_.tps_params_;

		float lambda = tps_params.lambda_;
		float kappa = tps_params.kappa_;
		int max_iter = tps_params.max_iter_;
		float iter_rate = tps_params.iter_rate_;
		float frac = 1.0f / max_iter;
		int iter = 0;
		while(iter < max_iter) {

			lineartps_tf->transform(source_data->pts_, src_indices, source_pts, source_data->norms_, source_norms);
			Data_Pack::Ptr compact_source_data(new Data_Pack(source_pts, 0, source_norms, 0));

			LinearICP::Ptr linear_icp(new LinearICP(target_data_));
			linear_icp->linear_params_ = linear_params;
			linear_icp->corrs_params_ = corrs_params_;
			Transform::Ptr tf = linear_icp->do_icp(compact_source_data);
			corrs->compute_correspondences(compact_source_data, target_data_, tf, corrs_params_);

			//Extract large error src_ids
			pcl::PointIndices::Ptr extracted_src_ids(new pcl::PointIndices);
			corrs->extract_larger_than_deviation(extracted_src_ids, 1.0f);

			//Cluster according to euclidean distance
			std::vector<pcl::PointIndices::Ptr> extracted_ids_cluster;
			float tolerance = lineartps_params_.ec_tolerance_;
			int min_size = 1;
			int max_size = 10000000;
			extract_euclidean_cluster(source_pts, extracted_src_ids, extracted_ids_cluster, tolerance, min_size, max_size);

			std::vector<Correspondences::Ptr> corrs_cluster;
			pcl::search::KdTree<pcl::PointXYZ>::Ptr source_kd(new pcl::search::KdTree<pcl::PointXYZ>);
			source_kd->setInputCloud(source_pts);
			float extend_width = tolerance / (sqrtf(3.0f) - 0.1f);
			for (int i = 0; i < extracted_ids_cluster.size(); i++) {
				pcl::PointIndices::Ptr part_indices = extracted_ids_cluster[i];
				pcl::PointIndices::Ptr extend_indices(new pcl::PointIndices);
				extend_pointcloud(source_pts, source_kd, part_indices, extend_width, extend_indices);
				if (extend_indices->indices.size() < 20) continue;
				Transform::Ptr part_tf = linear_icp->do_icp(compact_source_data, extend_indices);
				Correspondences::Ptr part_corrs(new Correspondences);
				part_corrs->compute_correspondences(compact_source_data, extend_indices, target_data_, part_tf, corrs_params_);
				corrs_cluster.push_back(part_corrs);
			}
			
			corrs_cluster.push_back(corrs);
			Correspondences::Ptr combined_corrs(new Correspondences);
			combine_corrs_cluster(corrs_cluster, combined_corrs);

			if(visualizer) {

				int real_size = combined_corrs->src_pts_->size();

				pcl::Correspondences s_t_corrs(real_size);
				for (int i = 0; i < real_size; i++) {
					s_t_corrs[i].index_query = i;
					s_t_corrs[i].index_match = i;
				}
				pcl::Correspondences s_s_corrs(real_size);
				for (int i = 0; i < real_size; i++) {
					s_s_corrs[i].index_query = i;
					s_s_corrs[i].index_match = src_indices->indices[ combined_corrs->src_ids_->indices[i] ];
				}
				visualizer->addPointCloud(combined_corrs->src_pts_, "src_pts_");
				visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "src_pts_");
				visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "src_pts_");

				visualizer->addPointCloud(combined_corrs->tgt_pts_, "tgt_pts_");
				visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0, "tgt_pts_");
				visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "tgt_pts_");

				visualizer->addCorrespondences<pcl::PointXYZ>(combined_corrs->src_pts_, combined_corrs->tgt_pts_, s_t_corrs, "s_t_correspondences");
				visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_t_correspondences");
				visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.0, "s_t_correspondences");

				visualizer->addCorrespondences<pcl::PointXYZ>(combined_corrs->src_pts_, source_data->pts_, s_s_corrs, "s_s_correspondences");
				visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_s_correspondences");
				visualizer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 1.0, "s_s_correspondences");

				visualizer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(combined_corrs->src_pts_, combined_corrs->src_norms_, 1, 0.003f, "src_norms_");
				visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "src_norms_");

				visualizer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(combined_corrs->tgt_pts_, combined_corrs->tgt_norms_, 1, 0.003f, "tgt_norms_");
				visualizer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "tgt_norms_");

				visualizer->spinOnce(spin_time+200000, true);

				visualizer->removePointCloud("src_pts_");
				visualizer->removePointCloud("tgt_pts_");
				visualizer->removeCorrespondences("s_t_correspondences");
				visualizer->removeCorrespondences("s_s_correspondences");
				visualizer->removePointCloud("src_norms_");
				visualizer->removePointCloud("tgt_norms_");
			}

			LinearTransform::Ptr linear_tf;

			////Approximate solution
			//{
			//	linear_tf = boost::dynamic_pointer_cast<LinearTransform>(tf);
			//	linear_tf->mat_ = linear_tf->mat_ * lineartps_tf->mat_;
			//	linear_tf->transform(source_data->pts_, src_indices, source_pts);

			//	std::vector<int> &src_ids =  corrs->src_ids_;
			//	pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = corrs->src_pts_;
			//	for (int i = 0; i < src_ids.size(); i++) src_pts->points[i] = source_pts->points[ src_indices->indices[src_ids[i]] ];
			//}

			//Theoretical solution
			{
				std::vector<int> &src_ids = combined_corrs->src_ids_->indices;
				pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = combined_corrs->src_pts_;
				pcl::PointCloud<pcl::Normal>::Ptr src_norms = combined_corrs->src_norms_;
				for (int i = 0; i < src_ids.size(); i++) {
					src_pts->points[i] = source_data->pts_->points[ src_indices->indices[src_ids[i]] ];
					src_norms->points[i] = source_data->norms_->points[ src_indices->indices[src_ids[i]] ];
				}

				pcl::PointCloud<pcl::Normal>::Ptr tgt_norms_tmp = corrs->tgt_norms_;
				combined_corrs->tgt_norms_ = 0;
				linear_tf = LinearICP::compute_linear(corrs, LinearICP::RIGID);
				combined_corrs->tgt_norms_ = tgt_norms_tmp;
				linear_tf->transform(source_data->pts_, src_indices, source_pts);
				for (int i = 0; i < src_ids.size(); i++) src_pts->points[i] = source_pts->points[ src_ids[i] ];
			}

			TPSTransform::Ptr tps_tf = ThinPlateSplinesICP::compute_tps(combined_corrs, lambda, kappa, (int)(combined_corrs->src_pts_->size() * frac * (iter+1)));

			lineartps_tf->mat_ = linear_tf->mat_;
			lineartps_tf->tps_ = tps_tf->tps_;

			kappa *= iter_rate;
			iter++;
		}
		
		std::vector<int> &src_ids = corrs->src_ids_->indices;
		pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_pts = corrs->src_pts_;
		for (int i = 0; i < src_ids.size(); i++) {
			MultiCorrespondence &multicorr = multicorrs[src_ids[i]];
			multicorr.tgt_mesh_ids_.push_back(tgt_mesh_id);
			multicorr.tgt_pts_.push_back(tgt_pts->points[i]);
		}

		return lineartps_tf;
	}

};

int main(int argc, char** argv) {

	process_args(argc, argv);
	load_target_list();
	loadlinks();

	load_src_data();
	create_sub_samples();

	Data_Pack::Ptr source_data(new Data_Pack(src_sub_cloud, 0, src_sub_norms, 0));

	for (tgt_mesh_id = 0; tgt_mesh_id < tgt_fmeshs.size(); tgt_mesh_id++) {
		if(tgt_mesh_id == src_mesh_id || linked_tgt_ids.find(tgt_mesh_id) == linked_tgt_ids.end()) continue;

		tgt_fmesh = tgt_fmeshs[tgt_mesh_id];
		get_file_path_name_extension(tgt_fmesh, tgt_path, tgt_name, tgt_extension);

		load_tgt_data();

		Data_Pack::Ptr target_data(new Data_Pack(tgt_cloud, tgt_kdtree, tgt_norms, tgt_bd));

		//LinearICP::Ptr icp(new LinearThinPlateSplinesICP(target_data));
		//ThinPlateSplinesICP::Ptr icp(new ThinPlateSplinesICP(target_data));
		LinearThinPlateSplinesICP::Ptr icp(new LinearThinPlateSplinesICP(target_data));
		icp->corrs_params_ = Correspondences::Parameters(dist_thresh, angle_thresh);
		icp->lineartps_params_ = LinearThinPlateSplinesICP::Parameters(linear_type, linear_max_iter, lambda, kappa, tps_max_iter, iter_rate, ec_tolerance);
		Transform::Ptr tf = icp->do_icp(source_data);

		if (visualizer) {
			visualizer->spinOnce(spin_time);
			visualizer->removePolygonMesh("target_mesh");
		}

		if (write_warped_src) {
			std::string output_name = output_directory + "/" + src_name + "-" + tgt_name + "_" + icp->get_classname() + ".ply";
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZ>);
			tf->transform(src_cloud, cloud_temp);
			pcl::copyPointCloud(*cloud_temp, *src_original);
			pcl::toPCLPointCloud2(*src_original, src_mesh->cloud);

			if(pcl::io::savePLYFile(output_name, *src_mesh) != 0) {
				PCL_ERROR ("Couldn't write file %s \n", output_name);
				exit (-1);
			}
		}
	}

	std::cout << src_mesh_id << " " << multicorrs.size() << " " << src_fmesh << std::endl;

	for (int i = 0; i < multicorrs.size(); i++) {
		MultiCorrespondence &multicorr = multicorrs[i];
		int sm_id = multicorr.src_mesh_id_;
		pcl::PointXYZ src_pt = multicorr.src_pt_;
		std::vector<int> &tm_ids = multicorr.tgt_mesh_ids_;
		std::vector<pcl::PointXYZ> &tgt_pts = multicorr.tgt_pts_;

		std::cout << src_pt.x << " " << src_pt.y << " " << src_pt.z << " " << tm_ids.size() << std::endl;

		for (int j = 0; j < tm_ids.size(); j++) {
			int tm_id = tm_ids[j];
			pcl::PointXYZ tgt_pt = tgt_pts[j];
			std::cout << tm_id << " " << tgt_pt.x << " " << tgt_pt.y << " " << tgt_pt.z << std::endl;
		}
	}

	return 0;
}




