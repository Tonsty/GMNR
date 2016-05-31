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

std::string src_fmesh_g, tgt_fmesh_g;
int src_mesh_id_g, tgt_mesh_id_g;

std::string flinks_g;
std::set<int> linked_tgt_ids_g;

std::vector<std::string> tgt_fmeshs_g;
pcl::PolygonMesh::Ptr src_mesh_g, tgt_mesh_g;

std::string src_path_g, src_name_g, src_extension_g;
std::string tgt_path_g, tgt_name_g, tgt_extension_g;

std::string tgt_fbd_g;
pcl::PointCloud<pcl::Boundary>::Ptr tgt_bd_g;

pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr src_original_g, tgt_original_g;
pcl::PointCloud<pcl::Normal>::Ptr src_norms_g, tgt_norms_g;
pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts_g, tgt_pts_g;

pcl::search::KdTree<pcl::PointXYZ>::Ptr tgt_kdtree_g;

pcl::PointCloud<pcl::PointXYZ>::Ptr src_sub_pts_g;
pcl::PointCloud<pcl::Normal>::Ptr src_sub_norms_g;

pcl::PLYReader reader_ply_g;
pcl::PLYWriter writer_ply_g;
pcl::PCDReader reader_pcd_g;

pcl::visualization::PCLVisualizer::Ptr visualizer_g;

std::string icp_type_g = "linear";

float sample_rate_g = 0.02f;
float dist_thresh_g = 0.01f;
float angle_thresh_g = 90;
LinearICP::LinearType linear_type_g = LinearICP::RIGID;
int linear_max_iter_g = 3;
float lambda_g = 1e-5f;
float kappa_g = 10.0f;
int tps_max_iter_g = 10;
float iter_rate_g = 0.5f;
float ec_tolerance_g = 0.08f;

std::string output_directory_g = ".";
bool write_sub_src_g = false;
bool write_warped_src_g = false;

int spin_time_g = 100;

struct MultiCorrespondence{
	pcl::PointXYZ src_pt_;
	int src_mesh_id_;
	std::vector<pcl::PointXYZ> tgt_pts_;
	std::vector<int> tgt_mesh_ids_;
};
std::vector<MultiCorrespondence> multicorrs_g;

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

	if (!tgt_fmesh_g.empty()) {
		src_mesh_id_g = 0;
		tgt_fmeshs_g.push_back(src_fmesh_g);
		tgt_fmeshs_g.push_back(tgt_fmesh_g);
	}else {
		std::string file_name;
		while(std::getline(std::cin, file_name)) {
			tgt_fmeshs_g.push_back(file_name);
		}
		if(tgt_fmeshs_g.empty()) {
			PCL_ERROR ("No input tgt PLY file\n");
			exit (-1);
		}
		src_mesh_id_g = -1;
		for (int i = 0; i < tgt_fmeshs_g.size(); i++) {
			if(src_fmesh_g == tgt_fmeshs_g[i]) {
				src_mesh_id_g = i;
				break;
			}
		}
		if(src_mesh_id_g < 0) {
			PCL_ERROR ("src PLY file NOT IN tgt list\n");
			exit (-1);
		}
	}
}

void loadlinks() {

	std::ifstream links_ifs(flinks_g, std::ios::in);
	if(links_ifs) {
		std::string line;
		while(std::getline(links_ifs, line)) {
			std::stringstream ss(line);
			int a, b;
			ss >> a >> b;
			if(a == src_mesh_id_g) linked_tgt_ids_g.insert(b);
			else if(b == src_mesh_id_g) linked_tgt_ids_g.insert(a);
		}
		links_ifs.close();
	}else {
		for (int i = 0; i < tgt_fmeshs_g.size(); i++) {
			linked_tgt_ids_g.insert(i);
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
	src_mesh_g.reset(new pcl::PolygonMesh);
	if (reader_ply_g.read(src_fmesh_g, *src_mesh_g) == 0) {
		std::cerr << "height: " << src_mesh_g->cloud.height << " , " << "width: " <<  src_mesh_g->cloud.width << std::endl;
		std::cerr << "polygons: " << src_mesh_g->polygons.size() << std::endl;
	} else {
		PCL_ERROR ("Couldn't read file %s \n", src_fmesh_g);
		exit(-1);
	}
	if (visualizer_g) {
		visualizer_g->addPolygonMesh(*src_mesh_g, "source_mesh");
		visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.3, "source_mesh");
		visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.8, 0.0, 0.0, "source_mesh");
		visualizer_g->spinOnce(spin_time_g);
	}

	src_original_g.reset(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	pcl::fromPCLPointCloud2(src_mesh_g->cloud, *src_original_g);
	src_pts_g.reset(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::copyPointCloud(*src_original_g, *src_pts_g);
	src_norms_g.reset(new pcl::PointCloud<pcl::Normal>);
	pcl::copyPointCloud(*src_original_g, *src_norms_g);
	normalize_norms(src_norms_g);
}

void load_tgt_data() {
	tgt_mesh_g.reset(new pcl::PolygonMesh);
	if (reader_ply_g.read(tgt_fmesh_g, *tgt_mesh_g) == 0) {
		std::cerr << "height: " << tgt_mesh_g->cloud.height << " , " << "width: " <<  tgt_mesh_g->cloud.width << std::endl;
		std::cerr << "polygons: " << tgt_mesh_g->polygons.size() << std::endl;
	} else {
		PCL_ERROR ("Couldn't read file %s \n", tgt_fmesh_g);
		exit (-1);
	}
	if (visualizer_g) {
		visualizer_g->addPolygonMesh(*tgt_mesh_g, "target_mesh");
		visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.3, "target_mesh");
		visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.8, "target_mesh");
		visualizer_g->spinOnce(spin_time_g);
	}

	tgt_fbd_g = replace_file_extension(tgt_fmesh_g, ".bd");
	tgt_bd_g.reset(new pcl::PointCloud<pcl::Boundary>);
    if (reader_pcd_g.read(tgt_fbd_g, *tgt_bd_g) != 0) { //* load the file
       PCL_ERROR ("Couldn't read file %s \n", tgt_fbd_g);
       exit (-1);
    }

	tgt_original_g.reset(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	pcl::fromPCLPointCloud2(tgt_mesh_g->cloud, *tgt_original_g);
	tgt_pts_g.reset(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::copyPointCloud(*tgt_original_g, *tgt_pts_g);
	tgt_norms_g.reset(new pcl::PointCloud<pcl::Normal>);
	pcl::copyPointCloud(*tgt_original_g, *tgt_norms_g);
	normalize_norms(tgt_norms_g);
	tgt_kdtree_g.reset(new pcl::search::KdTree<pcl::PointXYZ>);
	tgt_kdtree_g->setInputCloud(tgt_pts_g);
}

void process_args(int argc, char**argv) {

	std::vector<int> argv_ply_indices = pcl::console::parse_file_extension_argument(argc, argv, ".ply");
	if(argv_ply_indices.size() > 0) {
		src_fmesh_g = argv[argv_ply_indices[0]];
		get_file_path_name_extension(src_fmesh_g, src_path_g, src_name_g, src_extension_g);
	}
	else {
		PCL_ERROR ("No input src PLY file\n");
		exit (-1);
	}
	if(argv_ply_indices.size() > 1) {
		tgt_fmesh_g = argv[argv_ply_indices[1]];
		get_file_path_name_extension(tgt_fmesh_g, tgt_path_g, tgt_name_g, tgt_extension_g);
	}

	std::vector<int> argv_links_indices = pcl::console::parse_file_extension_argument(argc, argv, ".links");
	if(argv_links_indices.size() > 0) flinks_g = argv[argv_links_indices[0]];

	pcl::console::parse_argument(argc, argv, "--icp_type", icp_type_g);
	pcl::console::parse_argument(argc, argv, "--sample_rate", sample_rate_g);
	pcl::console::parse_argument(argc, argv, "--dist_thresh", dist_thresh_g);
	pcl::console::parse_argument(argc, argv, "--angle_thresh", angle_thresh_g);
	std::string linear_type_str;
	pcl::console::parse_argument(argc, argv, "--linear_type", linear_type_str);
	if(linear_type_str == "rigid") linear_type_g = LinearICP::RIGID;
	if(linear_type_str == "scale_rigid") linear_type_g = LinearICP::SCALE_RIGID;
	if(linear_type_str == "affine") linear_type_g = LinearICP::AFFINE;
	pcl::console::parse_argument(argc, argv, "--linear_max_iter", linear_max_iter_g);
	pcl::console::parse_argument(argc, argv, "--lambda", lambda_g);
	pcl::console::parse_argument(argc, argv, "--kappa", kappa_g);
	pcl::console::parse_argument(argc, argv, "--tps_max_iter", tps_max_iter_g);
	pcl::console::parse_argument(argc, argv, "--iter_rate", iter_rate_g);
	pcl::console::parse_argument(argc, argv, "--ec_tolerance", ec_tolerance_g);

	pcl::console::parse_argument(argc, argv, "--output_directory", output_directory_g);

	write_sub_src_g = pcl::console::find_switch(argc, argv, "--write_sub_src");
	write_warped_src_g = pcl::console::find_switch(argc, argv, "--write_warped_src");

	if(pcl::console::find_switch(argc, argv, "--visualize")) visualizer_g.reset(new pcl::visualization::PCLVisualizer());
	if(visualizer_g) visualizer_g->setBackgroundColor(0.5, 0.5, 0.5);
	pcl::console::parse_argument(argc, argv, "--spin_time", spin_time_g);

}

void create_sub_samples() {
	src_sub_pts_g.reset(new pcl::PointCloud<pcl::PointXYZ>);
	src_sub_norms_g.reset(new pcl::PointCloud<pcl::Normal>);

	//Sampler::Ptr sampler(new RandomSampler(sample_rate_g));
	Sampler::Ptr sampler(new GridSampler(0.02f, 0.02f, 0.02f));
	pcl::PointIndicesPtr sampled_indices(new pcl::PointIndices);
	sampler->sample(src_pts_g, src_sub_pts_g, sampled_indices);

	pcl::ExtractIndices<pcl::Normal> norms_extractor;
	norms_extractor.setIndices(sampled_indices);
	norms_extractor.setInputCloud(src_norms_g);
	norms_extractor.filter(*src_sub_norms_g);

	if (visualizer_g) {
		visualizer_g->addPointCloud(src_sub_pts_g, "source_sub");
		visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "source_sub");

		visualizer_g->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(src_sub_pts_g, src_sub_norms_g, 1, 0.003f, "source_sub_norms_");
		visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "source_sub_norms_");

		visualizer_g->spinOnce(200000 + spin_time_g);
	}
	if(write_sub_src_g) writer_ply_g.write(output_directory_g + "/" + src_name_g + "_sub.ply", *src_sub_pts_g);

	multicorrs_g.resize(src_sub_pts_g->size());
	for (int i = 0; i < src_sub_pts_g->size(); i++){
		multicorrs_g[i].src_mesh_id_ = src_mesh_id_g;
		multicorrs_g[i].src_pt_ = src_sub_pts_g->points[i];
	}
}

namespace gmnr {

	void Sampler::sample(pcl::PointCloud<pcl::PointXYZ>::Ptr in_pts, pcl::PointCloud<pcl::PointXYZ>::Ptr out_pts, pcl::PointIndices::Ptr out_indices) {
		pcl::PointIndices::Ptr in_indices(new pcl::PointIndices);
		for (int i = 0; i < in_pts->size(); i++) {
			in_indices->indices.push_back(i);
		}
		sample(in_pts, in_indices, out_pts, out_indices);
	}

	void RandomSampler::sample(pcl::PointCloud<pcl::PointXYZ>::Ptr in_pts, pcl::PointIndices::Ptr in_indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out_pts, pcl::PointIndices::Ptr out_indices) {
		
		int sample_size = sample_rate_ * in_indices->indices.size();
		if(out_pts) out_pts->resize(sample_size);
		if(out_indices) out_indices->indices.resize(sample_size);
		
		//srand((int)time(0));
		srand(0);
		//random sampling
		std::unordered_set<int> hashed_ids;
		int count = 0;
		while (count < sample_size) {
			int id = rand() % in_indices->indices.size();
			if(hashed_ids.find(id) == hashed_ids.end()) {
				hashed_ids.insert(id);
				if(out_pts) out_pts->points[count] = in_pts->points[ in_indices->indices[id] ];
				if(out_indices) out_indices->indices[count] = in_indices->indices[id];
				count++;
			}
		}
	}

	void BoundingBox::compute(pcl::PointCloud<pcl::PointXYZ>::Ptr pts) {

		min_pts_.x = pts->points[0].x;
		min_pts_.y = pts->points[0].y;
		min_pts_.z = pts->points[0].z;

		max_pts_ = min_pts_;

		for (int i = 1; i < pts->size(); i++) {
			min_pts_.x = min_pts_.x < pts->points[i].x ? min_pts_.x : pts->points[i].x;
			min_pts_.y = min_pts_.y < pts->points[i].y ? min_pts_.y : pts->points[i].y;
			min_pts_.z = min_pts_.z < pts->points[i].z ? min_pts_.z : pts->points[i].z;

			max_pts_.x = max_pts_.x > pts->points[i].x ? max_pts_.x : pts->points[i].x;
			max_pts_.y = max_pts_.y > pts->points[i].y ? max_pts_.y : pts->points[i].y;
			max_pts_.z = max_pts_.z > pts->points[i].z ? max_pts_.z : pts->points[i].z;
		}

		x_len_ = max_pts_.x - min_pts_.x;
		y_len_ = max_pts_.y - min_pts_.y;
		z_len_ = max_pts_.z - min_pts_.z;

		diag_len_ = x_len_ * x_len_ + y_len_ * y_len_ + z_len_ * z_len_;
		diag_len_ = sqrtf(diag_len_);
	}

	void BoundingBox::scale(float scale_factor) {
		float delta_x = x_len_ * (scale_factor - 1.0f) * 0.5f;
		float delta_y = y_len_ * (scale_factor - 1.0f) * 0.5f;
		float delta_z = z_len_ * (scale_factor - 1.0f) * 0.5f;

		min_pts_.x -= delta_x;
		min_pts_.y -= delta_y;
		min_pts_.z -= delta_z;

		max_pts_.x += delta_x;
		max_pts_.y += delta_y;
		max_pts_.z += delta_z;

		x_len_ *= scale_factor;
		y_len_ *= scale_factor;
		z_len_ *= scale_factor;
		diag_len_ *= scale_factor;
	}

	void GridSampler::sample(pcl::PointCloud<pcl::PointXYZ>::Ptr in_pts, pcl::PointIndices::Ptr in_indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out_pts, pcl::PointIndices::Ptr out_indices) {

		BoundingBox bbox;
		bbox.compute(in_pts);
		bbox.scale(1.05f);

		float min_x = bbox.min_pts_.x;
		float min_y = bbox.min_pts_.y;
		float min_z = bbox.min_pts_.z;

		std::stringstream ss;

		std::unordered_map<Voxel_Key, Voxel_Value> hashed_voxels;
		for (int i = 0; i < in_indices->indices.size(); i++) {
			int index_pt = in_indices->indices[i];
			pcl::PointXYZ pt = in_pts->points[index_pt];
			int id_x = (pt.x - min_x) / voxel_size_x_;
			int id_y = (pt.y - min_y) / voxel_size_y_;
			int id_z = (pt.z - min_z) / voxel_size_z_;

			ss.str("");
			ss << id_x << " " << id_y << " " << id_z;
			Voxel_Key key = ss.str();
			float vc_distance = ( pt.getVector3fMap() - voxel_center(min_x, min_y, min_z, id_x, id_y, id_z) ).norm();
			Voxel_Value value(i, vc_distance);
			std::unordered_map<Voxel_Key, Voxel_Value>::iterator it = hashed_voxels.find(key);
			if (it == hashed_voxels.end() || (*it).second.vc_distance_ > vc_distance) hashed_voxels[key] = value;
		}

		std::unordered_map<Voxel_Key, Voxel_Value>::iterator it = hashed_voxels.begin();
		while(it != hashed_voxels.end()) {
			int index_pt = in_indices->indices[(*it).second.id_];
			out_pts->push_back(in_pts->points[index_pt]);
			out_indices->indices.push_back(index_pt);
			it++;
		}
	}

	Eigen::Vector3f GridSampler::voxel_center(float min_x, float min_y, float min_z, int id_x, int id_y, int id_z) {
		Eigen::Vector3f vc;
		vc(0) = min_x + voxel_size_x_ * ( id_x + 0.5f );
		vc(1) = min_y + voxel_size_y_ * ( id_y + 0.5f );
		vc(2) = min_z + voxel_size_z_ * ( id_z + 0.5f );
		return vc;
	}

	void Correspondences::clear() {
		src_pts_->clear();
		tgt_pts_->clear();

		src_ids_->indices.clear();
		tgt_ids_->indices.clear();

		src_norms_->clear();
		tgt_norms_->clear();

		distance2s_.clear();
		stabiblites_.clear();
		sub_levels_.clear();
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
		bool accept_nonoverlap_corrs = params.accept_nonoverlap_corrs_;

		pcl::PointCloud<pcl::PointXYZ>::Ptr source_pts(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::Normal>::Ptr source_norms(new pcl::PointCloud<pcl::Normal>);

		tf->transform(source_data->pts_, src_indices, source_pts, source_data->norms_, source_norms);

		if (visualizer_g) {
			visualizer_g->addPointCloud(source_pts, "source_pts");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "source_pts");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "source_pts");

			visualizer_g->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(source_pts, source_norms, 1, 0.003f, "source_norms");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "source_norms");

			visualizer_g->spinOnce(spin_time_g, true);

			visualizer_g->removePointCloud("source_pts");
			visualizer_g->removePointCloud("source_norms");

		}

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
				stabiblites_.push_back(1.0f);
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

		if(visualizer_g) {

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
			visualizer_g->addPointCloud(src_pts_, "src_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "src_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "src_pts_");

			visualizer_g->addPointCloud(tgt_pts_, "tgt_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0, "tgt_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "tgt_pts_");

			visualizer_g->addCorrespondences<pcl::PointXYZ>(src_pts_, tgt_pts_, s_t_corrs, "s_t_correspondences");
			visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_t_correspondences");
			visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.0, "s_t_correspondences");

			visualizer_g->addCorrespondences<pcl::PointXYZ>(src_pts_, source_data->pts_, s_s_corrs, "s_s_correspondences");
			visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_s_correspondences");
			visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 1.0, "s_s_correspondences");

			visualizer_g->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(src_pts_, src_norms_, 1, 0.003f, "src_norms_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "src_norms_");

			visualizer_g->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(tgt_pts_, tgt_norms_, 1, 0.003f, "tgt_norms_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "tgt_norms_");

			visualizer_g->spinOnce(spin_time_g, true);

			visualizer_g->removePointCloud("src_pts_");
			visualizer_g->removePointCloud("tgt_pts_");
			visualizer_g->removeCorrespondences("s_t_correspondences");
			visualizer_g->removeCorrespondences("s_s_correspondences");
			visualizer_g->removePointCloud("src_norms_");
			visualizer_g->removePointCloud("tgt_norms_");
		}
	}

	void Correspondences::extract_larger_than_deviation(pcl::PointIndices::Ptr extracted_src_ids, float deviation_scale) {

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

	void Correspondences::remove_larger_than_deviation(float deviation_scale) {

		float sum_distance2 = 0.0f;
		for (int i = 0; i < distance2s_.size(); i++) {
			sum_distance2 += distance2s_[i];
		}
		float deviation_distance2 = sum_distance2 / distance2s_.size();
		float scaled_deviation_distance2 = deviation_distance2 * deviation_scale;

		pcl::PointCloud<pcl::PointXYZ>::iterator src_pts_it = src_pts_->begin();
		pcl::PointCloud<pcl::PointXYZ>::iterator tgt_pts_it = tgt_pts_->begin();
		std::vector<int>::iterator src_ids_it = src_ids_->indices.begin();
		std::vector<int>::iterator tgt_ids_it = tgt_ids_->indices.begin();
		pcl::PointCloud<pcl::Normal>::iterator src_norms_it = src_norms_->begin();
		pcl::PointCloud<pcl::Normal>::iterator tgt_norms_it = tgt_norms_->begin();
		std::vector<float>::iterator distance2s_it = distance2s_.begin();
		std::vector<float>::iterator stabiblites_it = stabiblites_.begin();
		std::vector<int>::iterator sub_levels_it = sub_levels_.begin();

		while(distance2s_it != distance2s_.end()) {
			if (*distance2s_it > scaled_deviation_distance2) {
				src_ids_it = src_ids_->indices.erase(src_ids_it);
				tgt_ids_it = tgt_ids_->indices.erase(tgt_ids_it);		
				src_pts_it = src_pts_->erase(src_pts_it);
				tgt_pts_it = tgt_pts_->erase(tgt_pts_it);
				src_norms_it = src_norms_->erase(src_norms_it);
				tgt_norms_it = tgt_norms_->erase(tgt_norms_it);
				distance2s_it = distance2s_.erase(distance2s_it);
				stabiblites_it = stabiblites_.erase(stabiblites_it);
				sub_levels_it = sub_levels_.erase(sub_levels_it);
				continue;
			}
			src_ids_it++;
			tgt_ids_it++;
			src_pts_it++;
			tgt_pts_it++;
			src_norms_it++;
			tgt_norms_it++;
			distance2s_it++;
			stabiblites_it++;
			sub_levels_it++;
		}

	}

	void Transform::transform(pcl::PointCloud<pcl::PointXYZ>::Ptr in, pcl::PointIndices::Ptr indices, pcl::PointCloud<pcl::PointXYZ>::Ptr out, 
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

	Transform::Ptr ICP::do_icp(Data_Pack::Ptr source_data, Transform::Ptr init_tf) {
		pcl::PointIndices::Ptr src_indices(new pcl::PointIndices);
		for (int i = 0; i < source_data->pts_->size(); i++) {
			src_indices->indices.push_back(i);
		}
		return do_icp(source_data, src_indices, init_tf);
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
				//if (tgt_norms) {
				//	NormalSet3D tnorms(tgt_norms->size(), 3);
				//	for (int i = 0; i < tgt_norms->size(); i++) tnorms.row(i) = tgt_norms->points[i].getNormalVector3fMap().cast<gmnr::Scalar>();
				//	PointToPlaneLinear point_to_plane(src, tgt, tnorms, true);
				//	return LinearTransform::Ptr(new LinearTransform(point_to_plane.transfomation().cast<float>()));
				//}
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
			corrs->remove_larger_than_deviation(2.0f);
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
			static float delta = 3e-3f;
			out_norms->resize(in_norms->size());
			for (int i = 0; i < in_norms->size(); i++) {
				pcl::Normal &in_norm = in_norms->points[i];
				pcl::Normal &out_norm = out_norms->points[i];

				out_norm.getNormalVector3fMap() = in_norm.getNormalVector3fMap();

				//Eigen::Vector3f pt_delta = in->points[i].getVector3fMap() + in_norm.getNormalVector3fMap() * delta;
				//Eigen::Vector3f pt_delta_tps = (tps_.evaluate(pt_delta.transpose().cast<gmnr::Scalar>())).transpose().cast<float>();

				//out_norm.getNormalVector3fMap() = (pt_delta_tps - out->points[i].getVector3fMap()).normalized(); 
				//if (out_norm.getNormalVector3fMap().dot(in_norm.getNormalVector3fMap()) < 0 ) out_norm.getNormalVector3fMap() *= -1.0f;
			}	
		}
	}

	TPSTransform::Ptr  ThinPlateSplinesICP::compute_tps(Correspondences::Ptr corrs, float lambda, float kappa, int n) {

		pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = corrs->src_pts_;
		pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_pts = corrs->tgt_pts_;

		PointSet3D src(src_pts->size(), 3), tgt(tgt_pts->size(), 3);
		for (int i = 0; i < src_pts->size(); i++) src.row(i) = src_pts->points[i].getVector3fMap().cast<gmnr::Scalar>();
		for (int i = 0; i < tgt_pts->size(); i++) tgt.row(i) = tgt_pts->points[i].getVector3fMap().cast<gmnr::Scalar>();

		std::cerr << "lambda = " << lambda << " kappa = " << kappa << std::endl;


		if (visualizer_g) {
			int corrs_size = n > 0 ? n : corrs->src_pts_->size();
			pcl::Correspondences s_t_corrs(corrs_size);
			for (int i = 0; i < corrs_size; i++) {
				s_t_corrs[i].index_query = i;
				s_t_corrs[i].index_match = i;
			}
			visualizer_g->addCorrespondences<pcl::PointXYZ>(src_pts, tgt_pts, s_t_corrs, "s_t_correspondences");
			visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_t_correspondences");
			visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.0, "s_t_correspondences");
			visualizer_g->spinOnce(spin_time_g, true);
			visualizer_g->removeCorrespondences("s_t_correspondences");
		}

		if(n > 0) return TPSTransform::Ptr(new TPSTransform(ApproxiTPSFunction(src, tgt, lambda, kappa, n)));
		return TPSTransform::Ptr(new TPSTransform(TPSFunction(src, tgt, lambda, kappa)));
	}

	Transform::Ptr ThinPlateSplinesICP::do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf) {

		Correspondences::Ptr corrs(new Correspondences);

		TPSTransform::Ptr tps_tf(new TPSTransform);
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
			corrs->remove_larger_than_deviation(2.0f);

			std::vector<int> &src_ids = corrs->src_ids_->indices;
			pcl::PointCloud<pcl::PointXYZ>::Ptr src_pts = corrs->src_pts_;
			for (int i = 0; i < src_ids.size(); i++) src_pts->points[i] = source_data->pts_->points[ src_ids[i] ];

			//tps_tf = compute_tps(corrs, lambda, kappa, (int)(corrs->src_pts_->size() * frac * (iter+1)));
			tps_tf = compute_tps(corrs, lambda, kappa);
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

		if(visualizer_g) {
			pcl::PointCloud<pcl::PointXYZ>::Ptr extracted_pts(new pcl::PointCloud<pcl::PointXYZ>);
			pcl::ExtractIndices<pcl::PointXYZ> pts_extractor;
			pts_extractor.setIndices(seed_indices);
			pts_extractor.setInputCloud(pts);
			pts_extractor.filter(*extracted_pts);

			visualizer_g->addPointCloud(extracted_pts, "large_deviation_src_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "large_deviation_src_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "large_deviation_src_pts_");	
			visualizer_g->spinOnce(spin_time_g);	
			visualizer_g->removePointCloud("large_deviation_src_pts_");

			pts_extractor.setIndices(extend_indices);
			pts_extractor.filter(*extracted_pts);

			visualizer_g->addPointCloud(extracted_pts, "extended_large_deviation_src_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "extended_large_deviation_src_pts_");
			visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "extended_large_deviation_src_pts_");	
			visualizer_g->spinOnce(spin_time_g);	
			visualizer_g->removePointCloud("extended_large_deviation_src_pts_");
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
			combined_corrs->stabiblites_.push_back(corrs_cluster[i]->stabiblites_[j]);
			combined_corrs->sub_levels_.push_back(corrs_cluster[i]->stabiblites_[j]);
		}
	}

	Transform::Ptr LinearThinPlateSplinesICP::do_icp(Data_Pack::Ptr source_data, pcl::PointIndices::Ptr src_indices, Transform::Ptr init_tf) {

		Correspondences::Ptr corrs(new Correspondences);
		pcl::PointCloud<pcl::PointXYZ>::Ptr source_pts(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::Normal>::Ptr source_norms(new pcl::PointCloud<pcl::Normal>);

		Transform::Ptr lineartps_tf(new LinearTPSTransform);
		if (init_tf) {
			boost::dynamic_pointer_cast<LinearTransform>(lineartps_tf)->mat_ = boost::dynamic_pointer_cast<LinearTransform>(init_tf)->mat_;
			boost::dynamic_pointer_cast<TPSTransform>(lineartps_tf)->tps_ = boost::dynamic_pointer_cast<TPSTransform>(init_tf)->tps_;
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

			ICP::Ptr linear_icp(new LinearICP(target_data_));
			boost::dynamic_pointer_cast<LinearICP>(linear_icp)->linear_params_ = linear_params;
			boost::dynamic_pointer_cast<LinearICP>(linear_icp)->corrs_params_ = corrs_params_;
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
				if(part_indices->indices.size() < 5) continue;
				pcl::PointIndices::Ptr extend_indices(new pcl::PointIndices);
				extend_pointcloud(source_pts, source_kd, part_indices, extend_width, extend_indices);
				if (extend_indices->indices.size() < 20) continue;
				Transform::Ptr part_tf = linear_icp->do_icp(compact_source_data, extend_indices);
				Correspondences::Ptr part_corrs(new Correspondences);
				part_corrs->compute_correspondences(compact_source_data, extend_indices, target_data_, part_tf, corrs_params_);
				part_corrs->remove_larger_than_deviation(2.0f);
				corrs_cluster.push_back(part_corrs);
			}
			
			corrs_cluster.push_back(corrs);
			Correspondences::Ptr combined_corrs(new Correspondences);
			combine_corrs_cluster(corrs_cluster, combined_corrs);

			combined_corrs->remove_larger_than_deviation(4.0f);

			if(visualizer_g) {

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
				visualizer_g->addPointCloud(combined_corrs->src_pts_, "src_pts_");
				visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "src_pts_");
				visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "src_pts_");

				visualizer_g->addPointCloud(combined_corrs->tgt_pts_, "tgt_pts_");
				visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0, "tgt_pts_");
				visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "tgt_pts_");

				visualizer_g->addCorrespondences<pcl::PointXYZ>(combined_corrs->src_pts_, combined_corrs->tgt_pts_, s_t_corrs, "s_t_correspondences");
				visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_t_correspondences");
				visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.0, "s_t_correspondences");

				visualizer_g->addCorrespondences<pcl::PointXYZ>(combined_corrs->src_pts_, source_data->pts_, s_s_corrs, "s_s_correspondences");
				visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, "s_s_correspondences");
				visualizer_g->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 1.0, "s_s_correspondences");

				visualizer_g->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(combined_corrs->src_pts_, combined_corrs->src_norms_, 1, 0.003f, "src_norms_");
				visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "src_norms_");

				visualizer_g->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(combined_corrs->tgt_pts_, combined_corrs->tgt_norms_, 1, 0.003f, "tgt_norms_");
				visualizer_g->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "tgt_norms_");

				visualizer_g->spinOnce(spin_time_g + 200000, true);

				visualizer_g->removePointCloud("src_pts_");
				visualizer_g->removePointCloud("tgt_pts_");
				visualizer_g->removeCorrespondences("s_t_correspondences");
				visualizer_g->removeCorrespondences("s_s_correspondences");
				visualizer_g->removePointCloud("src_norms_");
				visualizer_g->removePointCloud("tgt_norms_");
			}

			Transform::Ptr linear_tf;

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
				boost::dynamic_pointer_cast<Transform>(linear_tf)->transform(source_data->pts_, src_indices, source_pts);
				for (int i = 0; i < src_ids.size(); i++) src_pts->points[i] = source_pts->points[ src_ids[i] ];
			}

			//Transform::Ptr tps_tf = ThinPlateSplinesICP::compute_tps(combined_corrs, lambda, kappa, (int)(combined_corrs->src_pts_->size() * frac * (iter+1)));
			Transform::Ptr tps_tf = ThinPlateSplinesICP::compute_tps(combined_corrs, lambda, kappa);

			boost::dynamic_pointer_cast<LinearTransform>(lineartps_tf)->mat_ = boost::dynamic_pointer_cast<LinearTransform>(linear_tf)->mat_;
			boost::dynamic_pointer_cast<TPSTransform>(lineartps_tf)->tps_ = boost::dynamic_pointer_cast<TPSTransform>(tps_tf)->tps_;

			kappa *= iter_rate;
			iter++;
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

	Data_Pack::Ptr source_data(new Data_Pack(src_sub_pts_g, 0, src_sub_norms_g, 0));

	for (tgt_mesh_id_g = 0; tgt_mesh_id_g < tgt_fmeshs_g.size(); tgt_mesh_id_g++) {
		if(tgt_mesh_id_g == src_mesh_id_g || linked_tgt_ids_g.find(tgt_mesh_id_g) == linked_tgt_ids_g.end()) continue;

		tgt_fmesh_g = tgt_fmeshs_g[tgt_mesh_id_g];
		get_file_path_name_extension(tgt_fmesh_g, tgt_path_g, tgt_name_g, tgt_extension_g);

		load_tgt_data();

		Data_Pack::Ptr target_data(new Data_Pack(tgt_pts_g, tgt_kdtree_g, tgt_norms_g, tgt_bd_g));

		ICP::Ptr icp;
		if(icp_type_g == "linear") {
			LinearICP::Ptr linear_icp(new LinearICP(target_data));
			linear_icp->linear_params_ = LinearICP::Parameters(linear_type_g, linear_max_iter_g);
			linear_icp->corrs_params_ = Correspondences::Parameters(dist_thresh_g, angle_thresh_g);
			icp = linear_icp;
		} else if(icp_type_g == "tps") {
			ThinPlateSplinesICP::Ptr tps_icp(new ThinPlateSplinesICP(target_data));
			tps_icp->tps_params_ = ThinPlateSplinesICP::Parameters(lambda_g, kappa_g, tps_max_iter_g, iter_rate_g);
			tps_icp->corrs_params_ = Correspondences::Parameters(dist_thresh_g, angle_thresh_g);
			icp = tps_icp;
		} else if(icp_type_g == "lineartps") {
			LinearThinPlateSplinesICP::Ptr lineartps_icp(new LinearThinPlateSplinesICP(target_data));
			lineartps_icp->lineartps_params_ = LinearThinPlateSplinesICP::Parameters(linear_type_g, linear_max_iter_g, lambda_g, kappa_g, tps_max_iter_g, iter_rate_g, ec_tolerance_g);
			lineartps_icp->corrs_params_ = Correspondences::Parameters(dist_thresh_g, angle_thresh_g);
			icp = lineartps_icp;
		} else {
			std::cerr << "unknown icp_type" << std::endl;
			return -1;
		}

		Transform::Ptr tf = icp->do_icp(source_data);

		if (visualizer_g) {
			visualizer_g->spinOnce(spin_time_g);
			visualizer_g->removePolygonMesh("target_mesh");
		}

		if (write_warped_src_g) {
			std::string output_name = output_directory_g + "/" + src_name_g + "-" + tgt_name_g + "_" + icp->get_classname() + ".ply";
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZ>);
			tf->transform(src_pts_g, cloud_temp);
			pcl::copyPointCloud(*cloud_temp, *src_original_g);
			pcl::toPCLPointCloud2(*src_original_g, src_mesh_g->cloud);

			if(pcl::io::savePLYFile(output_name, *src_mesh_g) != 0) {
				PCL_ERROR ("Couldn't write file %s \n", output_name);
				exit (-1);
			}
		}
	}

	std::cout << src_mesh_id_g << " " << multicorrs_g.size() << " " << src_fmesh_g << std::endl;

	for (int i = 0; i < multicorrs_g.size(); i++) {
		MultiCorrespondence &multicorr = multicorrs_g[i];
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




