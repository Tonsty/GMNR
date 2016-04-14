/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2012, Willow Garage, Inc.
 * 
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id: $
 *
 * @author: Koen Buys - KU Leuven
 */

#include "manual_registration.h"

//PCL
#include <pcl/registration/icp.h>
#include <pcl/registration/transformation_estimation_svd_scale.h>
#include <pcl/registration/correspondence_rejection_surface_normal.h>
#include <pcl/registration/correspondence_rejection_one_to_one.h>

//QT4
#include <QtGui/QApplication>
#include <QtCore/QMutexLocker>
#include <QtCore/QEvent>
#include <QtCore/QObject>

// VTK
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkCamera.h>

using namespace pcl;
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
ManualRegistration::ManualRegistration (QWidget *parent): QMainWindow(parent)
{
  //Create a timer
  // vis_timer_ = new QTimer (this);
  // vis_timer_->start (5);//5ms

  // connect (vis_timer_, SIGNAL (timeout ()), this, SLOT (timeoutSlot()));

  ui_ = new Ui_Manual_Registration_MainWindow;
  ui_->setupUi (this);

  
  this->setWindowTitle ("PCL Manual Registration");

  // Set up the source window
  vis_src_.reset (new pcl::visualization::PCLVisualizer ("", false));
  //vis_src_->setSize ( 500, 500 );
  //vis_src_->resetCamera();
  //vis_src_->resetCameraViewpoint( "cloud_src_" );
  //vis_src_->addCoordinateSystem( 0.2, "reference_src" );
  
  ui_->qvtk_widget_src->SetRenderWindow (vis_src_->getRenderWindow ());
  vis_src_->setupInteractor (ui_->qvtk_widget_src->GetInteractor (), ui_->qvtk_widget_src->GetRenderWindow ());
  vis_src_->getInteractorStyle ()->setKeyboardModifier (pcl::visualization::INTERACTOR_KB_MOD_SHIFT);
  ui_->qvtk_widget_src->update ();

  vis_src_->registerPointPickingCallback (&ManualRegistration::SourcePointPickCallback, *this);
  //vis_src_->setBackgroundColor(1.0f, 1.0f, 1.0f);

  // Set up the destination window
  vis_dst_.reset (new pcl::visualization::PCLVisualizer ("", false));
  //vis_dst_->setSize( 500, 500 );
  //vis_dst_->resetCamera();
  //vis_dst_->resetCameraViewpoint( "cloud_dst_" );
  //vis_dst_->addCoordinateSystem( 0.2, "reference_dst" );
  
  ui_->qvtk_widget_dst->SetRenderWindow (vis_dst_->getRenderWindow ());
  vis_dst_->setupInteractor (ui_->qvtk_widget_dst->GetInteractor (), ui_->qvtk_widget_dst->GetRenderWindow ());
  vis_dst_->getInteractorStyle ()->setKeyboardModifier (pcl::visualization::INTERACTOR_KB_MOD_SHIFT);
  ui_->qvtk_widget_dst->update ();

  vis_dst_->registerPointPickingCallback (&ManualRegistration::DstPointPickCallback, *this);
  //vis_dst_->setBackgroundColor(1.0f, 1.0f, 1.0f);

  // Set up the window
  vis_.reset (new pcl::visualization::PCLVisualizer ("", false));
  //vis_->setSize ( 500, 500 );
  //vis_->resetCamera();
  //vis_->resetCameraViewpoint( "cloud_" );
  //vis_->addCoordinateSystem( 0.2, "reference" );

  ui_->qvtk_widget->SetRenderWindow (vis_->getRenderWindow ());
  vis_->setupInteractor (ui_->qvtk_widget->GetInteractor (), ui_->qvtk_widget->GetRenderWindow ());
  vis_->getInteractorStyle ()->setKeyboardModifier (pcl::visualization::INTERACTOR_KB_MOD_SHIFT);
  ui_->qvtk_widget->update ();

  //vis_src_->setBackgroundColor(1.0f, 1.0f, 1.0f);

  // Connect all buttons
  connect (ui_->confirmSrcPointButton, SIGNAL(clicked()), this, SLOT(confirmSrcPointPressed()));
  connect (ui_->confirmDstPointButton, SIGNAL(clicked()), this, SLOT(confirmDstPointPressed()));
  connect (ui_->calculateButton, SIGNAL(clicked()), this, SLOT(calculatePressed()));
  connect (ui_->clearButton, SIGNAL(clicked()), this, SLOT(clearPressed()));
  connect (ui_->orthoButton, SIGNAL(stateChanged(int)), this, SLOT(orthoChanged(int)));
  connect (ui_->applyTransformButton, SIGNAL(clicked()), this, SLOT(applyTransformPressed()));
  connect (ui_->refineButton, SIGNAL(clicked()), this, SLOT(refinePressed()));
  connect (ui_->undoButton, SIGNAL(clicked()), this, SLOT(undoPressed()));
  connect (ui_->safeButton, SIGNAL(clicked()), this, SLOT(safePressed()));

  cloud_src_modified_ = true; // first iteration is always a new pointcloud
  cloud_dst_modified_ = true;

  ball_radius_src_ = 1.0f;
  ball_radius_dst_ = 1.0f;

  transform_ = Eigen::Matrix4f::Identity();
}

void ManualRegistration::SourcePointPickCallback (const pcl::visualization::PointPickingEvent& event, void*)
{
  // Check to see if we got a valid point. Early exit.
  int idx = event.getPointIndex ();
  if (idx == -1)
    return;

  // Get the point that was picked
  event.getPoint (src_point_.x, src_point_.y, src_point_.z);
  PCL_INFO ("Src Window: Clicked point %d with X:%f Y:%f Z:%f\n", idx, src_point_.x, src_point_.y, src_point_.z);
  src_point_selected_ = true;

  vis_src_->removeShape( "point_src");
  vis_src_->addSphere( src_point_, ball_radius_src_, 1.0, 0.0, 0.0, "point_src");
}

void ManualRegistration::DstPointPickCallback (const pcl::visualization::PointPickingEvent& event, void*)
{
  // Check to see if we got a valid point. Early exit.
  int idx = event.getPointIndex ();
  if (idx == -1)
    return;

  // Get the point that was picked
  event.getPoint (dst_point_.x, dst_point_.y, dst_point_.z);
  PCL_INFO ("Dst Window: Clicked point %d with X:%f Y:%f Z:%f\n", idx, dst_point_.x, dst_point_.y, dst_point_.z);
  dst_point_selected_ = true;

  vis_dst_->removeShape( "point_dst");
  vis_dst_->addSphere( dst_point_, ball_radius_dst_, 0.0, 0.0, 1.0, "point_dst");
}

void ManualRegistration::confirmSrcPointPressed()
{
  if(src_point_selected_)
  {
    src_pc_.push_back(src_point_);
    PCL_INFO ("Selected %d source points\n", src_pc_.size());
    src_point_selected_ = false;
    src_pc_.width = src_pc_.points.size();
    stringstream sstream;
    sstream << src_pc_.size();
    vis_src_->removeShape( "point_src");
    vis_src_->addSphere( src_point_, ball_radius_src_, 1.0, 0.0, 0.0, "point_" + sstream.str());
  }
  else
  {
    PCL_INFO ("Please select a point in the source window first\n");
  }
}

void ManualRegistration::confirmDstPointPressed()
{
  if(dst_point_selected_)
  {
    dst_pc_.push_back(dst_point_);
    PCL_INFO ("Selected %d destination points\n", dst_pc_.size());
    dst_point_selected_ = false;
    dst_pc_.width = dst_pc_.points.size();
    stringstream sstream;
    sstream << dst_pc_.size();
    vis_dst_->removeShape( "point_dst");
    vis_dst_->addSphere( dst_point_, ball_radius_dst_, 0.0, 0.0, 1.0, "point_" + sstream.str());
  }
  else
  {
    PCL_INFO ("Please select a point in the destination window first\n");
  }
}

void ManualRegistration::calculatePressed()
{
  if(dst_pc_.size() != src_pc_.size())
  {
    PCL_INFO ("You haven't selected an equal amount of points, please do so\n");
    return;
  }
  if(dst_pc_.size() < 3)
  {
	  PCL_INFO ("You haven't selected an enough amount of points, please do so\n");
	  return;
  }
  pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ> tfe;
  tfe.estimateRigidTransformation(src_pc_, dst_pc_, transform_);

  std::cout << "Transform : " << std::endl << transform_ << std::endl;
}

void ManualRegistration::clearPressed()
{
	dst_point_selected_ = false;
	src_point_selected_ = false;
	src_pc_.clear();
	dst_pc_.clear();
	src_pc_.height = 1; src_pc_.width = 0;
	dst_pc_.height = 1; dst_pc_.width = 0;

	vis_src_->removeAllShapes();
	vis_dst_->removeAllShapes();

	ui_->qvtk_widget_src->update();
	ui_->qvtk_widget_dst->update();
}

void ManualRegistration::orthoChanged (int state)
{
  PCL_INFO ("Ortho state %d\n", state);
  if(state == 0) // Not selected
  {
    vis_src_->getRenderWindow ()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetParallelProjection(0);
    vis_dst_->getRenderWindow ()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetParallelProjection(0);
	vis_->getRenderWindow ()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetParallelProjection(0);
  }
  if(state == 2) // Selected
  {
    vis_src_->getRenderWindow ()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetParallelProjection(1);
    vis_dst_->getRenderWindow ()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetParallelProjection(1);
	vis_->getRenderWindow ()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetParallelProjection(1);
  }
  ui_->qvtk_widget_src->update();
  ui_->qvtk_widget_dst->update();
  ui_->qvtk_widget->update();
}

//TODO
void ManualRegistration::applyTransformPressed()
{
  pcl::PointCloud<PointT>::Ptr cloud_src_dup = cloud_src_->makeShared();
  pcl::transformPointCloudWithNormals(*cloud_src_, *cloud_src_dup, transform_);

  pcl::visualization::PointCloudColorHandlerRGBField<PointT> rgb_handler(cloud_src_);
  vis_->updatePointCloud(cloud_src_dup, rgb_handler, "cloud_src_");  
  ui_->qvtk_widget->update();
}

void ManualRegistration::refinePressed()
{
  pcl::IterativeClosestPoint<PointT, PointT> icp;

  icp.setInputSource(cloud_src_);
  icp.setInputTarget(cloud_dst_);

  pcl::registration::TransformationEstimationSVDScale<PointT, PointT>::Ptr est(new pcl::registration::TransformationEstimationSVDScale<PointT, PointT>);
  icp.setTransformationEstimation(est);

  icp.setMaxCorrespondenceDistance(0.002);

  //pcl::registration::CorrespondenceRejectorSurfaceNormal::Ptr normal_rejector(new pcl::registration::CorrespondenceRejectorSurfaceNormal);
  //icp.addCorrespondenceRejector(normal_rejector);

  pcl::registration::CorrespondenceRejectorOneToOne::Ptr one_to_one_rejector(new pcl::registration::CorrespondenceRejectorOneToOne);
  icp.addCorrespondenceRejector(one_to_one_rejector);

  //icp.getCorrespondenceRejectors()[0]->getClassName()
 
  pcl::PointCloud<PointT> cloud_final;
  icp.align(cloud_final, transform_);
  transform_ = icp.getFinalTransformation();
}

void ManualRegistration::undoPressed()
{
}

void ManualRegistration::safePressed()
{
}

void ManualRegistration::clearSrcVis()
{
  src_point_selected_ = false;
  src_pc_.clear();
  src_pc_.height = 1; src_pc_.width = 0;
  vis_src_->removeAllPointClouds();
  vis_src_->removeAllShapes();
  ui_->qvtk_widget_src->update();
}

void ManualRegistration::clearDstVis()
{
  dst_point_selected_ = false;
  dst_pc_.clear();
  dst_pc_.height = 1; dst_pc_.width = 0;
  vis_dst_->removeAllPointClouds();
  vis_dst_->removeAllShapes();
  ui_->qvtk_widget_dst->update();
}

void ManualRegistration::showSrcCloud()
{
  //vis_src_->addPointCloud<PointT>(cloud_src_, "cloud_src_");
  pcl::visualization::PointCloudColorHandlerRGBField<PointT> rgb_handler(cloud_src_);
  vis_src_->addPointCloud<PointT>(cloud_src_, rgb_handler, "cloud_src_");
  Eigen::Vector4f centroid;
  pcl::compute3DCentroid(*cloud_src_, centroid);
  PointT min_pt, max_pt;
  pcl::getMinMax3D(*cloud_src_, min_pt, max_pt);
  ball_radius_src_ = ( max_pt.getVector3fMap() - min_pt.getVector3fMap() ).norm() / 200;
  vis_src_->setCameraPosition(0, 0, 0, centroid(0), centroid(1), centroid(2), 0, 1, 0);  
}

void ManualRegistration::showDstCloud()
{
  //vis_dst_->addPointCloud<PointT>(cloud_dst_, "cloud_dst_");
  pcl::visualization::PointCloudColorHandlerRGBField<PointT> rgb_handler(cloud_dst_);
  vis_dst_->addPointCloud<PointT>(cloud_dst_, rgb_handler, "cloud_dst_");  
  Eigen::Vector4f centroid;
  pcl::compute3DCentroid(*cloud_dst_, centroid);
  PointT min_pt, max_pt;
  pcl::getMinMax3D(*cloud_dst_, min_pt, max_pt);
  ball_radius_dst_ = ( max_pt.getVector3fMap() - min_pt.getVector3fMap() ).norm() / 200;
  vis_dst_->setCameraPosition(0, 0, 0, centroid(0), centroid(1), centroid(2), 0, 1, 0);  
}

void ManualRegistration::showCloud()
{
	pcl::visualization::PointCloudColorHandlerRGBField<PointT> rgb_handler_src(cloud_src_), rgb_handler_dst(cloud_dst_);
	vis_->addPointCloud<PointT>(cloud_src_, rgb_handler_src, "cloud_src_");  
	vis_->addPointCloud<PointT>(cloud_dst_, rgb_handler_dst, "cloud_dst_");
	Eigen::Vector4f centroid_src;
	pcl::compute3DCentroid(*cloud_src_, centroid_src);
	Eigen::Vector4f centroid_dst;
	pcl::compute3DCentroid(*cloud_dst_, centroid_dst);
	Eigen::Vector4f centroid = (centroid_src + centroid_dst) * 0.5f;

	vis_->setCameraPosition(0, 0, 0, centroid(0), centroid(1), centroid(2), 0, 1, 0);  
}



