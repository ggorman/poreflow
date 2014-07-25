/*  Copyright (C) 2010 Imperial College London and others.
 *
 *  Please see the AUTHORS file in the main source directory for a
 *  full list of copyright holders.
 *
 *  Gerard Gorman
 *  Applied Modelling and Computation Group
 *  Department of Earth Science and Engineering
 *  Imperial College London
 *
 *  g.gorman@imperial.ac.uk
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials provided
 *  with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 *  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 *  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 *  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 */

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkTIFFReader.h>
#include <vtkImageData.h>
#include <vtkTriangleFilter.h>
#include <vtkImageToPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>

#include <algorithm>
#include <vector>
#include <set>

#include <cassert>
#include <cstdlib>

#include "CTImage.h"

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

CTImage::CTImage(){
  verbose = false;
  for(int i=0;i<3;i++)
    dims[i] = -1;
  resolution=1.0;
  image = NULL;
  raw_image = NULL;
  domain = NULL;
}

CTImage::~CTImage(){
  if(image!=NULL)
    delete image;
  if(raw_image!=NULL)
    delete raw_image;
  if(domain!=NULL)
    delete domain;
}

void CTImage::verbose_on(){
  verbose = true;
}

size_t CTImage::get_NNodes(){
  size_t NNodes = xyz.size()/3;
  if(verbose)
    std::cout<<"size_t get_NNodes() = "<<NNodes<<std::endl;
  return NNodes;
}

size_t CTImage::get_NElements(){
  size_t NElements = tets.size()/4;
  if(verbose)
    std::cout<<"size_t get_NElements() = "<<NElements<<std::endl;
  return NElements;
}

size_t CTImage::get_NFacets(){
  size_t NFacets = facets.size()/3;
  if(verbose)
    std::cout<<"size_t get_NFacets() = "<<NFacets<<std::endl;
  return NFacets;
}

int CTImage::read_raw_ese_image(const char *name, int slab_size){
  if(verbose)
    std::cout<<"int read_raw_ese_image(char *name, int slab_size)"<<std::endl;

  boost::filesystem::path image_dir(name);  
  stem = image_dir.stem();

  bool found_metadata=false;
  {
    boost::filesystem::path image_info = image_dir/std::string(stem.string()+"_info.txt");
    if(boost::filesystem::exists(image_info)){
      found_metadata=true;

      // Get the metadata.
      std::ifstream info;
      info.open(image_info.string().c_str());
      while(!info.eof()){
        std::string buffer;
        info>>buffer;
        if(buffer==std::string("Size:")){
          for(int i=0;i<3;i++)
            info>>dims[i];
        }else if(buffer==std::string("Resolution:")){
          info>>resolution;
          resolution*=1.0e-6;
        }
        if(slab_size>dims[0]){
          std::cerr<<"WARNING: slab_size larger than original image\n";
        }
      }
      info.close();
    }
  }

  if(!found_metadata){
    boost::filesystem::path image_info = image_dir/std::string(stem.string()+".dat");
    if(boost::filesystem::exists(image_info)){
      found_metadata=true;

      // Get the metadata.
      std::ifstream info;
      info.open(image_info.string().c_str());
      std::string buffer;
      info>>buffer;
      info>>buffer;
      for(int i=0;i<3;i++)
        info>>dims[i];

      double l, r;
      info>>l; info>>r;
      resolution = (r-l)/dims[0]*1.0e-6;

      info.close();
    }
    if(slab_size>dims[0]){
      std::cerr<<"WARNING: slab_size larger than original image\n";
    }
  }

  if(!found_metadata){
    std::cerr<<"ERROR: Could not find meta-data for CT-Image.";
    exit(-1);
  }

  // Read image.
  boost::filesystem::path image_filename = image_dir/std::string(stem.string()+".raw");

  image_size = dims[0]*dims[1]*dims[2];
  raw_image = new unsigned char[image_size];

  std::ifstream image_file;
  image_file.open(image_filename.string().c_str(), std::ios::binary);
  image_file>>raw_image;
  image_file.close();

  if(slab_size>0){
    int image_size_new = slab_size*slab_size*slab_size;
    unsigned char *raw_image_new = new unsigned char[image_size_new];
    int slice_size_new = slab_size*slab_size;
    int slice_size = dims[0]*dims[0];
    for(int i=0;i<slab_size;i++){
      for(int j=0;j<slab_size;j++){
	memcpy(raw_image_new+i*slice_size_new+j*slab_size, raw_image+i*slice_size+j*dims[0], slab_size);
      }
    }
    delete [] raw_image;
    raw_image = raw_image_new;
    for(int i=0;i<3;i++)
      dims[i] = slab_size;

    image_size = image_size_new;
  }

  for(int i=0;i<image_size;i++)
    if(raw_image[i]==0)
      raw_image[i] = 1;
    else
      raw_image[i] = 0;

  double spacing[] = {1,1,1};

  image = new CGAL::Image_3(_createImage(dims[0], dims[1], dims[2], 1,
	spacing[0], spacing[1], spacing[2],
	1, WK_FIXED, SGN_UNSIGNED)); 
  ImageIO_free(image->data());
  image->set_data((void*)raw_image); 

  // Domain
  domain = new Mesh_domain(*image);
}

int CTImage::create_hourglass(int size, int throat_width){
  image_size = size*size*size;
  raw_image = new unsigned char[image_size];
  for(int i=0;i<3;i++)
    dims[i] = size;
  resolution=1.0/size;

  double dx = 1.0/size;
  double A = (0.5 - throat_width*0.5*dx)*0.5;
  double offset = A + dx*throat_width*0.5;
  double two_pi = 2*3.14159265359;
  size_t pos = 0;
  for(size_t i=0;i<size;i++){
    for(size_t j=0;j<size;j++){
      for(size_t k=0;k<size;k++){
        double hourglass = offset + A*cos(k*dx*two_pi);
        double y = i*dx-0.5;
        double z = j*dx-0.5;
        double r = sqrt(y*y+z*z);
        if(r<=hourglass)
          raw_image[pos++] = 0;
        else
          raw_image[pos++] = 1;
      }
    }
  }
  return 0;
}

void CTImage::mesh(){
  if(verbose)
    std::cout<<"void mesh()\n";

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=25.0, 
      facet_size=1.0,
      cell_size=2.0,
      facet_distance=0.1);

  // Mesh_criteria criteria(facet_angle=25, facet_size=subsample*resolution,
  //                        cell_radius_edge_ratio=3, cell_size=subsample*resolution);

  //Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
  //                       cell_radius_edge_ratio=2, cell_size=5);

  // Mesh generation and optimization in one call
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(*domain, criteria,
      lloyd(time_limit=60),
      // odt(time_limit=60),
      exude(time_limit=60, sliver_bound=10));

  // C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(*domain, criteria);

  // Mesh generation and optimization in several call
  // C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
  //                                     no_perturb(), no_exude());

  // Output
  //

  // Export points
  Tr t = c3t3.triangulation();  // get triangulation (needed below)

  size_t index = 0;
  std::map<Point, size_t> coordToId;
  for(Tr::Point_iterator it=t.points_begin();it!=t.points_end();it++, index++){
    xyz.push_back(resolution*CGAL::to_double(it->x()));
    xyz.push_back(resolution*CGAL::to_double(it->y()));
    xyz.push_back(resolution*CGAL::to_double(it->z()));
    coordToId[*it] = index;
  }

  // Export elements
  for (C3t3::Cells_in_complex_iterator it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it){
    tets.push_back(coordToId[it->vertex(0)->point()]);
    tets.push_back(coordToId[it->vertex(1)->point()]);
    tets.push_back(coordToId[it->vertex(2)->point()]);
    tets.push_back(coordToId[it->vertex(3)->point()]);
  }

  // Find boundary facets
  size_t NNodes = get_NNodes();
  std::vector< std::set<int> > NEList(NNodes);

  size_t NElements = get_NElements();
  std::vector<int> EEList(NElements*4, -1);

  for(int i=0;i<NElements;i++){
    for(int j=0;j<4;j++){
      NEList[tets[i*4+j]].insert(i);
    }
  }
  for(int i=0;i<NElements;i++){
    for(int j=0;j<4;j++){
      int n1 = tets[i*4+(j+1)%4];
      int n2 = tets[i*4+(j+2)%4];
      int n3 = tets[i*4+(j+3)%4];

      for(std::set<int>::const_iterator it=NEList[n1].begin();it!=NEList[n1].end();++it){
	if(*it==i)
	  continue;
	if(NEList[n2].find(*it)!=NEList[n2].end() && NEList[n3].find(*it)!=NEList[n3].end()){
	  EEList[i*4+j] = *it;
	  break;
	}
      }
    }
  }
  for(int i=0;i<NElements;i++){
    if(EEList[i*4]==-1){
      facets.push_back(tets[i*4+1]);
      facets.push_back(tets[i*4+2]);
      facets.push_back(tets[i*4+3]);
    }
    if(EEList[i*4+1]==-1){
      facets.push_back(tets[i*4]);
      facets.push_back(tets[i*4+3]);
      facets.push_back(tets[i*4+2]);
    }
    if(EEList[i*4+2]==-1){
      facets.push_back(tets[i*4]);
      facets.push_back(tets[i*4+1]);
      facets.push_back(tets[i*4+3]);
    }
    if(EEList[i*4+3]==-1){
      facets.push_back(tets[i*4]);
      facets.push_back(tets[i*4+2]);
      facets.push_back(tets[i*4+1]);
    }
  }

  size_t NFacets = get_NFacets();

  // Label the boundary
  facet_ids = std::vector<int>(NFacets, 7);
  double dx=0.05*resolution;
  double dy=0.05*resolution;
  double dz=0.05*resolution;
  for(int i=0;i<NFacets;i++){
    double meanx = (xyz[facets[i*3]*3  ]+xyz[facets[i*3+1]*3  ]+xyz[facets[i*3+2]*3  ])/3;
    double meany = (xyz[facets[i*3]*3+1]+xyz[facets[i*3+1]*3+1]+xyz[facets[i*3+2]*3+1])/3;
    double meanz = (xyz[facets[i*3]*3+2]+xyz[facets[i*3+1]*3+2]+xyz[facets[i*3+2]*3+2])/3;

    if(meanx<dx)
      facet_ids[i] = 1;
    else if(((dims[0]-1)*resolution-meanx)<dx)
      facet_ids[i] = 2;
    else if(meany<dy)
      facet_ids[i] = 3;
    else if(((dims[1]-1)*resolution-meany)<dy)
      facet_ids[i] = 4;
    else if(meanz<dz)
      facet_ids[i] = 5;
    else if(((dims[2]-1)*resolution-meanz)<dz)
      facet_ids[i] = 6;
  }
}

void CTImage::trim_channels(int in_boundary, int out_boundary){
  if(verbose)
    std::cout<<"void trim_channels(int in_boundary, int out_boundary)"<<std::endl;

  // Create node-element adjancy list - delete invested elements as we go.
  size_t NNodes = get_NNodes();
  size_t NElements = get_NElements();

  std::vector< std::set<int> > NEList(NNodes);
  int count_positive=0, count_negative=0;
  for(int i=0;i<NElements;i++){
    if(tets[i*4]==-1)
      continue;

    double v = volume(&xyz[3*tets[i*4]], &xyz[3*tets[i*4+1]], &xyz[3*tets[i*4+2]], &xyz[3*tets[i*4+3]]);
    if(v<0){
      tets[i*4] = -1;
      count_negative++;
      continue;
    }else{
      count_positive++; 
    }

    for(int j=0;j<4;j++){
      NEList[tets[i*4+j]].insert(i);
    }
  }
  if(verbose)
    std::cout<<"Count of positive and negative volumes = "<<count_positive<<", "<<count_negative<<std::endl;

  // Create element-element adjancy list
  std::vector<int> EEList(NElements*4, -1);
  for(int i=0;i<NElements;i++){
    if(tets[i*4]==-1)
      continue;

    for(int j=0;j<4;j++){  
      std::set<int> edge_neighbours;
      set_intersection(NEList[tets[i*4+(j+1)%4]].begin(), NEList[tets[i*4+(j+1)%4]].end(),
	  NEList[tets[i*4+(j+2)%4]].begin(), NEList[tets[i*4+(j+2)%4]].end(),
	  inserter(edge_neighbours, edge_neighbours.begin()));

      std::set<int> neighbours;
      set_intersection(NEList[tets[i*4+(j+3)%4]].begin(), NEList[tets[i*4+(j+3)%4]].end(),
	  edge_neighbours.begin(), edge_neighbours.end(),
	  inserter(neighbours, neighbours.begin()));

      if(neighbours.size()==2){
	if(*neighbours.begin()==i)
	  EEList[i*4+j] = *neighbours.rbegin();
	else
	  EEList[i*4+j] = *neighbours.begin();
      }
#ifndef NDEBUG
      else{
	assert(neighbours.size()==1);
	assert(*neighbours.begin()==i);
      }
#endif
    }
  }

  // Create full facet ID list. Also, create the initial fronts for
  // the active region detection.
  std::map< std::set<int>, int> facet_id_lut;
  int NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    std::set<int> facet;
    for(int j=0;j<3;j++){
      facet.insert(facets[i*3+j]);
    }
    assert(facet.size()==3);
    assert(facet_id_lut.find(facet)==facet_id_lut.end());
    facet_id_lut[facet] = facet_ids[i];
  }

  std::set<int> front0, front1;
  std::map< std::set<int>, int> facet_element_lut;
  for(int i=0;i<NElements;i++){
    if(tets[i*4]==-1)
      continue;

    for(int j=0;j<4;j++){
      if(EEList[i*4+j]==-1){
	std::set<int> facet;
	for(int k=1;k<4;k++)
	  facet.insert(tets[i*4+(j+k)%4]);
	assert(facet.size()==3);
	assert(facet_element_lut.find(facet)==facet_element_lut.end());
	facet_element_lut[facet] = i;

	std::map< std::set<int>, int>::iterator facet_id_pair = facet_id_lut.find(facet);
	if(facet_id_pair==facet_id_lut.end()){
	  facet_ids.push_back(7);

	  if(j==0){
	    facets.push_back(tets[i*4+1]); facets.push_back(tets[i*4+3]); facets.push_back(tets[i*4+2]); 
	  }else if(j==1){
	    facets.push_back(tets[i*4]); facets.push_back(tets[i*4+3]); facets.push_back(tets[i*4+2]);
	  }else if(j==2){
	    facets.push_back(tets[i*4]); facets.push_back(tets[i*4+3]); facets.push_back(tets[i*4+1]);
	  }else if(j==3){
	    facets.push_back(tets[i*4]); facets.push_back(tets[i*4+2]); facets.push_back(tets[i*4+1]);
	  }
	}else{
	  if(facet_id_pair->second==in_boundary)
	    front0.insert(i);
	  else if(facet_id_pair->second==out_boundary)
	    front1.insert(i);
	}
      }
    }
  }

  // Advance front0
  std::vector<int> label(NElements, 0);
  while(!front0.empty()){
    // Get the next unprocessed element in the set.
    int seed = *front0.begin();
    front0.erase(front0.begin());
    if(label[seed]==1)
      continue;
    label[seed] = 1;

    for(int i=0;i<4;i++){
      int eid = EEList[seed*4+i];
      if(eid!=-1 && label[eid]!=1){
	front0.insert(eid);
      }
    }
  }

  // Advance back sweep using front1.
  while(!front1.empty()){
    // Get the next unprocessed element in the set.
    int seed = *front1.begin();
    front1.erase(front1.begin());
    if(label[seed]!=1) // ie was either never of interest or has been processed in the backsweep.
      continue;
    label[seed] = 2;

    for(int i=0;i<4;i++){
      int eid = EEList[seed*4+i];
      if(eid!=-1 && label[eid]==1){
	front1.insert(eid);
      }
    }
  }

  // Find active vertex set and create renumbering.
  std::map<int, int> renumbering;
  for(int i=0;i<NElements;i++){
    if(tets[i*4]==-1)
      continue;

    if(label[i]==2){
      for(int j=0;j<4;j++)
	renumbering.insert(std::pair<int, int>(tets[i*4+j], -1));
    }
  }

  // Create new compressed mesh.
  std::vector<double> xyz_new;
  std::vector<int> tets_new;
  std::vector<int> facets_new;
  std::vector<int> facet_ids_new;
  int cnt=0;
  for(std::map<int, int>::iterator it=renumbering.begin();it!=renumbering.end();++it){
    it->second = cnt++;

    xyz_new.push_back(xyz[(it->first)*3]);
    xyz_new.push_back(xyz[(it->first)*3+1]);
    xyz_new.push_back(xyz[(it->first)*3+2]);
  }
  for(int i=0;i<NElements;i++){
    if(tets[i*4]==-1)
      continue;

    if(label[i]==2){
      for(int j=0;j<4;j++){
	tets_new.push_back(renumbering[tets[i*4+j]]);
      }
    }
  }
  NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    std::set<int> facet;
    for(int j=0;j<3;j++)
      facet.insert(facets[i*3+j]);

    // Check if this is an orphaned facet from previous purge.
    if(facet_element_lut.find(facet)==facet_element_lut.end())
      continue;

    // Check if this is a newly orphaned facet.
    if(label[facet_element_lut[facet]]!=2){
      continue;
    }

    for(int j=0;j<3;j++){
      std::map<int, int>::iterator it=renumbering.find(facets[i*3+j]);
      assert(it!=renumbering.end());
      facets_new.push_back(it->second);
    }

    facet_ids_new.push_back(facet_ids[i]);
  }
  xyz.swap(xyz_new);
  tets.swap(tets_new);
  facets.swap(facets_new);
  facet_ids.swap(facet_ids_new);
}

// Write INR file.
void CTImage::write_inr(const char *filename){
  if(verbose)
    std::cout<<"void write_inr()"<<std::endl;

  if(filename==NULL)
    _writeImage(image->image(), std::string(stem.string()+".inr").c_str()); 
  else
    _writeImage(image->image(), filename);
}

// Write NRRD file.
void CTImage::write_nrrd(const char *filename){
  if(verbose)
    std::cout<<"void write_nrrd()"<<std::endl;

  std::ofstream file;
  if(filename==NULL)
    file.open(std::string(stem.string()+".nrrd").c_str());
  else
    file.open(filename);

  file<<"NRRD0004"<<std::endl
    <<"# Complete NRRD file format specification at:"<<std::endl
    <<"# http://teem.sourceforge.net/nrrd/format.html"<<std::endl
    <<"Content: Micro-CT scan of rock sample"<<std::endl
    <<"type: unsigned char"<<std::endl
    <<"dimension: 3"<<std::endl
    <<"space: 3D-right-handed"<<std::endl
    <<"sizes: "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<std::endl
    <<"spacings: "<<resolution<<" "<<resolution<<" "<<resolution<<std::endl
    <<"sample units: meters"<<std::endl
    <<"centerings: cell cell cell"<<std::endl
    <<"kinds: space space space"<<std::endl
    <<"endian: little"<<std::endl
    <<"encoding: raw"<<std::endl
    <<"space origin: (0,0,0)"<<std::endl
    <<"measurement frame: (1,0,0) (0,1,0) (0,0,1)"<<std::endl
    <<"data file: "<<stem.string()+".raw"<<std::endl;
  file.close();

  std::ofstream image_file;
  image_file.open(std::string(stem.string()+".raw").c_str(), std::ios::binary);
  image_file.write((const char *)raw_image, image_size);
  image_file.close();
}

// Write vox file.
void CTImage::write_vox(const char *filename){
  if(verbose)
    std::cout<<"void write_vox()"<<std::endl;

  std::ofstream file;
  if(filename==NULL)
    file.open(std::string(stem.string()+".vox").c_str());
  else
    file.open(filename);
  file<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<std::endl;
  file<<resolution<<" "<<resolution<<" "<<resolution<<std::endl;

  for(size_t i=0;i<image_size;i++)
    file<<(int)raw_image[i]<<" ";
  file<<std::endl;
  file.close();
}

// Write VTK unstructured grid file (*.vtu)
void CTImage::write_vtu(){
  if(verbose)
    std::cout<<"void write_vtu()"<<std::endl;

  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tets = vtkUnstructuredGrid::New();

  // Write out points
  size_t NNodes = get_NNodes();
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(NNodes);
  for(int i=0;i<NNodes;i++){
    pts->SetPoint(i, &(xyz[i*3]));
  }
  ug_tets->SetPoints(pts);

  size_t NElements = get_NElements();
  for(int i=0;i<NElements;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<4;j++)
      idlist->InsertNextId(tets[i*4+j]);
    ug_tets->InsertNextCell(10, idlist);
    idlist->Delete();
  }

  vtkXMLUnstructuredGridWriter *tet_writer = vtkXMLUnstructuredGridWriter::New();
  tet_writer->SetFileName(std::string(stem.string()+".vtu").c_str());
  tet_writer->SetInput(ug_tets);
  tet_writer->Write();

  ug_tets->Delete();
  tet_writer->Delete();

  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tris = vtkUnstructuredGrid::New();
  ug_tris->SetPoints(pts);

  // Write facet
  size_t NFacets = get_NFacets();
  for(int i=0;i<NFacets;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<3;j++)
      idlist->InsertNextId(facets[i*3+j]);
    ug_tris->InsertNextCell(5, idlist);
    idlist->Delete();
  }

  // Write facet_ids
  vtkIntArray *vtk_facet_ids = vtkIntArray::New();
  vtk_facet_ids->SetName("Boundary label");
  vtk_facet_ids->SetNumberOfComponents(1);
  vtk_facet_ids->SetNumberOfTuples(NFacets);
  for(int i=0;i<NFacets;i++)
    vtk_facet_ids->SetValue(i, facet_ids[i]);
  ug_tris->GetCellData()->AddArray(vtk_facet_ids);
  vtk_facet_ids->Delete();

  vtkXMLUnstructuredGridWriter *tri_writer = vtkXMLUnstructuredGridWriter::New();
  tri_writer->SetFileName(std::string(stem.string()+"_facets.vtu").c_str());
  tri_writer->SetInput(ug_tris);
  tri_writer->Write();

  ug_tris->Delete();
  tri_writer->Delete();  
}

int CTImage::write_gmsh(){
  if(verbose)
    std::cout<<"int write_gmsh()"<<std::endl;

  int NNodes = get_NNodes();
  int NElements = get_NElements();
  int NFacets = get_NFacets();

  ofstream file;
  file.open(std::string(stem.string()+".msh").c_str());
  file<<"$MeshFormat"<<std::endl
    <<"2.2 0 8"<<std::endl
    <<"$EndMeshFormat"<<std::endl
    <<"$Nodes"<<std::endl
    <<NNodes<<std::endl;
  file<<std::setprecision(std::numeric_limits<double>::digits10+1);
  for(size_t i=0;i<NNodes;i++){
    file<<i+1<<" "<<xyz[i*3]<<" "<<xyz[i*3+1]<<" "<<xyz[i*3+2]<<std::endl;
  }
  file<<"$EndNodes"<<std::endl
    <<"$Elements"<<std::endl
    <<NElements+NFacets<<std::endl;
  for(size_t i=0;i<NElements;i++){
    file<<i+1<<" 4 1 1 "<<tets[i*4]+1<<" "<<tets[i*4+1]+1<<" "<<tets[i*4+2]+1<<" "<<tets[i*4+3]+1<<std::endl;
  }
  for(size_t i=0;i<NFacets;i++){
    file<<i+NElements+1<<" 2 1 "<<facet_ids[i]<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<std::endl;
  }
  file<<"$EndElements"<<std::endl;
  file.close();

  return 0;
}

double CTImage::volume(const double *x0, const double *x1, const double *x2, const double *x3) const{

  double x01 = (x0[0] - x1[0]);
  double x02 = (x0[0] - x2[0]);
  double x03 = (x0[0] - x3[0]);

  double y01 = (x0[1] - x1[1]);
  double y02 = (x0[1] - x2[1]);
  double y03 = (x0[1] - x3[1]);

  double z01 = (x0[2] - x1[2]);
  double z02 = (x0[2] - x2[2]);
  double z03 = (x0[2] - x3[2]);

  return (-x03*(z02*y01 - z01*y02) + x02*(z03*y01 - z01*y03) - x01*(z03*y02 - z02*y03))/6;
}

