#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

typedef Tr::Point Point;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

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

#include <boost/filesystem.hpp>
#include "boost/filesystem/fstream.hpp"
#include "boost/filesystem/operations.hpp"

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int write_triangle_files_3d(std::string basename, std::vector<double> &xyz, std::vector<int> &cells, std::vector<int> &region_id, std::vector<int> &facets, std::vector<int> &boundary_id){
  // Write node file.
  std::ofstream node_file;
  node_file.open(std::string(basename+".node").c_str());

  int NNodes = xyz.size()/3;
  node_file<<NNodes<<" 3 0 0\n";
  for(int i=0;i<NNodes;i++){
    node_file<<i<<" "<<xyz[i*3]<<" "<<xyz[i*3+1]<<" "<<xyz[i*3+2]<<std::endl;
  }
  node_file.close();

  // Write ele file
  std::ofstream ele_file;
  ele_file.open(std::string(basename+".ele").c_str());

  int NCells = cells.size()/4;
  ele_file<<NCells<<" 4 1\n";
  for(int i=0;i<NCells;i++){
    ele_file<<i+1<<" "<<cells[i*4]+1<<" "<<cells[i*4+1]+1<<" "<<cells[i*4+2]+1<<" "<<cells[i*4+3]+1<<" "<<region_id[i]<<std::endl;
  }
  ele_file.close();

  // Write face file
  std::ofstream face_file;
  face_file.open(std::string(basename+".face").c_str());

  int NFacets = facets.size()/3;
  face_file<<NFacets<<" 1\n";
  for(int i=0;i<NFacets;i++){
    face_file<<i+1<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<" "<<boundary_id[i]<<std::endl;
  }
  face_file.close();

  return(0);
}

int main(int argc, char **argv){
  if(argc==1){
    std::cout<<"Usage: "<<argv[0]<< "CT-image [subsample]\n";
    return 0;
  }

  int subsample=1;
  if(argc==3){
    subsample = atoi(argv[2]);
  }

  boost::filesystem::path image_dir(argv[1]);  
  boost::filesystem::path stem = image_dir.stem();

  boost::filesystem::path image_info = image_dir/std::string(stem.string()+"_info.txt");
  boost::filesystem::path image_filename = image_dir/std::string(stem.string()+".raw");

  if(!boost::filesystem::exists(image_info)){
    std::cerr<<"ERROR: "<<image_info.string()<<" does not exist.";
    exit(-1);
  }
  
  if(!boost::filesystem::exists(image_filename)){
    std::cerr<<"ERROR: "<<image_filename.string()<<" does not exist.";
    exit(-1);
  }
  
  // Get the metadata.
  std::ifstream info;
  info.open(image_info.string().c_str());
  int dims[] = {-1, -1, -1};
  double resolution=1.0;
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
  }
  info.close();

  // Read image.
  int image_size = dims[0]*dims[1]*dims[2];
  unsigned char *raw_image = new unsigned char[image_size];

  std::ifstream image_file;
  image_file.open(image_filename.string().c_str(), std::ios::binary);
  image_file>>raw_image;
  image_file.close();

  double spacing[3];
  for(int i=0;i<3;i++)
    spacing[i]=resolution;
  CGAL::Image_3 image(_createImage(dims[0], dims[1], dims[2], 1,
                                   spacing[0], spacing[1], spacing[2],
                                   1, WK_FIXED, SGN_UNSIGNED)); 
  ImageIO_free(image.data());
  image.set_data((void*)(&raw_image[0])); 


  // Domain
  Mesh_domain domain(image);

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=25, facet_size=subsample*resolution,
                         cell_radius_edge_ratio=3, cell_size=subsample*resolution);

  // Mesh generation and optimization in one call
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      lloyd(time_limit=30),
                                      no_perturb(),
                                      exude(time_limit=10, sliver_bound=10));

  // C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Mesh generation and optimization in several call
  //C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
  //                                    no_perturb(), no_exude());

  // Output
  //

  // Export points
  Tr t = c3t3.triangulation();  // get triangulation (needed below)

  std::vector<double> xyz;
  std::map<Point, size_t> coordToId;
  size_t index = 0;
  for(Tr::Point_iterator it=t.points_begin();it!=t.points_end();it++, index++){
    xyz.push_back(CGAL::to_double(it->x()));
    xyz.push_back(CGAL::to_double(it->y()));
    xyz.push_back(CGAL::to_double(it->z()));
    coordToId[*it] = index;
  }
  int NNodes = xyz.size()/3;
  double bounds[6];
  bounds[0] = xyz[0];
  bounds[1] = xyz[0];
  bounds[2] = xyz[1];
  bounds[3] = xyz[1];
  bounds[4] = xyz[2];
  bounds[5] = xyz[2];
  for(int i=0;i<NNodes;i++){
    bounds[0] = std::min(bounds[0], xyz[i*3]);
    bounds[1] = std::max(bounds[1], xyz[i*3]);
    bounds[2] = std::min(bounds[2], xyz[i*3+1]);
    bounds[3] = std::max(bounds[3], xyz[i*3+1]);
    bounds[4] = std::min(bounds[4], xyz[i*3+2]);
    bounds[5] = std::max(bounds[5], xyz[i*3+2]);
  }

  // Export elements
  std::vector<int> enlist, subdomain_index;
  // for(Tr::Finite_cells_iterator it=t.finite_cells_begin();it!=t.finite_cells_end();it++){
  for (C3t3::Cells_in_complex_iterator it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it){
    enlist.push_back(coordToId[it->vertex(0)->point()]);
    enlist.push_back(coordToId[it->vertex(1)->point()]);
    enlist.push_back(coordToId[it->vertex(2)->point()]);
    enlist.push_back(coordToId[it->vertex(3)->point()]);
    
    subdomain_index.push_back(it->subdomain_index());
  }
  size_t NElements = enlist.size()/4;

  // Write out points
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(NNodes);
  for(int i=0;i<NNodes;i++){
    pts->SetPoint(i, &(xyz[i*3]));
  }

  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tets = vtkUnstructuredGrid::New();
  ug_tets->SetPoints(pts);
  
  for(int i=0;i<NElements;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<4;j++)
      idlist->InsertNextId(enlist[i*4+j]);
    ug_tets->InsertNextCell(10, idlist);
    idlist->Delete();
  }

  // Write subdomain index
  vtkIntArray *vtk_subdomain_index = vtkIntArray::New();
  vtk_subdomain_index->SetName("Subdomain_index");
  vtk_subdomain_index->SetNumberOfComponents(1);
  vtk_subdomain_index->SetNumberOfTuples(NElements);
  for(int i=0;i<NElements;i++)
    vtk_subdomain_index->SetValue(i, subdomain_index[i]);
  ug_tets->GetCellData()->AddArray(vtk_subdomain_index);
  vtk_subdomain_index->Delete();

  vtkXMLUnstructuredGridWriter *tet_writer = vtkXMLUnstructuredGridWriter::New();
  tet_writer->SetFileName(std::string(stem.string()+".vtu").c_str());
  tet_writer->SetInput(ug_tets);
  tet_writer->Write();
  
  ug_tets->Delete();
  tet_writer->Delete();

  // Find boundary facets
  std::vector<int> facets;
  std::vector< std::set<int> > NEList(NNodes);
  std::vector<int> EEList(NElements*4, -1);
  for(int i=0;i<NElements;i++){
    for(int j=0;j<4;j++){
      NEList[enlist[i*4+j]].insert(i);
    }
  }
  for(int i=0;i<NElements;i++){
    for(int j=0;j<4;j++){
      int n1 = enlist[i*4+(j+1)%4];
      int n2 = enlist[i*4+(j+2)%4];
      int n3 = enlist[i*4+(j+3)%4];

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
      facets.push_back(enlist[i*4+1]);
      facets.push_back(enlist[i*4+2]);
      facets.push_back(enlist[i*4+3]);
    }
    if(EEList[i*4+1]==-1){
      facets.push_back(enlist[i*4]);
      facets.push_back(enlist[i*4+3]);
      facets.push_back(enlist[i*4+2]);
    }
    if(EEList[i*4+2]==-1){
      facets.push_back(enlist[i*4]);
      facets.push_back(enlist[i*4+1]);
      facets.push_back(enlist[i*4+3]);
    }
    if(EEList[i*4+3]==-1){
      facets.push_back(enlist[i*4]);
      facets.push_back(enlist[i*4+2]);
      facets.push_back(enlist[i*4+1]);
    }
  }
  int NFacets = facets.size()/3;

  // Characterise the boundary
  std::vector< std::set<int> > fNEList(NNodes);
  std::vector<int> fEEList(NFacets*3, -1);
  for(int i=0;i<NFacets;i++){
    for(int j=0;j<3;j++){
      fNEList[facets[i*3+j]].insert(i);
    }
  }
  for(int i=0;i<NFacets;i++){
    for(int j=0;j<3;j++){
      int n1 = facets[i*3+(j+1)%3];
      int n2 = facets[i*3+(j+2)%3];
      
      for(std::set<int>::const_iterator it=fNEList[n1].begin();it!=fNEList[n1].end();++it){
        if(*it==i)
          continue;
        if(fNEList[n2].find(*it)!=fNEList[n2].end()){
          fEEList[i*3+j] = *it;
          break;
        }
      }
    }
  }

  std::vector<int> boundary_id(NFacets, -1);
  double dx=0.5*spacing[0];
  double dy=0.5*spacing[1];
  double dz=0.5*spacing[2];
  for(int i=0;i<NFacets;i++){
    double meanx = (xyz[facets[i*3]*3  ]+xyz[facets[i*3+1]*3  ]+xyz[facets[i*3+2]*3  ])/3;
    double meany = (xyz[facets[i*3]*3+1]+xyz[facets[i*3+1]*3+1]+xyz[facets[i*3+2]*3+1])/3;
    double meanz = (xyz[facets[i*3]*3+2]+xyz[facets[i*3+1]*3+2]+xyz[facets[i*3+2]*3+2])/3;

    std::multimap<double, int> boundary_proximity;
    boundary_proximity.insert(std::pair<double, int>(meanx-bounds[0], 1));
    boundary_proximity.insert(std::pair<double, int>(bounds[1]-meanx, 2));
    boundary_proximity.insert(std::pair<double, int>(meany-bounds[2], 3));
    boundary_proximity.insert(std::pair<double, int>(bounds[3]-meany, 4));
    boundary_proximity.insert(std::pair<double, int>(meanz-bounds[4], 5));
    boundary_proximity.insert(std::pair<double, int>(bounds[5]-meanz, 6));
    
    boundary_id[i] = boundary_proximity.begin()->second;
  }

  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tris = vtkUnstructuredGrid::New();
  ug_tris->SetPoints(pts);
  
  for(int i=0;i<NFacets;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<3;j++)
      idlist->InsertNextId(facets[i*3+j]);
    ug_tris->InsertNextCell(5, idlist);
    idlist->Delete();
  }

  // Write boundary_id
  vtkIntArray *vtk_boundary_id = vtkIntArray::New();
  vtk_boundary_id->SetName("Boundary id");
  vtk_boundary_id->SetNumberOfComponents(1);
  vtk_boundary_id->SetNumberOfTuples(NFacets);
  for(int i=0;i<NFacets;i++)
    vtk_boundary_id->SetValue(i, boundary_id[i]);
  ug_tris->GetCellData()->AddArray(vtk_boundary_id);
  vtk_boundary_id->Delete();

  vtkXMLUnstructuredGridWriter *tri_writer = vtkXMLUnstructuredGridWriter::New();
  tri_writer->SetFileName(std::string(stem.string()+"_facets.vtu").c_str());
  tri_writer->SetInput(ug_tris);
  tri_writer->Write();
  
  ug_tris->Delete();
  tri_writer->Delete();
  
  // write_triangle_files_3d(basename, xyz, enlist, subdomain_index, facets, boundary_id);

  return 0;
}
