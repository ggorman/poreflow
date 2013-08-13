/*
  This file is part of ImageMesher.

  ImageMesher is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.

  ImageMesher is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ImageMesher.  If not, see <http://www.gnu.org/licenses/>.

  Copyright holder: Gerard Gorman <gerard.j.gorman@gmail.com>

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
#include <vtkCell.h>
#include <vtkCellData.h>

#include <vector>
#include <map>
#include <set>
#include <string>

#include <iostream>
#include <fstream>

#include <cassert>

int write_triangle_files_2d(std::string basename, std::vector<double> &xy, std::vector<int> &cells, std::vector<int> &region_id, std::vector<int> &facets, std::vector<int> &boundary_id){
  // Write node file.
  std::ofstream node_file;
  node_file.open(std::string(basename+".node").c_str());

  int NNodes = xy.size()/2;
  node_file<<NNodes<<" 2 0 0\n";
  for(int i=0;i<NNodes;i++){
    node_file<<i<<" "<<xy[i*2]<<" "<<xy[i*2+1]<<std::endl;
  }
  node_file.close();

  // Write ele file
  std::ofstream ele_file;
  ele_file.open(std::string(basename+".ele").c_str());

  int NCells = cells.size()/3;
  ele_file<<NCells<<" 3 1\n";
  for(int i=0;i<NCells;i++){
    ele_file<<i<<" "<<cells[i*3]<<" "<<cells[i*3+1]<<" "<<cells[i*3+2]<<" "<<region_id[i]<<std::endl;
  }
  ele_file.close();

  // Write edge file
  std::ofstream edge_file;
  edge_file.open(std::string(basename+".edge").c_str());

  int NFacets = facets.size()/2;
  edge_file<<NFacets<<" 1\n";
  for(int i=0;i<NFacets;i++){
    edge_file<<i<<" "<<facets[i*2]<<" "<<facets[i*2+1]<<" "<<boundary_id[i]<<std::endl;
  }
  edge_file.close();

  return(0);

}

int main(int argc, char **argv){
  if(argc==1){
    std::cout<<"Usage: "<<argv[0]<< "image.tiff\n";
    return 0;
  }

  // Load image.
  vtkTIFFReader *tiffreader = vtkTIFFReader::New();
  tiffreader->SetFileName(argv[1]);
  tiffreader->Update();
  vtkImageData *vtk_image_data = tiffreader->GetOutput();

  int dims[3];
  vtk_image_data->GetDimensions(dims);

  double bounds[6];
  vtk_image_data->GetBounds(bounds);

  double spacing[3];
  vtk_image_data->GetSpacing(spacing);

  bounds[1]+=spacing[0];
  bounds[3]+=spacing[1];

  assert(dims[2]==1);

  // Convert data to polydata.
  vtkImageToPolyDataFilter *image2polydata = vtkImageToPolyDataFilter::New();
  image2polydata->SetOutputStyleToPixelize();
  image2polydata->SetInput(vtk_image_data);
  image2polydata->Update();

  // Convert polydata to triangles
  vtkTriangleFilter *triangle = vtkTriangleFilter::New();
  triangle->SetInput(image2polydata->GetOutput());
  triangle->Update();

  vtkPolyData *pd = vtkPolyData::New();
  pd->DeepCopy(triangle->GetOutput());
  pd->Update();

  int NCells = pd->GetNumberOfCells();
  vtkIntArray *subdomain_id = vtkIntArray::New();
  subdomain_id->SetName("subdomain id");
  subdomain_id->SetNumberOfTuples(NCells);
  subdomain_id->SetNumberOfComponents(1);
  std::map< std::vector<int>, int > lut;
  for(int i=0;i<NCells;i++){
    std::vector<int> tuple3;
    tuple3.push_back(pd->GetCellData()->GetArray(0)->GetTuple(i)[0]);
    tuple3.push_back(pd->GetCellData()->GetArray(0)->GetTuple(i)[1]);
    tuple3.push_back(pd->GetCellData()->GetArray(0)->GetTuple(i)[2]);
    lut[tuple3] = 0;
  }
  {
    int colour = 1;
    for(std::map< std::vector<int>, int >::iterator it=lut.begin();it!=lut.end();++it){
      it->second = colour++;
    }
  }
  for(int i=0;i<NCells;i++){
    std::vector<int> tuple3;
    tuple3.push_back(pd->GetCellData()->GetArray(0)->GetTuple(i)[0]);
    tuple3.push_back(pd->GetCellData()->GetArray(0)->GetTuple(i)[1]);
    tuple3.push_back(pd->GetCellData()->GetArray(0)->GetTuple(i)[2]);

    subdomain_id->SetValue(i, lut[tuple3]);
  }

  pd->GetCellData()->RemoveArray(0);
  pd->GetCellData()->AddArray(subdomain_id);
  subdomain_id->Delete();

  std::string basename(argv[1], std::string(argv[1]).find(".tif"));
  std::string filename(std::string(basename+".vtp"));

  vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(pd);
  writer->Write();

  int NNodes = pd->GetNumberOfPoints();
  std::vector<double> xy(NNodes*2);
  std::vector<int> cells(NCells*3), region_id(NCells), facets, boundary_id;

  for(int i=0;i<NNodes;i++){
    pd->GetPoints()->GetPoint(i, &(xy[i*2]));
  }
  std::vector< std::set<int> > NEList(NNodes);
  std::vector<int> EEList(NCells*3, -1);
  for(int i=0;i<NCells;i++){
    for(int j=0;j<3;j++){
      cells[i*3+j] = pd->GetCell(i)->GetPointId(j);
      NEList[cells[i*3+j]].insert(i);
    }
    region_id[i] = pd->GetCellData()->GetArray(0)->GetTuple1(i);
  }
  for(int i=0;i<NCells;i++){
    for(int j=0;j<3;j++){
      int n1 = cells[i*3+(j+1)%3];
      int n2 = cells[i*3+(j+2)%3];
      for(std::set<int>::const_iterator it=NEList[n1].begin();it!=NEList[n1].end();++it){
        if(*it==i)
          continue;
        if(NEList[n2].find(*it)!=NEList[n2].end()){
          EEList[i*3+j] = *it;
          break;
        }
      }
    }
  }
  for(int i=0;i<NCells;i++){
    if(EEList[i*3]==-1){
      facets.push_back(cells[i*3+1]);
      facets.push_back(cells[i*3+2]);
    }
    if(EEList[i*3+1]==-1){
      facets.push_back(cells[i*3+2]);
      facets.push_back(cells[i*3]);
    }
    if(EEList[i*3+2]==-1){
      facets.push_back(cells[i*3]);
      facets.push_back(cells[i*3+1]);
    }
  }
  int NFacets = facets.size()/2;
  double tol=spacing[0]*0.1;
  for(int i=0;i<NFacets;i++){
    if((fabs(xy[facets[i*2]*2]-bounds[0])<tol)&&(fabs(xy[facets[i*2+1]*2]-bounds[0])<tol))
      boundary_id.push_back(1);
    else if((fabs(xy[facets[i*2]*2]-bounds[1])<tol)&&(fabs(xy[facets[i*2+1]*2]-bounds[1])<tol))
      boundary_id.push_back(2);
    else if((fabs(xy[facets[i*2]*2+1]-bounds[2])<tol)&&(fabs(xy[facets[i*2+1]*2+1]-bounds[2])<tol))
      boundary_id.push_back(3);
    else if((fabs(xy[facets[i*2]*2+1]-bounds[3])<tol)&&(fabs(xy[facets[i*2+1]*2+1]-bounds[3])<tol))
      boundary_id.push_back(4);
    else{
      std::cerr<<"something foobar\n"
               <<bounds[0]<<", "<<bounds[1]<<", "<<bounds[2]<<", "<<bounds[3]<<"\n"
               <<facets[i*2]<<", "<<facets[i*2+1]<<"\n"
               <<xy[facets[i*2]*2]<<", "<<xy[facets[i*2]*2+1]<<"\n";
    }
  }

  write_triangle_files_2d(basename, xy, cells, region_id, facets, boundary_id);

  return 0;
}
