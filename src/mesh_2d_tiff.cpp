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

#include <vector>
#include <map>
#include <string>

#include <cassert>

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

  return 0;
}
