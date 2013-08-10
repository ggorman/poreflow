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

  vtkPolyData *ug = triangle->GetOutput();
  vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName("ug.vtp");
  writer->SetInput(ug);
  writer->Write();

  return 0;
}
