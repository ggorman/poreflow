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
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>

#include <fstream>
#include <cassert>
#include <limits>

#include "writers.h"

int write_vtk_file(std::string filename,
                   std::vector<double> &xyz,
                   std::vector<int> &tets, 
                   std::vector<int> &facets,
                   std::vector<int> &facet_ids){
  
  // Write out points
  int NNodes = xyz.size()/3;
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(NNodes);
  for(int i=0;i<NNodes;i++)
    pts->SetPoint(i, &(xyz[i*3]));
  
  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tets = vtkUnstructuredGrid::New();
  ug_tets->SetPoints(pts);

  int NTetra = tets.size()/4;
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<4;j++)
      idlist->InsertNextId(tets[i*4+j]);
    ug_tets->InsertNextCell(10, idlist);
    idlist->Delete();
  }
  
  vtkXMLUnstructuredGridWriter *tet_writer = vtkXMLUnstructuredGridWriter::New();
  tet_writer->SetFileName(std::string(filename+".vtu").c_str());
  tet_writer->SetInput(ug_tets);
  tet_writer->Write();
  
  ug_tets->Delete();
  tet_writer->Delete();

  if(facets.empty())
    return 0;

  // Write out facets
  vtkUnstructuredGrid *ug_facets = vtkUnstructuredGrid::New();
  ug_facets->SetPoints(pts);
  int NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<3;j++){
      idlist->InsertNextId(facets[i*3+j]);
    }
    ug_facets->InsertNextCell(5, idlist);
    idlist->Delete();
  }

  vtkIntArray *vtk_facet_ids = vtkIntArray::New();
  vtk_facet_ids->SetNumberOfTuples(NFacets);
  vtk_facet_ids->SetNumberOfComponents(1);
  vtk_facet_ids->SetName("Facet IDs");
  for(int i=0;i<NFacets;i++){
    vtk_facet_ids->SetValue(i, facet_ids[i]);
  }
  ug_facets->GetCellData()->AddArray(vtk_facet_ids);
  vtk_facet_ids->Delete();
  
  vtkXMLUnstructuredGridWriter *tri_writer = vtkXMLUnstructuredGridWriter::New();
  tri_writer->SetFileName(std::string(filename+"_facets.vtu").c_str());
  tri_writer->SetInput(ug_facets);
  tri_writer->Write();

  pts->Delete();
  ug_facets->Delete();
  tri_writer->Delete();

  return 0;
}

int write_triangle_file(std::string basename,
                        std::vector<double> &xyz,
                        std::vector<int> &tets, 
                        std::vector<int> &facets,
                        std::vector<int> &facet_ids){
  std::string filename_node = basename+".node";
  std::string filename_face = basename+".face";
  std::string filename_ele = basename+".ele";
  
  int NNodes = xyz.size()/3;
  int NTetra = tets.size()/4;
  int NFacets = facet_ids.size();
  assert(NFacets==facets.size()/3);

  ofstream nodefile;
  nodefile.open(std::string(basename+".node").c_str());
  nodefile<<NNodes<<" "<<3<<" "<<0<<" "<<0<<std::endl;
  nodefile<<std::setprecision(std::numeric_limits<double>::digits10+1);
  
  for(int i=0;i<NNodes;i++){
    nodefile<<i+1<<" "<<xyz[i*3]<<" "<<xyz[i*3+1]<<" "<<xyz[i*3+2]<<std::endl;
  }

  ofstream elefile;
  elefile.open(std::string(basename+".ele").c_str());
  elefile<<NTetra<<" "<<4<<" "<<1<<std::endl;

  for(int i=0;i<NTetra;i++){
    elefile<<i+1<<" "<<tets[i*4]+1<<" "<<tets[i*4+1]+1<<" "<<tets[i*4+2]+1<<" "<<tets[i*4+3]+1<<" 1"<<std::endl;
  }

  ofstream facefile;
  facefile.open(std::string(basename+".face").c_str());
  facefile<<NFacets<<" "<<1<<std::endl;
  for(int i=0;i<NFacets;i++){
    facefile<<i+1<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<" "<<facet_ids[i]<<std::endl;
  }
  
  return 0;
}

int write_gmsh_file(std::string basename,
		    std::vector<double> &xyz,
		    std::vector<int> &tets, 
		    std::vector<int> &facets,
		    std::vector<int> &facet_ids){
  
  int NNodes = xyz.size()/3;
  int NTetra = tets.size()/4;
  int NFacets = facet_ids.size();
  assert(NFacets==facets.size()/3);

  ofstream file;
  file.open(std::string(basename+".msh").c_str());
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
      <<NTetra+NFacets<<std::endl;
  for(size_t i=0;i<NTetra;i++){
    file<<i+1<<" 4 1 1 "<<tets[i*4]+1<<" "<<tets[i*4+1]+1<<" "<<tets[i*4+2]+1<<" "<<tets[i*4+3]+1<<std::endl;
  }
  for(size_t i=0;i<NFacets;i++){
    file<<i+NTetra+1<<" 2 1 "<<facet_ids[i]<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<std::endl;
  }
  file<<"$EndElements"<<std::endl;
  file.close();

  return 0;
}

