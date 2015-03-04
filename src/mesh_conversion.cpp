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

#include <algorithm>
#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <string>
#include <vector>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCellType.h>
#include <vtkPoints.h>
#include <vtkCell.h>

#include "writers.h"
#include "mesh_conversion.h"

int create_domain(int axis,
		  std::vector<double> &xyz,
                  std::vector<int> &tets, 
                  std::vector<int> &facets,
                  std::vector<int> &facet_ids){
  
  // Create node-element adjancy list.
  size_t NNodes = xyz.size()/3;

  std::vector< std::unordered_set<int> > NEList(NNodes);
  int NTetra = tets.size()/4;

  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    double v = volume(xyz.data()+3*tets[i*4],
     		      xyz.data()+3*tets[i*4+1], 
		      xyz.data()+3*tets[i*4+2],
		      xyz.data()+3*tets[i*4+3]);
    if(v<0)
      std::swap(tets[i*4+2], tets[i*4+3]);
    
    
    for(int j=0;j<4;j++)
      NEList[tets[i*4+j]].insert(i);
  }

  std::vector<int> EEList(NTetra*4);
#pragma omp parallel
  {
    // Initialise and ensure 1st touch placement.
#pragma omp for
    for(int i=0;i<4*NTetra;i++)
      EEList[i] = -1;

#pragma omp for
    for(int i=0;i<NTetra;i++){
      if(tets[i*4]==-1)
        continue;

      int n0 = tets[i*4];
      int n1 = tets[i*4+1];
      int n2 = tets[i*4+2];
      int n3 = tets[i*4+3];

      for(auto& e : NEList[n0]){
        if(e!=i){
          if(NEList[n1].find(e)!=NEList[n1].end()){
            if(EEList[i*4+3]==-1 && NEList[n2].find(e)!=NEList[n2].end()){
              EEList[i*4+3] = e;
            }else if(EEList[i*4+2]==-1 && NEList[n3].find(e)!=NEList[n3].end()){
              EEList[i*4+2] = e;
            }
          }
        }
        if(EEList[i*4+3]!=-1 && EEList[i*4+2]!=-1)
          break;
      }

      for(auto& e : NEList[n2]){
        if(e!=i){
          if(NEList[n3].find(e)!=NEList[n3].end()){
            if(EEList[i*4+1]==-1 && NEList[n0].find(e)!=NEList[n0].end()){
              EEList[i*4+1] = e;
            }else if(EEList[i*4]==-1 && NEList[n1].find(e)!=NEList[n1].end()){
              EEList[i*4] = e;
            }
          }
        }
        if(EEList[i*4+1]!=-1 && EEList[i*4]!=-1)
          break;
      }
    }
  }
  NEList.clear();

  // Calculate the bounding box.
  double bbox[] = {xyz[0], xyz[0],
		   xyz[1], xyz[1],
		   xyz[2], xyz[2]};
  for(size_t i=1;i<NNodes;i++){
    for(size_t j=0;j<3;j++){
      bbox[j*2  ] = std::min(bbox[j*2  ], xyz[i*3+j]);
      bbox[j*2+1] = std::max(bbox[j*2+1], xyz[i*3+j]);
    }
  }

  // Calculate the a element size - use the l-infinity norm.
  size_t livecnt=0;
  double eta=0.0;
  for(size_t i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;

    livecnt++;

    int vid = tets[i*4];
    double lbbox[] = {xyz[vid*3],   xyz[vid*3],
		      xyz[vid*3+1], xyz[vid*3+1],
		      xyz[vid*3+2], xyz[vid*3+2]};
    for(size_t j=1;j<4;j++){
      vid = tets[i*4+j];
      for(size_t k=0;k<3;k++){
	lbbox[k*2  ] = std::min(lbbox[k*2  ], xyz[vid*3+k]);
	lbbox[k*2+1] = std::max(lbbox[k*2+1], xyz[vid*3+k]);
      }
    }
    eta += ((lbbox[1]-lbbox[0])+
	    (lbbox[3]-lbbox[2])+
	    (lbbox[5]-lbbox[4]));
  }
  eta/=(livecnt*3);   // i.e. the mean element size
  
  // Define what we mean by a "small" distance.
  eta*=0.1;
  
  // Calculate the facet list, facet id's and the initial forward and
  // backward fronts.
  std::set<int> front0, front1;
  for(size_t i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    for(size_t j=0;j<4;j++){
      bool is_facet = false;
      int facet[3];
      if(EEList[i*4+j]==-1){
	is_facet=true;
	switch(j){
	case 0:
	  facet[0] = tets[i*4+1];
	  facet[1] = tets[i*4+3];
	  facet[2] = tets[i*4+2];
	  break;
	case 1:
	  facet[0] = tets[i*4+0];
	  facet[1] = tets[i*4+2];
	  facet[2] = tets[i*4+3];
	  break;
	case 2:
	  facet[0] = tets[i*4+0];
	  facet[1] = tets[i*4+3];
	  facet[2] = tets[i*4+1];
	  break;
	case 3:
	  facet[0] = tets[i*4+0];
	  facet[1] = tets[i*4+1];
	  facet[2] = tets[i*4+2];
	  break;
	}
      }
      
      if(is_facet){
	// Decide boundary id.
	double mean_x = (xyz[facet[0]*3+axis]+
			 xyz[facet[1]*3+axis]+
			 xyz[facet[2]*3+axis])/3.0;
	
	if(fabs(mean_x-bbox[axis*2])<eta){
	  front0.insert(front0.end(), i);
	}else if(fabs(mean_x-bbox[axis*2+1])<eta){
	  front1.insert(front1.end(), i);
	}
      }
    }
  }
  
  // Advance front0
  std::vector<int> label(NTetra, 0);
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
    
    if(label[seed]!=1) // i.e. was either never of interest or has been processed in the backsweep.
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
  for(int i=0;i<NTetra;i++){
    if(label[i]==2){
      for(int j=0;j<4;j++)
        renumbering.insert(std::pair<int, int>(tets[i*4+j], -1));
    }
  }

  // Create new compressed mesh.
  std::vector<double> xyz_new;
  std::vector<int> tets_new;
  int cnt=0;
  for(auto& it : renumbering){
    it.second = cnt++;
    
    xyz_new.push_back(xyz[(it.first)*3]);
    xyz_new.push_back(xyz[(it.first)*3+1]);
    xyz_new.push_back(xyz[(it.first)*3+2]);
  }
  for(int i=0;i<NTetra;i++){
    if(label[i]==2){
      for(int j=0;j<4;j++){
        tets_new.push_back(renumbering[tets[i*4+j]]);
      }
    }
  }

  // Re-create facets.
  facets.clear();
  facet_ids.clear();
  for(size_t i=0;i<NTetra;i++){
    if(label[i]!=2)
      continue;
    
    for(size_t j=0;j<4;j++){
      bool is_facet = false;
      int facet[3];
      if(EEList[i*4+j]==-1){
	is_facet=true;
	switch(j){
	case 0:
	  facet[0] = tets[i*4+1];
	  facet[1] = tets[i*4+3];
	  facet[2] = tets[i*4+2];
	  break;
	case 1:
	  facet[0] = tets[i*4+0];
	  facet[1] = tets[i*4+2];
	  facet[2] = tets[i*4+3];
	  break;
	case 2:
	  facet[0] = tets[i*4+0];
	  facet[1] = tets[i*4+3];
	  facet[2] = tets[i*4+1];
	  break;
	case 3:
	  facet[0] = tets[i*4+0];
	  facet[1] = tets[i*4+1];
	  facet[2] = tets[i*4+2];
	  break;
	}
      }
      
      if(is_facet){
	// Insert new facet.
	for(int k=0;k<3;k++)
	  facets.push_back(renumbering[facet[k]]);
	
	// Decide boundary id.
	double mean_xyz[3];
	for(int k=0;k<3;k++)
	  mean_xyz[k] = (xyz[facet[0]*3+k]+xyz[facet[1]*3+k]+xyz[facet[2]*3+k])/3.0;
	
	if(fabs(mean_xyz[0]-bbox[0])<eta){
	  facet_ids.push_back(1);
	}else if(fabs(mean_xyz[0]-bbox[1])<eta){
	  facet_ids.push_back(2);
	}else if(fabs(mean_xyz[1]-bbox[2])<eta){
	  facet_ids.push_back(3);
	}else if(fabs(mean_xyz[1]-bbox[3])<eta){
	  facet_ids.push_back(4);
	}else if(fabs(mean_xyz[2]-bbox[4])<eta){
	  facet_ids.push_back(5);
	}else if(fabs(mean_xyz[2]-bbox[5])<eta){
	  facet_ids.push_back(6);
	}else{
	  facet_ids.push_back(7);
	}
      }
    }
  }

  xyz.swap(xyz_new);
  tets.swap(tets_new);

  return 0;
}


double read_resolution_from_nhdr(std::string filename){
  
  double resolution = 1.0;

  // Open file.
  std::ifstream nhdr(filename);
  if(!nhdr.is_open()){
    std::cerr<<"ERROR: Cannot open NHDR file: "<<filename<<std::endl;
    return resolution;
  }

  // Parse NHDR file
  std::string line;
  std::getline(nhdr, line);
  if(line.compare(0, 8, "NRRD0004")!=0){
    std::cerr<<"ERROR: Unexpected first line in NHDR file. Expected NRRD0004, got "<<line<<std::endl;
  }

  std::map<std::string, std::string> nrrd_dict;
  while(std::getline(nhdr, line)){
    // Skip comments.
    if(line.compare(0, 1, "#")==0){
      continue;
    }
      
    int delimiter_pos = line.find(":");
    std::string key = line.substr(0, delimiter_pos);

    int offset = line.substr(delimiter_pos+1).find_first_not_of(" \t");
    std::string value = line.substr(delimiter_pos+1+offset);
      
    if(key=="space directions"){	
      int loc = value.find("(");
      int length = value.substr(loc+1).find(",");
      resolution *= atof(value.substr(loc+1, length).c_str());
      break;
    }else if(key=="space units"){
      int loc = value.find('"');
      int length = value.substr(loc+1).find('"');
      if(value.substr(loc+1, length)=="µm"){
        resolution *= 1.0e-6;
      }else{
        std::cerr<<"WARNING: units not handled. Set manually."<<std::endl;
      }
      break;
    }
  }

  nhdr.close();

  return resolution;
}

   
void read_tarantula_mesh_file(std::string filename, std::string nhdr_filename,
                              bool toggle_material,
                              std::vector<double> &xyz,
                              std::vector<int> &tets){

  double resolution = 1.0;
  
  // Get meta-data from NHDR file if specified.
  if(!nhdr_filename.empty()){
    resolution = read_resolution_from_nhdr(nhdr_filename);
  }

  // Read Tarantuala file
  std::ifstream infile;
  infile.open(filename.c_str());
  
  // Read header
  std::string throwaway;
  std::getline(infile, throwaway); // line 1
  std::getline(infile, throwaway); // line 2
  int NNodes;
  infile>>NNodes;
  
  // Read vertices
  xyz.resize(NNodes*3);
  for(int i=0;i<NNodes;i++){
    infile>>xyz[i*3];
    infile>>xyz[i*3+1]; 
    infile>>xyz[i*3+2];
  }

  // Rescale if necessary.
  if(resolution!=1.0){
    for(int i=0;i<NNodes*3;i++){
      xyz[i]*=resolution;
    }
  }

  // throwaway trash
  std::getline(infile, throwaway); 
  std::getline(infile, throwaway); 
  std::getline(infile, throwaway);

  // Read elements
  int NTetra, nloc;
  infile>>NTetra;
  tets.resize(NTetra*4);
  for(int i=0;i<NTetra;i++){
    infile>>nloc;
    assert(nloc==4);
    infile>>tets[i*4];
    infile>>tets[i*4+1];
    infile>>tets[i*4+2];
    infile>>tets[i*4+3];
  }

  // Read materials
  std::vector< std::vector<size_t> > materials;
  while(infile.good()){
    // Stream through file until we find material data.
    std::getline(infile, throwaway); 
    if(throwaway.substr(0, 4)!="mat1" && throwaway.substr(0, 4)!="mat2")
      continue;

    // Junk next line.
    std::getline(infile, throwaway);
    
    // Get number of cells of this material
    size_t cnt;
    infile>>cnt;
    std::vector<size_t> cells(cnt);
    for(int i=0;i<cnt;i++)
      infile>>cells[i];
    materials.push_back(cells);
  }

  // Junking the rest of the file.
  infile.close();

  if(materials.size()>1){
    assert(materials.size()==2);
    size_t select=0;
    if(toggle_material)
      select = 1;

    // Create the mask.
    std::vector<bool> mask(NTetra, false);
    for(std::vector<size_t>::const_iterator it=materials[select].begin();it!=materials[select].end();++it){
      mask[*it] = true;
    }

    // Turn off masked tets.
    for(size_t i=0;i<NTetra;i++){
      if(mask[i])
        tets[i*4] = -1;
    }
  }
}

void read_vtk_mesh_file(std::string filename, std::string nhdr_filename,
                       std::vector<double> &xyz,
                       std::vector<int> &tets){

  double resolution = 1.0;
  
  // Get meta-data from NHDR file if specified.
  if(!nhdr_filename.empty()){
    resolution = read_resolution_from_nhdr(nhdr_filename);
  } 

  // Read VTK file
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->DeepCopy(reader->GetOutput());

  // Read vertices
  int NNodes = pd->GetNumberOfPoints();
  xyz.resize(NNodes*3);
  for(int i=0;i<NNodes;i++){
    pd->GetPoints()->GetPoint(i, xyz.data()+i*3);
  }

  // Rescale if necessary.
  if(resolution!=1.0){
    for(int i=0;i<NNodes*3;i++){
      xyz[i]*=resolution;
    }
  }

  // Read facets - cleaver writes out tetrahedra by writing 4 facets..
  int nfacets = pd->GetNumberOfCells();
  for(int i=0;i<nfacets;i+=4){
    for(int j=0;j<4;j++){
      int cell_type = pd->GetCell(i+j)->GetCellType();
      if(cell_type!=VTK_TRIANGLE){
        std::cerr<<"ERROR("<<__FILE__<<"): unsupported element type.\n";
        exit(-1);
      }
    }

    vtkCell *cell0 = pd->GetCell(i);
    for(int j=0;j<3;j++){
      tets.push_back(cell0->GetPointId(j));
    }
     
    vtkCell *cell1 = pd->GetCell(i+1);
    tets.push_back(cell1->GetPointId(2));
  }
}

double volume(const double *x0, const double *x1, const double *x2, const double *x3){

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


