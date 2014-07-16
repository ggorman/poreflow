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
#include <map>
#include <string>
#include <vector>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <getopt.h>

#include "writers.h"

void usage(char *cmd){
  std::cerr<<"\nUsage: "<<cmd<<" [options ...] [Tarantula mesh file]\n"
           <<"\nOptions:\n"
           <<" -h, --help\n\tHelp! Prints this message.\n"
           <<" -v, --verbose\n\tVerbose output.\n"
           <<" -t, --toggle\n\tToggle the material selection for the mesh.\n";
  return;
}

int parse_arguments(int argc, char **argv,
                    std::string &filename, bool &verbose, bool &toggle_material){

  // Set defaults
  verbose = false;
  toggle_material = false;

  if(argc==1){
    usage(argv[0]);
    exit(0);
  }

  struct option longOptions[] = {
    {"help",    0, 0, 'h'},
    {"verbose", 0, 0, 'v'},
    {"toggle",  0, 0, 't'},
    {0, 0, 0, 0}
  };

  int optionIndex = 0;
  int verbosity = 0;
  int c;
  const char *shortopts = "hvt";

  // Set opterr to nonzero to make getopt print error messages
  opterr=1;
  while (true){
    c = getopt_long(argc, argv, shortopts, longOptions, &optionIndex);
    
    if (c == -1) break;
    
    switch (c){
    case 'h':
      usage(argv[0]);
      break;
    case 'v':
      verbose = true;
      break;
    case 't':
      toggle_material = true;
      break;
    case '?':
      // missing argument only returns ':' if the option string starts with ':'
      // but this seems to stop the printing of error messages by getopt?
      std::cerr<<"ERROR: unknown option or missing argument\n";
      usage(argv[0]);
      exit(-1);
    case ':':
      std::cerr<<"ERROR: missing argument\n";
      usage(argv[0]);
      exit(-1);
    default:
      // unexpected:
      std::cerr<<"ERROR: getopt returned unrecognized character code\n";
      exit(-1);
    }
  }

  filename = std::string(argv[argc-1]);
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.good()){
    std::cerr<<"ERROR: Cannot read file: "<<filename<<std::endl;
    usage(argv[0]);
    exit(0);
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

int read_tarantula_mesh_file(std::string filename, bool toggle_material,
                             std::vector<double> &xyz,
                             std::vector<int> &tets){
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
    if(throwaway.substr(0, 3)!="mat")
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

  // I have no idea what mat0 is.
  // mat1 appears to be the zero valued voxels
  // mat2 appears to be the one valued voxels
  assert(materials.size()==3);

  // Junking the rest of the file.
  infile.close();

  size_t select=2;
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

int create_domain(std::vector<double> &xyz,
                  std::vector<int> &tets, 
                  std::vector<int> &facets,
                  std::vector<int> &facet_ids){
  
  // Create node-element adjancy list.
  size_t NNodes = xyz.size()/3;

  std::vector< std::set<int> > NEList(NNodes);
  int NTetra = tets.size()/4;
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    double v = volume(&xyz[3*tets[i*4]], &xyz[3*tets[i*4+1]], &xyz[3*tets[i*4+2]], &xyz[3*tets[i*4+3]]);
    if(v<0){
      std::cerr<<"WARNING: Inverted element encountered. Deleting and continuing.\n";
      for(int j=0;j<4;j++)
	std::cerr<<i<<", "<<tets[i*4+j]<<": "<<xyz[3*tets[i*4+j]]<<" "<<xyz[3*tets[i*4+j]+1]<<" "<<xyz[3*tets[i*4+j]+2]<<std::endl;

      tets[i*4] = -1;
    }else{
      for(int j=0;j<4;j++){
	NEList[tets[i*4+j]].insert(i);
      }
    }
  }

  // Create element-element adjancy list
  std::vector<int> EEList(NTetra*4, -1);
  for(int i=0;i<NTetra;i++){
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
	double mean_x = (xyz[facet[0]*3]+
			 xyz[facet[1]*3]+
			 xyz[facet[2]*3])/3.0;
	
	if(fabs(mean_x-bbox[0])<eta){
	  front0.insert(front0.end(), i);
	}else if(fabs(mean_x-bbox[1])<eta){
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
  for(std::map<int, int>::iterator it=renumbering.begin();it!=renumbering.end();++it){
    it->second = cnt++;
    
    xyz_new.push_back(xyz[(it->first)*3]);
    xyz_new.push_back(xyz[(it->first)*3+1]);
    xyz_new.push_back(xyz[(it->first)*3+2]);
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

int main(int argc, char **argv){
  std::string filename;
  bool verbose, toggle_material;
  parse_arguments(argc, argv, filename, verbose, toggle_material);

  std::vector<double> xyz;
  std::vector<int> tets;
  std::string basename = filename.substr(0, filename.size()-4);
  
  if(verbose)
    std::cout<<"INFO: Reading "<<filename<<std::endl;
  read_tarantula_mesh_file(filename, toggle_material, xyz, tets);
  if(verbose)
    std::cout<<"INFO: Finished reading "<<filename<<std::endl;
  
  // Generate facets and trim disconnnected parts of the domain.
  std::vector<int> facets, facet_ids;
  if(verbose){
    std::cout<<"INFO: Create the active domain."<<filename<<std::endl;
    write_vtk_file(basename+"_original", xyz, tets, facets, facet_ids);
  }  

  create_domain(xyz, tets, facets, facet_ids);
  
  if(verbose) 
    std::cout<<"INFO: Active domain created."<<filename<<std::endl;
  
  if(verbose){
    std::cout<<"INFO: Writing out mesh."<<filename<<std::endl;
    write_vtk_file(basename, xyz, tets, facets, facet_ids);
  }

  write_gmsh_file(basename, xyz, tets, facets, facet_ids);
  if(verbose)
    std::cout<<"INFO: Finished."<<filename<<std::endl;

  return 0;
}
